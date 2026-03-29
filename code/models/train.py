"""
train.py — Training pipeline for the CACE model.

Supports two tasks:
  - "prognostic" : Cox proportional-hazards survival modelling
  - "icb"        : binary ICB response classification (BCE)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Tuple

import numpy as np
import torch
import torch.nn as nn
from torch import Tensor
from torch.utils.data import DataLoader

# ---------------------------------------------------------------------------
# Path bootstrap — allow running as a top-level script
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from code.config import CACE_PARAMS, CHECKPOINT_DIR, PROCESSED_DIR, SEED
from code.models.cace_model import CACE


# ---------------------------------------------------------------------------
# Loss functions
# ---------------------------------------------------------------------------

def cox_ph_loss(risk_scores: Tensor, times: Tensor, events: Tensor) -> Tensor:
    """Negative partial log-likelihood for the Cox proportional-hazards model.

    Parameters
    ----------
    risk_scores : (N,)  — raw (log) risk scores from the prognostic head
    times       : (N,)  — observed survival / censoring times
    events      : (N,)  — event indicators (1 = event occurred, 0 = censored)

    Returns
    -------
    Tensor : scalar loss (mean over uncensored observations)

    Notes
    -----
    Breslow approximation for ties.  Sorting by *descending* time means the
    risk set for observation i is {j : j >= i in the sorted order}.
    """
    # Sort by descending observed time
    order = torch.argsort(times, descending=True)
    risk_scores = risk_scores[order]
    events = events[order]

    # log-sum-exp over the risk set (cumulative from the top)
    log_cumsum_exp = torch.logcumsumexp(risk_scores, dim=0)

    # Partial log-likelihood contributions for event cases
    uncensored = events.bool()
    loss = -(risk_scores[uncensored] - log_cumsum_exp[uncensored]).mean()
    return loss


# ---------------------------------------------------------------------------
# Training loop helpers
# ---------------------------------------------------------------------------

def train_epoch(
    model: CACE,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    task: str,
    device: torch.device,
) -> float:
    """Run one full training epoch.

    Parameters
    ----------
    model     : CACE model
    loader    : DataLoader yielding batches appropriate for *task*
                  "prognostic" batches → (x, times, events)
                  "icb"        batches → (x, labels)
    optimizer : optimiser instance
    task      : "prognostic" or "icb"
    device    : torch device

    Returns
    -------
    float : mean batch loss over the epoch
    """
    model.train()
    total_loss = 0.0
    n_batches = 0

    for batch in loader:
        optimizer.zero_grad()

        if task == "prognostic":
            x, times, events = batch
            x = x.to(device)
            times = times.to(device)
            events = events.to(device)
            out = model(x, task="prognostic")
            loss = cox_ph_loss(out["risk_score"], times, events)

        elif task == "icb":
            x, labels = batch
            x = x.to(device)
            labels = labels.float().to(device)
            out = model(x, task="icb")
            loss = nn.functional.binary_cross_entropy_with_logits(
                out["icb_logit"], labels
            )

        else:
            raise ValueError(f"Unknown task: {task!r}. Choose 'prognostic' or 'icb'.")

        loss.backward()
        # Gradient clipping
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        total_loss += loss.item()
        n_batches += 1

    return total_loss / max(n_batches, 1)


def evaluate(
    model: CACE,
    loader: DataLoader,
    task: str,
    device: torch.device,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Evaluate the model on a data loader without gradient tracking.

    Parameters
    ----------
    model  : CACE model
    loader : DataLoader (same batch format as for train_epoch)
    task   : "prognostic" or "icb"
    device : torch device

    Returns
    -------
    loss        : float — mean loss over all batches
    predictions : np.ndarray — concatenated model outputs
                    "prognostic" → risk scores
                    "icb"        → ICB logits
    labels      : np.ndarray — concatenated ground-truth labels
                    "prognostic" → event indicators (1-D)
                    "icb"        → binary response labels
    """
    model.eval()
    total_loss = 0.0
    n_batches = 0
    all_preds: list[np.ndarray] = []
    all_labels: list[np.ndarray] = []

    with torch.no_grad():
        for batch in loader:
            if task == "prognostic":
                x, times, events = batch
                x = x.to(device)
                times = times.to(device)
                events = events.to(device)
                out = model(x, task="prognostic")
                loss = cox_ph_loss(out["risk_score"], times, events)
                all_preds.append(out["risk_score"].cpu().numpy())
                all_labels.append(events.cpu().numpy())

            elif task == "icb":
                x, labels = batch
                x = x.to(device)
                labels = labels.float().to(device)
                out = model(x, task="icb")
                loss = nn.functional.binary_cross_entropy_with_logits(
                    out["icb_logit"], labels
                )
                all_preds.append(out["icb_logit"].cpu().numpy())
                all_labels.append(labels.cpu().numpy())

            else:
                raise ValueError(f"Unknown task: {task!r}.")

            total_loss += loss.item()
            n_batches += 1

    predictions = np.concatenate(all_preds, axis=0) if all_preds else np.array([])
    labels_arr = np.concatenate(all_labels, axis=0) if all_labels else np.array([])
    mean_loss = total_loss / max(n_batches, 1)

    return mean_loss, predictions, labels_arr


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[train] device: {device}")

    # ------------------------------------------------------------------
    # Load communication tensor and metadata
    # ------------------------------------------------------------------
    comm_dir = PROCESSED_DIR / "communication"

    tensor_path = comm_dir / "communication_tensor.pt"
    meta_path = comm_dir / "communication_meta.json"

    if not tensor_path.exists() or not meta_path.exists():
        print(
            f"[train] Communication tensor not found at {comm_dir}.\n"
            "        Awaiting TCGA deconvolution data — nothing to train on yet."
        )
        _init_model_stub(device)
        return

    comm_tensor = torch.load(tensor_path, map_location="cpu")
    with open(meta_path) as fh:
        meta = json.load(fh)

    print(f"[train] Loaded communication tensor: {tuple(comm_tensor.shape)}")
    print(f"[train] Metadata: {meta}")

    n_tex = meta["n_tex"]
    n_tam = meta["n_tam"]
    n_lr = meta["n_lr"]

    # ------------------------------------------------------------------
    # Initialise CACE model
    # ------------------------------------------------------------------
    model = CACE(
        n_tex=n_tex,
        n_tam=n_tam,
        n_lr=n_lr,
        d_model=CACE_PARAMS["d_model"],
        n_heads=CACE_PARAMS["n_heads"],
        n_layers=CACE_PARAMS["n_layers"],
        dropout=CACE_PARAMS["dropout"],
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[train] CACE trainable parameters: {n_params:,}")

    # ------------------------------------------------------------------
    # Optimiser
    # ------------------------------------------------------------------
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=CACE_PARAMS["lr"],
        weight_decay=CACE_PARAMS["weight_decay"],
    )

    # ------------------------------------------------------------------
    # Save model config
    # ------------------------------------------------------------------
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    config_path = CHECKPOINT_DIR / "cace_config.json"
    model_config = {
        "n_tex": n_tex,
        "n_tam": n_tam,
        "n_lr": n_lr,
        **CACE_PARAMS,
    }
    with open(config_path, "w") as fh:
        json.dump(model_config, fh, indent=2)
    print(f"[train] Model config saved to {config_path}")

    # ------------------------------------------------------------------
    # Actual training data not yet available
    # ------------------------------------------------------------------
    print(
        "[train] Awaiting TCGA deconvolution data — "
        "training loop will run once bulk deconvolution is complete."
    )


def _init_model_stub(device: torch.device) -> None:
    """Initialise a CACE model with placeholder dimensions and report parameters."""
    # Placeholder dimensions matching typical HCC Tex/TAM landscape
    n_tex, n_tam, n_lr = 3, 6, 50

    model = CACE(
        n_tex=n_tex,
        n_tam=n_tam,
        n_lr=n_lr,
        d_model=CACE_PARAMS["d_model"],
        n_heads=CACE_PARAMS["n_heads"],
        n_layers=CACE_PARAMS["n_layers"],
        dropout=CACE_PARAMS["dropout"],
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[train] Stub CACE model — trainable parameters: {n_params:,}")

    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    config_path = CHECKPOINT_DIR / "cace_config.json"
    model_config = {
        "n_tex": n_tex,
        "n_tam": n_tam,
        "n_lr": n_lr,
        **CACE_PARAMS,
    }
    with open(config_path, "w") as fh:
        json.dump(model_config, fh, indent=2)
    print(f"[train] Stub model config saved to {config_path}")


if __name__ == "__main__":
    main()
