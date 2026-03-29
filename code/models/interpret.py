"""
interpret.py — Interpretation utilities for the CACE model.

Provides three public functions:

    extract_attention_weights   — summarise per-layer, per-head attention
                                  into a (tex_state × tam_type) DataFrame
    rank_lr_pairs_shap          — rank L-R pairs by SHAP KernelExplainer
                                  importance
    get_top_communication_axes  — combine attention and SHAP rankings into
                                  a single interpretable DataFrame
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd
import shap
import torch


# ---------------------------------------------------------------------------
# 1. Attention-weight extraction
# ---------------------------------------------------------------------------

def extract_attention_weights(
    model: "torch.nn.Module",
    tensor: np.ndarray,
    metadata: Dict[str, List[str]],
) -> pd.DataFrame:
    """Summarise cross-attention weights as a (tex_state × tam_type) DataFrame.

    Parameters
    ----------
    model    : trained CACE model
    tensor   : np.ndarray, shape (n_samples, n_tex*n_tam, n_lr)
               Communication tensor – one row per patient/sample.
    metadata : dict with keys
                   "tex_states" : list of str, length n_tex
                   "tam_types"  : list of str, length n_tam

    Returns
    -------
    pd.DataFrame, shape (n_tex, n_tam)
        Index   = tex_states, columns = tam_types.
        Each cell holds the average attention weight that the corresponding
        Tex state places on the corresponding TAM type, averaged across all
        samples, all attention heads, and all cross-attention layers.
    """
    tex_states: List[str] = metadata["tex_states"]
    tam_types: List[str] = metadata["tam_types"]

    device = next(model.parameters()).device
    x = torch.tensor(tensor, dtype=torch.float32, device=device)

    model.eval()
    with torch.no_grad():
        out = model(x)

    # attn_weights: list of (batch, n_heads, n_tex, n_tam), one per layer
    attn_weights_list = out["attn_weights"]

    # Stack layers → (n_layers, batch, n_heads, n_tex, n_tam)
    stacked = torch.stack(attn_weights_list, dim=0)

    # Average over layers, samples, and heads → (n_tex, n_tam)
    avg = stacked.mean(dim=(0, 1, 2))  # mean over layers, batch, heads

    avg_np = avg.cpu().numpy()  # (n_tex, n_tam)

    df = pd.DataFrame(avg_np, index=tex_states, columns=tam_types)
    df.index.name = "tex_state"
    df.columns.name = "tam_type"

    return df


# ---------------------------------------------------------------------------
# 2. L-R pair ranking via SHAP
# ---------------------------------------------------------------------------

def rank_lr_pairs_shap(
    model: "torch.nn.Module",
    tensor: np.ndarray,
    metadata: Dict[str, List[str]],
    n_background: int = 50,
) -> pd.DataFrame:
    """Rank L-R pairs by their SHAP importance for the prognostic risk score.

    Parameters
    ----------
    model        : trained CACE model
    tensor       : np.ndarray, shape (n_samples, n_tex*n_tam, n_lr)
    metadata     : dict with keys "tex_states", "tam_types";
                   also optionally "lr_pairs" (list of str, length n_lr)
    n_background : number of background samples used by KernelExplainer

    Returns
    -------
    pd.DataFrame with columns ["lr_pair", "shap_importance"]
        Sorted in descending order of mean absolute SHAP value across all
        samples and (tex, tam) state-pair positions.
    """
    n_samples, n_pairs, n_lr = tensor.shape

    # Derive L-R pair names
    if "lr_pairs" in metadata:
        lr_names: List[str] = list(metadata["lr_pairs"])
    else:
        lr_names = [f"lr_{i}" for i in range(n_lr)]

    # Flatten tensor to 2D for SHAP: (n_samples, n_pairs * n_lr)
    flat = tensor.reshape(n_samples, n_pairs * n_lr).astype(np.float32)

    device = next(model.parameters()).device

    def predict_fn(x_flat: np.ndarray) -> np.ndarray:
        """Map (n, n_pairs*n_lr) → (n,) risk scores."""
        x_tensor = torch.tensor(
            x_flat.reshape(-1, n_pairs, n_lr),
            dtype=torch.float32,
            device=device,
        )
        model.eval()
        with torch.no_grad():
            out = model(x_tensor, task="prognostic")
        return out["risk_score"].cpu().numpy()

    # Select background samples (subsample if necessary)
    rng = np.random.default_rng(42)
    n_bg = min(n_background, n_samples)
    bg_idx = rng.choice(n_samples, size=n_bg, replace=False)
    background = flat[bg_idx]

    explainer = shap.KernelExplainer(predict_fn, background)
    # shap_values: (n_samples, n_pairs * n_lr)
    shap_values = explainer.shap_values(flat, silent=True)
    shap_values = np.array(shap_values)  # ensure ndarray

    # Reshape back to 3D: (n_samples, n_pairs, n_lr)
    shap_3d = shap_values.reshape(n_samples, n_pairs, n_lr)

    # Average |SHAP| over samples and state pairs → (n_lr,)
    importance = np.abs(shap_3d).mean(axis=(0, 1))

    df = pd.DataFrame({
        "lr_pair": lr_names,
        "shap_importance": importance,
    })
    df = df.sort_values("shap_importance", ascending=False).reset_index(drop=True)

    return df


# ---------------------------------------------------------------------------
# 3. Top communication axes
# ---------------------------------------------------------------------------

def get_top_communication_axes(
    attn_df: pd.DataFrame,
    lr_ranking: pd.DataFrame,
    top_k: int = 10,
) -> pd.DataFrame:
    """Identify the most important (L-R pair, Tex state, TAM type) triplets.

    For each of the top-k L-R pairs (by SHAP importance), this function
    identifies which (tex_state, tam_type) combination exhibits the highest
    average attention weight and returns a combined summary.

    Parameters
    ----------
    attn_df    : pd.DataFrame (n_tex × n_tam) from extract_attention_weights
                 Index   = tex_states, columns = tam_types.
    lr_ranking : pd.DataFrame with columns ["lr_pair", "shap_importance"]
                 from rank_lr_pairs_shap, already sorted descending.
    top_k      : number of top L-R pairs to include (default 10)

    Returns
    -------
    pd.DataFrame with columns:
        ["lr_pair", "tex_state", "tam_type", "attention", "shap_importance"]
    One row per top L-R pair, with the (tex_state, tam_type) pair that has
    the highest attention weight for that L-R pair.
    """
    top_lr = lr_ranking.head(top_k).reset_index(drop=True)

    rows = []
    for _, row in top_lr.iterrows():
        # Find the (tex_state, tam_type) with maximum attention
        # attn_df has tex_states as index and tam_types as columns
        max_idx = attn_df.values.argmax()
        tex_idx, tam_idx = np.unravel_index(max_idx, attn_df.shape)
        best_tex = attn_df.index[tex_idx]
        best_tam = attn_df.columns[tam_idx]
        best_attn = float(attn_df.iloc[tex_idx, tam_idx])

        rows.append({
            "lr_pair": row["lr_pair"],
            "tex_state": best_tex,
            "tam_type": best_tam,
            "attention": best_attn,
            "shap_importance": float(row["shap_importance"]),
        })

    result = pd.DataFrame(rows, columns=[
        "lr_pair", "tex_state", "tam_type", "attention", "shap_importance"
    ])
    return result
