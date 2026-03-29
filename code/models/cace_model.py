"""
cace_model.py — Cross-Attention Communication Encoder (CACE)

Models ligand-receptor-mediated communication between T cell exhaustion states
(Tex) and tumour-associated macrophage (TAM) subtypes.  Given a communication
tensor it produces:
  1. A d_model-dimensional embedding vector
  2. A prognostic risk score for Cox PH survival modelling
  3. An ICB response logit for binary classification

Input shape
-----------
x : (batch, n_tex * n_tam, n_lr)
    Flattened cell-state-pair axis carrying n_lr ligand-receptor activity
    scores for each (Tex, TAM) combination.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import torch
import torch.nn as nn
from torch import Tensor


# ---------------------------------------------------------------------------
# Sub-modules
# ---------------------------------------------------------------------------

class CommunicationEncoder(nn.Module):
    """Project each cell-state-pair's L-R vector into d_model space.

    Architecture:
        Linear(n_lr, d_model) → LayerNorm → GELU → Dropout
        → Linear(d_model, d_model) → LayerNorm
    """

    def __init__(self, n_lr: int, d_model: int, dropout: float = 0.3) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_lr, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, d_model),
            nn.LayerNorm(d_model),
        )

    def forward(self, x: Tensor) -> Tensor:
        """
        Parameters
        ----------
        x : (..., n_lr)

        Returns
        -------
        Tensor : (..., d_model)
        """
        return self.net(x)


class CrossAttentionBlock(nn.Module):
    """One layer of multi-head cross-attention where Tex states attend to TAMs.

    Uses *pre-norm* style residual connections.

    Architecture
    ------------
    Attention sub-layer:
        LN(query) + LN(key/value) → MultiheadAttention → residual on query
    FFN sub-layer:
        LN(x) → Linear(d_model, 4*d_model) → GELU → Dropout
               → Linear(4*d_model, d_model) → Dropout → residual
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        dropout: float = 0.3,
    ) -> None:
        super().__init__()
        self.norm_q = nn.LayerNorm(d_model)
        self.norm_kv = nn.LayerNorm(d_model)
        self.attn = nn.MultiheadAttention(
            embed_dim=d_model,
            num_heads=n_heads,
            dropout=dropout,
            batch_first=True,
        )
        self.norm_ffn = nn.LayerNorm(d_model)
        self.ffn = nn.Sequential(
            nn.Linear(d_model, d_model * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 4, d_model),
            nn.Dropout(dropout),
        )

    def forward(
        self, query: Tensor, key_value: Tensor
    ) -> Tuple[Tensor, Tensor]:
        """
        Parameters
        ----------
        query     : (batch, n_tex, d_model)  — Tex representations
        key_value : (batch, n_tam, d_model)  — TAM representations

        Returns
        -------
        out         : (batch, n_tex, d_model)  — updated Tex representations
        attn_weights: (batch, n_heads, n_tex, n_tam)
        """
        # Pre-norm cross-attention with residual on query
        q_norm = self.norm_q(query)
        kv_norm = self.norm_kv(key_value)
        attn_out, attn_weights = self.attn(
            q_norm, kv_norm, kv_norm, average_attn_weights=False
        )
        x = query + attn_out

        # Pre-norm FFN with residual
        x = x + self.ffn(self.norm_ffn(x))

        return x, attn_weights  # attn_weights: (batch, n_heads, n_tex, n_tam)


# ---------------------------------------------------------------------------
# Main model
# ---------------------------------------------------------------------------

class CACE(nn.Module):
    """Cross-Attention Communication Encoder.

    Parameters
    ----------
    n_tex    : number of T cell exhaustion states
    n_tam    : number of TAM subtypes
    n_lr     : number of ligand-receptor pairs
    d_model  : hidden dimension (default 128)
    n_heads  : attention heads (default 8)
    n_layers : number of cross-attention blocks (default 2)
    dropout  : dropout rate (default 0.3)
    """

    def __init__(
        self,
        n_tex: int,
        n_tam: int,
        n_lr: int,
        d_model: int = 128,
        n_heads: int = 8,
        n_layers: int = 2,
        dropout: float = 0.3,
    ) -> None:
        super().__init__()
        self.n_tex = n_tex
        self.n_tam = n_tam
        self.n_lr = n_lr
        self.d_model = d_model

        # L-R pair encoder shared across all cell-state pairs
        self.lr_encoder = CommunicationEncoder(n_lr, d_model, dropout)

        # Learnable positional embeddings for Tex and TAM identities
        self.tex_pos = nn.Embedding(n_tex, d_model)
        self.tam_pos = nn.Embedding(n_tam, d_model)

        # Cross-attention layers (Tex attends to TAM)
        self.cross_attn_layers = nn.ModuleList(
            [CrossAttentionBlock(d_model, n_heads, dropout) for _ in range(n_layers)]
        )

        # Pool: flatten n_tex representations into single embedding vector
        self.pool = nn.Sequential(
            nn.Linear(d_model * n_tex, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
        )

        # Task heads
        # Prognostic head: scalar risk score for Cox PH loss
        self.prognostic_head = nn.Linear(d_model, 1)

        # ICB response head: binary classification logit
        self.icb_head = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

        self._init_weights()

    # ------------------------------------------------------------------
    # Weight initialisation
    # ------------------------------------------------------------------

    def _init_weights(self) -> None:
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Embedding):
                nn.init.normal_(module.weight, mean=0.0, std=0.02)

    # ------------------------------------------------------------------
    # Core encode method
    # ------------------------------------------------------------------

    def encode(self, x: Tensor) -> Tuple[Tensor, List[Tensor]]:
        """Encode communication tensor into an embedding.

        Parameters
        ----------
        x : (batch, n_tex * n_tam, n_lr)

        Returns
        -------
        embedding    : (batch, d_model)
        attn_weights : list of (batch, n_heads, n_tex, n_tam), one per layer
        """
        batch = x.shape[0]

        # Reshape: (batch, n_tex, n_tam, n_lr)
        x = x.view(batch, self.n_tex, self.n_tam, self.n_lr)

        # Encode each (Tex, TAM) pair's L-R vector → (batch, n_tex, n_tam, d_model)
        encoded = self.lr_encoder(x)

        # Build Tex query representation: mean over TAM dim + positional embedding
        # Shape: (batch, n_tex, d_model)
        tex_repr = encoded.mean(dim=2)
        tex_idx = torch.arange(self.n_tex, device=x.device)
        tex_repr = tex_repr + self.tex_pos(tex_idx).unsqueeze(0)

        # Build TAM key/value representation: mean over Tex dim + positional embedding
        # Shape: (batch, n_tam, d_model)
        tam_repr = encoded.mean(dim=1)
        tam_idx = torch.arange(self.n_tam, device=x.device)
        tam_repr = tam_repr + self.tam_pos(tam_idx).unsqueeze(0)

        # Apply cross-attention layers: Tex attends to TAM
        attn_weights_list: List[Tensor] = []
        query = tex_repr
        for layer in self.cross_attn_layers:
            query, attn_w = layer(query, tam_repr)
            attn_weights_list.append(attn_w)  # (batch, n_heads, n_tex, n_tam)

        # Pool: flatten (n_tex, d_model) → d_model
        # (batch, n_tex * d_model)
        flat = query.reshape(batch, self.n_tex * self.d_model)
        embedding = self.pool(flat)  # (batch, d_model)

        return embedding, attn_weights_list

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    def forward(
        self,
        x: Tensor,
        task: str = "both",
    ) -> dict:
        """Run the model for one or both tasks.

        Parameters
        ----------
        x    : (batch, n_tex * n_tam, n_lr)
        task : "both" | "prognostic" | "icb"

        Returns
        -------
        dict with keys:
            "embedding"   : (batch, d_model)        always present
            "attn_weights": list[(batch, n_heads, n_tex, n_tam)]  always present
            "risk_score"  : (batch,)                if task in {"both","prognostic"}
            "icb_logit"   : (batch,)                if task in {"both","icb"}
        """
        embedding, attn_weights = self.encode(x)

        out: dict = {
            "embedding": embedding,
            "attn_weights": attn_weights,
        }

        if task in ("both", "prognostic"):
            risk = self.prognostic_head(embedding).squeeze(-1)  # (batch,)
            out["risk_score"] = risk

        if task in ("both", "icb"):
            icb = self.icb_head(embedding).squeeze(-1)  # (batch,)
            out["icb_logit"] = icb

        return out

    # ------------------------------------------------------------------
    # Convenience method
    # ------------------------------------------------------------------

    def get_communication_score(self, x: Tensor) -> Tensor:
        """Return a scalar prognostic risk score per sample.

        Parameters
        ----------
        x : (batch, n_tex * n_tam, n_lr)

        Returns
        -------
        Tensor : (batch,)
        """
        out = self.forward(x, task="prognostic")
        return out["risk_score"]
