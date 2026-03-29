"""
deconvolution.py — Bulk-to-communication transfer via NNLS deconvolution.

This module bridges scRNA-seq–derived communication references to bulk
RNA-seq cohorts (TCGA-LIHC, GSE285963).  The pipeline is:

    1. build_signature_matrix  — derive a gene × cell-state signature from
                                  the scRNA-seq atlas (CD8+ T cells + TAMs)
    2. estimate_fractions_nnls — project bulk samples onto the signature via
                                  non-negative least squares to obtain cell-
                                  state fraction estimates
    3. build_pseudo_communication — weight the scRNA-seq communication
                                     reference by the estimated fractions to
                                     produce a per-sample communication tensor
                                     that CACE can consume

A helper function build_comm_reference converts a scRNA-seq communication
tensor + metadata dict into the (tex, tam, lr_pair) → strength dict that
build_pseudo_communication expects.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.optimize import nnls


# ---------------------------------------------------------------------------
# 1. Signature matrix
# ---------------------------------------------------------------------------

def build_signature_matrix(
    cd8_adata: AnnData,
    tam_adata: AnnData,
    n_markers: int = 200,
) -> pd.DataFrame:
    """Build a gene × cell-state signature matrix from scRNA-seq AnnData objects.

    The function concatenates the CD8+ T cell and TAM AnnData objects, runs
    Wilcoxon rank-sum differential expression per cell state, and returns the
    mean expression of the top ``n_markers`` marker genes for each state.

    Parameters
    ----------
    cd8_adata : AnnData
        CD8+ T cell AnnData.  Must have ``obs["exhaustion_state"]`` with values
        such as "progenitor", "transitional", "terminal".
    tam_adata : AnnData
        TAM AnnData.  Must have ``obs["tam_subtype"]`` with values such as
        "SPP1", "TREM2", "FOLR2", etc.
    n_markers : int
        Maximum number of marker genes to include per cell state.  The actual
        number may be smaller if fewer significant markers are found.

    Returns
    -------
    pd.DataFrame
        Shape ``(n_genes, n_states)``.  Index: gene names.  Columns: cell-state
        labels (exhaustion states then TAM subtypes).
    """
    # --- build a unified "deconv_label" column ---
    cd8_copy = cd8_adata.copy()
    tam_copy = tam_adata.copy()

    cd8_copy.obs["deconv_label"] = cd8_copy.obs["exhaustion_state"].astype(str)
    tam_copy.obs["deconv_label"] = tam_copy.obs["tam_subtype"].astype(str)

    # Concatenate; batch_key prevents obs-name collisions
    combined = sc.concat(
        [cd8_copy, tam_copy],
        join="inner",   # only genes present in both
        label="source",
        keys=["cd8", "tam"],
    )

    # Ensure raw counts are normalised for marker detection
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    # Rank genes per cell state (Wilcoxon)
    sc.tl.rank_genes_groups(
        combined,
        groupby="deconv_label",
        method="wilcoxon",
        n_genes=n_markers,
        use_raw=False,
    )

    # Collect the union of marker genes across all states
    rank_result = combined.uns["rank_genes_groups"]
    cell_states: List[str] = list(rank_result["names"].dtype.names)

    marker_genes: List[str] = []
    for state in cell_states:
        genes = list(rank_result["names"][state])
        marker_genes.extend(genes)

    # Deduplicate while preserving order
    seen: set = set()
    unique_markers: List[str] = []
    for g in marker_genes:
        if g not in seen:
            seen.add(g)
            unique_markers.append(g)

    # Filter to genes actually in the AnnData
    available = set(combined.var_names)
    unique_markers = [g for g in unique_markers if g in available]

    # Build signature: mean expression per state for marker genes
    sig_dict: Dict[str, pd.Series] = {}
    expr_df = pd.DataFrame(
        combined[:, unique_markers].X.toarray()
        if hasattr(combined[:, unique_markers].X, "toarray")
        else combined[:, unique_markers].X,
        index=combined.obs_names,
        columns=unique_markers,
    )
    labels = combined.obs["deconv_label"]

    for state in cell_states:
        mask = labels == state
        if mask.sum() == 0:
            sig_dict[state] = pd.Series(0.0, index=unique_markers)
        else:
            sig_dict[state] = expr_df.loc[mask].mean(axis=0)

    signature = pd.DataFrame(sig_dict)  # genes × states
    return signature


# ---------------------------------------------------------------------------
# 2. NNLS fraction estimation
# ---------------------------------------------------------------------------

def estimate_fractions_nnls(
    bulk_expr: pd.DataFrame,
    signature: pd.DataFrame,
) -> pd.DataFrame:
    """Estimate cell-state fractions in bulk RNA-seq samples via NNLS.

    For each bulk sample the function solves::

        min_{f >= 0}  || S f - b ||_2

    where *S* is the signature matrix and *b* is the bulk expression vector,
    then normalises *f* to sum to 1.

    Parameters
    ----------
    bulk_expr : pd.DataFrame
        Shape ``(n_genes, n_samples)``.  Index: gene names.  Columns: sample IDs.
    signature : pd.DataFrame
        Shape ``(n_genes, n_states)``.  Index: gene names.  Columns: cell-state
        labels (as returned by :func:`build_signature_matrix`).

    Returns
    -------
    pd.DataFrame
        Shape ``(n_samples, n_states)``.  Index: sample IDs.  Columns: cell-state
        labels.  Each row sums to 1 (or is all-zeros if NNLS returned a zero
        solution).
    """
    # Align genes: keep the intersection present in both matrices
    common_genes = signature.index.intersection(bulk_expr.index)
    if len(common_genes) == 0:
        raise ValueError(
            "No genes in common between bulk_expr and signature.  "
            "Check that both use the same gene identifier (e.g. HGNC symbol)."
        )

    S = signature.loc[common_genes].values.astype(float)  # (n_genes, n_states)
    B = bulk_expr.loc[common_genes].values.astype(float)   # (n_genes, n_samples)

    n_samples = B.shape[1]
    n_states = S.shape[1]
    fracs = np.zeros((n_samples, n_states), dtype=float)

    for i in range(n_samples):
        b = B[:, i]
        coef, _ = nnls(S, b)
        total = coef.sum()
        if total > 0:
            fracs[i] = coef / total
        # else: leave as zeros (degenerate sample)

    return pd.DataFrame(
        fracs,
        index=bulk_expr.columns,
        columns=signature.columns,
    )


# ---------------------------------------------------------------------------
# 3. Pseudo-communication tensor
# ---------------------------------------------------------------------------

def build_pseudo_communication(
    fractions: pd.DataFrame,
    comm_reference: Dict[Tuple[str, str, str], float],
    tex_states: List[str],
    tam_types: List[str],
    lr_pairs: List[str],
) -> np.ndarray:
    """Compute a per-sample pseudo-communication tensor from deconvolved fractions.

    For each sample *s*, each (Tex, TAM) pair index *p*, and each L-R pair *l*::

        pseudo_comm[s, p, l] = frac(s, tex) × frac(s, tam) × ref(tex, tam, lr)

    where the cell-pair axis is ordered as::

        p = tex_idx * n_tam + tam_idx

    Parameters
    ----------
    fractions : pd.DataFrame
        Shape ``(n_samples, n_states)``.  As returned by
        :func:`estimate_fractions_nnls`.  Must contain columns for all states
        listed in ``tex_states`` and ``tam_types``.
    comm_reference : dict
        Mapping ``(tex_state, tam_type, lr_pair) → float`` giving the reference
        communication strength derived from scRNA-seq data (see
        :func:`build_comm_reference`).
    tex_states : list of str
        Ordered list of T cell exhaustion state names.
    tam_types : list of str
        Ordered list of TAM subtype names.
    lr_pairs : list of str
        Ordered list of ligand-receptor pair identifiers.

    Returns
    -------
    np.ndarray
        Shape ``(n_samples, n_tex * n_tam, n_lr)``.
    """
    n_samples = len(fractions)
    n_tex = len(tex_states)
    n_tam = len(tam_types)
    n_lr = len(lr_pairs)
    n_pairs = n_tex * n_tam

    output = np.zeros((n_samples, n_pairs, n_lr), dtype=float)

    # Pre-build reference array: (n_tex, n_tam, n_lr)
    ref = np.zeros((n_tex, n_tam, n_lr), dtype=float)
    for ti, tex in enumerate(tex_states):
        for mi, tam in enumerate(tam_types):
            for li, lr in enumerate(lr_pairs):
                ref[ti, mi, li] = comm_reference.get((tex, tam, lr), 0.0)

    # Fraction arrays: (n_samples, n_tex) and (n_samples, n_tam)
    # Missing columns are silently treated as zero.
    tex_fracs = np.zeros((n_samples, n_tex), dtype=float)
    for ti, tex in enumerate(tex_states):
        if tex in fractions.columns:
            tex_fracs[:, ti] = fractions[tex].values

    tam_fracs = np.zeros((n_samples, n_tam), dtype=float)
    for mi, tam in enumerate(tam_types):
        if tam in fractions.columns:
            tam_fracs[:, mi] = fractions[tam].values

    # Compute outer product per sample and weight by reference
    # tex_fracs[:, ti] * tam_fracs[:, mi] → (n_samples,)
    for ti in range(n_tex):
        for mi in range(n_tam):
            pair_idx = ti * n_tam + mi
            # (n_samples,) outer broadcast with (n_lr,) → (n_samples, n_lr)
            weight = tex_fracs[:, ti] * tam_fracs[:, mi]  # (n_samples,)
            output[:, pair_idx, :] = np.outer(weight, ref[ti, mi])

    return output


# ---------------------------------------------------------------------------
# 4. Helper: build comm_reference from scRNA-seq tensor
# ---------------------------------------------------------------------------

def build_comm_reference(
    comm_tensor: np.ndarray,
    metadata: Dict,
) -> Dict[Tuple[str, str, str], float]:
    """Build a communication reference dict by averaging a scRNA-seq tensor.

    Parameters
    ----------
    comm_tensor : np.ndarray
        Shape ``(n_patients, n_tex * n_tam, n_lr)``.  Per-patient communication
        strengths from the scRNA-seq cohort.
    metadata : dict
        Must contain:

        ``"tex_states"`` : list of str — Tex state names in index order.
        ``"tam_types"``  : list of str — TAM subtype names in index order.
        ``"lr_pairs"``   : list of str — L-R pair names in index order.

    Returns
    -------
    dict
        Mapping ``(tex_state, tam_type, lr_pair) → mean_strength``.
    """
    tex_states: List[str] = metadata["tex_states"]
    tam_types: List[str] = metadata["tam_types"]
    lr_pairs: List[str] = metadata["lr_pairs"]

    n_tex = len(tex_states)
    n_tam = len(tam_types)
    n_lr = len(lr_pairs)

    expected_pairs = n_tex * n_tam
    if comm_tensor.shape[1] != expected_pairs:
        raise ValueError(
            f"comm_tensor has {comm_tensor.shape[1]} cell-pair entries but "
            f"n_tex ({n_tex}) × n_tam ({n_tam}) = {expected_pairs}."
        )
    if comm_tensor.shape[2] != n_lr:
        raise ValueError(
            f"comm_tensor has {comm_tensor.shape[2]} L-R entries but "
            f"len(lr_pairs) = {n_lr}."
        )

    # Average across patients (axis 0)
    mean_tensor = comm_tensor.mean(axis=0)  # (n_tex * n_tam, n_lr)

    ref: Dict[Tuple[str, str, str], float] = {}
    for ti, tex in enumerate(tex_states):
        for mi, tam in enumerate(tam_types):
            pair_idx = ti * n_tam + mi
            for li, lr in enumerate(lr_pairs):
                ref[(tex, tam, lr)] = float(mean_tensor[pair_idx, li])

    return ref
