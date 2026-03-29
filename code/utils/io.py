"""
utils/io.py — Data I/O helpers for HCC scRNA-seq project.

Each loader returns an AnnData object with a ``dataset`` column added to
``adata.obs``.  Call ``save_h5ad`` to persist any AnnData to disk.
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Optional

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def decompress_gz(src: Path, dest: Optional[Path] = None) -> Path:
    """Decompress a .gz file in-place or to an explicit destination.

    Parameters
    ----------
    src:
        Path to the ``.gz`` file.
    dest:
        Destination path for the decompressed file.  Defaults to *src* with
        the ``.gz`` extension stripped.

    Returns
    -------
    Path
        Path to the decompressed file.
    """
    if dest is None:
        dest = src.with_suffix("")  # strip .gz

    if dest.exists():
        print(f"[skip] already decompressed: {dest.name}")
        return dest

    print(f"[decompress] {src.name} -> {dest.name}")
    with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return dest


def load_10x_mtx(
    matrix_path: Path,
    barcodes_path: Path,
    genes_path: Path,
    dataset_name: str,
) -> ad.AnnData:
    """Load a 10x-style (matrix, barcodes, genes) triplet into AnnData.

    Parameters
    ----------
    matrix_path:
        Path to the ``.mtx`` (or ``.mtx.gz``) file.
    barcodes_path:
        Path to the barcodes ``.tsv`` (or ``.tsv.gz``) file.
    genes_path:
        Path to the genes/features ``.tsv`` (or ``.tsv.gz``) file.
    dataset_name:
        Value written to ``adata.obs["dataset"]``.

    Returns
    -------
    ad.AnnData
        Cells x genes AnnData with sparse count matrix.
    """
    # Read barcodes
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].tolist()

    # Read genes — first column is Ensembl ID, second (if present) is symbol
    genes_df = pd.read_csv(genes_path, header=None, sep="\t")
    gene_ids = genes_df.iloc[:, 0].tolist()
    gene_names = (
        genes_df.iloc[:, 1].tolist() if genes_df.shape[1] > 1 else gene_ids
    )

    # Read sparse matrix (genes x cells in MTX format → transpose)
    mat = scipy.io.mmread(str(matrix_path))
    mat = scipy.sparse.csr_matrix(mat).T  # cells x genes

    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame({"gene_ids": gene_ids}, index=gene_names),
    )
    adata.obs["dataset"] = dataset_name
    return adata


# ---------------------------------------------------------------------------
# Dataset-specific loaders
# ---------------------------------------------------------------------------

def load_gse149614(raw_dir: Optional[Path] = None) -> ad.AnnData:
    """Load GSE149614 (dense count matrix + metadata).

    The supplemental file is a genes x cells tab-separated count matrix.

    Parameters
    ----------
    raw_dir:
        Directory containing the downloaded supplemental files.
        Defaults to the value in ``config.DATASETS["GSE149614"]["raw_dir"]``.

    Returns
    -------
    ad.AnnData
    """
    if raw_dir is None:
        from code.config import DATASETS
        raw_dir = DATASETS["GSE149614"]["raw_dir"]

    raw_dir = Path(raw_dir)

    count_gz   = raw_dir / "GSE149614_HCC.scRNAseq.S71271.count.txt.gz"
    meta_gz    = raw_dir / "GSE149614_HCC.scRNAseq.S71271.metadata.txt.gz"

    print(f"[load_gse149614] reading count matrix …")
    counts = pd.read_csv(count_gz, sep="\t", index_col=0, compression="gzip")
    # counts: genes x cells → transpose to cells x genes
    counts = counts.T

    print(f"[load_gse149614] reading metadata …")
    meta = pd.read_csv(meta_gz, sep="\t", index_col=0, compression="gzip")

    # Align metadata to cells present in count matrix
    shared = counts.index.intersection(meta.index)
    counts = counts.loc[shared]
    meta   = meta.loc[shared]

    adata = ad.AnnData(
        X=scipy.sparse.csr_matrix(counts.values),
        obs=meta.copy(),
        var=pd.DataFrame(index=counts.columns),
    )
    adata.obs["dataset"] = "GSE149614"
    print(f"[load_gse149614] {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def load_gse140228(raw_dir: Optional[Path] = None) -> ad.AnnData:
    """Load GSE140228 Droplet (10x CD45-sorted) data only.

    Smart-seq2 data from this accession is intentionally excluded.
    Adds ``cd45_sorted=True`` to ``adata.obs``.

    Parameters
    ----------
    raw_dir:
        Directory containing the downloaded supplemental files.

    Returns
    -------
    ad.AnnData
    """
    if raw_dir is None:
        from code.config import DATASETS
        raw_dir = DATASETS["GSE140228"]["raw_dir"]

    raw_dir = Path(raw_dir)

    matrix_gz   = raw_dir / "GSE140228_Droplet_raw_data_matrix.mtx.gz"
    barcodes_gz = raw_dir / "GSE140228_Droplet_raw_data_barcodes.tsv.gz"
    genes_gz    = raw_dir / "GSE140228_Droplet_raw_data_genes.tsv.gz"
    cellinfo_gz = raw_dir / "GSE140228_Droplet_cell_info.txt.gz"

    # Decompress to temp files if needed
    matrix_f   = decompress_gz(matrix_gz)
    barcodes_f = decompress_gz(barcodes_gz)
    genes_f    = decompress_gz(genes_gz)

    print("[load_gse140228] loading Droplet MTX …")
    adata = load_10x_mtx(matrix_f, barcodes_f, genes_f, dataset_name="GSE140228")

    # Attach cell metadata
    print("[load_gse140228] attaching cell info …")
    cell_info = pd.read_csv(cellinfo_gz, sep="\t", index_col=0, compression="gzip")
    shared = adata.obs_names.intersection(cell_info.index)
    adata  = adata[shared].copy()
    adata.obs = adata.obs.join(cell_info, how="left")

    adata.obs["cd45_sorted"] = True
    adata.obs["dataset"]     = "GSE140228"

    print(f"[load_gse140228] {adata.n_obs} cells × {adata.n_vars} genes (Droplet only)")
    return adata


def load_gse156625(raw_dir: Optional[Path] = None) -> ad.AnnData:
    """Load GSE156625 pre-processed h5ad object.

    Parameters
    ----------
    raw_dir:
        Directory containing ``GSE156625_HCCscanpyobj.h5ad.gz``.

    Returns
    -------
    ad.AnnData
    """
    if raw_dir is None:
        from code.config import DATASETS
        raw_dir = DATASETS["GSE156625"]["raw_dir"]

    raw_dir = Path(raw_dir)
    h5ad_gz = raw_dir / "GSE156625_HCCscanpyobj.h5ad.gz"

    # Decompress .h5ad.gz → .h5ad
    h5ad_path = decompress_gz(h5ad_gz)

    print("[load_gse156625] reading h5ad …")
    adata = ad.read_h5ad(h5ad_path)
    adata.obs["dataset"] = "GSE156625"

    print(f"[load_gse156625] {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def load_gse151530(raw_dir: Optional[Path] = None) -> ad.AnnData:
    """Load GSE151530 and filter to HCC samples only.

    Uses the ``Cancer_type`` column in ``CellInfo.txt`` to retain only cells
    annotated as ``HCC`` (hepatocellular carcinoma).

    Parameters
    ----------
    raw_dir:
        Directory containing the downloaded supplemental files.

    Returns
    -------
    ad.AnnData
        HCC cells only.
    """
    if raw_dir is None:
        from code.config import DATASETS
        raw_dir = DATASETS["GSE151530"]["raw_dir"]

    raw_dir = Path(raw_dir)

    matrix_gz   = raw_dir / "GSE151530_matrix.mtx.gz"
    barcodes_gz = raw_dir / "GSE151530_barcodes.tsv.gz"
    genes_gz    = raw_dir / "GSE151530_genes.tsv.gz"
    cellinfo_gz = raw_dir / "GSE151530_CellInfo.txt.gz"

    matrix_f   = decompress_gz(matrix_gz)
    barcodes_f = decompress_gz(barcodes_gz)
    genes_f    = decompress_gz(genes_gz)

    print("[load_gse151530] loading MTX …")
    adata = load_10x_mtx(matrix_f, barcodes_f, genes_f, dataset_name="GSE151530")

    print("[load_gse151530] reading CellInfo …")
    cell_info = pd.read_csv(cellinfo_gz, sep="\t", index_col=0, compression="gzip")

    # Attach metadata
    shared = adata.obs_names.intersection(cell_info.index)
    adata  = adata[shared].copy()
    adata.obs = adata.obs.join(cell_info, how="left")

    # Filter to HCC only
    if "Cancer_type" in adata.obs.columns:
        n_before = adata.n_obs
        adata = adata[adata.obs["Cancer_type"] == "HCC"].copy()
        print(
            f"[load_gse151530] kept {adata.n_obs}/{n_before} HCC cells "
            f"(filtered by Cancer_type=='HCC')"
        )
    else:
        print("[load_gse151530] WARNING: 'Cancer_type' column not found; returning all cells")

    adata.obs["dataset"] = "GSE151530"
    print(f"[load_gse151530] {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------

def save_h5ad(adata: ad.AnnData, path: Path, compression: str = "gzip") -> None:
    """Save AnnData to an h5ad file, creating parent directories as needed.

    Parameters
    ----------
    adata:
        AnnData object to save.
    path:
        Destination file path (should end in ``.h5ad``).
    compression:
        Compression algorithm passed to ``anndata.write_h5ad``.
        Use ``None`` for no compression (faster writes, larger files).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[save_h5ad] writing {adata.n_obs} cells × {adata.n_vars} genes → {path}")
    adata.write_h5ad(path, compression=compression)

    size_mb = path.stat().st_size / 1_048_576
    print(f"[save_h5ad] done — file size: {size_mb:.1f} MB")
