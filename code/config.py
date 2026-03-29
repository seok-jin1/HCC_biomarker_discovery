"""
config.py — Central configuration for HCC scRNA-seq analysis project.

Defines path constants, dataset metadata, cell-type marker gene lists,
QC/model hyperparameters, and a helper to verify the data mount.
"""

from pathlib import Path
import os

# ---------------------------------------------------------------------------
# Random seed
# ---------------------------------------------------------------------------
SEED = 42

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
# Large data files live on an external drive mounted at /mnt/e
DATA_ROOT = Path("/mnt/e/no_exp_paper/data")

RAW_DIR        = DATA_ROOT / "raw"
PROCESSED_DIR  = DATA_ROOT / "processed"
EXTERNAL_DIR   = DATA_ROOT / "external"
ICB_DIR        = EXTERNAL_DIR / "ICB_cohorts"
CHECKPOINT_DIR = DATA_ROOT / "model_checkpoints"

# Code and results live in the WSL home directory
PROJECT_ROOT = Path("/home/laugh/no_exp_paper")

CODE_DIR    = PROJECT_ROOT / "code"
RESULTS_DIR = PROJECT_ROOT / "results"
FIG_DIR     = RESULTS_DIR / "figures"
TABLE_DIR   = RESULTS_DIR / "tables"

# ---------------------------------------------------------------------------
# Primary scRNA-seq datasets
# ---------------------------------------------------------------------------
DATASETS = {
    "GSE149614": {
        "raw_dir": RAW_DIR / "GSE149614",
        "platform": "10x Chromium",
        "sorted": False,
        "urls": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz",
        ],
    },
    "GSE140228": {
        "raw_dir": RAW_DIR / "GSE140228",
        "platform": "10x Chromium (Droplet) + Smart-seq2",
        "sorted": True,   # CD45-sorted
        "urls": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/GSE140228_UMI_counts_Droplet.mtx.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/GSE140228_UMI_counts_Droplet_barcodes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/GSE140228_UMI_counts_Droplet_genes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz",
        ],
    },
    "GSE156625": {
        "raw_dir": RAW_DIR / "GSE156625",
        "platform": "10x Chromium",
        "sorted": False,
        "urls": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/GSE156625_HCCscanpyobj.h5ad.gz",
        ],
    },
    "GSE151530": {
        "raw_dir": RAW_DIR / "GSE151530",
        "platform": "10x Chromium",
        "sorted": False,
        "urls": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_matrix.mtx.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_barcodes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_genes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_Info.txt.gz",
        ],
    },
}

# ---------------------------------------------------------------------------
# Validation / bulk datasets
# ---------------------------------------------------------------------------
VALIDATION = {
    "TCGA_LIHC": {
        "expression": EXTERNAL_DIR / "TCGA_LIHC" / "TCGA-LIHC.htseq_fpkm.tsv.gz",
        "clinical":   EXTERNAL_DIR / "TCGA_LIHC" / "TCGA-LIHC.survival.tsv",
    },
    "GSE14520": {
        "expression":    EXTERNAL_DIR / "GSE14520" / "GSE14520_Extra_Supplement.txt.gz",
        "raw_dir":       EXTERNAL_DIR / "GSE14520",
    },
    "GSE285963": {
        "expression": ICB_DIR / "GSE285963" / "GSE285963_CU24_n84_readcount.gct.gz",
        "raw_dir":    ICB_DIR / "GSE285963",
    },
    "GSE206325": {
        "raw_dir": ICB_DIR / "GSE206325",
    },
    "GSE140901": {
        "raw_tar":    ICB_DIR / "GSE140901" / "GSE140901_RAW.tar",
        "processed":  ICB_DIR / "GSE140901" / "GSE140901_processed_data.txt.gz",
        "raw_dir":    ICB_DIR / "GSE140901",
    },
}

# ---------------------------------------------------------------------------
# Cell-type marker genes
# ---------------------------------------------------------------------------
MAJOR_MARKERS = {
    "T_cell":      ["CD3D", "CD3E", "CD3G", "CD247", "TRAC"],
    "Myeloid":     ["CD68", "CD14", "LYZ", "CSF1R", "ITGAM"],
    "B_cell":      ["CD19", "MS4A1", "CD79A", "CD79B", "IGHM"],
    "NK":          ["NCAM1", "KLRB1", "KLRD1", "NKG7", "GNLY"],
    "Fibroblast":  ["COL1A1", "COL1A2", "COL3A1", "ACTA2", "FAP"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "ENG", "CLDN5"],
    "Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19", "ALB"],
}

# ---------------------------------------------------------------------------
# T-cell exhaustion state markers
# ---------------------------------------------------------------------------
EXHAUSTION_MARKERS = {
    "progenitor":   ["TCF7", "SELL", "LEF1", "CCR7"],
    "transitional": ["PDCD1", "TOX", "CXCL13", "CD38"],
    "terminal":     ["HAVCR2", "LAG3", "ENTPD1", "LAYN", "TIGIT", "CTLA4"],
}

# ---------------------------------------------------------------------------
# Tumour-associated macrophage (TAM) subtype markers
# ---------------------------------------------------------------------------
TAM_MARKERS = {
    "SPP1":  ["SPP1", "APOC1", "TREM2", "CD63"],
    "TREM2": ["TREM2", "GPNMB", "CD9", "LGALS3"],
    "FOLR2": ["FOLR2", "LYVE1", "TIMD4", "MRC1"],
    "FCN1":  ["FCN1", "S100A8", "S100A9", "VCAN"],
    "C1Q":   ["C1QA", "C1QB", "C1QC", "SEPP1"],
    "ISG":   ["IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1"],
}

# ---------------------------------------------------------------------------
# Quality-control parameters
# ---------------------------------------------------------------------------
QC_PARAMS = {
    "min_genes":    200,
    "max_genes":    6000,
    "max_pct_mito": 20,
    "min_cells":    3,
}

# ---------------------------------------------------------------------------
# scVI integration parameters
# ---------------------------------------------------------------------------
SCVI_PARAMS = {
    "n_latent":   30,
    "n_layers":   2,
    "max_epochs": 200,
    "batch_key":  "dataset",
}

# ---------------------------------------------------------------------------
# CACE (cross-attention cell encoder) parameters
# ---------------------------------------------------------------------------
CACE_PARAMS = {
    "d_model":      128,
    "n_heads":      8,
    "n_layers":     2,
    "dropout":      0.3,
    "lr":           1e-4,
    "weight_decay": 1e-5,
    "max_epochs":   100,
    "patience":     15,
}

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def check_data_mount() -> bool:
    """Verify that the external data drive is mounted and DATA_ROOT exists.

    Returns
    -------
    bool
        True if the mount is accessible, False otherwise (with a warning).
    """
    if DATA_ROOT.exists():
        print(f"[OK] DATA_ROOT is accessible: {DATA_ROOT}")
        return True
    else:
        mount_point = Path("/mnt/e")
        if not mount_point.exists():
            print(
                f"[WARNING] /mnt/e is not mounted. "
                "Mount the external drive with:\n"
                "  sudo mount -t drvfs E: /mnt/e"
            )
        else:
            print(
                f"[WARNING] DATA_ROOT does not exist: {DATA_ROOT}\n"
                f"  /mnt/e is present but {DATA_ROOT} has not been created yet.\n"
                f"  Run: mkdir -p {DATA_ROOT}"
            )
        return False


if __name__ == "__main__":
    check_data_mount()
