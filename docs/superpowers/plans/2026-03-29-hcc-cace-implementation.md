# HCC CACE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a cross-attention deep learning framework (CACE) that models T cell exhaustion–TAM communication networks from integrated HCC scRNA-seq data and predicts ICB response.

**Architecture:** Integrate 4 public HCC scRNA-seq datasets (~350K cells) → discover T cell exhaustion states and TAM subtypes → quantify cell-cell communication via CellChat → train Cross-Attention Communication Encoder (CACE) with dual prognostic/ICB heads → validate on TCGA-LIHC (n=371), GSE14520 (n=488), GSE285963 (n=73 ICB), GSE206325 (n=18 ICB scRNA-seq).

**Tech Stack:** Python (Scanpy, scVI, PyTorch, lifelines, SHAP), R (CellChat v2, NicheNet, pySCENIC via rpy2 or standalone)

**Data locations:**
- Raw/processed data: `/mnt/e/no_exp_paper/data/`
- Code/results: `/home/laugh/no_exp_paper/`

---

## File Structure

```
code/
├── config.py                          # 경로 상수, random seed, 공통 파라미터
├── requirements.txt                   # Python dependencies
├── environment.yml                    # Conda environment (Python + R)
├── preprocessing/
│   ├── 00_download.sh                 # GEO 데이터 다운로드 스크립트
│   ├── 01_load_datasets.py            # 4개 데이터셋 로드 → 개별 h5ad
│   ├── 02_qc_filtering.py             # QC, doublet removal, HCC 필터링
│   ├── 03_integration.py              # scVI batch integration
│   └── 04_annotation.py               # Major cell type annotation
├── analysis/
│   ├── 01_tcell_exhaustion.py         # CD8+ T cell subsetting, pseudotime, states
│   ├── 02_tam_subtyping.py            # TAM subsetting, NMF, subtype annotation
│   ├── 03_cellchat.R                  # CellChat v2 L-R interaction
│   ├── 04_nichenet.R                  # NicheNet upstream ligand analysis
│   ├── 05_scenic.py                   # pySCENIC TF activity
│   └── 06_communication_tensor.py     # 3D tensor 구축, differential communication
├── models/
│   ├── cace_model.py                  # CACE PyTorch architecture
│   ├── deconvolution.py               # Bulk → cell state fraction → pseudo-comm matrix
│   ├── train.py                       # Training loop (multi-task Cox + BCE)
│   ├── evaluate.py                    # Evaluation metrics, comparison vs baselines
│   └── interpret.py                   # Attention weights, SHAP, integrated gradients
├── validation/
│   ├── 01_prognostic.py               # TCGA-LIHC + GSE14520 KM, Cox, tdAUC
│   ├── 02_icb_response.py             # GSE285963 + GSE140901 ROC, comparison
│   ├── 03_icb_scrna.py                # GSE206325 scRNA-seq mechanistic validation
│   └── 04_ablation.py                 # Ablation study
├── utils/
│   ├── io.py                          # 데이터 I/O helpers
│   ├── plotting.py                    # 공통 시각화 함수
│   └── metrics.py                     # Survival, AUC, communication metrics
└── notebooks/
    ├── EDA_01_per_dataset.ipynb        # 데이터셋별 탐색
    ├── EDA_02_integrated.ipynb         # 통합 후 탐색
    └── Figure_01_to_08.ipynb           # Main figure 생성
```

---

## Task 1: Project Setup — config, environment, download

**Files:**
- Create: `code/config.py`
- Create: `code/requirements.txt`
- Create: `code/preprocessing/00_download.sh`
- Create: `code/utils/io.py`
- Create: `.gitignore`

- [ ] **Step 1: Create config.py**

```python
# code/config.py
from pathlib import Path
import os

SEED = 42

# === Data paths (external drive) ===
DATA_ROOT = Path("/mnt/e/no_exp_paper/data")
RAW_DIR = DATA_ROOT / "raw"
PROCESSED_DIR = DATA_ROOT / "processed"
EXTERNAL_DIR = DATA_ROOT / "external"
ICB_DIR = EXTERNAL_DIR / "ICB_cohorts"
CHECKPOINT_DIR = DATA_ROOT / "model_checkpoints"

# === Project paths (local) ===
PROJECT_ROOT = Path("/home/laugh/no_exp_paper")
CODE_DIR = PROJECT_ROOT / "code"
RESULTS_DIR = PROJECT_ROOT / "results"
FIG_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"

# === Dataset-specific ===
DATASETS = {
    "GSE149614": {
        "raw_dir": RAW_DIR / "GSE149614",
        "url_count": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
        "url_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz",
        "platform": "10x_v2",
        "sorted": False,
        "species": "human",
    },
    "GSE140228": {
        "raw_dir": RAW_DIR / "GSE140228",
        "url_prefix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/",
        "files": [
            "GSE140228_UMI_counts_Droplet.mtx.gz",
            "GSE140228_UMI_counts_Droplet_barcodes.tsv.gz",
            "GSE140228_UMI_counts_Droplet_genes.tsv.gz",
            "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz",
        ],
        "platform": "10x_v2",
        "sorted": True,  # CD45+ sorted
        "species": "human",
    },
    "GSE156625": {
        "raw_dir": RAW_DIR / "GSE156625",
        "url_h5ad": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/GSE156625_HCCscanpyobj.h5ad.gz",
        "platform": "10x_v2v3",
        "sorted": False,
        "species": "human",  # HCC-only h5ad
    },
    "GSE151530": {
        "raw_dir": RAW_DIR / "GSE151530",
        "url_prefix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/",
        "files": [
            "GSE151530_matrix.mtx.gz",
            "GSE151530_barcodes.tsv.gz",
            "GSE151530_genes.tsv.gz",
            "GSE151530_Info.txt.gz",
        ],
        "platform": "10x_v2",
        "sorted": False,
        "species": "human",  # HCC+iCCA mixed → filter needed
    },
}

# === Validation datasets ===
VALIDATION = {
    "TCGA_LIHC": EXTERNAL_DIR / "TCGA_LIHC",
    "GSE14520": EXTERNAL_DIR / "GSE14520",
    "GSE285963": ICB_DIR / "GSE285963",
    "GSE206325": ICB_DIR / "GSE206325",
    "GSE140901": ICB_DIR / "GSE140901",
}

# === Cell type markers ===
MAJOR_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "CD2"],
    "Myeloid": ["CD68", "CD14", "LYZ", "CSF1R"],
    "B_cell": ["CD79A", "MS4A1", "CD19"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "Fibroblast": ["COL1A1", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Epithelial": ["EPCAM", "KRT19", "ALB"],
}

EXHAUSTION_MARKERS = {
    "progenitor": ["TCF7", "SELL", "LEF1", "CCR7"],
    "transitional": ["PDCD1", "TOX", "CXCL13", "CD38"],
    "terminal": ["HAVCR2", "LAG3", "ENTPD1", "LAYN", "TIGIT", "CTLA4"],
}

TAM_MARKERS = {
    "SPP1": ["SPP1", "MARCO", "VEGFA"],
    "TREM2": ["TREM2", "CD9", "GPNMB"],
    "FOLR2": ["FOLR2", "SEPP1", "SLC40A1"],
    "FCN1": ["FCN1", "S100A8", "S100A9", "VCAN"],
    "C1Q": ["C1QA", "C1QB", "C1QC"],
    "ISG": ["ISG15", "IFIT1", "MX1", "IFI6"],
}

# === Analysis parameters ===
QC_PARAMS = {
    "min_genes": 200,
    "max_genes": 6000,
    "max_pct_mito": 20,
    "min_cells": 3,
}

SCVI_PARAMS = {
    "n_latent": 30,
    "n_layers": 2,
    "max_epochs": 200,
    "batch_key": "dataset",
}

CACE_PARAMS = {
    "d_model": 128,
    "n_heads": 8,
    "n_layers": 2,
    "dropout": 0.3,
    "lr": 1e-4,
    "weight_decay": 1e-5,
    "max_epochs": 100,
    "patience": 15,
}

def check_data_mount():
    """외장하드 마운트 확인"""
    if not DATA_ROOT.exists():
        raise RuntimeError(
            f"External drive not mounted at {DATA_ROOT}. "
            "Run: ls /mnt/e/ to check mount status."
        )
    print(f"✓ Data root accessible: {DATA_ROOT}")
```

- [ ] **Step 2: Create requirements.txt**

```
# code/requirements.txt
# === Core scRNA-seq ===
scanpy>=1.10
anndata>=0.10
scvi-tools>=1.1
scrublet>=0.2.3
leidenalg>=0.10
pyscenic>=0.12

# === Trajectory ===
scvelo>=0.3
cellrank>=2.0

# === Communication (Python side) ===
liana>=1.2

# === ML / DL ===
torch>=2.2
scikit-learn>=1.4
xgboost>=2.0
lifelines>=0.28
shap>=0.44
captum>=0.7

# === Visualization ===
matplotlib>=3.8
seaborn>=0.13

# === Stats / Utils ===
scipy>=1.12
statsmodels>=0.14
pandas>=2.1
numpy>=1.26
tqdm
```

- [ ] **Step 3: Create download script**

```bash
#!/bin/bash
# code/preprocessing/00_download.sh
# Usage: bash code/preprocessing/00_download.sh
set -euo pipefail

DATA_ROOT="/mnt/e/no_exp_paper/data"

echo "=== Checking mount ==="
if [ ! -d "$DATA_ROOT" ]; then
    echo "ERROR: External drive not mounted at $DATA_ROOT"
    exit 1
fi

echo "=== GSE149614 (Lu et al. 2022) ==="
cd "$DATA_ROOT/raw/GSE149614"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz"

echo "=== GSE140228 (Zhang et al. 2019) ==="
cd "$DATA_ROOT/raw/GSE140228"
for f in UMI_counts_Droplet.mtx.gz UMI_counts_Droplet_barcodes.tsv.gz \
         UMI_counts_Droplet_genes.tsv.gz UMI_counts_Droplet_cellinfo.tsv.gz; do
    wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/GSE140228_${f}"
done

echo "=== GSE156625 (Sharma et al. 2020) — HCC-only h5ad ==="
cd "$DATA_ROOT/raw/GSE156625"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/GSE156625_HCCscanpyobj.h5ad.gz"

echo "=== GSE151530 (Ma et al. 2021) ==="
cd "$DATA_ROOT/raw/GSE151530"
for f in matrix.mtx.gz barcodes.tsv.gz genes.tsv.gz Info.txt.gz; do
    wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_${f}"
done

echo "=== Validation: GSE285963 (ICB RNA-seq) ==="
cd "$DATA_ROOT/external/ICB_cohorts/GSE285963"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285963/suppl/GSE285963_CU24_n84_readcount.gct.gz"

echo "=== Validation: GSE140901 (ICB NanoString) ==="
cd "$DATA_ROOT/external/ICB_cohorts/GSE140901"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140901/suppl/GSE140901_RAW.tar"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140901/suppl/GSE140901_processed_data.txt.gz"

echo "=== Validation: GSE14520 (Bulk microarray) ==="
cd "$DATA_ROOT/external/GSE14520"
wget -nc "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz"

echo "=== DONE ==="
du -sh "$DATA_ROOT/raw/"*
du -sh "$DATA_ROOT/external/"*
```

- [ ] **Step 4: Create io.py utilities**

```python
# code/utils/io.py
"""Data I/O helpers for HCC CACE project."""
import gzip
import shutil
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix


def decompress_gz(gz_path: Path, out_path: Path | None = None) -> Path:
    """Decompress .gz file if not already done."""
    if out_path is None:
        out_path = gz_path.with_suffix("")
    if out_path.exists():
        print(f"  Already decompressed: {out_path.name}")
        return out_path
    print(f"  Decompressing: {gz_path.name}")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out_path


def load_10x_mtx(data_dir: Path, prefix: str = "") -> ad.AnnData:
    """Load 10x-style matrix.mtx + barcodes.tsv + genes.tsv."""
    mtx_file = list(data_dir.glob(f"{prefix}*matrix.mtx*"))[0]
    bc_file = list(data_dir.glob(f"{prefix}*barcodes.tsv*"))[0]
    gene_file = list(data_dir.glob(f"{prefix}*genes.tsv*"))[0]

    # Decompress if needed
    if mtx_file.suffix == ".gz":
        mtx_file = decompress_gz(mtx_file)
    if bc_file.suffix == ".gz":
        bc_file = decompress_gz(bc_file)
    if gene_file.suffix == ".gz":
        gene_file = decompress_gz(gene_file)

    mat = csr_matrix(mmread(mtx_file).T)
    barcodes = pd.read_csv(bc_file, header=None, sep="\t")[0].values
    genes = pd.read_csv(gene_file, header=None, sep="\t")
    gene_names = genes.iloc[:, -1].values  # last column = gene symbol

    adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_names))
    adata.var_names_make_unique()
    return adata


def load_gse149614(raw_dir: Path) -> ad.AnnData:
    """Load GSE149614 count matrix + metadata."""
    count_gz = raw_dir / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    meta_gz = raw_dir / "GSE149614_HCC.metadata.updated.txt.gz"

    print("Loading GSE149614 count matrix...")
    counts = pd.read_csv(count_gz, sep="\t", index_col=0)
    meta = pd.read_csv(meta_gz, sep="\t", index_col=0)

    adata = ad.AnnData(X=csr_matrix(counts.T.values), obs=meta.loc[counts.columns], var=pd.DataFrame(index=counts.index))
    adata.var_names_make_unique()
    adata.obs["dataset"] = "GSE149614"
    return adata


def load_gse140228(raw_dir: Path) -> ad.AnnData:
    """Load GSE140228 droplet data (CD45+ sorted immune cells)."""
    print("Loading GSE140228 droplet...")
    adata = load_10x_mtx(raw_dir, prefix="GSE140228_UMI_counts_Droplet")

    cellinfo_gz = raw_dir / "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz"
    cellinfo = pd.read_csv(cellinfo_gz, sep="\t", index_col=0)
    adata.obs = adata.obs.join(cellinfo, how="left")
    adata.obs["dataset"] = "GSE140228"
    adata.obs["cd45_sorted"] = True
    return adata


def load_gse156625(raw_dir: Path) -> ad.AnnData:
    """Load GSE156625 HCC-only h5ad."""
    h5ad_gz = raw_dir / "GSE156625_HCCscanpyobj.h5ad.gz"
    h5ad_file = h5ad_gz.with_suffix("")  # remove .gz

    if not h5ad_file.exists():
        decompress_gz(h5ad_gz, h5ad_file)

    print("Loading GSE156625 HCC h5ad...")
    adata = sc.read_h5ad(h5ad_file)
    adata.obs["dataset"] = "GSE156625"
    return adata


def load_gse151530(raw_dir: Path) -> ad.AnnData:
    """Load GSE151530 and filter to HCC only."""
    print("Loading GSE151530...")
    adata = load_10x_mtx(raw_dir, prefix="GSE151530")

    info_gz = raw_dir / "GSE151530_Info.txt.gz"
    info = pd.read_csv(info_gz, sep="\t", index_col=0)
    adata.obs = adata.obs.join(info, how="left")
    adata.obs["dataset"] = "GSE151530"

    # Filter to HCC only (exclude iCCA)
    if "Cancer_type" in adata.obs.columns:
        n_before = adata.n_obs
        adata = adata[adata.obs["Cancer_type"].str.contains("HCC", case=False, na=False)].copy()
        print(f"  Filtered HCC: {n_before} → {adata.n_obs} cells")
    return adata


def save_h5ad(adata: ad.AnnData, path: Path, compression: str = "gzip"):
    """Save AnnData with directory creation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression=compression)
    size_mb = path.stat().st_size / 1e6
    print(f"Saved: {path} ({size_mb:.1f} MB, {adata.n_obs} cells × {adata.n_vars} genes)")
```

- [ ] **Step 5: Create .gitignore**

```
# .gitignore
*.h5ad
*.h5
*.rds
*.RData
*.gct
*.gct.gz
*.mtx
*.mtx.gz
*.tar
*.CEL
*.cel
__pycache__/
.ipynb_checkpoints/
*.pyc
.DS_Store
data/
/mnt/
*.log
wandb/
lightning_logs/
```

- [ ] **Step 6: Verify and commit**

Run:
```bash
cd /home/laugh/no_exp_paper
python -c "from code.config import check_data_mount; check_data_mount()"
```
Expected: `✓ Data root accessible: /mnt/e/no_exp_paper/data`

```bash
git init
git add code/config.py code/requirements.txt code/preprocessing/00_download.sh code/utils/io.py .gitignore
git commit -m "feat: project setup — config, download script, io utilities"
```

---

## Task 2: Data Download and Per-Dataset Loading

**Files:**
- Create: `code/preprocessing/01_load_datasets.py`
- Use: `code/preprocessing/00_download.sh`

- [ ] **Step 1: Run download script**

```bash
bash /home/laugh/no_exp_paper/code/preprocessing/00_download.sh
```
Expected: All files downloaded to `/mnt/e/no_exp_paper/data/raw/` and `/mnt/e/.../external/`. Check with `du -sh`.

- [ ] **Step 2: Write 01_load_datasets.py**

```python
# code/preprocessing/01_load_datasets.py
"""Load all 4 scRNA-seq datasets into individual h5ad files."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

from config import DATASETS, PROCESSED_DIR, check_data_mount
from utils.io import (
    load_gse149614, load_gse140228, load_gse156625, load_gse151530, save_h5ad
)

def main():
    check_data_mount()
    out_dir = PROCESSED_DIR / "per_dataset"
    out_dir.mkdir(parents=True, exist_ok=True)

    loaders = {
        "GSE149614": load_gse149614,
        "GSE140228": load_gse140228,
        "GSE156625": load_gse156625,
        "GSE151530": load_gse151530,
    }

    for name, loader in loaders.items():
        print(f"\n{'='*60}")
        print(f"Loading {name}")
        print(f"{'='*60}")
        raw_dir = DATASETS[name]["raw_dir"]
        adata = loader(raw_dir)
        print(f"  Shape: {adata.shape}")
        print(f"  obs columns: {list(adata.obs.columns)}")
        save_h5ad(adata, out_dir / f"{name}.h5ad")

    print("\n✓ All datasets loaded successfully.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run and validate**

```bash
cd /home/laugh/no_exp_paper
python code/preprocessing/01_load_datasets.py
```
Expected: 4 h5ad files in `/mnt/e/no_exp_paper/data/processed/per_dataset/`.
Validate: Each file has >5000 cells, `dataset` column in obs.

- [ ] **Step 4: Commit**

```bash
git add code/preprocessing/01_load_datasets.py
git commit -m "feat: load all 4 HCC scRNA-seq datasets from GEO"
```

---

## Task 3: QC and Filtering

**Files:**
- Create: `code/preprocessing/02_qc_filtering.py`

- [ ] **Step 1: Write QC script**

```python
# code/preprocessing/02_qc_filtering.py
"""QC filtering: gene/cell counts, mito%, doublets."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import scrublet as scr
import numpy as np
import pandas as pd
from config import PROCESSED_DIR, QC_PARAMS, SEED, check_data_mount
from utils.io import save_h5ad

def run_qc(adata, dataset_name: str):
    """Run QC on a single dataset."""
    print(f"\n--- QC: {dataset_name} (n={adata.n_obs}) ---")
    n_start = adata.n_obs

    # Basic metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=QC_PARAMS["min_genes"])
    adata = adata[adata.obs["n_genes_by_counts"] < QC_PARAMS["max_genes"]].copy()
    adata = adata[adata.obs["pct_counts_mt"] < QC_PARAMS["max_pct_mito"]].copy()

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=QC_PARAMS["min_cells"])

    # Doublet detection (Scrublet)
    scrub = scr.Scrublet(adata.X, random_state=SEED)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets
    adata = adata[~adata.obs["predicted_doublet"]].copy()

    n_end = adata.n_obs
    print(f"  {dataset_name}: {n_start} → {n_end} cells ({n_start - n_end} removed, {(n_start-n_end)/n_start*100:.1f}%)")
    return adata

def main():
    check_data_mount()
    in_dir = PROCESSED_DIR / "per_dataset"
    out_dir = PROCESSED_DIR / "qc_filtered"
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = ["GSE149614", "GSE140228", "GSE156625", "GSE151530"]
    summary = []

    for name in datasets:
        adata = sc.read_h5ad(in_dir / f"{name}.h5ad")
        n_before = adata.n_obs
        adata = run_qc(adata, name)
        summary.append({"dataset": name, "before": n_before, "after": adata.n_obs})
        save_h5ad(adata, out_dir / f"{name}_qc.h5ad")

    # Summary table
    df = pd.DataFrame(summary)
    df["removed_pct"] = ((df["before"] - df["after"]) / df["before"] * 100).round(1)
    print(f"\n{'='*60}")
    print("QC Summary:")
    print(df.to_string(index=False))
    print(f"Total cells after QC: {df['after'].sum()}")
    df.to_csv(out_dir / "qc_summary.csv", index=False)

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run and validate**

```bash
python code/preprocessing/02_qc_filtering.py
```
Expected: 4 QC-filtered h5ad files + `qc_summary.csv`. Total cells should be >250K. Each dataset retains >70% of original cells.

- [ ] **Step 3: Commit**

```bash
git add code/preprocessing/02_qc_filtering.py
git commit -m "feat: QC filtering — gene counts, mito%, doublet removal"
```

---

## Task 4: Batch Integration (scVI)

**Files:**
- Create: `code/preprocessing/03_integration.py`

- [ ] **Step 1: Write integration script**

```python
# code/preprocessing/03_integration.py
"""Batch integration using scVI across 4 datasets."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import scvi
import anndata as ad
import numpy as np
from config import PROCESSED_DIR, SCVI_PARAMS, SEED, check_data_mount
from utils.io import save_h5ad

def main():
    check_data_mount()
    scvi.settings.seed = SEED
    in_dir = PROCESSED_DIR / "qc_filtered"

    # 1. Load and concatenate
    print("Loading QC-filtered datasets...")
    datasets = ["GSE149614", "GSE140228", "GSE156625", "GSE151530"]
    adatas = []
    for name in datasets:
        a = sc.read_h5ad(in_dir / f"{name}_qc.h5ad")
        adatas.append(a)
        print(f"  {name}: {a.n_obs} cells")

    adata = ad.concat(adatas, join="inner", label="dataset", keys=datasets)
    print(f"\nConcatenated: {adata.n_obs} cells × {adata.n_vars} genes")

    # 2. Preprocessing for scVI
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # HVG selection
    sc.pp.highly_variable_genes(
        adata, n_top_genes=3000, flavor="seurat_v3",
        layer="counts", batch_key="dataset"
    )
    print(f"HVGs selected: {adata.var['highly_variable'].sum()}")

    # 3. scVI integration
    print("\nSetting up scVI model...")
    scvi.model.SCVI.setup_anndata(
        adata, layer="counts", batch_key=SCVI_PARAMS["batch_key"]
    )
    model = scvi.model.SCVI(
        adata,
        n_latent=SCVI_PARAMS["n_latent"],
        n_layers=SCVI_PARAMS["n_layers"],
    )
    print("Training scVI...")
    model.train(max_epochs=SCVI_PARAMS["max_epochs"], early_stopping=True)

    # Save latent representation
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # 4. Neighbors + UMAP on scVI latent
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=1.0, random_state=SEED)

    # 5. Save
    out_path = PROCESSED_DIR / "integrated.h5ad"
    save_h5ad(adata, out_path)

    # Save scVI model
    model_dir = PROCESSED_DIR / "scvi_model"
    model.save(str(model_dir), overwrite=True)
    print(f"scVI model saved: {model_dir}")

    # 6. Quick stats
    print(f"\n{'='*60}")
    print(f"Integration complete:")
    print(f"  Total cells: {adata.n_obs}")
    print(f"  Total genes: {adata.n_vars}")
    print(f"  Datasets: {adata.obs['dataset'].value_counts().to_dict()}")
    print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run** (may take 30-60 min depending on GPU)

```bash
python code/preprocessing/03_integration.py
```
Expected: `integrated.h5ad` in processed dir. UMAP should show mixing of datasets (check in notebook).

- [ ] **Step 3: Visual QC in notebook**

Open `code/notebooks/EDA_02_integrated.ipynb` and verify:
- Datasets mix well on UMAP (no dataset-specific islands)
- Leiden clusters contain cells from multiple datasets
- Known marker genes show expected patterns

- [ ] **Step 4: Commit**

```bash
git add code/preprocessing/03_integration.py
git commit -m "feat: scVI batch integration across 4 HCC datasets"
```

---

## Task 5: Major Cell Type Annotation

**Files:**
- Create: `code/preprocessing/04_annotation.py`

- [ ] **Step 1: Write annotation script**

```python
# code/preprocessing/04_annotation.py
"""Annotate major cell types in integrated dataset."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import numpy as np
import pandas as pd
from config import PROCESSED_DIR, MAJOR_MARKERS, SEED, check_data_mount
from utils.io import save_h5ad

def score_cell_types(adata):
    """Score each cell for major cell type markers."""
    for ct, markers in MAJOR_MARKERS.items():
        present = [m for m in markers if m in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ct}", random_state=SEED)
        else:
            print(f"  Warning: {ct} — only {len(present)} markers found: {present}")

def assign_cell_types(adata):
    """Assign cell types based on highest marker score per cell."""
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    ct_names = [c.replace("score_", "") for c in score_cols]
    scores = adata.obs[score_cols].values
    best_idx = np.argmax(scores, axis=1)
    adata.obs["cell_type_auto"] = [ct_names[i] for i in best_idx]
    return adata

def main():
    check_data_mount()
    adata = sc.read_h5ad(PROCESSED_DIR / "integrated.h5ad")
    print(f"Loaded: {adata.n_obs} cells, {adata.obs['leiden'].nunique()} clusters")

    # 1. Score cell types
    print("\nScoring cell types...")
    score_cell_types(adata)

    # 2. Auto-assign
    adata = assign_cell_types(adata)
    print("\nAuto cell type distribution:")
    print(adata.obs["cell_type_auto"].value_counts())

    # 3. Cluster-level annotation (majority vote per leiden cluster)
    cluster_ct = adata.obs.groupby("leiden")["cell_type_auto"].agg(
        lambda x: x.value_counts().index[0]
    )
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_ct)
    print("\nCluster-level cell type distribution:")
    print(adata.obs["cell_type"].value_counts())

    # 4. Save subsets for downstream
    save_h5ad(adata, PROCESSED_DIR / "integrated_annotated.h5ad")

    # T cell subset
    tcell = adata[adata.obs["cell_type"] == "T_cell"].copy()
    save_h5ad(tcell, PROCESSED_DIR / "tcell_subset.h5ad")

    # Myeloid subset
    myeloid = adata[adata.obs["cell_type"] == "Myeloid"].copy()
    save_h5ad(myeloid, PROCESSED_DIR / "myeloid_subset.h5ad")

    print(f"\n✓ Annotation complete.")
    print(f"  T cells: {tcell.n_obs}")
    print(f"  Myeloid: {myeloid.n_obs}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run and validate**

```bash
python code/preprocessing/04_annotation.py
```
Expected: T cell >40K, Myeloid >20K. Saved `integrated_annotated.h5ad`, `tcell_subset.h5ad`, `myeloid_subset.h5ad`.

- [ ] **Step 3: Commit**

```bash
git add code/preprocessing/04_annotation.py
git commit -m "feat: major cell type annotation — T cell and myeloid subsets"
```

---

## Task 6: T Cell Exhaustion State Analysis

**Files:**
- Create: `code/analysis/01_tcell_exhaustion.py`
- Create: `code/utils/metrics.py`

- [ ] **Step 1: Create metrics utilities**

```python
# code/utils/metrics.py
"""Analysis metrics and scoring utilities."""
import numpy as np
import pandas as pd
import scanpy as sc


def score_gene_modules(adata, modules: dict, seed: int = 42):
    """Score cells for multiple gene modules."""
    for name, genes in modules.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"module_{name}", random_state=seed)
        else:
            print(f"  Warning: module '{name}' — only {len(present)}/{len(genes)} genes found")
    return adata


def assign_exhaustion_state(adata, progenitor_col="module_progenitor",
                            transitional_col="module_transitional",
                            terminal_col="module_terminal"):
    """Assign exhaustion state based on module scores."""
    scores = adata.obs[[progenitor_col, transitional_col, terminal_col]].values
    states = ["Progenitor_Tex", "Transitional_Tex", "Terminal_Tex"]
    best_idx = np.argmax(scores, axis=1)
    adata.obs["exhaustion_state"] = [states[i] for i in best_idx]
    return adata
```

- [ ] **Step 2: Write T cell exhaustion analysis**

```python
# code/analysis/01_tcell_exhaustion.py
"""CD8+ T cell exhaustion trajectory and state assignment."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import scvelo as scv
import numpy as np
from config import PROCESSED_DIR, EXHAUSTION_MARKERS, SEED, check_data_mount
from utils.io import save_h5ad
from utils.metrics import score_gene_modules, assign_exhaustion_state

def main():
    check_data_mount()
    adata = sc.read_h5ad(PROCESSED_DIR / "tcell_subset.h5ad")
    print(f"T cells loaded: {adata.n_obs}")

    # 1. Re-cluster T cells
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.tl.pca(adata, n_comps=30, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.8, random_state=SEED, key_added="leiden_tcell")

    # 2. Separate CD4 vs CD8
    cd8_genes = ["CD8A", "CD8B"]
    cd4_genes = ["CD4"]
    for g in cd8_genes + cd4_genes:
        if g in adata.var_names:
            sc.tl.score_genes(adata, [g], score_name=f"score_{g}", random_state=SEED)
    adata.obs["is_cd8"] = adata.obs.get("score_CD8A", 0) > adata.obs.get("score_CD4", 0)
    cd8 = adata[adata.obs["is_cd8"]].copy()
    print(f"CD8+ T cells: {cd8.n_obs}")

    # 3. Re-cluster CD8
    sc.pp.neighbors(cd8, n_neighbors=20, n_pcs=30, random_state=SEED)
    sc.tl.umap(cd8, random_state=SEED)
    sc.tl.leiden(cd8, resolution=0.6, random_state=SEED, key_added="leiden_cd8")

    # 4. Score exhaustion modules
    cd8 = score_gene_modules(cd8, EXHAUSTION_MARKERS, seed=SEED)

    # 5. Assign exhaustion states
    cd8 = assign_exhaustion_state(cd8)
    print("\nExhaustion state distribution:")
    print(cd8.obs["exhaustion_state"].value_counts())

    # 6. Pseudotime (diffusion pseudotime as fallback for scVelo)
    sc.tl.diffmap(cd8, random_state=SEED)
    # Set root as progenitor-enriched cell
    progenitor_scores = cd8.obs["module_progenitor"].values
    root_cell = np.argmax(progenitor_scores)
    cd8.uns["iroot"] = root_cell
    sc.tl.dpt(cd8)
    print(f"Pseudotime range: {cd8.obs['dpt_pseudotime'].min():.3f} – {cd8.obs['dpt_pseudotime'].max():.3f}")

    # 7. DEG per exhaustion state
    sc.tl.rank_genes_groups(cd8, groupby="exhaustion_state", method="wilcoxon")

    # 8. Save
    save_h5ad(cd8, PROCESSED_DIR / "cd8_exhaustion.h5ad")
    save_h5ad(adata, PROCESSED_DIR / "tcell_annotated.h5ad")  # all T cells with CD4/CD8 label

    print("\n✓ T cell exhaustion analysis complete.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run and validate**

```bash
python code/analysis/01_tcell_exhaustion.py
```
Expected: 3 exhaustion states identified in CD8+ T cells. `cd8_exhaustion.h5ad` saved. Each state should have >1000 cells.

- [ ] **Step 4: Commit**

```bash
git add code/analysis/01_tcell_exhaustion.py code/utils/metrics.py
git commit -m "feat: CD8+ T cell exhaustion trajectory and state assignment"
```

---

## Task 7: TAM Subtype Analysis

**Files:**
- Create: `code/analysis/02_tam_subtyping.py`

- [ ] **Step 1: Write TAM subtyping script**

```python
# code/analysis/02_tam_subtyping.py
"""TAM subtype discovery and annotation."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import numpy as np
from config import PROCESSED_DIR, TAM_MARKERS, SEED, check_data_mount
from utils.io import save_h5ad
from utils.metrics import score_gene_modules

def main():
    check_data_mount()
    adata = sc.read_h5ad(PROCESSED_DIR / "myeloid_subset.h5ad")
    print(f"Myeloid cells loaded: {adata.n_obs}")

    # 1. Re-cluster myeloid
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.tl.pca(adata, n_comps=30, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.8, random_state=SEED, key_added="leiden_myeloid")

    # 2. Filter to macrophages (remove DCs, monocytes if clearly separate)
    # Score for macrophage identity
    macro_markers = ["CD68", "CD163", "MRC1", "MSR1"]
    present = [m for m in macro_markers if m in adata.var_names]
    sc.tl.score_genes(adata, present, score_name="macro_score", random_state=SEED)

    # 3. Score TAM subtypes
    adata = score_gene_modules(adata, TAM_MARKERS, seed=SEED)

    # 4. Assign TAM subtype (highest module score)
    module_cols = [c for c in adata.obs.columns if c.startswith("module_")]
    subtype_names = [c.replace("module_", "") for c in module_cols]
    scores = adata.obs[module_cols].values
    best_idx = np.argmax(scores, axis=1)
    adata.obs["tam_subtype"] = [f"{subtype_names[i]}+_TAM" for i in best_idx]

    print("\nTAM subtype distribution:")
    print(adata.obs["tam_subtype"].value_counts())

    # 5. DEG per TAM subtype
    sc.tl.rank_genes_groups(adata, groupby="tam_subtype", method="wilcoxon")

    # 6. Save
    save_h5ad(adata, PROCESSED_DIR / "tam_subtyped.h5ad")
    print("\n✓ TAM subtyping complete.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run and validate**

```bash
python code/analysis/02_tam_subtyping.py
```
Expected: 5-6 TAM subtypes. Each subtype >500 cells.

- [ ] **Step 3: Commit**

```bash
git add code/analysis/02_tam_subtyping.py
git commit -m "feat: TAM subtype annotation — SPP1, TREM2, FOLR2, FCN1, C1Q, ISG"
```

---

## Task 8: Cell-Cell Communication (CellChat)

**Files:**
- Create: `code/analysis/03_cellchat.R`
- Create: `code/analysis/06_communication_tensor.py`

- [ ] **Step 1: Write CellChat R script**

```r
# code/analysis/03_cellchat.R
# CellChat v2: T cell exhaustion state ↔ TAM subtype L-R interactions
# Run: Rscript code/analysis/03_cellchat.R

library(CellChat)
library(anndata)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
processed_dir <- ifelse(length(args) > 0, args[1], "/mnt/e/no_exp_paper/data/processed")

cat("=== Loading data ===\n")
# Load CD8 exhaustion and TAM data, merge for CellChat
cd8 <- read_h5ad(file.path(processed_dir, "cd8_exhaustion.h5ad"))
tam <- read_h5ad(file.path(processed_dir, "tam_subtyped.h5ad"))

# Combine into one matrix
# (In practice, extract counts and metadata, build CellChat object)
# This is a template — exact implementation depends on anndata R compatibility

# Create CellChat object per patient for patient-level communication
patients <- unique(c(cd8$obs$patient_id, tam$obs$patient_id))

results_list <- list()

for (pt in patients) {
  cat(sprintf("\nProcessing patient: %s\n", pt))

  # Subset to this patient
  cd8_pt <- cd8[cd8$obs$patient_id == pt, ]
  tam_pt <- tam[tam$obs$patient_id == pt, ]

  # Skip if too few cells
  if (nrow(cd8_pt) < 20 || nrow(tam_pt) < 20) {
    cat(sprintf("  Skipping %s: too few cells (CD8=%d, TAM=%d)\n",
                pt, nrow(cd8_pt), nrow(tam_pt)))
    next
  }

  # Combine and create CellChat
  # ... (merge count matrices, set cell labels)
  # cellchat <- createCellChat(combined_data, group.by = "cell_state")
  # cellchat@DB <- CellChatDB.human
  # cellchat <- subsetData(cellchat)
  # cellchat <- identifyOverExpressedGenes(cellchat)
  # cellchat <- identifyOverExpressedInteractions(cellchat)
  # cellchat <- computeCommunProb(cellchat, type = "triMean")
  # cellchat <- computeCommunProbPathway(cellchat)

  # Extract communication probability matrix
  # results_list[[pt]] <- cellchat@net
}

# Save results
save(results_list, file = file.path(processed_dir, "cellchat_results.RData"))
cat("\n✓ CellChat analysis complete.\n")
```

Note: CellChat R script needs adaptation based on actual data format. The exact implementation depends on how anndata loads in R. Alternative: use LIANA+ in Python (see step 2).

- [ ] **Step 2: Write Python alternative using LIANA+**

```python
# code/analysis/03_communication_liana.py
"""Cell-cell communication using LIANA+ (Python alternative to CellChat)."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import liana as li
import pandas as pd
import numpy as np
from config import PROCESSED_DIR, SEED, check_data_mount
from utils.io import save_h5ad

def main():
    check_data_mount()

    # 1. Load and merge CD8 + TAM
    cd8 = sc.read_h5ad(PROCESSED_DIR / "cd8_exhaustion.h5ad")
    tam = sc.read_h5ad(PROCESSED_DIR / "tam_subtyped.h5ad")

    cd8.obs["cell_state"] = cd8.obs["exhaustion_state"]
    tam.obs["cell_state"] = tam.obs["tam_subtype"]

    import anndata as ad
    combined = ad.concat([cd8, tam], join="inner")
    combined.raw = combined
    print(f"Combined: {combined.n_obs} cells, states: {combined.obs['cell_state'].nunique()}")

    # 2. Run LIANA (multi-method consensus)
    li.mt.rank_aggregate(
        combined,
        groupby="cell_state",
        resource_name="consensus",
        use_raw=False,
        verbose=True,
    )
    liana_res = combined.uns["liana_res"]
    print(f"LIANA results: {len(liana_res)} interactions")

    # 3. Filter to T cell ↔ TAM interactions
    tex_states = cd8.obs["exhaustion_state"].unique().tolist()
    tam_types = tam.obs["tam_subtype"].unique().tolist()

    mask = (
        (liana_res["source"].isin(tex_states) & liana_res["target"].isin(tam_types)) |
        (liana_res["source"].isin(tam_types) & liana_res["target"].isin(tex_states))
    )
    tex_tam_interactions = liana_res[mask].copy()
    print(f"T cell–TAM interactions: {len(tex_tam_interactions)}")

    # 4. Save
    tex_tam_interactions.to_csv(PROCESSED_DIR / "liana_tex_tam_interactions.csv", index=False)
    liana_res.to_csv(PROCESSED_DIR / "liana_all_interactions.csv", index=False)

    # 5. Per-patient communication (for tensor construction)
    if "patient_id" in combined.obs.columns:
        patient_results = {}
        for pt in combined.obs["patient_id"].unique():
            pt_data = combined[combined.obs["patient_id"] == pt].copy()
            if pt_data.n_obs < 50:
                continue
            li.mt.rank_aggregate(pt_data, groupby="cell_state", resource_name="consensus", use_raw=False, verbose=False)
            patient_results[pt] = pt_data.uns["liana_res"]
            print(f"  Patient {pt}: {len(pt_data.uns['liana_res'])} interactions")

        # Save per-patient results
        import pickle
        with open(PROCESSED_DIR / "liana_per_patient.pkl", "wb") as f:
            pickle.dump(patient_results, f)

    print("\n✓ Communication analysis complete.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Write communication tensor builder**

```python
# code/analysis/06_communication_tensor.py
"""Build 3D communication tensor from LIANA/CellChat results."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import pickle
import numpy as np
import pandas as pd
from config import PROCESSED_DIR, check_data_mount

def build_tensor(patient_results: dict, tex_states: list, tam_types: list, top_n_lr: int = 100):
    """
    Build [patients × (T_states × TAM_types) × LR_pairs] tensor.
    Returns tensor and metadata.
    """
    # Identify top L-R pairs across all patients
    all_lr = []
    for pt, df in patient_results.items():
        lr_col = df["ligand_complex"].astype(str) + "_" + df["receptor_complex"].astype(str)
        all_lr.extend(lr_col.tolist())
    lr_counts = pd.Series(all_lr).value_counts()
    top_lr_pairs = lr_counts.head(top_n_lr).index.tolist()

    # Build tensor
    patients = sorted(patient_results.keys())
    n_patients = len(patients)
    n_tex = len(tex_states)
    n_tam = len(tam_types)
    n_lr = len(top_lr_pairs)

    tensor = np.zeros((n_patients, n_tex * n_tam, n_lr))

    for i, pt in enumerate(patients):
        df = patient_results[pt]
        df["lr_pair"] = df["ligand_complex"].astype(str) + "_" + df["receptor_complex"].astype(str)

        for j_tex, tex in enumerate(tex_states):
            for j_tam, tam in enumerate(tam_types):
                flat_idx = j_tex * n_tam + j_tam
                # Both directions
                mask = (
                    ((df["source"] == tex) & (df["target"] == tam)) |
                    ((df["source"] == tam) & (df["target"] == tex))
                )
                sub = df[mask]
                for k, lr in enumerate(top_lr_pairs):
                    lr_mask = sub["lr_pair"] == lr
                    if lr_mask.any():
                        # Use aggregate_rank (lower = stronger interaction)
                        rank_val = sub.loc[lr_mask, "magnitude_rank"].values[0]
                        tensor[i, flat_idx, k] = 1.0 / (rank_val + 1)  # invert rank to strength

    metadata = {
        "patients": patients,
        "tex_states": tex_states,
        "tam_types": tam_types,
        "lr_pairs": top_lr_pairs,
        "shape_desc": "(patients, tex*tam_pairs, lr_pairs)",
    }
    return tensor, metadata

def main():
    check_data_mount()

    with open(PROCESSED_DIR / "liana_per_patient.pkl", "rb") as f:
        patient_results = pickle.load(f)
    print(f"Patients with communication data: {len(patient_results)}")

    # Get state lists from data
    import scanpy as sc
    cd8 = sc.read_h5ad(PROCESSED_DIR / "cd8_exhaustion.h5ad")
    tam = sc.read_h5ad(PROCESSED_DIR / "tam_subtyped.h5ad")
    tex_states = sorted(cd8.obs["exhaustion_state"].unique().tolist())
    tam_types = sorted(tam.obs["tam_subtype"].unique().tolist())

    tensor, metadata = build_tensor(patient_results, tex_states, tam_types, top_n_lr=100)
    print(f"Communication tensor shape: {tensor.shape}")
    print(f"  Patients: {tensor.shape[0]}")
    print(f"  Cell state pairs: {tensor.shape[1]} ({len(tex_states)} Tex × {len(tam_types)} TAM)")
    print(f"  L-R pairs: {tensor.shape[2]}")

    # Save
    np.save(PROCESSED_DIR / "communication" / "comm_tensor.npy", tensor)
    pd.to_pickle(metadata, PROCESSED_DIR / "communication" / "comm_metadata.pkl")
    print("\n✓ Communication tensor saved.")

if __name__ == "__main__":
    import os
    os.makedirs(PROCESSED_DIR / "communication", exist_ok=True)
    main()
```

- [ ] **Step 4: Run sequentially and commit**

```bash
python code/analysis/03_communication_liana.py
python code/analysis/06_communication_tensor.py
```

```bash
git add code/analysis/03_cellchat.R code/analysis/03_communication_liana.py code/analysis/06_communication_tensor.py
git commit -m "feat: cell-cell communication analysis and tensor construction"
```

---

## Task 9: CACE Model Architecture

**Files:**
- Create: `code/models/cace_model.py`

- [ ] **Step 1: Write CACE model**

```python
# code/models/cace_model.py
"""Cross-Attention Communication Encoder (CACE)."""
import torch
import torch.nn as nn
import torch.nn.functional as F
import math


class CommunicationEncoder(nn.Module):
    """Encode L-R interaction features for a cell state pair."""
    def __init__(self, n_lr: int, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(n_lr, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, d_model),
            nn.LayerNorm(d_model),
        )

    def forward(self, x):
        # x: (batch, n_lr)
        return self.proj(x)


class CrossAttentionBlock(nn.Module):
    """Multi-head cross-attention: T cell states attend to TAM subtypes."""
    def __init__(self, d_model: int, n_heads: int, dropout: float = 0.1):
        super().__init__()
        self.attn = nn.MultiheadAttention(d_model, n_heads, dropout=dropout, batch_first=True)
        self.norm1 = nn.LayerNorm(d_model)
        self.ffn = nn.Sequential(
            nn.Linear(d_model, d_model * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 4, d_model),
            nn.Dropout(dropout),
        )
        self.norm2 = nn.LayerNorm(d_model)

    def forward(self, query, key_value):
        # query: (batch, n_tex, d_model) — T cell states
        # key_value: (batch, n_tam, d_model) — TAM subtypes
        attn_out, attn_weights = self.attn(query, key_value, key_value)
        x = self.norm1(query + attn_out)
        x = self.norm2(x + self.ffn(x))
        return x, attn_weights


class CACE(nn.Module):
    """
    Cross-Attention Communication Encoder.

    Input: communication tensor (batch, n_tex * n_tam, n_lr)
    Output: communication embedding + task-specific predictions
    """
    def __init__(self, n_tex: int, n_tam: int, n_lr: int,
                 d_model: int = 128, n_heads: int = 8, n_layers: int = 2,
                 dropout: float = 0.3):
        super().__init__()
        self.n_tex = n_tex
        self.n_tam = n_tam
        self.n_lr = n_lr
        self.d_model = d_model

        # Encode each cell state pair's L-R vector
        self.lr_encoder = CommunicationEncoder(n_lr, d_model, dropout)

        # Learnable positional embeddings for T cell states and TAM subtypes
        self.tex_pos = nn.Embedding(n_tex, d_model)
        self.tam_pos = nn.Embedding(n_tam, d_model)

        # Cross-attention layers
        self.cross_attn_layers = nn.ModuleList([
            CrossAttentionBlock(d_model, n_heads, dropout)
            for _ in range(n_layers)
        ])

        # Pool to single embedding
        self.pool = nn.Sequential(
            nn.Linear(d_model * n_tex, d_model),
            nn.LayerNorm(d_model),
            nn.GELU(),
        )

        # Task heads
        self.prognostic_head = nn.Linear(d_model, 1)  # Cox PH risk score
        self.icb_head = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),  # BCEWithLogits
        )

    def encode(self, x):
        """
        x: (batch, n_tex * n_tam, n_lr)
        Returns: communication embedding (batch, d_model), attention weights
        """
        batch_size = x.shape[0]

        # Reshape to (batch, n_tex, n_tam, n_lr)
        x = x.view(batch_size, self.n_tex, self.n_tam, self.n_lr)

        # Encode each (tex, tam) pair's L-R vector
        # → (batch, n_tex, n_tam, d_model)
        encoded = self.lr_encoder(x)

        # Build T cell query: aggregate over TAM dim for each T state
        # → (batch, n_tex, d_model)
        tex_repr = encoded.mean(dim=2)  # avg over TAM types
        tex_repr = tex_repr + self.tex_pos(torch.arange(self.n_tex, device=x.device))

        # Build TAM key/value: aggregate over Tex dim for each TAM type
        # → (batch, n_tam, d_model)
        tam_repr = encoded.mean(dim=1)  # avg over Tex states
        tam_repr = tam_repr + self.tam_pos(torch.arange(self.n_tam, device=x.device))

        # Cross-attention: T cell states attend to TAM subtypes
        all_attn_weights = []
        for layer in self.cross_attn_layers:
            tex_repr, attn_w = layer(tex_repr, tam_repr)
            all_attn_weights.append(attn_w)

        # Pool: flatten T cell states → single embedding
        pooled = tex_repr.view(batch_size, -1)  # (batch, n_tex * d_model)
        embedding = self.pool(pooled)  # (batch, d_model)

        return embedding, all_attn_weights

    def forward(self, x, task="both"):
        embedding, attn_weights = self.encode(x)

        out = {"embedding": embedding, "attn_weights": attn_weights}

        if task in ("prognostic", "both"):
            out["risk_score"] = self.prognostic_head(embedding).squeeze(-1)
        if task in ("icb", "both"):
            out["icb_logit"] = self.icb_head(embedding).squeeze(-1)

        return out

    def get_communication_score(self, x):
        """Single scalar summary of communication pattern."""
        embedding, _ = self.encode(x)
        return self.prognostic_head(embedding).squeeze(-1)
```

- [ ] **Step 2: Quick test**

```bash
cd /home/laugh/no_exp_paper
python -c "
import torch
from code.models.cace_model import CACE
model = CACE(n_tex=3, n_tam=6, n_lr=100, d_model=128, n_heads=8, n_layers=2)
x = torch.randn(4, 18, 100)  # batch=4, 3*6=18 pairs, 100 LR
out = model(x)
print(f'Embedding: {out[\"embedding\"].shape}')
print(f'Risk score: {out[\"risk_score\"].shape}')
print(f'ICB logit: {out[\"icb_logit\"].shape}')
print(f'Attention weights: {len(out[\"attn_weights\"])} layers, shape={out[\"attn_weights\"][0].shape}')
print('✓ CACE model works.')
"
```
Expected: Embedding (4, 128), Risk score (4,), ICB logit (4,), Attention weights shape (4, 8, 3, 6).

- [ ] **Step 3: Commit**

```bash
git add code/models/cace_model.py
git commit -m "feat: CACE model — cross-attention communication encoder architecture"
```

---

## Task 10: Deconvolution — Bulk to Communication Score Transfer

**Files:**
- Create: `code/models/deconvolution.py`

- [ ] **Step 1: Write deconvolution module**

```python
# code/models/deconvolution.py
"""Transfer scRNA-seq communication reference to bulk RNA-seq."""
import numpy as np
import pandas as pd
from pathlib import Path


def build_signature_matrix(cd8_adata, tam_adata, n_markers: int = 200):
    """Build CIBERSORTx-style signature matrix from scRNA-seq reference."""
    import scanpy as sc
    import anndata as ad

    combined = ad.concat([cd8_adata, tam_adata], join="inner")
    combined.obs["deconv_label"] = pd.Categorical(
        combined.obs.get("exhaustion_state", combined.obs.get("tam_subtype", "unknown"))
    )

    # Find marker genes per state
    sc.tl.rank_genes_groups(combined, groupby="deconv_label", method="wilcoxon", n_genes=n_markers)
    markers = pd.DataFrame(combined.uns["rank_genes_groups"]["names"])

    # Build signature: mean expression per cell state
    sig = pd.DataFrame(index=combined.var_names)
    for label in combined.obs["deconv_label"].cat.categories:
        subset = combined[combined.obs["deconv_label"] == label]
        sig[label] = np.array(subset.X.mean(axis=0)).flatten()

    # Keep only marker genes
    all_markers = set()
    for col in markers.columns:
        all_markers.update(markers[col].tolist())
    sig = sig.loc[sig.index.isin(all_markers)]
    return sig


def estimate_fractions_nnls(bulk_expr: pd.DataFrame, signature: pd.DataFrame) -> pd.DataFrame:
    """Estimate cell state fractions using NNLS (no CIBERSORTx dependency)."""
    from scipy.optimize import nnls

    common_genes = bulk_expr.index.intersection(signature.index)
    bulk_sub = bulk_expr.loc[common_genes].values
    sig_sub = signature.loc[common_genes].values

    fractions = np.zeros((bulk_sub.shape[1], sig_sub.shape[1]))
    for i in range(bulk_sub.shape[1]):
        coef, _ = nnls(sig_sub, bulk_sub[:, i])
        fractions[i] = coef / (coef.sum() + 1e-8)

    return pd.DataFrame(
        fractions,
        index=bulk_expr.columns,
        columns=signature.columns,
    )


def build_pseudo_communication(fractions: pd.DataFrame, comm_reference: dict,
                               tex_states: list, tam_types: list, lr_pairs: list) -> np.ndarray:
    """
    Build pseudo-communication matrix for bulk samples.

    comm_reference: mean communication strength per (tex, tam, lr) from scRNA-seq
    fractions: estimated cell state proportions per bulk sample
    """
    n_samples = len(fractions)
    n_tex = len(tex_states)
    n_tam = len(tam_types)
    n_lr = len(lr_pairs)

    tensor = np.zeros((n_samples, n_tex * n_tam, n_lr))

    for i in range(n_samples):
        for j_tex, tex in enumerate(tex_states):
            for j_tam, tam in enumerate(tam_types):
                flat_idx = j_tex * n_tam + j_tam
                tex_frac = fractions.iloc[i].get(tex, 0)
                tam_frac = fractions.iloc[i].get(tam, 0)

                for k, lr in enumerate(lr_pairs):
                    ref_strength = comm_reference.get((tex, tam, lr), 0)
                    # Pseudo-communication = fraction product × reference strength
                    tensor[i, flat_idx, k] = tex_frac * tam_frac * ref_strength

    return tensor
```

- [ ] **Step 2: Commit**

```bash
git add code/models/deconvolution.py
git commit -m "feat: bulk deconvolution — NNLS fractions and pseudo-communication tensor"
```

---

## Task 11: Training Pipeline

**Files:**
- Create: `code/models/train.py`

- [ ] **Step 1: Write training script**

```python
# code/models/train.py
"""CACE training: multi-task (Cox PH + ICB response)."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import KFold
from config import PROCESSED_DIR, CHECKPOINT_DIR, CACE_PARAMS, SEED, check_data_mount
from models.cace_model import CACE

torch.manual_seed(SEED)
np.random.seed(SEED)


def cox_ph_loss(risk_scores, times, events):
    """Negative partial log-likelihood for Cox PH model."""
    sorted_idx = torch.argsort(times, descending=True)
    risk_scores = risk_scores[sorted_idx]
    events = events[sorted_idx]

    log_cumsum_exp = torch.logcumsumexp(risk_scores, dim=0)
    loss = -torch.mean((risk_scores - log_cumsum_exp) * events)
    return loss


def train_epoch(model, loader, optimizer, task="prognostic", device="cpu"):
    model.train()
    total_loss = 0
    for batch in loader:
        if task == "prognostic":
            x, times, events = [b.to(device) for b in batch]
            out = model(x, task="prognostic")
            loss = cox_ph_loss(out["risk_score"], times, events)
        elif task == "icb":
            x, labels = [b.to(device) for b in batch]
            out = model(x, task="icb")
            loss = nn.functional.binary_cross_entropy_with_logits(out["icb_logit"], labels)

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)


@torch.no_grad()
def evaluate(model, loader, task="prognostic", device="cpu"):
    model.eval()
    total_loss = 0
    preds, labels_all = [], []
    for batch in loader:
        if task == "prognostic":
            x, times, events = [b.to(device) for b in batch]
            out = model(x, task="prognostic")
            loss = cox_ph_loss(out["risk_score"], times, events)
            preds.append(out["risk_score"].cpu())
        elif task == "icb":
            x, labels = [b.to(device) for b in batch]
            out = model(x, task="icb")
            loss = nn.functional.binary_cross_entropy_with_logits(out["icb_logit"], labels)
            preds.append(torch.sigmoid(out["icb_logit"]).cpu())
            labels_all.append(labels.cpu())
        total_loss += loss.item()
    return total_loss / len(loader), torch.cat(preds), torch.cat(labels_all) if labels_all else None


def main():
    check_data_mount()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    # Load communication tensor and clinical data
    comm_tensor = np.load(PROCESSED_DIR / "communication" / "comm_tensor.npy")
    metadata = pd.read_pickle(PROCESSED_DIR / "communication" / "comm_metadata.pkl")

    n_tex = len(metadata["tex_states"])
    n_tam = len(metadata["tam_types"])
    n_lr = len(metadata["lr_pairs"])

    # TODO: Load TCGA-LIHC pseudo-communication tensor and survival data
    # TODO: Load GSE285963 pseudo-communication tensor and response labels
    # (These will be created by deconvolution.py applied to bulk data)

    # Placeholder: demonstrate training loop structure
    print(f"Tensor shape: {comm_tensor.shape}")
    print(f"n_tex={n_tex}, n_tam={n_tam}, n_lr={n_lr}")

    model = CACE(
        n_tex=n_tex, n_tam=n_tam, n_lr=n_lr,
        d_model=CACE_PARAMS["d_model"],
        n_heads=CACE_PARAMS["n_heads"],
        n_layers=CACE_PARAMS["n_layers"],
        dropout=CACE_PARAMS["dropout"],
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=CACE_PARAMS["lr"],
        weight_decay=CACE_PARAMS["weight_decay"],
    )

    print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")
    print("✓ Training pipeline ready. Awaiting bulk deconvolution data.")

    # Save model config
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    torch.save({
        "model_config": {
            "n_tex": n_tex, "n_tam": n_tam, "n_lr": n_lr,
            **CACE_PARAMS,
        },
    }, CHECKPOINT_DIR / "model_config.pt")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add code/models/train.py
git commit -m "feat: CACE training pipeline — Cox PH + ICB multi-task"
```

---

## Task 12: Validation Scripts

**Files:**
- Create: `code/validation/01_prognostic.py`
- Create: `code/validation/02_icb_response.py`
- Create: `code/validation/04_ablation.py`
- Create: `code/utils/plotting.py`

- [ ] **Step 1: Write plotting utilities**

```python
# code/utils/plotting.py
"""Common plotting functions."""
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42  # editable text in PDF
matplotlib.rcParams["font.size"] = 10
import numpy as np
import pandas as pd


def km_plot(times, events, groups, group_labels=None, ax=None, title="", colors=None):
    """Kaplan-Meier survival plot with log-rank p-value."""
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))

    unique_groups = sorted(set(groups))
    if colors is None:
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_groups)))

    for i, g in enumerate(unique_groups):
        mask = np.array(groups) == g
        kmf = KaplanMeierFitter()
        label = group_labels[i] if group_labels else str(g)
        kmf.fit(times[mask], events[mask], label=f"{label} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[i])

    # Log-rank test (if 2 groups)
    if len(unique_groups) == 2:
        m0 = np.array(groups) == unique_groups[0]
        m1 = np.array(groups) == unique_groups[1]
        result = logrank_test(times[m0], times[m1], events[m0], events[m1])
        ax.text(0.95, 0.95, f"p = {result.p_value:.2e}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9)

    ax.set_title(title)
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    return ax


def roc_plot(y_true, y_score, ax=None, label="", color=None):
    """ROC curve with AUC."""
    from sklearn.metrics import roc_curve, auc
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(fpr, tpr, color=color, label=f"{label} (AUC={roc_auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend(loc="lower right")
    return ax, roc_auc
```

- [ ] **Step 2: Write prognostic validation**

```python
# code/validation/01_prognostic.py
"""Prognostic validation on TCGA-LIHC and GSE14520."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import torch
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.metrics import roc_auc_score
from config import PROCESSED_DIR, EXTERNAL_DIR, CHECKPOINT_DIR, RESULTS_DIR, check_data_mount
from utils.plotting import km_plot
import matplotlib.pyplot as plt

def main():
    check_data_mount()
    # 1. Load trained CACE model
    # 2. Load TCGA-LIHC pseudo-communication tensor
    # 3. Compute communication scores
    # 4. KM analysis (median split)
    # 5. Multivariate Cox (score + stage + AFP + etiology)
    # 6. Time-dependent AUC (1yr, 3yr, 5yr)
    # 7. Compare vs TIDE, ESTIMATE, GEP
    # 8. Repeat on GSE14520
    # 9. Save figures to results/figures/main/
    print("Prognostic validation — awaiting trained model and deconvolution data.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Write ICB response validation**

```python
# code/validation/02_icb_response.py
"""ICB response prediction validation on GSE285963 + GSE140901."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from config import PROCESSED_DIR, ICB_DIR, RESULTS_DIR, check_data_mount

def main():
    check_data_mount()
    # 1. Load GSE285963 (RNA-seq) → deconvolution → pseudo-communication
    # 2. Compute communication scores via CACE
    # 3. ROC-AUC for response prediction
    # 4. Compare vs TIDE, GEP, TMB
    # 5. Waterfall plot: communication score vs response
    # 6. Attention weight visualization
    # 7. Repeat on GSE140901 (NanoString subset)
    # 8. Save figures
    print("ICB validation — awaiting trained model and deconvolution data.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Write ablation study**

```python
# code/validation/04_ablation.py
"""Ablation study for CACE model."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

from config import check_data_mount

def main():
    check_data_mount()
    # Ablation experiments:
    # 1. CACE (full) vs CACE without cross-attention (simple concat + MLP)
    # 2. Communication features vs cell abundance features only
    # 3. Vary n_tex (2, 3, 4, 5) and n_tam (3, 4, 5, 6)
    # 4. Random L-R permutation test (shuffle L-R pair labels)
    # 5. Compare vs XGBoost, Random Survival Forest on same features
    # 6. Report: bar chart of AUC across ablation conditions
    print("Ablation study — awaiting trained model.")

if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Commit**

```bash
git add code/validation/ code/utils/plotting.py
git commit -m "feat: validation scripts — prognostic, ICB, ablation (skeleton)"
```

---

## Task 13: Interpretation Module

**Files:**
- Create: `code/models/interpret.py`

- [ ] **Step 1: Write interpretation script**

```python
# code/models/interpret.py
"""CACE interpretation: attention weights, SHAP, communication axis ranking."""
import torch
import numpy as np
import pandas as pd
import shap


def extract_attention_weights(model, tensor, metadata):
    """Extract and average attention weights across samples."""
    model.eval()
    with torch.no_grad():
        x = torch.tensor(tensor, dtype=torch.float32)
        out = model(x)
        attn_weights = out["attn_weights"]

    # Average over samples and heads → (n_tex, n_tam)
    avg_attn = torch.stack(attn_weights).mean(dim=(0, 1, 2))  # avg layers, batch, heads
    attn_df = pd.DataFrame(
        avg_attn.numpy(),
        index=metadata["tex_states"],
        columns=metadata["tam_types"],
    )
    return attn_df


def rank_lr_pairs_shap(model, tensor, metadata, n_background=50):
    """Rank L-R pairs by SHAP importance."""
    model.eval()
    x_np = tensor.reshape(tensor.shape[0], -1)  # flatten to 2D

    background = x_np[np.random.choice(len(x_np), min(n_background, len(x_np)), replace=False)]

    def predict_fn(x_flat):
        x_t = torch.tensor(x_flat.reshape(-1, tensor.shape[1], tensor.shape[2]), dtype=torch.float32)
        with torch.no_grad():
            return model(x_t, task="prognostic")["risk_score"].numpy()

    explainer = shap.KernelExplainer(predict_fn, background)
    shap_values = explainer.shap_values(x_np)

    # Reshape SHAP values to (samples, n_tex*n_tam, n_lr) and sum over state pairs
    shap_3d = shap_values.reshape(tensor.shape)
    lr_importance = np.abs(shap_3d).mean(axis=(0, 1))  # avg over samples and state pairs

    lr_ranking = pd.DataFrame({
        "lr_pair": metadata["lr_pairs"],
        "shap_importance": lr_importance,
    }).sort_values("shap_importance", ascending=False)

    return lr_ranking
```

- [ ] **Step 2: Commit**

```bash
git add code/models/interpret.py
git commit -m "feat: CACE interpretation — attention weights and SHAP ranking"
```

---

## Task 14: pySCENIC Analysis (TF Activity)

**Files:**
- Create: `code/analysis/05_scenic.py`

- [ ] **Step 1: Write SCENIC script**

```python
# code/analysis/05_scenic.py
"""pySCENIC: TF activity for exhaustion states and TAM subtypes."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

import scanpy as sc
import pandas as pd
import numpy as np
from config import PROCESSED_DIR, SEED, check_data_mount

def run_scenic_pipeline(adata, label_key: str, out_prefix: str):
    """Run pySCENIC on a subset. Requires pre-downloaded databases."""
    # pySCENIC is resource-intensive; this is a template
    # In practice, run via CLI:
    # 1. pyscenic grn {loom} {tfs} -o adj.csv
    # 2. pyscenic ctx adj.csv {db} --annotations_fname {motifs} -o reg.csv
    # 3. pyscenic aucell {loom} reg.csv -o auc.loom
    print(f"SCENIC pipeline for {out_prefix} — run via CLI (see comments in script)")
    print("Required files:")
    print("  - TF list: https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt")
    print("  - Ranking DB: https://resources.aertslab.org/cistarget/databases/")
    print("  - Motif annotations: https://resources.aertslab.org/cistarget/motif2tf/")

def main():
    check_data_mount()
    cd8 = sc.read_h5ad(PROCESSED_DIR / "cd8_exhaustion.h5ad")
    tam = sc.read_h5ad(PROCESSED_DIR / "tam_subtyped.h5ad")

    run_scenic_pipeline(cd8, "exhaustion_state", "cd8_tex")
    run_scenic_pipeline(tam, "tam_subtype", "tam")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add code/analysis/05_scenic.py
git commit -m "feat: pySCENIC TF activity analysis template"
```

---

## Task 15: NicheNet Analysis

**Files:**
- Create: `code/analysis/04_nichenet.R`

- [ ] **Step 1: Write NicheNet R script**

```r
# code/analysis/04_nichenet.R
# NicheNet: TAM ligand → T cell exhaustion gene regulatory potential
# Run: Rscript code/analysis/04_nichenet.R

library(nichenetr)
library(tidyverse)

cat("=== NicheNet Analysis ===\n")
cat("This script identifies TAM-derived ligands that regulate\n")
cat("T cell exhaustion gene programs.\n\n")

# Load NicheNet networks
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

cat("NicheNet reference networks loaded.\n")
cat("Proceed with TAM sender → T cell receiver analysis.\n")
cat("(Requires marker gene lists from Python analysis — export from h5ad)\n")
```

- [ ] **Step 2: Commit**

```bash
git add code/analysis/04_nichenet.R
git commit -m "feat: NicheNet TAM→T cell ligand analysis template"
```

---

## Execution Dependencies

```
Task 1 (setup) → Task 2 (download/load) → Task 3 (QC) → Task 4 (integration) → Task 5 (annotation)
                                                                                      ↓
                                                              Task 6 (T cell) + Task 7 (TAM) [parallel]
                                                                           ↓
                                                              Task 8 (communication) → Task 10 (deconvolution)
                                                                           ↓                    ↓
                                                              Task 14 (SCENIC)      Task 9 (CACE model)
                                                              Task 15 (NicheNet)         ↓
                                                                                    Task 11 (training)
                                                                                         ↓
                                                                           Task 12 (validation) + Task 13 (interpretation) [parallel]
```

## Checkpoint Reviews

**Checkpoint 1 (after Task 5):** Review integrated UMAP, cell type proportions, batch mixing. Go/no-go for downstream analysis.

**Checkpoint 2 (after Tasks 6+7):** Review exhaustion states and TAM subtypes. Confirm biologically meaningful separation. Adjust markers/resolution if needed.

**Checkpoint 3 (after Task 8):** Review communication tensor. Confirm T cell–TAM interactions are non-trivial. Check patient-level heterogeneity.

**Checkpoint 4 (after Task 11):** Review CACE training. Check overfitting, ablation results. Decide if model architecture needs changes.

**Checkpoint 5 (after Task 12):** Review all validation results. Go/no-go for manuscript writing.
