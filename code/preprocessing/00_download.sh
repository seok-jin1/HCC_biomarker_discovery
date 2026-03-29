#!/usr/bin/env bash
# =============================================================================
# 00_download.sh — Download all HCC scRNA-seq datasets from GEO FTP
#
# Usage:
#   bash code/preprocessing/00_download.sh
#
# Prerequisites:
#   - External drive mounted at /mnt/e
#   - wget available in PATH
#
# Notes:
#   - Uses wget --no-clobber (-nc) so already-downloaded files are skipped.
#   - Run from project root: /home/laugh/no_exp_paper
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 1. Verify /mnt/e is mounted
# ---------------------------------------------------------------------------
if [[ ! -d /mnt/e ]]; then
    echo "[ERROR] /mnt/e is not mounted."
    echo "        Mount the external drive first:"
    echo "          sudo mount -t drvfs E: /mnt/e"
    exit 1
fi

DATA_ROOT="/mnt/e/no_exp_paper/data"
RAW_DIR="${DATA_ROOT}/raw"
EXTERNAL_DIR="${DATA_ROOT}/external"

echo "[INFO] DATA_ROOT = ${DATA_ROOT}"

# ---------------------------------------------------------------------------
# Helper: download a single URL into a target directory (no-clobber)
# ---------------------------------------------------------------------------
download() {
    local url="$1"
    local dest_dir="$2"
    mkdir -p "${dest_dir}"
    echo "  -> $(basename "${url}")"
    wget --no-clobber --progress=bar:force:noscroll \
         --directory-prefix="${dest_dir}" \
         "${url}" || true
}

# =============================================================================
# PRIMARY scRNA-seq DATASETS
# =============================================================================

# ---------------------------------------------------------------------------
# GSE149614 — HCC tumour microenvironment (10x Chromium)
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE149614 ..."
GSE149614_DIR="${RAW_DIR}/GSE149614"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl"

download "${BASE}/GSE149614_HCC.scRNAseq.S71271.count.txt.gz"    "${GSE149614_DIR}"
download "${BASE}/GSE149614_HCC.scRNAseq.S71271.metadata.txt.gz" "${GSE149614_DIR}"

# ---------------------------------------------------------------------------
# GSE140228 — CD45-sorted TILs, Droplet + Smart-seq2 (10x Chromium)
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE140228 (Droplet data only) ..."
GSE140228_DIR="${RAW_DIR}/GSE140228"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl"

download "${BASE}/GSE140228_Droplet_raw_data_matrix.mtx.gz"   "${GSE140228_DIR}"
download "${BASE}/GSE140228_Droplet_raw_data_barcodes.tsv.gz" "${GSE140228_DIR}"
download "${BASE}/GSE140228_Droplet_raw_data_genes.tsv.gz"    "${GSE140228_DIR}"
download "${BASE}/GSE140228_Droplet_cell_info.txt.gz"         "${GSE140228_DIR}"

# ---------------------------------------------------------------------------
# GSE156625 — Pre-processed HCC h5ad object
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE156625 ..."
GSE156625_DIR="${RAW_DIR}/GSE156625"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl"

download "${BASE}/GSE156625_HCCscanpyobj.h5ad.gz" "${GSE156625_DIR}"

# ---------------------------------------------------------------------------
# GSE151530 — Multi-region HCC (10x Chromium)
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE151530 ..."
GSE151530_DIR="${RAW_DIR}/GSE151530"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl"

download "${BASE}/GSE151530_matrix.mtx.gz"    "${GSE151530_DIR}"
download "${BASE}/GSE151530_barcodes.tsv.gz"  "${GSE151530_DIR}"
download "${BASE}/GSE151530_genes.tsv.gz"     "${GSE151530_DIR}"
download "${BASE}/GSE151530_CellInfo.txt.gz"  "${GSE151530_DIR}"

# =============================================================================
# VALIDATION / BULK DATASETS
# =============================================================================

# ---------------------------------------------------------------------------
# GSE285963 — Bulk RNA-seq (read count matrix)
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE285963 ..."
GSE285963_DIR="${EXTERNAL_DIR}/GSE285963"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285963/suppl"

download "${BASE}/GSE285963_readcount.gct.gz" "${GSE285963_DIR}"

# ---------------------------------------------------------------------------
# GSE140901 — Immunotherapy cohort (RAW + processed)
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE140901 ..."
GSE140901_DIR="${EXTERNAL_DIR}/GSE140901"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140901/suppl"

download "${BASE}/GSE140901_RAW.tar"                   "${GSE140901_DIR}"
download "${BASE}/GSE140901_processed_data.txt.gz"     "${GSE140901_DIR}"

# ---------------------------------------------------------------------------
# GSE14520 — HCC microarray with survival data
# ---------------------------------------------------------------------------
echo ""
echo "[INFO] Downloading GSE14520 ..."
GSE14520_DIR="${EXTERNAL_DIR}/GSE14520"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl"

download "${BASE}/GSE14520_Extra_Supplement.txt.gz" "${GSE14520_DIR}"

# =============================================================================
# Done
# =============================================================================
echo ""
echo "[INFO] All downloads complete."
echo "       Data root: ${DATA_ROOT}"
echo ""
echo "       Next step: run code/preprocessing/01_qc_filter.py"
