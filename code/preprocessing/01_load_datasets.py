"""
01_load_datasets.py — Load all 4 HCC scRNA-seq datasets into individual h5ad files.

Reads raw data from the external drive (DATA_ROOT) and writes one .h5ad per
dataset under PROCESSED_DIR / "per_dataset".

Usage
-----
    python code/preprocessing/01_load_datasets.py

NOTE: These datasets are large (total ~1.5 GB compressed).  Loading all four
in sequence requires substantial RAM.  Run on a machine with ≥32 GB RAM or
load datasets one at a time if memory is constrained.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the project root is on sys.path so "code.*" imports resolve correctly
# regardless of where the script is invoked from.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))

from code.config import check_data_mount, DATASETS, PROCESSED_DIR
from code.utils.io import (
    load_gse149614,
    load_gse140228,
    load_gse156625,
    load_gse151530,
    save_h5ad,
)


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Verify external drive is mounted
    # ------------------------------------------------------------------
    mounted = check_data_mount()
    if not mounted:
        sys.exit(1)

    # ------------------------------------------------------------------
    # 2. Create output directory
    # ------------------------------------------------------------------
    out_dir = PROCESSED_DIR / "per_dataset"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[01_load_datasets] output directory: {out_dir}\n")

    # ------------------------------------------------------------------
    # 3. Define datasets to process (name → loader, raw_dir, output path)
    # ------------------------------------------------------------------
    tasks = [
        ("GSE149614", load_gse149614, DATASETS["GSE149614"]["raw_dir"]),
        ("GSE140228", load_gse140228, DATASETS["GSE140228"]["raw_dir"]),
        ("GSE156625", load_gse156625, DATASETS["GSE156625"]["raw_dir"]),
        ("GSE151530", load_gse151530, DATASETS["GSE151530"]["raw_dir"]),
    ]

    summary: list[dict] = []

    # ------------------------------------------------------------------
    # 4. Load, inspect, and save each dataset
    # ------------------------------------------------------------------
    for name, loader_fn, raw_dir in tasks:
        print(f"{'=' * 60}")
        print(f"[01_load_datasets] Loading {name} from {raw_dir}")
        print(f"{'=' * 60}")

        adata = loader_fn(raw_dir)

        print(f"\n  Shape      : {adata.shape[0]} cells × {adata.shape[1]} genes")
        print(f"  obs columns: {list(adata.obs.columns)}\n")

        out_path = out_dir / f"{name}.h5ad"
        save_h5ad(adata, out_path)

        summary.append(
            {
                "dataset": name,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "obs_cols": list(adata.obs.columns),
                "out_path": str(out_path),
            }
        )

        # Free memory before loading the next dataset
        del adata
        print()

    # ------------------------------------------------------------------
    # 5. Print summary
    # ------------------------------------------------------------------
    print(f"{'=' * 60}")
    print("[01_load_datasets] SUMMARY")
    print(f"{'=' * 60}")
    total_cells = 0
    for entry in summary:
        print(
            f"  {entry['dataset']:<12} "
            f"{entry['n_cells']:>7,} cells × {entry['n_genes']:>6,} genes  "
            f"→ {entry['out_path']}"
        )
        total_cells += entry["n_cells"]
    print(f"{'─' * 60}")
    print(f"  {'TOTAL':<12} {total_cells:>7,} cells across {len(summary)} datasets")
    print(f"{'=' * 60}")
    print("[01_load_datasets] All datasets loaded and saved successfully.")


if __name__ == "__main__":
    main()
