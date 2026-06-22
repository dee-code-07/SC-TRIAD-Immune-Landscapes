#!/usr/bin/env python3
# sc-triad project
# script: 01_grn.py
# purpose: GRNBoost2 co-expression network inference (pyscenic step 1)
#
# approach: use pyscenic CLI subprocess instead of arboreto Python API directly
#
# why CLI instead of Python API:
#   arboreto 0.1.6 is incompatible with dask >= 2023.x due to the migration
#   from dask.dataframe to dask_expr. the breakage is inside arboreto/core.py
#   at the from_delayed() call — no client-side fix resolves this.
#   the pyscenic CLI uses the same arboreto internally but with a pinned
#   internal dask call path that works correctly when invoked via subprocess.
#   alternatively: pin dask==2022.10.0 before running (see setup notes below).
#
# setup (run once before submitting):
#   conda activate sctriad
#   pip install "dask[distributed]==2022.10.0" "distributed==2022.10.0"
#
# GRNBoost2 is a tree-based (gradient boosting) co-expression method.
# it does NOT infer directionality — that comes from cisTarget in step 2.
# reference: moerman et al., bioinformatics 2019
#
# runtime: 4-8 hours on 16 cores for 155k cells
# memory:  ~64-128 GB
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

import os
import sys
import logging
import subprocess
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import loompy

# ── logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "[%(asctime)s] %(levelname)s %(message)s",
    datefmt = "%H:%M:%S"
)
log = logging.getLogger("sc-triad-grn")

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    = Path.home() / "sc-triad" / "03_pyscenic"
LOOM    = BASE / "input"   / "triad_pbmc_raw.loom"
TF_LIST = BASE / "databases" / "allTFs_hg38.txt"
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ADJ_OUT = OUT_DIR / "01_adjacencies.tsv"

# ── workers ────────────────────────────────────────────────────────────────────
N_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
log.info(f"workers: {N_WORKERS}")


def verify_inputs():
    """Check all required input files exist before starting."""
    missing = []
    for f in [LOOM, TF_LIST]:
        if not f.exists():
            missing.append(str(f))
    if missing:
        for m in missing:
            log.error(f"  missing: {m}")
        raise FileNotFoundError(f"Missing input files: {missing}")
    log.info(f"  loom: {LOOM} ({LOOM.stat().st_size / 1e9:.2f} GB)")
    log.info(f"  TF list: {TF_LIST} ({sum(1 for _ in open(TF_LIST))} TFs)")


def verify_loom_structure():
    """
    Verify loom has correct structure for pyscenic.
    pyscenic CLI requires:
      - row attribute 'Gene' containing gene symbols
      - col attribute 'CellID' containing cell barcodes
    """
    log.info(f"verifying loom structure: {LOOM}")
    with loompy.connect(str(LOOM), mode="r") as ds:
        log.info(f"  shape (genes x cells): {ds.shape}")
        log.info(f"  row attrs: {list(ds.ra.keys())}")
        log.info(f"  col attrs: {list(ds.ca.keys())}")

        if "Gene" not in ds.ra:
            raise ValueError(
                "Loom missing 'Gene' row attribute. "
                "pyscenic CLI requires row attribute named exactly 'Gene'."
            )
        if "CellID" not in ds.ca:
            raise ValueError(
                "Loom missing 'CellID' col attribute. "
                "pyscenic CLI requires col attribute named exactly 'CellID'."
            )

        n_genes = ds.shape[0]
        n_cells = ds.shape[1]
        log.info(f"  genes: {n_genes} | cells: {n_cells}")
        log.info(f"  sample genes: {list(ds.ra['Gene'][:5])}")
        log.info(f"  sample cells: {list(ds.ca['CellID'][:3])}")

    return n_genes, n_cells


def run_pyscenic_grn():
    """
    Run GRNBoost2 via pyscenic CLI.

    pyscenic grn command:
      pyscenic grn <loom> <tf_list> -o <output> --num_workers <n>

    the CLI handles all dask setup internally and is tested against
    the installed pyscenic version — bypasses arboreto API breakage.
    """
    log.info("running pyscenic grn (GRNBoost2 via CLI)...")
    log.info(f"  loom:     {LOOM}")
    log.info(f"  TF list:  {TF_LIST}")
    log.info(f"  output:   {ADJ_OUT}")
    log.info(f"  workers:  {N_WORKERS}")
    log.info(f"  expected runtime: 4-8 hours")

    cmd = [
        "pyscenic", "grn",
        str(LOOM),
        str(TF_LIST),
        "--output",       str(ADJ_OUT),
        "--num_workers",  str(N_WORKERS),
        "--seed",         "42",
        "--method",       "grnboost2",
    ]

    log.info(f"command: {' '.join(cmd)}")

    # run with output streaming to log in real time
    # stdout and stderr both captured so progress is visible in .out file
    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.STDOUT,
        text   = True,
        bufsize = 1,
    )

    # stream output line by line
    for line in process.stdout:
        line = line.rstrip()
        if line:
            log.info(f"  [pyscenic] {line}")

    process.wait()

    if process.returncode != 0:
        raise RuntimeError(
            f"pyscenic grn failed with exit code {process.returncode}. "
            "Check the log output above for details."
        )

    log.info("pyscenic grn completed successfully")


def summarise_adjacencies():
    """Load and summarise the output adjacency file."""
    if not ADJ_OUT.exists():
        raise FileNotFoundError(
            f"Adjacency file not created: {ADJ_OUT}\n"
            "pyscenic grn may have failed silently — check logs above."
        )

    log.info(f"loading adjacencies for summary: {ADJ_OUT}")
    adj = pd.read_csv(str(ADJ_OUT), sep="\t")
    log.info(f"  total TF-target links: {len(adj):,}")
    log.info(f"  unique TFs: {adj['TF'].nunique()}")
    log.info(f"  unique targets: {adj['target'].nunique()}")
    log.info(f"  importance score range: {adj['importance'].min():.4f} - {adj['importance'].max():.4f}")
    log.info(f"  file size: {ADJ_OUT.stat().st_size / 1e6:.1f} MB")

    # top 20 TFs by mean importance
    top_tfs = (
        adj.groupby("TF")["importance"]
        .mean()
        .sort_values(ascending=False)
        .head(20)
    )
    log.info("top 20 TFs by mean importance:")
    for tf, imp in top_tfs.items():
        log.info(f"  {tf}: {imp:.4f}")

    # check for biologically expected TFs
    expected_tfs = ["TBX21", "EOMES", "IRF4", "STAT3", "IRF1",
                    "SPI1", "CEBPB", "IRF8", "RELA", "TCF7"]
    found = [tf for tf in expected_tfs if tf in adj["TF"].values]
    missing = [tf for tf in expected_tfs if tf not in adj["TF"].values]
    log.info(f"expected TFs found: {found}")
    if missing:
        log.warning(f"expected TFs NOT found: {missing}")
        log.warning("  this may indicate gene name issues or low expression")


def main():
    log.info("=== SC-TRIAD GRNBoost2 pipeline (pyscenic CLI) ===")
    log.info(f"started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # verify environment
    try:
        import pyscenic
        log.info(f"pyscenic version: {pyscenic.__version__}")
    except ImportError:
        raise ImportError("pyscenic not installed. Run: pip install pyscenic==0.12.1")

    import dask
    log.info(f"dask version: {dask.__version__}")
    # warn if dask version may be incompatible with arboreto Python API
    # (CLI path is unaffected but good to document)
    dask_major = int(dask.__version__.split(".")[0])
    dask_minor = int(dask.__version__.split(".")[1]) if len(dask.__version__.split(".")) > 1 else 0
    if dask_major >= 2023 or (dask_major == 2022 and dask_minor >= 11):
        log.warning(
            f"dask {dask.__version__} detected. arboreto Python API is incompatible "
            "with dask >= 2022.11. Using pyscenic CLI (unaffected by this issue)."
        )

    # run pipeline
    verify_inputs()
    n_genes, n_cells = verify_loom_structure()
    log.info(f"proceeding with {n_genes} genes x {n_cells} cells")

    run_pyscenic_grn()
    summarise_adjacencies()

    log.info(f"finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("next step: submit 02_ctx.sh (cisTarget pruning)")


if __name__ == "__main__":
    main()