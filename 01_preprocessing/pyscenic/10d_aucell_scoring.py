#!/usr/bin/env python3
# sc-triad project
# script: 03_aucell.py
# purpose: AUCell per-cell regulon activity scoring (pyscenic step 3)
#
# what this does:
#   loads the pruned regulons from step 2
#   scores each cell for each regulon's activity using the AUCell algorithm
#   AUCell ranks each cell's expressed genes and tests whether regulon genes
#   are enriched at the top of that ranking (like a per-cell GSEA)
#   output: cells x regulons activity score matrix
#
# the AUC score for each cell x regulon pair represents the fraction of
# the top-ranked genes in that cell that belong to the regulon
# a high AUC = the TF's target genes are highly expressed in that cell
# = the TF is active in that cell
#
# reference: aibar et al., nature methods 2017
#
# runtime: 1-3 hours for 155k cells x typical regulon count
# memory:  ~32 GB
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

import os
import pickle
import logging
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import loompy
from pyscenic.aucell import aucell

# ── logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "[%(asctime)s] %(levelname)s %(message)s",
    datefmt = "%H:%M:%S"
)
log = logging.getLogger("sc-triad-aucell")

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    = Path.home() / "sc-triad" / "03_pyscenic"
LOOM    = BASE / "input"  / "triad_pbmc_raw.loom"
OUT_DIR = BASE / "output"

REGULONS_PKL = OUT_DIR / "02_regulons.pkl"
AUCCELL_OUT  = OUT_DIR / "03_auc_mtx.csv"
METADATA_OUT = OUT_DIR / "03_cell_metadata.csv"

# AUCell parameter: fraction of genes in the "high-expression" tail
# 0.05 means top 5% of expressed genes per cell define the active set
# standard pyscenic default — do not change without good reason
AUC_THRESHOLD = 0.05

N_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
log.info(f"workers: {N_WORKERS}")


def load_expression_and_metadata(loom_path: Path):
    """
    Load expression matrix and cell metadata from loom.
    Returns expression as DataFrame (cells x genes) and metadata DataFrame.
    """
    log.info(f"loading loom: {loom_path}")
    with loompy.connect(str(loom_path), mode="r") as ds:
        log.info(f"  shape (genes x cells): {ds.shape}")

        # gene names
        if "Gene" in ds.ra:
            gene_names = list(ds.ra["Gene"])
        elif "var_names" in ds.ra:
            gene_names = list(ds.ra["var_names"])
        else:
            gene_names = list(ds.ra[list(ds.ra.keys())[0]])

        # cell barcodes
        if "CellID" in ds.ca:
            cell_ids = list(ds.ca["CellID"])
        elif "obs_names" in ds.ca:
            cell_ids = list(ds.ca["obs_names"])
        else:
            cell_ids = [f"cell_{i}" for i in range(ds.shape[1])]

        # expression matrix — transpose to cells x genes
        log.info("  loading count matrix...")
        ex_matrix = ds[:, :].T.astype(np.float32)
        log.info(f"  matrix: {ex_matrix.shape} (cells x genes)")

        # collect available cell metadata
        meta_dict = {"cell_id": cell_ids}
        for attr in ["cell_type", "disease", "group", "sample_label",
                     "sample", "CellType", "Disease", "Group"]:
            if attr in ds.ca:
                # clean attribute name for output
                clean_name = attr.lower().replace("celltype", "cell_type")
                meta_dict[clean_name] = list(ds.ca[attr])

        meta_df = pd.DataFrame(meta_dict)
        log.info(f"  metadata columns: {list(meta_df.columns)}")

    ex_df = pd.DataFrame(ex_matrix, index=cell_ids, columns=gene_names)
    return ex_df, meta_df


def main():
    log.info("=== SC-TRIAD AUCell regulon activity scoring ===")
    log.info(f"started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # ── load regulons ──────────────────────────────────────────────────────────
    log.info(f"loading regulons: {REGULONS_PKL}")
    if not REGULONS_PKL.exists():
        raise FileNotFoundError(
            f"Regulon pickle not found: {REGULONS_PKL}\n"
            "Run 02_ctx.py first"
        )

    with open(str(REGULONS_PKL), "rb") as f:
        regulons = pickle.load(f)

    log.info(f"  regulons loaded: {len(regulons)}")

    # filter regulons with too few targets (< 5 genes)
    # small regulons are unreliable for AUCell scoring
    regulons_filt = [r for r in regulons if len(r.genes) >= 5]
    log.info(f"  regulons >= 5 targets: {len(regulons_filt)}")

    if len(regulons_filt) == 0:
        raise ValueError("No regulons with >= 5 target genes")

    # regulon name -> gene set dict for logging
    log.info("sample regulons (name: n_targets):")
    for reg in sorted(regulons_filt, key=lambda r: len(r.genes), reverse=True)[:10]:
        log.info(f"  {reg.name}: {len(reg.genes)}")

    # ── load expression matrix ─────────────────────────────────────────────────
    ex_df, meta_df = load_expression_and_metadata(LOOM)

    # ── run AUCell ─────────────────────────────────────────────────────────────
    # AUCell is fast even for 155k cells — typically 30-90 minutes
    log.info(f"running AUCell (auc_threshold={AUC_THRESHOLD})...")
    log.info(f"  cells: {ex_df.shape[0]} | regulons: {len(regulons_filt)}")
    log.info("  expected runtime: 30-90 minutes")

    auc_mtx = aucell(
        ex_df,
        regulons_filt,
        auc_threshold = AUC_THRESHOLD,
        num_workers   = N_WORKERS,
        seed          = 42,
    )

    log.info(f"AUCell complete: {auc_mtx.shape} (cells x regulons)")
    log.info(f"  AUC score range: {auc_mtx.values.min():.4f} - {auc_mtx.values.max():.4f}")
    log.info(f"  mean AUC across all cells/regulons: {auc_mtx.values.mean():.4f}")

    # ── save outputs ───────────────────────────────────────────────────────────
    log.info(f"saving AUC matrix: {AUCCELL_OUT}")
    # save as csv — will be loaded into R for downstream analysis
    # index = cell barcodes, columns = regulon names
    auc_mtx.to_csv(str(AUCCELL_OUT))
    log.info(f"  size: {AUCCELL_OUT.stat().st_size / 1e6:.1f} MB")

    # save metadata with matching cell order
    # ensures cell_id order matches auc_mtx row order
    meta_ordered = meta_df.set_index("cell_id").reindex(auc_mtx.index).reset_index()
    meta_ordered.to_csv(str(METADATA_OUT), index=False)
    log.info(f"  metadata saved: {METADATA_OUT}")

    # ── quick sanity check ─────────────────────────────────────────────────────
    # top regulons by mean AUC across all cells
    # high mean AUC = TF broadly active across cell types
    top_regulons = auc_mtx.mean(axis=0).sort_values(ascending=False).head(20)
    log.info("top 20 regulons by mean AUC (broadly active TFs):")
    for reg, auc_val in top_regulons.items():
        log.info(f"  {reg}: {auc_val:.4f}")

    # regulons with highest variance = most cell-type-specific TFs
    var_regulons = auc_mtx.var(axis=0).sort_values(ascending=False).head(20)
    log.info("top 20 regulons by AUC variance (cell-type-specific TFs):")
    for reg, var_val in var_regulons.items():
        log.info(f"  {reg}: {var_val:.6f}")

    log.info(f"finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("next step: run 04_scenic_analysis.R in R")
    log.info(f"  input for R: {AUCCELL_OUT}")
    log.info(f"  metadata:    {METADATA_OUT}")


if __name__ == "__main__":
    main()
