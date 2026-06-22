# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# sc-triad project
# script: 02_ctx.py
# purpose: cisTarget motif enrichment and regulon pruning (pyscenic step 2)
#
# fix log:
#   ex_mtx is a required positional argument in pyscenic 0.12.1
#   must pass the full expression matrix (cells x genes DataFrame)
#   loaded from the same loom file used in step 1
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
from pyscenic.prune import prune2df, df2regulons
from pyscenic.utils import modules_from_adjacencies
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

logging.basicConfig(
    level   = logging.INFO,
    format  = "[%(asctime)s] %(levelname)s %(message)s",
    datefmt = "%H:%M:%S"
)
log = logging.getLogger("sc-triad-ctx")

BASE     = Path.home() / "sc-triad" / "03_pyscenic"
DB_DIR   = BASE / "databases"
OUT_DIR  = BASE / "output"
LOOM     = BASE / "input" / "triad_pbmc_raw.loom"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ADJ_IN       = OUT_DIR / "01_adjacencies.tsv"
MOTIF_ANN    = DB_DIR  / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
DB_10KB      = DB_DIR  / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
DB_500BP     = DB_DIR  / "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
REGULONS_OUT     = OUT_DIR / "02_regulons.csv"
REGULONS_PKL_OUT = OUT_DIR / "02_regulons.pkl"

NES_THRESHOLD = 3.0
N_WORKERS     = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
log.info(f"workers: {N_WORKERS}")


def load_expression_matrix(loom_path):
    """
    Load expression matrix from loom as DataFrame (cells x genes).
    This is required by modules_from_adjacencies for rho dichotomization --
    it uses the expression matrix to compute TF-target correlations and
    distinguish activating from repressing relationships.
    """
    log.info(f"loading expression matrix from loom: {loom_path}")
    with loompy.connect(str(loom_path), mode="r") as ds:
        log.info(f"  loom shape (genes x cells): {ds.shape}")

        # gene names
        gene_names = list(ds.ra["Gene"])
        # cell barcodes
        cell_ids   = list(ds.ca["CellID"])

        log.info(f"  loading matrix into memory (cells x genes)...")
        # loom: genes x cells -- transpose to cells x genes for pyscenic
        ex_matrix = ds[:, :].T.astype(np.float32)
        log.info(f"  matrix shape: {ex_matrix.shape}")

    ex_df = pd.DataFrame(ex_matrix, index=cell_ids, columns=gene_names)
    log.info(f"  expression DataFrame: {ex_df.shape} (cells x genes)")
    del ex_matrix
    return ex_df


def load_databases():
    dbs = []
    for db_path in [DB_10KB, DB_500BP]:
        if not db_path.exists():
            log.warning(f"database not found, skipping: {db_path.name}")
            continue
        log.info(f"loading database: {db_path.name}")
        db = RankingDatabase(fname=str(db_path), name=db_path.stem)
        dbs.append(db)
        log.info(f"  loaded: {db_path.name}")
    if not dbs:
        raise FileNotFoundError("No ranking databases found in " + str(DB_DIR))
    log.info(f"total databases loaded: {len(dbs)}")
    return dbs


def main():
    log.info("=== SC-TRIAD cisTarget regulon pruning ===")
    log.info(f"started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    import pyscenic
    log.info(f"pyscenic: {pyscenic.__version__} | numpy: {np.__version__}")

    # load adjacencies
    log.info(f"loading adjacencies: {ADJ_IN}")
    if not ADJ_IN.exists():
        raise FileNotFoundError(f"Adjacency file not found: {ADJ_IN}")
    adjacencies = pd.read_csv(str(ADJ_IN), sep="\t")
    log.info(f"  links: {len(adjacencies):,} | TFs: {adjacencies['TF'].nunique()} | targets: {adjacencies['target'].nunique()}")

    # load expression matrix -- required for rho dichotomization
    ex_df = load_expression_matrix(LOOM)

    # convert adjacencies to modules
    # ex_mtx is positional and required in pyscenic 0.12.1
    # it is used to compute TF-target expression correlations (rho)
    # which separates activating (rho > 0.03) from repressing modules
    # keep_only_activating=True retains only positive regulators
    log.info("converting adjacencies to modules...")
    log.info("  thresholds=(0.75, 0.9) | top_n_targets=(50,) | min_genes=20")
    log.info("  keep_only_activating=True | rho_threshold=0.03")

    modules = list(
        modules_from_adjacencies(
            adjacencies,
            ex_df,
            thresholds        = (0.75, 0.9),
            top_n_targets     = (50,),
            top_n_regulators  = (5, 10, 50),
            min_genes         = 20,
            rho_dichotomize   = True,
            keep_only_activating = True,
            rho_threshold     = 0.03,
            rho_mask_dropouts = False,
        )
    )

    # free expression matrix -- no longer needed after module creation
    del ex_df
    log.info(f"  modules: {len(modules)}")

    if len(modules) == 0:
        raise ValueError(
            "No modules created. Possible causes:\n"
            "  1. all TF-target correlations below rho_threshold (0.03)\n"
            "  2. all modules below min_genes (20) after filtering\n"
            "  try lowering min_genes to 10 if this persists"
        )

    sizes = sorted([len(m) for m in modules], reverse=True)
    log.info(f"  size range: {min(sizes)} - {max(sizes)} | median: {sizes[len(sizes)//2]}")
    log.info(f"  unique TFs: {len(set(m.transcription_factor for m in modules))}")

    # load databases
    dbs = load_databases()

    # verify motif annotation
    if not MOTIF_ANN.exists():
        raise FileNotFoundError(f"Motif annotation not found: {MOTIF_ANN}")
    log.info(f"motif annotation: {MOTIF_ANN.stat().st_size / 1e6:.1f} MB")

    # run cisTarget
    log.info(f"running cisTarget (NES>={NES_THRESHOLD}, {N_WORKERS} workers)...")
    log.info("  expected runtime: 1-3 hours")

    df = prune2df(
        rnkdbs                  = dbs,
        modules                 = modules,
        motif_annotations_fname = str(MOTIF_ANN),
        nes_threshold           = NES_THRESHOLD,  # <-- explicitly applying your variable
        num_workers             = N_WORKERS,
    )
    log.info(f"cisTarget complete: {len(df)} motif enrichments")

    if len(df) == 0:
        log.error("cisTarget returned 0 enrichments")
        log.error("possible causes:")
        log.error("  1. gene name mismatch between adjacencies and database")
        log.error("  2. NES threshold too stringent (current: %.1f)" % NES_THRESHOLD)
        raise ValueError("cisTarget returned 0 enrichments")

    # convert to regulons
    log.info("converting to regulons...")
    regulons = df2regulons(df)
    log.info(f"  regulons: {len(regulons)}")

    if len(regulons) == 0:
        raise ValueError("df2regulons returned 0 regulons")

    reg_sizes = sorted([len(r.genes) for r in regulons], reverse=True)
    log.info(f"  size range: {min(reg_sizes)} - {max(reg_sizes)} | median: {reg_sizes[len(reg_sizes)//2]}")

    log.info("top 20 regulons by target count:")
    for reg in sorted(regulons, key=lambda r: len(r.genes), reverse=True)[:20]:
        log.info(f"  {reg.name}: {len(reg.genes)} targets")

    # check expected TFs
    reg_tfs  = set(r.transcription_factor for r in regulons)
    expected = ["TBX21", "EOMES", "IRF4", "STAT3", "IRF1", "SPI1", "CEBPB", "IRF8", "RELA", "TCF7"]
    log.info(f"expected TFs found: {[tf for tf in expected if tf in reg_tfs]}")
    log.info(f"expected TFs missing: {[tf for tf in expected if tf not in reg_tfs]}")

    # save
    log.info(f"saving: {REGULONS_OUT}")
    df.to_csv(str(REGULONS_OUT), index=False)

    log.info(f"saving: {REGULONS_PKL_OUT}")
    with open(str(REGULONS_PKL_OUT), "wb") as fh:
        pickle.dump(regulons, fh)

    log.info(f"  csv: {REGULONS_OUT.stat().st_size / 1e6:.1f} MB")
    log.info(f"  pkl: {REGULONS_PKL_OUT.stat().st_size / 1e6:.1f} MB")
    log.info(f"finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("next step: submit 03_aucell.sh")


if __name__ == "__main__":
    main()