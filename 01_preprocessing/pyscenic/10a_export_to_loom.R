# sc-triad project
# script: 00_export_to_loom.R
# purpose: export triad_integrated seurat object to loom format for pyscenic
#
# approach: bypass SeuratDisk and loomR entirely.
#   both have known HDF5 file-locking issues with Seurat v5 large objects.
#   instead, export three plain files from R, then assemble the loom in python
#   using loompy directly. this is more robust and fully reproducible.
#
# what this script produces (R side):
#   counts_filtered.mtx.gz  — sparse count matrix (genes x cells, market format)
#   barcodes.txt            — cell barcodes (colnames of matrix)
#   features.txt            — gene names (rownames of matrix)
#   cell_metadata.csv       — cell-level annotations for AUCell stratification
#   expressed_genes.txt     — filtered gene list for GRNBoost2
#
# a companion python script (00b_build_loom.py) then assembles these into
# a valid loom file. this two-step approach avoids all HDF5 locking issues.
#
# pyscenic requires raw counts — NOT normalized. SCT and lognorm must NOT be used.
#
# runtime: ~10 min for matrix export (155k cells)
# memory:  ~32 GB
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

set.seed(42)

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
RDS_IN  <- file.path(BASE, "02_scrna", "pbmc", "03_integration",
                     "triad_integrated.rds")
OUT_DIR <- file.path(BASE, "03_pyscenic", "input")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(BASE, "03_pyscenic", "00_export.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log("pyscenic matrix export started")
log(paste("R version:", R.version$version.string))
log(paste("seurat:", packageVersion("Seurat")))
log(paste("Matrix:", packageVersion("Matrix")))

# ── load object ────────────────────────────────────────────────────────────────
log("loading triad_integrated.rds...")
obj <- readRDS(RDS_IN)
log(paste("loaded:", ncol(obj), "cells |", nrow(obj), "genes"))

# verify patch applied
stopifnot(!any(obj$cell_type %in% c("Platelet", "CD16 NK")))
log("patch verified: no Platelet or CD16 NK labels")

log("cell type distribution:")
ct_tab <- sort(table(obj$cell_type), decreasing = TRUE)
for (ct in names(ct_tab)) log(paste0("  ", ct, ": ", ct_tab[ct]))

# ── extract raw RNA counts ─────────────────────────────────────────────────────
# seurat v5: layers may be split per sample after integration — join first
log("joining RNA layers...")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

log("extracting raw count matrix...")
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
log(paste("raw matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells"))
log(paste("matrix class:", class(counts_mat)[1]))
log(paste("total counts:", format(sum(counts_mat), big.mark = ",")))
log(paste("% non-zero entries:",
          round(100 * nnzero(counts_mat) / prod(dim(counts_mat)), 2)))

# ── filter expressed genes ─────────────────────────────────────────────────────
# >= 1% detection rate: balances signal quality vs GRNBoost2 runtime
# removing mito/ribo: these reflect cell stress/size, not TF-driven regulation
# reference: van de sande et al., nature protocols 2020

log("filtering expressed genes...")
pct_expressed   <- rowMeans(counts_mat > 0)
expressed_genes <- names(pct_expressed[pct_expressed >= 0.01])
n_before        <- length(expressed_genes)

expressed_genes <- expressed_genes[
  !grepl("^MT-|^RPL|^RPS|^MRPL|^MRPS", expressed_genes)
]

log(paste("genes before filter:", nrow(counts_mat)))
log(paste("genes >= 1% detection:", n_before))
log(paste("genes after mito/ribo removal:", length(expressed_genes)))

# subset to filtered genes
counts_filt <- counts_mat[expressed_genes, ]
log(paste("final matrix:", nrow(counts_filt), "genes x", ncol(counts_filt), "cells"))

# ── write sparse matrix in 10x market format ───────────────────────────────────
# writeMM writes the standard Matrix Market format
# barcodes and features written as plain text
# python loompy reads these natively via scipy.io.mmread

log("writing sparse matrix to disk...")

# delete any existing corrupted loom before proceeding
loom_path <- file.path(OUT_DIR, "triad_pbmc_raw.loom")
if (file.exists(loom_path)) {
  file.remove(loom_path)
  log("removed existing (potentially corrupted) loom file")
}

mtx_path      <- file.path(OUT_DIR, "counts_filtered.mtx")
barcodes_path <- file.path(OUT_DIR, "barcodes.txt")
features_path <- file.path(OUT_DIR, "features.txt")
gene_list_path <- file.path(OUT_DIR, "expressed_genes.txt")

# write matrix — this is the slow step (~5-10 min for 155k cells)
log("  writing counts_filtered.mtx (slow step, ~5-10 min)...")
writeMM(counts_filt, file = mtx_path)
log(paste("  matrix written:", round(file.size(mtx_path) / 1e9, 2), "GB"))

# compress to save disk space
log("  compressing matrix...")
system2("gzip", args = c("-f", mtx_path))
log(paste("  compressed:", round(file.size(paste0(mtx_path, ".gz")) / 1e6, 1), "MB"))

# write barcodes (cell IDs)
writeLines(colnames(counts_filt), barcodes_path)
log(paste("  barcodes written:", length(colnames(counts_filt))))

# write features (gene names)
writeLines(rownames(counts_filt), features_path)
log(paste("  features written:", length(rownames(counts_filt))))

# write expressed gene list (same as features, for GRNBoost2 reference)
writeLines(expressed_genes, gene_list_path)
log(paste("  gene list written:", gene_list_path))

# ── write cell metadata ────────────────────────────────────────────────────────
log("writing cell metadata...")

# carefully select only columns that definitely exist
# rownames of metadata = cell barcodes — these must match barcodes.txt exactly
available_cols <- colnames(obj@meta.data)
log(paste("  available metadata columns:", paste(available_cols, collapse = ", ")))

wanted_cols <- c("sample_label", "disease", "group", "cell_type",
                 "nCount_RNA", "nFeature_RNA", "percent_mt")
present_cols <- wanted_cols[wanted_cols %in% available_cols]
missing_cols <- wanted_cols[!wanted_cols %in% available_cols]

if (length(missing_cols) > 0) {
  log(paste("  WARNING: columns not found:", paste(missing_cols, collapse = ", ")))
}

meta_df <- obj@meta.data[colnames(counts_filt), present_cols, drop = FALSE]
meta_df$cell_id <- rownames(meta_df)
meta_df <- meta_df[, c("cell_id", present_cols)]

meta_path <- file.path(OUT_DIR, "cell_metadata.csv")
write.csv(meta_df, meta_path, row.names = FALSE)
log(paste("  metadata written:", nrow(meta_df), "cells x", ncol(meta_df), "columns"))

# verify barcode order matches between matrix and metadata
stopifnot(all(meta_df$cell_id == colnames(counts_filt)))
log("  barcode order verified: metadata matches matrix columns exactly")

# ── write python loom assembly script ─────────────────────────────────────────
# this python script is generated automatically from R so paths are correct
# it reads the mtx/txt files and builds a valid loom using loompy
# loompy is the reference implementation — no HDF5 locking issues

python_script <- sprintf('#!/usr/bin/env python3
# sc-triad project
# script: 00b_build_loom.py (auto-generated by 00_export_to_loom.R)
# purpose: assemble pyscenic loom from R-exported matrix files
#
# this script is called automatically by 00_export_loom.sh
# do not edit paths manually — they are injected by R
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe

import os
import gzip
import logging
import numpy as np
import scipy.io
import loompy

logging.basicConfig(
    level   = logging.INFO,
    format  = "[%%(asctime)s] %%(levelname)s %%(message)s",
    datefmt = "%%H:%%M:%%S"
)
log = logging.getLogger("sc-triad-loom")

# paths (injected by R)
MTX_GZ    = "%s"
BARCODES  = "%s"
FEATURES  = "%s"
METADATA  = "%s"
LOOM_OUT  = "%s"

def main():
    log.info("=== SC-TRIAD loom assembly ===")

    # delete corrupted loom if it exists
    if os.path.exists(LOOM_OUT):
        os.remove(LOOM_OUT)
        log.info(f"removed existing loom: {LOOM_OUT}")

    # load sparse matrix
    log.info(f"loading count matrix: {MTX_GZ}")
    with gzip.open(MTX_GZ, "rb") as f:
        matrix = scipy.io.mmread(f).tocsc()
    log.info(f"  matrix shape (genes x cells): {matrix.shape}")
    log.info(f"  matrix type: {type(matrix)}")
    log.info(f"  non-zero entries: {matrix.nnz:,}")

    # load barcodes and features
    with open(BARCODES) as f:
        barcodes = [line.strip() for line in f]
    with open(FEATURES) as f:
        gene_names = [line.strip() for line in f]

    log.info(f"  barcodes: {len(barcodes)}")
    log.info(f"  genes: {len(gene_names)}")

    # verify dimensions match
    assert matrix.shape[0] == len(gene_names), (
        f"gene count mismatch: matrix rows={matrix.shape[0]}, "
        f"features={len(gene_names)}"
    )
    assert matrix.shape[1] == len(barcodes), (
        f"cell count mismatch: matrix cols={matrix.shape[1]}, "
        f"barcodes={len(barcodes)}"
    )

    # load metadata
    import pandas as pd
    log.info(f"loading metadata: {METADATA}")
    meta = pd.read_csv(METADATA)
    log.info(f"  metadata shape: {meta.shape}")
    log.info(f"  metadata columns: {list(meta.columns)}")

    # verify barcode order matches metadata
    assert list(meta["cell_id"]) == barcodes, (
        "Barcode order mismatch between metadata and matrix. "
        "This should not happen — check 00_export_to_loom.R output."
    )

    # build column attributes (cell-level)
    # CellID is mandatory for pyscenic — it reads this as the cell identifier
    col_attrs = {"CellID": np.array(barcodes)}
    for col in ["cell_type", "disease", "group", "sample_label"]:
        if col in meta.columns:
            col_attrs[col] = np.array(meta[col].astype(str))
            log.info(f"  added col_attr: {col}")

    # build row attributes (gene-level)
    # Gene is mandatory for pyscenic — it reads this as the gene identifier
    row_attrs = {"Gene": np.array(gene_names)}

    # write loom
    # loompy.create expects: matrix as numpy array (genes x cells)
    # convert sparse to dense in chunks to avoid memory spike
    log.info(f"writing loom: {LOOM_OUT}")
    log.info("  this may take 10-30 minutes for 155k cells...")

    # loompy.create_append allows chunked writing — safer for large matrices
    # chunk size of 10000 cells keeps peak memory per chunk manageable
    CHUNK = 10000
    n_cells = matrix.shape[1]
    n_chunks = (n_cells + CHUNK - 1) // CHUNK

    log.info(f"  writing in {n_chunks} chunks of {CHUNK} cells...")

    for i in range(n_chunks):
        start = i * CHUNK
        end   = min(start + CHUNK, n_cells)

        chunk_mat  = matrix[:, start:end].toarray().astype(np.float32)

        # col_attrs for this chunk
        chunk_col = {k: v[start:end] for k, v in col_attrs.items()}

        if i == 0:
            # create loom on first chunk
            loompy.create(
                filename  = LOOM_OUT,
                layers    = chunk_mat,
                row_attrs = row_attrs,
                col_attrs = chunk_col
            )
            log.info(f"  created loom with first chunk ({end} cells)")
        else:
            # append subsequent chunks
            with loompy.connect(LOOM_OUT) as ds:
                ds.add_columns(chunk_mat, col_attrs=chunk_col)

        if (i + 1) %% 5 == 0 or (i + 1) == n_chunks:
            log.info(f"  progress: {end}/{n_cells} cells ({100*end/n_cells:.1f}%%)")

    # verify written loom
    with loompy.connect(LOOM_OUT, mode="r") as ds:
        log.info(f"loom written successfully:")
        log.info(f"  shape (genes x cells): {ds.shape}")
        log.info(f"  row attrs: {list(ds.ra.keys())}")
        log.info(f"  col attrs: {list(ds.ca.keys())}")

    size_gb = os.path.getsize(LOOM_OUT) / 1e9
    log.info(f"  file size: {size_gb:.2f} GB")
    log.info("loom assembly complete")
    log.info(f"next step: submit 01_grn.sh")


if __name__ == "__main__":
    main()
',
  paste0(mtx_path, ".gz"),
  barcodes_path,
  features_path,
  meta_path,
  loom_path
)

python_script_path <- file.path(OUT_DIR, "00b_build_loom.py")
writeLines(python_script, python_script_path)
log(paste("python loom assembly script written:", python_script_path))

# ── summary ────────────────────────────────────────────────────────────────────
log("")
log("=== R export summary ===")
log(paste("matrix (gz):     ", paste0(mtx_path, ".gz")))
log(paste("barcodes:        ", barcodes_path))
log(paste("features:        ", features_path))
log(paste("metadata:        ", meta_path))
log(paste("python script:   ", python_script_path))
log(paste("expressed genes: ", length(expressed_genes)))
log(paste("cells:           ", ncol(counts_filt)))
log("")
log("R export complete.")
log("the SLURM script (00_export_loom.sh) will now call python 00b_build_loom.py")
log("to assemble the final loom file.")
