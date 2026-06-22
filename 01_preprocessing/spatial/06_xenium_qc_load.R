# sc-triad project
# script: 01_xenium_qc_load.R
# purpose: load and QC Xenium airway spatial data (GSE269354)
#          GSM8313612 = sample 0013717 (asthma_healthy label)
#          GSM8313613 = sample 0013532 (resections_all label)
#
# xenium key differences from visium:
#   - already cell-segmented (no RCTD deconvolution needed)
#   - ~300-500 gene panel (not whole transcriptome)
#   - single-molecule resolution: transcripts.csv.gz available
#   - QC: nCounts, nFeatures, cell area (from cells.csv)
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
})

set.seed(42)
options(future.globals.maxSize = 16 * 1024^3)

# ── paths ──────────────────────────────────────────────────────────────────────
BASE     <- file.path(Sys.getenv("HOME"), "sc-triad")
RAW_DIR  <- file.path(BASE, "01_raw_data", "03_asthma", "gse269354_airway_spatial")
OUT_DIR  <- file.path(BASE, "04_spatial", "asthma_lung")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")
OBJ_DIR  <- file.path(OUT_DIR, "objects")
LOG_FILE <- file.path(OUT_DIR, "01_xenium_qc.log")

for (d in c(FIG_DIR, TAB_DIR, OBJ_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log("SC-TRIAD Xenium asthma lung analysis started")
log(paste("Seurat:", packageVersion("Seurat")))

# ── sample manifest ────────────────────────────────────────────────────────────
# NOTE: "asthma_healthy" label in GEO refers to sample 0013717
# "resections_all" refers to sample 0013532
# Biological grouping MUST be verified against the paper before differential analysis
# Until verified, treat as: sample_A and sample_B
SAMPLES <- list(
  sample_A = list(
    gsm     = "GSM8313612",
    prefix  = "GSM8313612_0013717_asthma_healthy_cell",
    label   = "0013717",
    note    = "asthma_healthy_label"
  ),
  sample_B = list(
    gsm     = "GSM8313613",
    prefix  = "GSM8313613_0013532_resections_all_cell",
    label   = "0013532",
    note    = "resections_all_label"
  )
)

save_tiff <- function(plot, fname, width = 12, height = 10) {
  path <- file.path(FIG_DIR, fname)
  tiff(path, width = width, height = height, units = "in",
       res = 300, compression = "lzw")
  print(plot)
  dev.off()
  log(paste("  saved:", fname))
}

# ── step 1: load Xenium data ───────────────────────────────────────────────────
# GEO Xenium datasets often lack the massive transcripts.csv file, causing 
# Seurat's ReadXenium() to crash. We use a robust manual loader instead.

load_xenium_sample <- function(raw_dir, prefix, label) {
  log(paste("loading:", label))

  file_map <- list(
    barcodes_src = file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz")),
    features_src = file.path(raw_dir, paste0(prefix, "_features.tsv.gz")),
    matrix_src   = file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz")),
    cells_src    = file.path(raw_dir, paste0(prefix, "_cells.csv.gz"))
  )

  # 1. Load Matrix
  mat <- ReadMtx(
    mtx      = file_map$matrix_src,
    cells    = file_map$barcodes_src,
    features = file_map$features_src
  )
  
  # 2. If ReadMtx returns a list (due to Gene Expression & Control probes), bind them
  if (inherits(mat, "list")) {
    mat <- do.call(rbind, mat)
  }

  # 3. Load cell centroids
  cells_df <- read.csv(gzfile(file_map$cells_src))
  log(paste("  cells.csv parsed successfully"))

  coords <- data.frame(
    x = cells_df$x_centroid,
    y = cells_df$y_centroid,
    row.names = cells_df$cell_id
  )

  # 4. Create Seurat Object (Explicitly assign "Xenium" assay for downstream QC)
  obj <- CreateSeuratObject(counts = mat, project = label, assay = "Xenium", min.cells = 1)
  
  # 5. Add spatial FOV (Centroids only, dropping boundary dependency)
  centroids <- CreateCentroids(coords[colnames(obj), , drop = FALSE])
  coords_obj <- CreateFOV(
    coords   = list(centroids = centroids),
    type     = "centroids",
    assay    = "Xenium"
  )
  obj[["fov"]] <- coords_obj
  
  # 6. Add metadata
  obj$sample_id   <- label
  obj$sample_gsm  <- prefix
  obj$tissue_type <- "lung_airway"
  
  # 7. Add cell area natively
  if ("cell_area" %in% colnames(cells_df)) {
    area_vec <- setNames(cells_df$cell_area, cells_df$cell_id)
    common <- intersect(names(area_vec), colnames(obj))
    obj$cell_area <- area_vec[common]
    log(paste("  cell_area added for", length(common), "cells"))
  }

  log(paste("  cells:", ncol(obj), "| genes:", nrow(obj)))
  return(obj)
}

# ── load both samples ──────────────────────────────────────────────────────────
xenium_list <- list()
for (samp_name in names(SAMPLES)) {
  samp <- SAMPLES[[samp_name]]
  xenium_list[[samp_name]] <- tryCatch(
    load_xenium_sample(RAW_DIR, samp$prefix, samp$label),
    error = function(e) {
      log(paste("FAILED:", samp_name, "-", conditionMessage(e)))
      NULL
    }
  )
}
xenium_list <- Filter(Negate(is.null), xenium_list)
log(paste("samples loaded:", length(xenium_list)))

# ── step 2: QC per sample ─────────────────────────────────────────────────────
# Xenium QC metrics differ from scRNA-seq:
#   - nCount_Xenium: transcripts per cell (already filtered by QV threshold)
#   - nFeature_Xenium: genes detected per cell
#   - cell_area: from Xenium segmentation (um^2)
#   - NO mitochondrial % filter (panel is curated, MT genes rarely included)
#   - Negative controls: "NegControlProbe" and "UnassignedCodeword" entries
#     in features — use as QC check, then remove from analysis genes

qc_results <- list()

for (samp_name in names(xenium_list)) {
  obj <- xenium_list[[samp_name]]
  label <- obj$sample_id[1]
  log(paste("QC:", label))

  # identify and flag control probes (not real genes)
  all_features  <- rownames(obj)
  ctrl_features <- grep("NegControl|Unassigned|BLANK", all_features,
                        value = TRUE, ignore.case = TRUE)
  real_features <- setdiff(all_features, ctrl_features)
  log(paste("  real genes:", length(real_features),
            "| control probes:", length(ctrl_features)))

  # control probe count as QC metric
  if (length(ctrl_features) > 0) {
    obj[["pct_neg_control"]] <- colSums(
      GetAssayData(obj, layer = "counts")[ctrl_features, , drop = FALSE]
    ) / (obj$nCount_Xenium + 1) * 100
  } else {
    obj[["pct_neg_control"]] <- 0
  }

  # raw QC summary
  raw_summary <- data.frame(
    sample         = label,
    n_cells_raw    = ncol(obj),
    median_counts  = round(median(obj$nCount_Xenium)),
    median_genes   = round(median(obj$nFeature_Xenium)),
    median_area_um2 = if ("cell_area" %in% colnames(obj@meta.data))
                        round(median(obj$cell_area, na.rm = TRUE)) else NA,
    pct_neg_control_median = round(median(obj$pct_neg_control), 3)
  )
  log(paste("  raw cells:", raw_summary$n_cells_raw,
            "| median counts:", raw_summary$median_counts,
            "| median genes:", raw_summary$median_genes))

  # QC violin pre-filter
  p_vio_pre <- VlnPlot(obj,
    features = c("nCount_Xenium", "nFeature_Xenium"),
    pt.size = 0, ncol = 2) &
    theme_classic(base_size = 9) &
    theme(axis.text.x = element_blank())

  save_tiff(p_vio_pre,
            paste0("01_", label, "_qc_violin_prefilter.tiff"),
            width = 10, height = 5)

  # adaptive MAD filter (same philosophy as PBMC QC)
  md <- obj@meta.data
  lo_count  <- median(md$nCount_Xenium)  - 3 * mad(md$nCount_Xenium)
  hi_count  <- median(md$nCount_Xenium)  + 3 * mad(md$nCount_Xenium)
  lo_feat   <- median(md$nFeature_Xenium) - 3 * mad(md$nFeature_Xenium)
  hi_feat   <- median(md$nFeature_Xenium) + 3 * mad(md$nFeature_Xenium)

  # also filter: cells with 0 real gene counts, high neg control signal
  pass <- md$nCount_Xenium  >= max(lo_count, 1) &
          md$nCount_Xenium  <= hi_count &
          md$nFeature_Xenium >= max(lo_feat, 1) &
          md$pct_neg_control  < 5   # <5% of counts from negative controls

  # area filter if available (remove segmentation artifacts)
  if ("cell_area" %in% colnames(md) && !all(is.na(md$cell_area))) {
    lo_area <- median(md$cell_area, na.rm = TRUE) - 3 * mad(md$cell_area, na.rm = TRUE)
    hi_area <- median(md$cell_area, na.rm = TRUE) + 3 * mad(md$cell_area, na.rm = TRUE)
    pass    <- pass &
               md$cell_area >= max(lo_area, 10) &
               md$cell_area <= hi_area
  }

  n_removed <- sum(!pass)
  log(paste("  MAD filter: removed", n_removed, "cells (",
            round(100 * n_removed / ncol(obj), 1), "%)"))

  obj_filt <- obj[, pass]

  # remove control probes from assay
  if (length(real_features) > 0 && length(ctrl_features) > 0) {
    obj_filt <- obj_filt[real_features, ]
    log(paste("  control probes removed. Final genes:", nrow(obj_filt)))
  }

  post_summary <- data.frame(
    sample          = label,
    n_cells_final   = ncol(obj_filt),
    pct_retained    = round(100 * ncol(obj_filt) / ncol(obj), 1),
    median_counts   = round(median(obj_filt$nCount_Xenium)),
    median_genes    = round(median(obj_filt$nFeature_Xenium)),
    count_lo_thresh = round(max(lo_count, 1)),
    count_hi_thresh = round(hi_count),
    feat_lo_thresh  = round(max(lo_feat, 1))
  )

  qc_results[[samp_name]] <- list(
    raw     = raw_summary,
    post    = post_summary,
    object  = obj_filt
  )

  xenium_list[[samp_name]] <- obj_filt
  log(paste("  final cells:", ncol(obj_filt)))
}

# QC summary table
qc_summary <- do.call(rbind, lapply(qc_results, function(r) {
  cbind(r$raw[, c("sample","n_cells_raw","median_counts","median_genes")],
        r$post[, c("n_cells_final","pct_retained")])
}))
write.csv(qc_summary, file.path(TAB_DIR, "xenium_qc_summary.csv"), row.names = FALSE)
log("QC summary saved")

# ── step 3: merge samples for joint analysis ───────────────────────────────────
# Tag cells with sample prefix before merging to avoid barcode collision
log("merging samples...")
for (samp_name in names(xenium_list)) {
  xenium_list[[samp_name]] <- RenameCells(
    xenium_list[[samp_name]],
    add.cell.id = xenium_list[[samp_name]]$sample_id[1]
  )
}

if (length(xenium_list) == 2) {
  xenium_merged <- merge(xenium_list[[1]], xenium_list[[2]])
} else if (length(xenium_list) == 1) {
  xenium_merged <- xenium_list[[1]]
  log("WARNING: only one sample loaded — no differential analysis possible")
} else {
  stop("No samples loaded. Check file paths.")
}

log(paste("merged object:", ncol(xenium_merged), "cells |",
          nrow(xenium_merged), "genes"))

# save checkpoint
saveRDS(xenium_merged, file.path(OBJ_DIR, "xenium_merged_qc.rds"))
log("QC checkpoint saved: xenium_merged_qc.rds")
log("next: 02_xenium_cluster_annotate.R")