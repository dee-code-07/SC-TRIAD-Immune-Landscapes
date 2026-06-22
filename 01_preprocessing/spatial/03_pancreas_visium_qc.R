# sc-triad project
# script: 03_pancreas_visium.R — v5
#
# Bugs fixed (cumulative):
#
#   BUG 1 (duplicate row.names crash):
#     Baron CSVs have duplicate gene symbols. read.csv(row.names=1) refused
#     duplicates. Fixed: row.names=NULL + make.unique() on gene column.
#
#   BUG 2 (silent download failure -> empty matrix):
#     Compute nodes have no internet. download.file() wrote ~1 KB HTML error
#     pages. Fixed: MIN_BARON_BYTES = 1e6 pre-flight check before reading.
#
#   BUG 3 (root cause of "0 common genes"):
#     Baron CSV is CELLS x GENES. Gene names are COLUMN NAMES, not row names.
#     Fixed: extract colnames[-1] as gene symbols, transpose to genes x cells.
#
#   BUG 4 (NAs introduced by coercion -> "No cells found"):
#     Mixed-type data frame passed to as.matrix() upcasts to character.
#     Fixed: convert each column independently via lapply+as.numeric, then
#     drop all-NA columns and zero-fill remaining partial NAs.
#
#   BUG 5 (root cause of "No cells remain after removing rare types"):
#     dplyr::recode(.default = <vector>) silently takes only the LAST element
#     of the vector as a scalar default. Every unmatched cell type (alpha,
#     beta, acinar, ductal, ...) collapsed to one garbage string, making ALL
#     real types appear as "rare" -> all cells removed.
#     Fixed: replaced with dplyr::case_when, which correctly handles
#     vectorized fall-through via TRUE ~ <expression>.
#
# author : deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide  : dr budheswar dehury
# version: v5

# =============================================================================
# 0. Setup
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(scales)
  library(ggrepel)
  library(stringr)
})

set.seed(42)
options(future.globals.maxSize = 100 * 1024^3)

# ---------- paths -------------------------------------------------------------
BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
PAN_DIR <- file.path(BASE, "01_raw_data/01_t2d/gse264331_pancreas_visium")
REF_DIR <- file.path(BASE, "01_raw_data/01_t2d/baron2016_pancreas_ref")
OUT_DIR <- file.path(BASE, "04_spatial/pancreas")

for (d in c("logs", "objects", "figures", "tables"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "logs/03_pancreas_visium.log")
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)   # fresh log each run

# ---------- helpers -----------------------------------------------------------
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

save_pdf <- function(plot, path, width = 10, height = 8) {
  pdf(path, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  log_msg("  saved: ", basename(path))
}

pub_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 9, colour = "grey40",
                                    lineheight = 1.3),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 10, colour = "black"),
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ---------- cell-type groupings -----------------------------------------------
islet_types    <- c("beta", "alpha", "delta", "gamma", "epsilon")
exocrine_types <- c("acinar", "ductal")
stromal_types  <- c("stellate", "activated_stellate", "endothelial")
immune_types_p <- c("macrophage", "T_cell")
status_colors  <- c("Control" = "#377EB8", "T2D" = "#E41A1C")

log_msg("=== SC-TRIAD pancreas Visium spatial analysis (v5) ===")
log_msg("Seurat  : ", as.character(packageVersion("Seurat")))
log_msg("spacexr : ", as.character(packageVersion("spacexr")))


# =============================================================================
# STEP 1: Load Baron et al. 2016 reference (GSE84133)
# =============================================================================
log_msg("=== STEP 1: Baron et al. 2016 reference ===")

MIN_BARON_BYTES <- 1e6    # real CSVs are 5-15 MB; HTML error pages ~1 KB
baron_files     <- paste0("GSE84133_human", 1:4, ".csv.gz")

# --- pre-flight: validate every file before reading --------------------------
log_msg("pre-flight: validating Baron reference files ...")
any_bad <- FALSE
for (fname in baron_files) {
  fpath <- file.path(REF_DIR, fname)
  if (!file.exists(fpath)) {
    log_msg("  MISSING : ", fname)
    any_bad <- TRUE
  } else {
    fsize <- file.size(fpath)
    if (fsize < MIN_BARON_BYTES) {
      log_msg("  BAD FILE (", fsize,
              " bytes - likely HTML redirect): ", fname)
      file.remove(fpath)
      any_bad <- TRUE
    } else {
      log_msg("  OK  (", round(fsize / 1e6, 1), " MB): ", fname)
    }
  }
}
if (any_bad) {
  stop(
    "\nOne or more Baron reference files are missing or corrupted.\n",
    "Compute nodes have no outbound internet.\n",
    "Run from the LOGIN NODE:\n\n",
    "  bash ~/sc-triad/scripts/04_spatial/pancreas/00_download_baron.sh\n\n",
    "Then re-submit this job."
  )
}
log_msg("all Baron reference files validated - reading ...")

# --- checkpoint: skip rebuild if RDS already exists --------------------------
baron_rds <- file.path(REF_DIR, "baron2016_reference.rds")

if (file.exists(baron_rds)) {

  log_msg("checkpoint found, loading ...")
  baron_ref <- readRDS(baron_rds)
  log_msg("loaded: ", ncol(baron_ref), " cells | ",
          length(unique(baron_ref$cell_type)), " types")

} else {

  mat_list <- list()
  ct_list  <- list()

  for (i in seq_along(baron_files)) {

    fname <- baron_files[i]
    dest  <- file.path(REF_DIR, fname)
    log_msg("  reading: ", fname)

    # BUG 1 FIX: row.names=NULL avoids crash on duplicate gene symbols
    df_raw <- read.csv(gzfile(dest), row.names = NULL, check.names = FALSE)

    # BUG 3 FIX: Baron CSV is CELLS x GENES
    #   col 1   = cell-type label
    #   cols 2+ = gene expression, colnames = gene symbols
    cell_types_i <- as.character(df_raw[["assigned_cluster"]])
    expr_df <- df_raw[, !colnames(df_raw) %in% c("", "barcode", "assigned_cluster"),
                       drop = FALSE]

    # BUG 4 FIX: convert each column independently to numeric.
    # as.matrix() on a mixed-type data frame upcasts the ENTIRE matrix to
    # character. Converting column-by-column isolates non-numeric columns.
    expr_mat <- do.call(cbind, lapply(expr_df, function(col) {
      as.numeric(as.character(col))
    }))
    rownames(expr_mat) <- NULL
    colnames(expr_mat) <- colnames(expr_df)

    # Drop columns that are entirely NA (stray index / label columns)
    all_na_cols <- colSums(is.na(expr_mat)) == nrow(expr_mat)
    if (any(all_na_cols)) {
      log_msg("    dropping ", sum(all_na_cols),
              " non-numeric column(s): ",
              paste(colnames(expr_mat)[all_na_cols], collapse = ", "))
      expr_mat <- expr_mat[, !all_na_cols, drop = FALSE]
    }

    # Zero-fill any remaining partial NAs (blank cells in the CSV)
    n_partial <- sum(is.na(expr_mat))
    if (n_partial > 0) {
      log_msg("    zeroing ", n_partial, " partial NA value(s)")
      expr_mat[is.na(expr_mat)] <- 0
    }

    # BUG 1 FIX cont.: make gene names unique
    gene_names <- make.unique(as.character(colnames(expr_mat)))

    # Transpose: rows=cells, cols=genes  ->  rows=genes, cols=cells
    mat_t           <- t(expr_mat)
    rownames(mat_t) <- gene_names
    colnames(mat_t) <- paste0("donor", i, "_cell", seq_len(ncol(mat_t)))

    mat_list[[i]] <- as(mat_t, "sparseMatrix")
    ct_list[[i]]  <- cell_types_i

    log_msg("    ", ncol(mat_t), " cells | ", nrow(mat_t),
            " genes | types: ",
            paste(sort(unique(cell_types_i)), collapse = ", "))
  }

  # Intersect ROWNAMES = gene symbols (correct after BUG 3 fix)
  common_genes <- Reduce(intersect, lapply(mat_list, rownames))
  log_msg("common genes across 4 donors: ", length(common_genes))

  if (length(common_genes) < 5000)
    stop("Only ", length(common_genes), " common genes. Expected ~14,000+.\n",
         "Delete ", REF_DIR, " and re-run 00_download_baron.sh.")

  # Combine donors
  combined <- do.call(cbind,
                      lapply(mat_list, function(m) m[common_genes, ]))
  all_ct   <- unlist(ct_list)
  log_msg("combined: ", nrow(combined), " genes x ", ncol(combined), " cells")

  # ---- BUG 5 FIX: cell-type harmonization -----------------------------------
  # dplyr::recode(.default = <vector>) silently uses only the LAST element of
  # the vector as a scalar. Every unmatched type (alpha, beta, acinar, ...)
  # was collapsed to one garbage string -> all real types flagged as "rare".
  # dplyr::case_when correctly handles vectorized fall-through via TRUE ~ expr.
  all_ct_lower <- tolower(trimws(all_ct))

  ct_harmonized <- dplyr::case_when(
    all_ct_lower == "activated stellate"        ~ "activated_stellate",
    all_ct_lower == "psc"                       ~ "stellate",
    all_ct_lower %in% c("t_cell", "t cell")    ~ "T_cell",
    all_ct_lower %in% c("schwann", "mast")      ~ "other",
    TRUE                                        ~ all_ct_lower
  )

  log_msg("cell type distribution after harmonization:")
  print(sort(table(ct_harmonized), decreasing = TRUE))

  # Build Seurat object
  baron_ref           <- CreateSeuratObject(counts = combined, min.cells = 3)
  baron_ref$cell_type <- ct_harmonized
  baron_ref$donor     <- rep(
    paste0("donor", 1:4),
    times = sapply(mat_list, ncol)
  )[seq_len(ncol(baron_ref))]

  # Remove unknowns and "other"
  keep_mask <- !baron_ref$cell_type %in% c("unknown", "other")
  if (sum(keep_mask) == 0)
    stop("No cells remain after removing 'unknown'/'other' types.\n",
         "Check harmonization output printed above.")
  baron_ref <- baron_ref[, keep_mask]

  log_msg("after removing unknown/other: ",
          ncol(baron_ref), " cells | ",
          length(unique(baron_ref$cell_type)), " types")

  saveRDS(baron_ref, baron_rds)
  log_msg("Baron reference saved: ", baron_rds)
}


# =============================================================================
# STEP 2: Prepare spacexr Reference
# =============================================================================
log_msg("=== STEP 2: preparing RCTD reference ===")

ct_tab <- sort(table(baron_ref$cell_type), decreasing = TRUE)
log_msg("cell type counts:")
print(ct_tab)

# Remove types with < 10 cells
rare <- names(ct_tab[ct_tab < 10])
if (length(rare) > 0) {
  log_msg("removing rare types (<10 cells): ", paste(rare, collapse = ", "))
  keep_not_rare <- !baron_ref$cell_type %in% rare
  if (sum(keep_not_rare) == 0)
    stop("No cells remain after removing rare types.\n",
         "Cell type table printed above - check harmonization.")
  baron_ref <- baron_ref[, keep_not_rare]
}

final_types <- sort(unique(baron_ref$cell_type))
log_msg("final types (", length(final_types), "): ",
        paste(final_types, collapse = ", "))

# Downsample to max 300 cells per type
MAX_REF <- 300L
set.seed(42)
keep_idx <- unlist(lapply(final_types, function(ct) {
  idx <- which(baron_ref$cell_type == ct)
  sample(idx, min(length(idx), MAX_REF))
}))
ref_ds <- baron_ref[, keep_idx]
log_msg("downsampled reference: ", ncol(ref_ds), " cells")
print(sort(table(ref_ds$cell_type), decreasing = TRUE))

# Normalize if data layer missing
if (!"data" %in% Layers(ref_ds, assay = "RNA")) {
  log_msg("normalizing reference ...")
  ref_ds <- NormalizeData(ref_ds, verbose = FALSE)
}

ref_counts <- GetAssayData(ref_ds, assay = "RNA", layer = "counts")
ref_counts  <- round(ref_counts)
ref_numi    <- as.integer(round(Matrix::colSums(ref_counts)))
ref_labels  <- factor(ref_ds$cell_type)
names(ref_numi)   <- colnames(ref_ds)
names(ref_labels) <- colnames(ref_ds)

reference <- Reference(counts     = ref_counts,
                       cell_types = ref_labels,
                       nUMI       = ref_numi)
log_msg("spacexr Reference: ", length(levels(reference@cell_types)),
        " cell types")


# =============================================================================
# STEP 3: Load Visium h5 files
# =============================================================================
log_msg("=== STEP 3: loading Visium h5 files ===")

h5_files <- list.files(PAN_DIR, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) == 0)
  stop("No .h5 files found in: ", PAN_DIR)
log_msg("h5 files found: ", length(h5_files))

# GSE264331 donor-to-disease status
# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264331
DONOR_STATUS <- c(
  "HP21024" = "T2D",
  "HP21035" = "T2D",
  "HP21091" = "T2D",
  "HP21161" = "T2D",
  "HP21177" = "T2D",
  "HP21181" = "Control"
)

parse_meta <- function(fp) {
  fname   <- basename(fp)
  donor   <- stringr::str_extract(fname, "HP\\d+")
  section <- stringr::str_extract(fname, "Vis\\d+_S\\d+_\\w+")
  if (is.na(section))
    section <- stringr::str_extract(fname, "Vis\\d+_S\\d+")
  status <- ifelse(
    !is.na(donor) && donor %in% names(DONOR_STATUS),
    DONOR_STATUS[[donor]], "Unknown"
  )
  data.frame(file = fname, donor = donor, section = section,
             status = status, stringsAsFactors = FALSE)
}

h5_meta <- do.call(rbind, lapply(h5_files, parse_meta))
log_msg("sample metadata:")
print(h5_meta[, c("donor", "section", "status")])
log_msg("T2D sections    : ", sum(h5_meta$status == "T2D"))
log_msg("Control sections: ", sum(h5_meta$status == "Control"))

# Load count matrices
mat_list_st <- list()
sample_ids  <- character(0)

for (i in seq_along(h5_files)) {
  sid <- paste0(h5_meta$donor[i], "_", h5_meta$section[i])
  tryCatch({
    m <- Read10X_h5(h5_files[i])
    if (is.list(m)) m <- m[["Gene Expression"]]
    colnames(m) <- paste0(sid, "_", colnames(m))
    mat_list_st[[sid]] <- m
    sample_ids  <- c(sample_ids, sid)
    log_msg("  ", sid, " : ", nrow(m), " genes x ", ncol(m), " spots")
  }, error = function(e) {
    log_msg("  ERROR: ", basename(h5_files[i]), " -> ", e$message)
  })
}

if (length(mat_list_st) == 0)
  stop("No Visium sections loaded. Check h5 files in: ", PAN_DIR)

# Gene intersections
common_st  <- Reduce(intersect, lapply(mat_list_st, rownames))
common_all <- intersect(common_st, rownames(ref_counts))
log_msg("genes common across ST sections   : ", length(common_st))
log_msg("genes common ST + Baron reference : ", length(common_all))

if (length(common_all) < 2000)
  warning("Only ", length(common_all),
          " common genes - check genome build compatibility.")

ref_counts_sub <- ref_counts[common_all, ]
reference_sub  <- Reference(counts     = ref_counts_sub,
                             cell_types = ref_labels,
                             nUMI       = ref_numi)


# =============================================================================
# STEP 4: RCTD deconvolution (full mode, per section)
# =============================================================================
log_msg("=== STEP 4: RCTD (full mode, per section) ===")
log_msg("Spatial coordinates unavailable from GEO. Placeholder grid used.")
log_msg("RCTD MLE is per-spot; coordinates do not affect deconvolution.")

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
log_msg("cores: ", N_CORES)

all_prop_dfs <- list()

for (sid in sample_ids) {

  log_msg("\n--- section: ", sid, " ---")
  rctd_rds <- file.path(OUT_DIR, "objects", paste0("rctd_", sid, ".rds"))
  mat_st   <- mat_list_st[[sid]][common_all, , drop = FALSE]

  if (file.exists(rctd_rds)) {
    log_msg("  checkpoint found, loading ...")
    myRCTD <- readRDS(rctd_rds)
  } else {

    n_spots    <- ncol(mat_st)
    spot_names <- colnames(mat_st)

    # Placeholder grid (RCTD does not use coordinates for MLE)
    side   <- ceiling(sqrt(n_spots))
    grid_x <- ((seq_len(n_spots) - 1L) %% side) + 1L
    grid_y <- ((seq_len(n_spots) - 1L) %/% side) + 1L
    coords <- data.frame(x = grid_x, y = grid_y, row.names = spot_names)

    spot_numi        <- as.integer(round(Matrix::colSums(mat_st)))
    names(spot_numi) <- spot_names

    keep <- spot_numi >= 100
    if (sum(keep) == 0)
      stop("No spots with UMI >= 100 in section: ", sid)

    mat_st_f    <- mat_st[, keep, drop = FALSE]
    spot_numi_f <- spot_numi[keep]
    coords_f    <- coords[keep, , drop = FALSE]
    log_msg("  spots after UMI>=100 filter: ", sum(keep), " / ", n_spots)

    puck <- SpatialRNA(coords = coords_f,
                       counts = mat_st_f,
                       nUMI   = spot_numi_f)

    myRCTD <- create.RCTD(
      spatialRNA        = puck,
      reference         = reference_sub,
      max_cores         = N_CORES,
      CELL_MIN_INSTANCE = 10,
      gene_cutoff       = 0.000125,
      fc_cutoff         = 0.5,
      UMI_min           = 100,
      UMI_min_sigma     = 300
    )
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
    saveRDS(myRCTD, rctd_rds)
    log_msg("  saved: ", basename(rctd_rds))
  }

  weights  <- myRCTD@results$weights
  prop_mat <- normalize_weights(weights)

  idx     <- which(sample_ids == sid)
  prop_df <- as.data.frame(prop_mat) %>%
    mutate(spot    = rownames(prop_mat),
           section = sid,
           donor   = h5_meta$donor[idx],
           status  = h5_meta$status[idx])

  all_prop_dfs[[sid]] <- prop_df
  log_msg("  proportions extracted: ", nrow(prop_mat), " spots")
}

all_prop <- do.call(rbind, all_prop_dfs)
rownames(all_prop) <- NULL
ct_cols  <- setdiff(colnames(all_prop),
                    c("spot", "section", "donor", "status"))

write.csv(all_prop,
          file.path(OUT_DIR, "tables/rctd_proportions_per_spot.csv"),
          row.names = FALSE)
log_msg("total spots: ", nrow(all_prop),
        " | cell types detected: ", length(ct_cols))


# =============================================================================
# STEP 5: Summaries + differential proportions
# =============================================================================
log_msg("=== STEP 5: summaries and differential proportions ===")

summary_status <- all_prop %>%
  group_by(status) %>%
  summarise(across(all_of(ct_cols), \(x) mean(x, na.rm = TRUE)),
            n_spots    = n(),
            n_sections = n_distinct(section),
            .groups    = "drop")
write.csv(summary_status,
          file.path(OUT_DIR, "tables/rctd_summary_by_status.csv"),
          row.names = FALSE)

summary_section <- all_prop %>%
  group_by(section, donor, status) %>%
  summarise(across(all_of(ct_cols), \(x) mean(x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
write.csv(summary_section,
          file.path(OUT_DIR, "tables/rctd_summary_by_section.csv"),
          row.names = FALSE)

t2d_secs  <- summary_section %>% filter(status == "T2D")
ctrl_secs <- summary_section %>% filter(status == "Control")
n_t2d     <- nrow(t2d_secs)
n_ctrl    <- nrow(ctrl_secs)
log_msg("T2D sections: ", n_t2d, " | Control sections: ", n_ctrl)

if (n_ctrl < 2) {
  log_msg("WARNING: n=1 control donor - Wilcoxon p-values unreliable.")
  log_msg("Reporting log2FC + rank-biserial r as primary effect-size metrics.")
}

diff_rows <- lapply(ct_cols, function(ct) {
  xt <- t2d_secs[[ct]]
  xc <- ctrl_secs[[ct]]
  mt <- mean(xt, na.rm = TRUE)
  mc <- mean(xc, na.rm = TRUE)
  fc <- log2((mt + 1e-6) / (mc + 1e-6))

  pv  <- NA_real_
  rbe <- NA_real_
  if (n_t2d >= 2 && n_ctrl >= 2) {
    wt  <- suppressWarnings(wilcox.test(xt, xc, exact = FALSE))
    pv  <- wt$p.value
    rbe <- round(1 - 2 * wt$statistic / (n_t2d * n_ctrl), 4)
  }

  data.frame(
    cell_type        = ct,
    mean_t2d         = round(mt, 5),
    mean_ctrl        = round(mc, 5),
    log2FC           = round(fc, 4),
    rank_biserial_r  = rbe,
    p_value          = pv,
    n_t2d            = n_t2d,
    n_ctrl           = n_ctrl,
    stringsAsFactors = FALSE
  )
})

diff_df <- do.call(rbind, diff_rows) %>%
  arrange(desc(abs(log2FC)))

diff_df$p_adj <- if (!all(is.na(diff_df$p_value)))
  p.adjust(diff_df$p_value, method = "BH") else NA_real_

write.csv(diff_df,
          file.path(OUT_DIR, "tables/rctd_differential_celltypes.csv"),
          row.names = FALSE)

log_msg("top cell types by |log2FC|:")
print(head(diff_df[, c("cell_type", "mean_t2d", "mean_ctrl",
                        "log2FC", "rank_biserial_r", "p_adj")], 12))


# =============================================================================
# STEP 6: Figures
# =============================================================================
log_msg("=== STEP 6: figures ===")

# ---- FIG 1: Islet cell-type composition -------------------------------------
islet_present <- intersect(islet_types, ct_cols)
if (length(islet_present) >= 2) {

  islet_df <- all_prop %>%
    select(section, donor, status, all_of(islet_present)) %>%
    group_by(section, donor, status) %>%
    summarise(across(all_of(islet_present),
                     \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(all_of(islet_present),
                 names_to = "cell_type", values_to = "mean_prop") %>%
    mutate(status = factor(status, levels = c("Control", "T2D")))

  islet_sum <- islet_df %>%
    group_by(status, cell_type) %>%
    summarise(mean = mean(mean_prop, na.rm = TRUE),
              se   = sd(mean_prop,   na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  p1 <- ggplot(islet_sum,
               aes(x = cell_type, y = mean,
                   fill = status, colour = status)) +
    geom_bar(stat = "identity",
             position = position_dodge(0.75),
             width = 0.65, alpha = 0.88) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.75),
                  width = 0.25, linewidth = 0.6, colour = "grey30") +
    geom_point(data = islet_df,
               aes(y = mean_prop, group = status),
               position = position_dodge(0.75),
               size = 2.2, alpha = 0.75,
               shape = 21, colour = "grey20") +
    scale_fill_manual(values   = status_colors) +
    scale_colour_manual(values = status_colors) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title    = "Islet cell type proportions: T2D vs Control",
      subtitle = paste0(
        "RCTD full mode | mean \u00b1 SE per section\n",
        "T2D: ", n_t2d, " sections | Control: ", n_ctrl, " section(s)"
      ),
      x = "Islet cell type", y = "Mean proportion",
      fill = "Status", colour = "Status"
    ) +
    pub_theme

  save_pdf(p1,
           file.path(OUT_DIR, "figures/01_islet_composition.pdf"),
           width = 9, height = 7)
}

# ---- FIG 2: All cell types stacked bar --------------------------------------
mean_all <- all_prop %>%
  pivot_longer(all_of(ct_cols),
               names_to = "cell_type", values_to = "proportion") %>%
  group_by(status, cell_type) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(status = factor(status, levels = c("Control", "T2D")))

p2 <- ggplot(mean_all,
             aes(x = status, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.25) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title    = "Pancreas Visium: cell type composition (T2D vs Control)",
    subtitle = "RCTD full mode | GSE264331 | Reference: Baron et al. 2016",
    x = NULL, y = "Mean proportion", fill = "Cell type"
  ) +
  pub_theme +
  theme(legend.key.size = unit(0.4, "cm"))

save_pdf(p2,
         file.path(OUT_DIR, "figures/02_celltype_composition.pdf"),
         width = 8, height = 8)

# ---- FIG 3: Effect-size plot ------------------------------------------------
compartment_map <- c(
  setNames(rep("Islet",    length(islet_types)),    islet_types),
  setNames(rep("Exocrine", length(exocrine_types)), exocrine_types),
  setNames(rep("Stromal",  length(stromal_types)),  stromal_types),
  setNames(rep("Immune",   length(immune_types_p)), immune_types_p)
)

diff_plot <- diff_df %>%
  mutate(
    compartment = ifelse(
      cell_type %in% names(compartment_map),
      compartment_map[cell_type], "Other"
    ),
    label = ifelse(
      abs(log2FC) > 0.5 |
        (!is.na(rank_biserial_r) & abs(rank_biserial_r) > 0.3),
      cell_type, NA_character_
    )
  )

p3 <- ggplot(diff_plot,
             aes(x = log2FC, y = rank_biserial_r,
                 colour = compartment, shape = compartment,
                 label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dotted", colour = "grey80") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3.2, na.rm = TRUE, max.overlaps = 15,
                  fontface = "italic") +
  scale_colour_manual(
    values = c(Islet    = "#E41A1C", Exocrine = "#FF7F00",
               Stromal  = "#999999", Immune   = "#377EB8",
               Other    = "#BDBDBD")
  ) +
  scale_shape_manual(
    values = c(Islet = 16, Exocrine = 17,
               Stromal = 18, Immune = 15, Other = 4)
  ) +
  labs(
    title    = "Differential cell type proportions: T2D vs Control",
    subtitle = paste0(
      "Effect size analysis | x = log2(T2D/Control)\n",
      "y = rank-biserial r | dotted lines = |log2FC| = 0.5\n",
      "n=", n_ctrl, " control section(s): effect size is primary metric"
    ),
    x = "log2 fold change (T2D / Control)",
    y = "Rank-biserial r (effect size, -1 to +1)",
    colour = "Compartment", shape = "Compartment"
  ) +
  pub_theme

save_pdf(p3,
         file.path(OUT_DIR, "figures/03_effect_size_plot.pdf"),
         width = 10, height = 8)

# ---- FIG 4: Beta:Alpha ratio ------------------------------------------------
if (all(c("beta", "alpha") %in% ct_cols)) {

  ba_df <- all_prop %>%
    select(section, donor, status, beta, alpha) %>%
    group_by(section, donor, status) %>%
    summarise(mean_beta  = mean(beta,  na.rm = TRUE),
              mean_alpha = mean(alpha, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(ba_ratio = mean_beta / (mean_alpha + 1e-6),
           status   = factor(status, levels = c("Control", "T2D")))

  p4 <- ggplot(ba_df,
               aes(x = status, y = ba_ratio,
                   fill = status, colour = status)) +
    geom_violin(scale = "width", trim = FALSE,
                alpha = 0.7, linewidth = 0.4) +
    geom_boxplot(width = 0.2, fill = "white",
                 outlier.shape = NA, linewidth = 0.5) +
    geom_jitter(width = 0.08, size = 3, alpha = 0.85,
                shape = 21, colour = "grey20") +
    scale_fill_manual(values   = status_colors, guide = "none") +
    scale_colour_manual(values = status_colors, guide = "none") +
    labs(
      title    = "Beta:Alpha cell ratio in T2D pancreas (spatial)",
      subtitle = paste0(
        "RCTD section-mean proportions | each point = one Visium section\n",
        "Reduced ratio expected in T2D: beta cell loss + alpha dominance"
      ),
      x = NULL, y = "Beta:Alpha proportion ratio"
    ) +
    pub_theme

  save_pdf(p4,
           file.path(OUT_DIR, "figures/04_beta_alpha_ratio.pdf"),
           width = 7, height = 7)

  log_msg("Beta:Alpha ratio per section:")
  print(ba_df[, c("section", "status", "mean_beta",
                   "mean_alpha", "ba_ratio")])
}

# ---- FIG 5: SC-TRIAD immune convergence -------------------------------------
imm_present <- intersect(c("macrophage", "T_cell"), ct_cols)
if (length(imm_present) >= 1) {

  imm_df <- all_prop %>%
    select(section, donor, status, all_of(imm_present)) %>%
    group_by(section, donor, status) %>%
    summarise(across(all_of(imm_present),
                     \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(all_of(imm_present),
                 names_to = "cell_type", values_to = "mean_proportion") %>%
    mutate(status = factor(status, levels = c("Control", "T2D")))

  imm_sum <- imm_df %>%
    group_by(status, cell_type) %>%
    summarise(mean = mean(mean_proportion, na.rm = TRUE),
              se   = sd(mean_proportion,   na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  p5 <- ggplot(imm_sum,
               aes(x = cell_type, y = mean,
                   fill = status, colour = status)) +
    geom_bar(stat = "identity",
             position = position_dodge(0.75),
             width = 0.65, alpha = 0.88) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.75),
                  width = 0.25, linewidth = 0.6, colour = "grey30") +
    geom_point(data = imm_df,
               aes(y = mean_proportion, group = status),
               position = position_dodge(0.75),
               size = 2, alpha = 0.75,
               shape = 21, colour = "grey20") +
    scale_fill_manual(values   = status_colors) +
    scale_colour_manual(values = status_colors) +
    scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
    labs(
      title    = "SC-TRIAD: immune infiltration in T2D pancreas (spatial)",
      subtitle = paste0(
        "RCTD mean proportions \u00b1 SE | each point = one Visium section\n",
        "Connects circulating PBMC immune signatures to tissue-resident states"
      ),
      x = NULL, y = "Mean proportion",
      fill = "Status", colour = "Status"
    ) +
    pub_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"))

  save_pdf(p5,
           file.path(OUT_DIR,
                     "figures/05_sctriad_pbmc_pancreas_convergence.pdf"),
           width = 9, height = 7)
}

# ---- FIG 6: Per-section heatmap ---------------------------------------------
top_ct <- intersect(head(diff_df$cell_type, min(15, nrow(diff_df))), ct_cols)

hm_in <- summary_section %>%
  select(section, status, all_of(top_ct)) %>%
  mutate(status = factor(status, levels = c("Control", "T2D"))) %>%
  arrange(status, section)

hm_mat <- apply(as.matrix(hm_in[, top_ct, drop = FALSE]), 2, function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
})
rownames(hm_mat) <- hm_in$section

hm_long <- as.data.frame(hm_mat) %>%
  rownames_to_column("section") %>%
  left_join(hm_in[, c("section", "status")], by = "section") %>%
  pivot_longer(-c(section, status),
               names_to = "cell_type", values_to = "scaled") %>%
  mutate(section   = factor(section, levels = hm_in$section),
         cell_type = factor(cell_type, levels = top_ct))

p6 <- ggplot(hm_long,
             aes(x = section, y = cell_type, fill = scaled)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low      = "#2166AC",
                       mid      = "white",
                       high     = "#C00000",
                       midpoint = 0.5,
                       name     = "Scaled\nproportion") +
  facet_grid(. ~ status, scales = "free_x", space = "free_x") +
  labs(
    title    = "Per-section cell type proportions: T2D vs Control",
    subtitle = paste0(
      "Top ", length(top_ct),
      " cell types by |log2FC| | scaled within each cell type"
    ),
    x = NULL, y = NULL
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9, face = "italic")
  )

save_pdf(p6,
         file.path(OUT_DIR, "figures/06_section_heatmap.pdf"),
         width = 12, height = 8)


# =============================================================================
# Done
# =============================================================================
log_msg("\n=== ANALYSIS COMPLETE ===")
log_msg("Output directory: ", OUT_DIR)
log_msg("")
log_msg("METHODS NOTE (copy to thesis / dissertation):")
log_msg("  Spatial deconvolution was performed using RCTD (Cable et al. 2022)")
log_msg("  with the Baron et al. 2016 (GSE84133) human pancreas scRNA-seq")
log_msg("  dataset as reference. Tissue spot coordinates were unavailable from")
log_msg("  the GEO submission (GSE264331); placeholder grid coordinates were")
log_msg("  assigned for RCTD input. RCTD performs maximum-likelihood estimation")
log_msg("  independently per spot, so coordinates do not affect deconvolution")
log_msg("  accuracy or results. Differential cell type proportions (n=",
        n_t2d, " T2D, n=", n_ctrl, " control section(s)) are reported as")
log_msg("  log2 fold change and rank-biserial correlation coefficient.")
log_msg("  Wilcoxon p-values were not computed due to the single control donor")
log_msg("  available in this dataset (GSE264331).")