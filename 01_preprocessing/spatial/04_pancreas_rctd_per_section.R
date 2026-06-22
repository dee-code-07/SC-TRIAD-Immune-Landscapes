# sc-triad project
# script: 03a_rctd_per_section.R  — v2 (spot subsampling fix)
#
# Root cause of time limit kills:
#   chooseSigma() inside create.RCTD() scales O(n_spots). With ~4800 spots
#   and 15424 genes, it takes >12 hours before run.RCTD() even starts.
#   All 16 array tasks were killed at this exact step.
#
# Fix: subsample MAX_SPOTS = 2000 spots per section (UMI-weighted random).
#   Methodologically valid because:
#     (a) We compute section-level MEAN proportions, not spot-level spatial maps
#     (b) 2000 spots is statistically representative of 4800
#     (c) UMI-weighted sampling preferentially selects information-rich spots
#   Expected runtime after fix: ~3-4 hours per section (within 12h wall time).
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(Matrix)
  library(stringr)
})

set.seed(42)
options(future.globals.maxSize = 100 * 1024^3)

MAX_SPOTS <- 500L   # subsample target — reduces chooseSigma from >12h to ~3-4h

BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
PAN_DIR <- file.path(BASE, "01_raw_data/01_t2d/gse264331_pancreas_visium")
REF_DIR <- file.path(BASE, "01_raw_data/01_t2d/baron2016_pancreas_ref")
OUT_DIR <- file.path(BASE, "04_spatial/pancreas")

for (d in c("logs", "objects"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

args        <- commandArgs(trailingOnly = TRUE)
section_idx <- as.integer(args[1])
if (is.na(section_idx) || section_idx < 1)
  stop("Usage: Rscript 03a_rctd_per_section.R <section_index 1..16>")

LOG_FILE <- file.path(OUT_DIR, "logs",
                      paste0("rctd_section_", section_idx, ".log"))
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log("=== RCTD per-section v2 | index: ", section_idx,
    " | MAX_SPOTS: ", MAX_SPOTS, " ===")

# metadata
h5_files <- list.files(PAN_DIR, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) == 0) stop("No .h5 files in: ", PAN_DIR)
if (section_idx > length(h5_files))
  stop("section_idx ", section_idx, " > n files (", length(h5_files), ")")

DONOR_STATUS <- c(
  "HP21024" = "T2D",  "HP21035" = "T2D",  "HP21091" = "T2D",
  "HP21161" = "T2D",  "HP21177" = "T2D",  "HP21181" = "Control"
)

parse_meta <- function(fp) {
  fname   <- basename(fp)
  donor   <- str_extract(fname, "HP\\d+")
  section <- str_extract(fname, "Vis\\d+_S\\d+_\\w+")
  if (is.na(section)) section <- str_extract(fname, "Vis\\d+_S\\d+")
  status  <- ifelse(!is.na(donor) && donor %in% names(DONOR_STATUS),
                    DONOR_STATUS[[donor]], "Unknown")
  data.frame(file = fname, donor = donor, section = section,
             status = status, stringsAsFactors = FALSE)
}

h5_meta  <- do.call(rbind, lapply(h5_files, parse_meta))
this_h5  <- h5_files[section_idx]
this_sid <- paste0(h5_meta$donor[section_idx], "_",
                   h5_meta$section[section_idx])

log("section ID : ", this_sid)
log("status     : ", h5_meta$status[section_idx])

# checkpoint
rctd_rds <- file.path(OUT_DIR, "objects", paste0("rctd_", this_sid, ".rds"))
if (file.exists(rctd_rds)) {
  log("checkpoint found — already done. exiting.")
  quit(status = 0)
}

# load Baron reference
baron_rds <- file.path(REF_DIR, "baron2016_reference.rds")
if (!file.exists(baron_rds))
  stop("Baron reference RDS not found: ", baron_rds)

log("loading Baron reference ...")
baron_ref <- readRDS(baron_rds)
log("loaded: ", ncol(baron_ref), " cells | ",
    length(unique(baron_ref$cell_type)), " types")

ct_tab <- table(baron_ref$cell_type)
rare   <- names(ct_tab[ct_tab < 10])
if (length(rare) > 0)
  baron_ref <- baron_ref[, !baron_ref$cell_type %in% rare]

final_types <- sort(unique(baron_ref$cell_type))
log("types (", length(final_types), "): ", paste(final_types, collapse = ", "))

set.seed(42)
keep_idx <- unlist(lapply(final_types, function(ct) {
  idx <- which(baron_ref$cell_type == ct)
  sample(idx, min(length(idx), 300L))
}))
ref_ds <- baron_ref[, keep_idx]

if (!"data" %in% Layers(ref_ds, assay = "RNA"))
  ref_ds <- NormalizeData(ref_ds, verbose = FALSE)

ref_counts <- round(GetAssayData(ref_ds, assay = "RNA", layer = "counts"))
ref_numi   <- as.integer(round(Matrix::colSums(ref_counts)))
ref_labels <- factor(ref_ds$cell_type)
names(ref_numi)   <- colnames(ref_ds)
names(ref_labels) <- colnames(ref_ds)

# load h5
log("loading: ", basename(this_h5))
mat_st <- Read10X_h5(this_h5)
if (is.list(mat_st)) mat_st <- mat_st[["Gene Expression"]]
colnames(mat_st) <- paste0(this_sid, "_", colnames(mat_st))
log("raw: ", nrow(mat_st), " genes x ", ncol(mat_st), " spots")

# gene overlap
common_genes   <- intersect(rownames(mat_st), rownames(ref_counts))
log("common genes: ", length(common_genes))
if (length(common_genes) < 2000)
  stop("Only ", length(common_genes), " common genes.")

mat_sub        <- mat_st[common_genes, , drop = FALSE]
ref_counts_sub <- ref_counts[common_genes, ]

# nUMI warning is expected and harmless
reference_sub <- Reference(counts     = ref_counts_sub,
                           cell_types = ref_labels,
                           nUMI       = ref_numi)

# filter low-UMI spots
spot_numi        <- as.integer(round(Matrix::colSums(mat_sub)))
names(spot_numi) <- colnames(mat_sub)
keep             <- spot_numi >= 100
log("spots after UMI>=100: ", sum(keep), " / ", ncol(mat_sub))
if (sum(keep) == 0) stop("No spots with UMI >= 100.")

mat_pass  <- mat_sub[, keep, drop = FALSE]
numi_pass <- spot_numi[keep]
n_pass    <- sum(keep)

# SPOT SUBSAMPLING — the fix for chooseSigma timeout
if (n_pass > MAX_SPOTS) {
  set.seed(section_idx)   # per-section seed for reproducibility
  probs      <- numi_pass / sum(numi_pass)   # UMI-weighted
  sel        <- sample(n_pass, MAX_SPOTS, replace = FALSE, prob = probs)
  mat_final  <- mat_pass[, sel, drop = FALSE]
  numi_final <- numi_pass[sel]
  log("subsampled: ", n_pass, " -> ", MAX_SPOTS, " spots (UMI-weighted)")
} else {
  mat_final  <- mat_pass
  numi_final <- numi_pass
  log("no subsampling: ", n_pass, " <= ", MAX_SPOTS)
}

n_spots    <- ncol(mat_final)
spot_names <- colnames(mat_final)

# placeholder grid coordinates
side   <- ceiling(sqrt(n_spots))
coords <- data.frame(
  x = ((seq_len(n_spots) - 1L) %% side) + 1L,
  y = ((seq_len(n_spots) - 1L) %/% side) + 1L,
  row.names = spot_names
)

# run RCTD
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
log("RCTD full mode | cores: ", N_CORES,
    " | spots: ", n_spots, " | genes: ", length(common_genes))

puck <- SpatialRNA(coords = coords, counts = mat_final, nUMI = numi_final)

log("create.RCTD ... (chooseSigma is the slow step, ~2-3h)")
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

log("run.RCTD ...")
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

saveRDS(myRCTD, rctd_rds)
log("saved: ", basename(rctd_rds))

weights  <- myRCTD@results$weights
prop_mat <- normalize_weights(weights)
log("spots deconvolved : ", nrow(prop_mat))
log("mean beta         : ", round(mean(prop_mat[, "beta"],  na.rm = TRUE), 4))
log("mean alpha        : ", round(mean(prop_mat[, "alpha"], na.rm = TRUE), 4))
log("=== COMPLETE: ", this_sid, " ===")