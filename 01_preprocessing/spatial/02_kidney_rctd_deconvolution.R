# sc-triad project
# script: 02_rctd_deconvolution.R
# purpose: RCTD cell type deconvolution of kidney Visium spatial data
#          using GSE211785 scRNA/snRNA as reference
#
# fixes vs corrupted version:
#   - st_sub created BEFORE spot_numi (ordering was broken by sed)
#   - ref_numi uses as.integer(round()) — spacexr requires integer nUMI
#   - spot_numi uses as.integer(round()) — same requirement
#   - ref_counts left as dgCMatrix with double @x — spacexr requires double
#   - no storage.mode() calls — incompatible with S4 sparse matrices
#   - all figures saved as PDF
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

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
})

set.seed(42)
options(future.globals.maxSize = 100 * 1024^3)

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
ST_DIR  <- file.path(BASE, "01_raw_data/02_htn/gse211785_kidney/visium")
OUT_DIR <- file.path(BASE, "04_spatial/kidney")

for (d in c("logs", "objects", "figures", "tables"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "logs/02_rctd.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

# ── publication PDF theme ──────────────────────────────────────────────────────
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

save_pdf <- function(plot, path, width = 10, height = 8) {
  pdf(path, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  log(paste("  PDF saved:", basename(path)))
}

# ── cell type group definitions ────────────────────────────────────────────────
status_colors <- c(Control = "#7CAE00", Disease = "#C77CFF")

immune_types  <- c("CD4T", "CD8T", "NK", "B_Naive", "B_memory",
                   "CD14_Mono", "CD16_Mono", "Mac", "Neutrophil",
                   "Baso_Mast", "pDC", "cDC", "Plasma_Cells")
tubular_types <- c("PT_S1", "PT_S2", "PT_S3", "iPT",
                   "M_TAL", "C_TAL", "DCT1", "DCT2",
                   "CNT", "PC", "IC_A", "IC_B",
                   "Des-Thin_Limb", "Ascending_Thin_LOH", "Macula_Densa")
stromal_types <- c("Fibroblast_1", "Fibroblast_2", "MyoFib_VSMC",
                   "GS_Stromal", "Mes", "Neural_Cells")
vascular_types <- c("Endo_GC", "Endo_Peritubular", "Endo_Lymphatic",
                    "Podo", "PEC")

log("=== RCTD deconvolution pipeline ===")
log(paste("spacexr:", packageVersion("spacexr")))
log(paste("Seurat:", packageVersion("Seurat")))

# ── step 1: load ST count matrix ───────────────────────────────────────────────
log("loading ST count matrix...")
st_counts <- readRDS(
  file.path(ST_DIR, "GSE211785_EXPORT_ST_counts.rds")
)
log(paste("ST matrix:", nrow(st_counts), "genes x", ncol(st_counts), "spots"))
log(paste("ST values are integers:",
          all(st_counts@x == floor(st_counts@x))))

# ── step 2: load ST metadata ───────────────────────────────────────────────────
log("loading ST metadata...")
st_meta <- read.table(
  gzfile(file.path(ST_DIR, "GSE211785_ST_metadata.txt.gz")),
  header    = TRUE,
  sep       = "\t",
  row.names = 1
)
log(paste("ST metadata spots:", nrow(st_meta)))
log("Status distribution:")
print(table(st_meta$Status))
log("Sample distribution:")
print(table(st_meta$orig.ident))

# ── step 3: align ST counts and metadata ──────────────────────────────────────
common_spots <- intersect(colnames(st_counts), rownames(st_meta))
log(paste("matching spots:", length(common_spots),
          "of", ncol(st_counts), "total"))

if (length(common_spots) < ncol(st_counts) * 0.95)
  stop("Less than 95% of ST spots found in metadata — check barcode format")

st_counts <- st_counts[, common_spots]
st_meta   <- st_meta[common_spots, , drop = FALSE]
log("ST counts and metadata aligned")

# ── step 4: load kidney reference ─────────────────────────────────────────────
log("loading kidney reference (HKD + Control)...")
ref_obj <- readRDS(
  file.path(OUT_DIR, "objects/kidney_ref_HKD_Control.rds")
)
log(paste("reference cells:", ncol(ref_obj)))
log("cell type counts:")
print(sort(table(ref_obj$cell_type), decreasing = TRUE))

# ── step 5: clean reference cell types ────────────────────────────────────────
log("cleaning reference cell types...")

# rename cell types containing "/" — prohibited by spacexr
ref_obj$cell_type[ref_obj$cell_type == "MyoFib/VSMC"] <- "MyoFib_VSMC"
ref_obj$cell_type[ref_obj$cell_type == "Baso/Mast"]   <- "Baso_Mast"
log("  MyoFib/VSMC renamed to MyoFib_VSMC")
log("  Baso/Mast renamed to Baso_Mast")

# remove RBC (artefact in SC_RNA, not present in Visium tissue sections)
n_before  <- ncol(ref_obj)
ref_obj   <- ref_obj[, ref_obj$cell_type != "RBC"]
log(paste("  RBC removed:", n_before - ncol(ref_obj), "cells"))

# check remaining counts
ct_counts <- sort(table(ref_obj$cell_type), decreasing = TRUE)
log("cell type counts after cleaning:")
print(ct_counts)

# remove types below RCTD minimum (25 cells)
rare_types <- names(ct_counts[ct_counts < 25])
if (length(rare_types) > 0) {
  log(paste("removing rare types (<25 cells):",
            paste(rare_types, collapse = ", ")))
  ref_obj <- ref_obj[, !ref_obj$cell_type %in% rare_types]
}

final_types <- sort(unique(ref_obj$cell_type))
log(paste("final cell types:", length(final_types)))

# ── step 6: downsample reference ──────────────────────────────────────────────
MAX_CELLS <- 500L
log(paste("downsampling to max", MAX_CELLS, "cells per type..."))

keep_cells <- unlist(lapply(final_types, function(ct) {
  ct_cells <- colnames(ref_obj)[ref_obj$cell_type == ct]
  sample(ct_cells, min(length(ct_cells), MAX_CELLS))
}))
ref_ds <- ref_obj[, keep_cells]
log(paste("downsampled reference:", ncol(ref_ds), "cells"))
print(sort(table(ref_ds$cell_type), decreasing = TRUE))

# ── step 7: extract raw counts for reference ───────────────────────────────────
# dgCMatrix @x MUST be double — do NOT coerce to integer
# only nUMI vectors need to be R integer type (spacexr check_UMI requirement)

log("extracting reference counts...")
ref_counts <- GetAssayData(ref_ds, assay = "RNA", layer = "counts")

# round to whole numbers; leave as double (dgCMatrix requirement)
ref_counts <- round(ref_counts)
log(paste("ref_counts class:", class(ref_counts)))
log(paste("ref_counts @x type:", typeof(ref_counts@x), "(must be double)"))
log(paste("ref_counts whole numbers:",
          all(ref_counts@x == floor(ref_counts@x))))

# nUMI for reference — must be R integer type
ref_numi <- as.integer(round(Matrix::colSums(ref_counts)))
names(ref_numi) <- colnames(ref_ds)
log(paste("ref_numi type:", typeof(ref_numi), "(must be integer)"))

# cell type labels
ref_labels <- factor(ref_ds$cell_type)
names(ref_labels) <- colnames(ref_ds)

# ── step 8: find gene overlap ──────────────────────────────────────────────────
common_genes <- intersect(rownames(ref_counts), rownames(st_counts))
log(paste("genes in reference:", nrow(ref_counts)))
log(paste("genes in ST:        ", nrow(st_counts)))
log(paste("common genes:       ", length(common_genes)))

if (length(common_genes) < 3000)
  warning(paste("Only", length(common_genes),
                "common genes — deconvolution accuracy may be reduced"))

# subset BOTH to common genes
ref_counts <- ref_counts[common_genes, ]

# st_sub: ST matrix subset to common genes
# IMPORTANT: st_sub must be created BEFORE spot_numi
st_sub <- st_counts[common_genes, ]
log(paste("st_sub created:", nrow(st_sub), "genes x", ncol(st_sub), "spots"))

# ── step 9: build spacexr Reference ───────────────────────────────────────────
log("building spacexr Reference object...")
reference <- Reference(
  counts     = ref_counts,
  cell_types = ref_labels,
  nUMI       = ref_numi
)
log(paste("Reference created:", length(levels(reference@cell_types)),
          "cell types"))

# ── step 10: build SpatialRNA object ──────────────────────────────────────────
# tissue_positions_list.csv not available from GEO
# placeholder grid coordinates assigned — deconvolution is spot-independent
# all biological results unaffected by coordinate values
# methods note: "spatial coordinates unavailable from GEO submission;
# RCTD performed on count data; results summarized by disease status"

log("building SpatialRNA object (placeholder grid coordinates)...")

n_spots    <- ncol(st_sub)
spot_names <- colnames(st_sub)

# assign sequential grid positions
side     <- ceiling(sqrt(n_spots))
grid_x   <- ((seq_len(n_spots) - 1) %% side) + 1
grid_y   <- ((seq_len(n_spots) - 1) %/% side) + 1
coords   <- data.frame(x = grid_x, y = grid_y,
                       row.names = spot_names)

# spot nUMI — must be R integer type (created AFTER st_sub)
spot_numi <- as.integer(round(Matrix::colSums(st_sub)))
names(spot_numi) <- spot_names
log(paste("spot_numi type:", typeof(spot_numi), "(must be integer)"))

puck <- SpatialRNA(
  coords = coords,
  counts = st_sub,
  nUMI   = spot_numi
)
log(paste("SpatialRNA object created:", n_spots, "spots"))

# ── step 11: run RCTD ─────────────────────────────────────────────────────────
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
log(paste("running RCTD full mode |", N_CORES, "cores"))
log("expected runtime: 4-8 hours for 37,143 spots")

myRCTD <- create.RCTD(
  spatialRNA        = puck,
  reference         = reference,
  max_cores         = N_CORES,
  CELL_MIN_INSTANCE = 25,
  gene_cutoff       = 0.000125,
  fc_cutoff         = 0.5,
  UMI_min           = 100,
  UMI_min_sigma     = 300,
  class_df          = NULL
)

myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
log("RCTD complete")

saveRDS(myRCTD, file.path(OUT_DIR, "objects/rctd_results_full.rds"))
log("RCTD object saved")

# ── step 12: extract and normalize proportions ────────────────────────────────
log("extracting RCTD results...")
weights     <- myRCTD@results$weights
prop_matrix <- normalize_weights(weights)
log(paste("proportion matrix:", nrow(prop_matrix), "spots x",
          ncol(prop_matrix), "cell types"))

# add metadata
prop_df <- as.data.frame(prop_matrix)
prop_df$spot              <- rownames(prop_matrix)
prop_df$Status            <- st_meta[rownames(prop_matrix), "Status"]
prop_df$orig.ident        <- st_meta[rownames(prop_matrix), "orig.ident"]
prop_df$published_celltype <- st_meta[rownames(prop_matrix), "celltype"]
prop_df$nCount            <- st_meta[rownames(prop_matrix), "nCount_Spatial"]

ct_cols <- colnames(prop_matrix)

write.csv(prop_df,
          file.path(OUT_DIR, "tables/rctd_proportions_per_spot.csv"),
          row.names = FALSE)
log(paste("proportions saved:", nrow(prop_df), "spots"))

# ── step 13: summarise by status and sample ───────────────────────────────────
log("summarising by disease status...")

summary_status <- prop_df %>%
  group_by(Status) %>%
  summarise(across(all_of(ct_cols), \(x) mean(x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")

write.csv(summary_status,
          file.path(OUT_DIR, "tables/rctd_summary_by_status.csv"),
          row.names = FALSE)

summary_sample <- prop_df %>%
  group_by(orig.ident, Status) %>%
  summarise(across(all_of(ct_cols), \(x) mean(x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")

write.csv(summary_sample,
          file.path(OUT_DIR, "tables/rctd_summary_by_sample.csv"),
          row.names = FALSE)

# ── step 14: differential cell type analysis ───────────────────────────────────
# sample-level means as input (pseudobulk) to avoid pseudoreplication
log("differential cell type proportions (Disease vs Control)...")

diff_results <- lapply(ct_cols, function(ct) {
  ctrl_vals <- summary_sample %>%
    filter(Status == "Control") %>% pull(!!sym(ct))
  dis_vals  <- summary_sample %>%
    filter(Status == "Disease") %>% pull(!!sym(ct))

  if (length(ctrl_vals) < 2 || length(dis_vals) < 2) return(NULL)

  wt <- tryCatch(
    wilcox.test(dis_vals, ctrl_vals, exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(wt)) return(NULL)

  data.frame(
    cell_type      = ct,
    mean_ctrl      = mean(ctrl_vals, na.rm = TRUE),
    mean_disease   = mean(dis_vals,  na.rm = TRUE),
    delta          = mean(dis_vals,  na.rm = TRUE) -
                     mean(ctrl_vals, na.rm = TRUE),
    p_value        = wt$p.value,
    n_ctrl_samples = length(ctrl_vals),
    n_dis_samples  = length(dis_vals)
  )
})

diff_df <- do.call(rbind, Filter(Negate(is.null), diff_results))
diff_df$p_adj  <- p.adjust(diff_df$p_value, method = "BH")
diff_df$log2FC <- log2((diff_df$mean_disease + 1e-6) /
                        (diff_df$mean_ctrl    + 1e-6))
diff_df <- diff_df[order(diff_df$p_adj), ]

write.csv(diff_df,
          file.path(OUT_DIR, "tables/rctd_differential_celltypes.csv"),
          row.names = FALSE)

log("top differential cell types:")
print(head(diff_df[, c("cell_type", "mean_ctrl", "mean_disease",
                        "delta", "log2FC", "p_adj")], 15))

# ── step 15: FIGURE 1 — reference composition ─────────────────────────────────
log("generating figures...")

ref_comp <- as.data.frame(table(ref_ds$cell_type)) %>%
  rename(cell_type = Var1, n_cells = Freq) %>%
  mutate(
    compartment = case_when(
      cell_type %in% tubular_types  ~ "Tubular",
      cell_type %in% vascular_types ~ "Vascular/Glomerular",
      cell_type %in% immune_types   ~ "Immune",
      cell_type %in% stromal_types  ~ "Stromal",
      TRUE ~ "Other"
    ),
    cell_type = reorder(cell_type, n_cells)
  )

p1 <- ggplot(ref_comp, aes(x = n_cells, y = cell_type,
                            fill = compartment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    Tubular                = "#2171B5",
    "Vascular/Glomerular"  = "#FB6A4A",
    Immune                 = "#238B45",
    Stromal                = "#FD8D3C",
    Other                  = "#BDBDBD"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title    = "RCTD reference: cell type composition",
       subtitle = paste0(ncol(ref_ds), " cells | ",
                         length(final_types), " types | ",
                         "max 500 per type | RBC excluded"),
       x = "Number of cells", y = NULL,
       fill = "Compartment") +
  pub_theme +
  theme(axis.text.y = element_text(size = 8))

save_pdf(p1,
         file.path(OUT_DIR, "figures/01_reference_composition.pdf"),
         width = 10, height = 10)

# ── step 16: FIGURE 2 — cell type proportions by status ───────────────────────
mean_prop <- prop_df %>%
  select(Status, all_of(ct_cols)) %>%
  pivot_longer(all_of(ct_cols),
               names_to  = "cell_type",
               values_to = "proportion") %>%
  group_by(Status, cell_type) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    Status    = factor(Status, levels = c("Control", "Disease")),
    cell_type = factor(cell_type,
                       levels = rev(diff_df$cell_type))
  )

p2 <- ggplot(mean_prop,
             aes(x = Status, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  labs(title    = "Kidney Visium: cell type composition by disease status",
       subtitle = paste0("RCTD full mode | ",
                         sum(prop_df$Status == "Control"),
                         " Control spots | ",
                         sum(prop_df$Status == "Disease"),
                         " Disease (HKD) spots"),
       x = NULL, y = "Mean proportion",
       fill = "Cell type") +
  pub_theme +
  theme(legend.key.size = unit(0.35, "cm"),
        legend.text     = element_text(size = 7))

save_pdf(p2,
         file.path(OUT_DIR, "figures/02_celltype_proportions_by_status.pdf"),
         width = 10, height = 9)

# ── step 17: FIGURE 3 — differential proportions volcano ──────────────────────
diff_plot <- diff_df %>%
  mutate(
    sig = case_when(
      p_adj < 0.05 & delta > 0  ~ "Enriched in HKD",
      p_adj < 0.05 & delta < 0  ~ "Depleted in HKD",
      TRUE ~ "NS"
    ),
    compartment = case_when(
      cell_type %in% tubular_types  ~ "Tubular",
      cell_type %in% vascular_types ~ "Vascular",
      cell_type %in% immune_types   ~ "Immune",
      cell_type %in% stromal_types  ~ "Stromal",
      TRUE ~ "Other"
    ),
    label = ifelse(p_adj < 0.1, as.character(cell_type), NA_character_)
  )

p3 <- ggplot(diff_plot,
             aes(x = log2FC, y = -log10(p_adj + 1e-10),
                 colour = sig, shape = compartment)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_text_repel(aes(label = label), size = 3.2, na.rm = TRUE,
                  max.overlaps = 20, colour = "black",
                  fontface = "italic") +
  scale_colour_manual(values = c(
    "Enriched in HKD"  = "#C00000",
    "Depleted in HKD"  = "#2166AC",
    "NS"               = "grey75"
  )) +
  scale_shape_manual(values = c(
    Tubular = 16, Vascular = 17,
    Immune  = 15, Stromal  = 18, Other = 4
  )) +
  labs(title    = "Differential cell type proportions: HKD vs Control",
       subtitle = "RCTD full mode | sample-level Wilcoxon | BH-adjusted p",
       x        = "log2 fold change (HKD / Control)",
       y        = expression(-log[10](p[adj])),
       colour   = "Direction", shape = "Compartment") +
  pub_theme

save_pdf(p3,
         file.path(OUT_DIR, "figures/03_differential_celltypes_volcano.pdf"),
         width = 10, height = 8)

# ── step 18: FIGURE 4 — immune infiltration ────────────────────────────────────
immune_present <- intersect(immune_types, ct_cols)

immune_df <- prop_df %>%
  select(spot, Status, all_of(immune_present)) %>%
  pivot_longer(all_of(immune_present),
               names_to  = "cell_type",
               values_to = "proportion") %>%
  filter(proportion > 0.001) %>%
  mutate(
    Status    = factor(Status, levels = c("Control", "Disease")),
    cell_type = factor(cell_type, levels = immune_present)
  )

p4 <- ggplot(immune_df,
             aes(x = Status, y = proportion,
                 fill = Status, colour = Status)) +
  geom_violin(scale = "width", trim = FALSE,
              alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white",
               outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values   = status_colors, guide = "none") +
  scale_colour_manual(values = status_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) +
  labs(title    = "Immune cell infiltration in HKD kidney (spatial)",
       subtitle = paste0(
         "RCTD spot-level proportions | spots with >0.1% immune content\n",
         "HKD = Hypertensive Kidney Disease"
       ),
       x = NULL, y = "Cell type proportion") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_pdf(p4,
         file.path(OUT_DIR, "figures/04_immune_infiltration_HKD.pdf"),
         width = 14, height = 12)

# ── step 19: FIGURE 5 — SC-TRIAD cross-tissue convergence ─────────────────────
# key manuscript figure: immune states elevated in PBMC (CellChat / DEG results)
# are also spatially enriched in HKD kidney tissue
# connects circulating findings to tissue-resident architecture

cross_types <- intersect(
  c("NK", "CD8T", "CD4T", "CD16_Mono", "CD14_Mono",
    "Mac", "Endo_Peritubular", "Fibroblast_1"),
  ct_cols
)

ct_sample <- prop_df %>%
  select(orig.ident, Status, all_of(cross_types)) %>%
  group_by(orig.ident, Status) %>%
  summarise(across(all_of(cross_types), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_longer(all_of(cross_types),
               names_to  = "cell_type",
               values_to = "mean_proportion") %>%
  mutate(
    Status    = factor(Status, levels = c("Control", "Disease")),
    cell_type = factor(cell_type, levels = cross_types)
  )

ct_summary <- ct_sample %>%
  group_by(Status, cell_type) %>%
  summarise(
    mean = mean(mean_proportion, na.rm = TRUE),
    se   = sd(mean_proportion,  na.rm = TRUE) /
           sqrt(sum(!is.na(mean_proportion))),
    .groups = "drop"
  )

p5 <- ggplot(ct_summary,
             aes(x = cell_type, y = mean,
                 fill = Status, colour = Status)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.25, linewidth = 0.6, colour = "grey30"
  ) +
  geom_point(
    data     = ct_sample,
    aes(x = cell_type, y = mean_proportion, group = Status),
    position = position_dodge(width = 0.8),
    size = 1.8, alpha = 0.6, shape = 21, colour = "grey30"
  ) +
  scale_fill_manual(values   = status_colors) +
  scale_colour_manual(values = status_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title    = "SC-TRIAD: tissue-resident immune landscape in HKD kidney",
    subtitle = paste0(
      "RCTD mean proportions ± SE per Visium section | points = individual sections\n",
      "Connects circulating PBMC signatures to tissue-resident spatial states"
    ),
    x = NULL, y = "Mean cell type proportion",
    fill = "Status", colour = "Status"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1,
                                    face = "italic", size = 10))

save_pdf(p5,
         file.path(OUT_DIR,
                   "figures/05_sctriad_pbmc_kidney_convergence.pdf"),
         width = 12, height = 7)

# ── step 20: FIGURE 6 — validation vs published labels ────────────────────────
prop_mat_only <- as.matrix(prop_df[, ct_cols])
rctd_dominant <- apply(prop_mat_only, 1, function(x) names(x)[which.max(x)])

broad_label <- function(x) {
  dplyr::case_when(
    x %in% tubular_types  ~ "Tubular",
    x %in% immune_types   ~ "Immune",
    x %in% stromal_types  ~ "Stromal",
    x %in% vascular_types ~ "Vascular",
    TRUE ~ "Other"
  )
}

# map published spot labels (different nomenclature) to broad categories
pub_to_broad <- function(x) {
  dplyr::case_when(
    grepl("PT|TAL|DCT|CNT|PC|IC|LOH|Macula|Thin", x) ~ "Tubular",
    grepl("Endo|Podo|PEC|GC|Peritubular|Lymphatic",  x) ~ "Vascular",
    grepl("Fib|Myo|VSMC|Strom|Mes|Neural",           x) ~ "Stromal",
    grepl("CD4|CD8|NK|Mac|Mono|Neutro|Baso|DC|Plasma|B_", x) ~ "Immune",
    TRUE ~ "Other"
  )
}

validation_df <- data.frame(
  spot          = prop_df$spot,
  rctd_dominant = rctd_dominant,
  rctd_broad    = broad_label(rctd_dominant),
  published     = prop_df$published_celltype,
  pub_broad     = pub_to_broad(prop_df$published_celltype),
  Status        = prop_df$Status,
  stringsAsFactors = FALSE
)

concordance <- mean(validation_df$rctd_broad == validation_df$pub_broad,
                    na.rm = TRUE)
log(paste("RCTD vs published broad concordance:",
          round(concordance * 100, 1), "%"))

conf_df <- as.data.frame(
  table(RCTD = validation_df$rctd_broad,
        Published = validation_df$pub_broad)
) %>%
  group_by(Published) %>%
  mutate(prop = Freq / sum(Freq)) %>%
  ungroup()

p6 <- ggplot(conf_df,
             aes(x = Published, y = RCTD, fill = prop)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            size = 3.5, colour = "black") +
  scale_fill_gradient2(
    low = "#f7fbff", mid = "#6baed6", high = "#08306b",
    midpoint = 0.5,
    labels   = percent_format(),
    name     = "Proportion"
  ) +
  labs(
    title    = "RCTD validation: dominant cell type vs published labels",
    subtitle = paste0("Broad category concordance: ",
                      round(concordance * 100, 1), "%\n",
                      "RCTD dominant = highest proportion per spot"),
    x = "Published label (broad category)",
    y = "RCTD dominant type (broad category)"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_pdf(p6,
         file.path(OUT_DIR,
                   "figures/06_rctd_validation_vs_published.pdf"),
         width = 9, height = 7)

write.csv(validation_df,
          file.path(OUT_DIR, "tables/rctd_validation_vs_published.csv"),
          row.names = FALSE)

# ── step 21: final summary ─────────────────────────────────────────────────────
log("")
log("=== RCTD pipeline complete ===")
log(paste("spots deconvolved:", nrow(prop_df)))
log(paste("cell types:", length(ct_cols)))
log(paste("RCTD vs published concordance:",
          round(concordance * 100, 1), "%"))

sig_diff <- diff_df[!is.na(diff_df$p_adj) & diff_df$p_adj < 0.05, ]
log(paste("significant cell type changes:", nrow(sig_diff)))

if (nrow(sig_diff) > 0) {
  enriched <- sig_diff$cell_type[sig_diff$delta > 0]
  depleted <- sig_diff$cell_type[sig_diff$delta < 0]
  if (length(enriched) > 0)
    log(paste("  enriched in HKD:", paste(enriched, collapse = ", ")))
  if (length(depleted) > 0)
    log(paste("  depleted in HKD:", paste(depleted, collapse = ", ")))
}

log("")
log("figures saved:")
for (f in list.files(file.path(OUT_DIR, "figures"), pattern = "\\.pdf$"))
  log(paste(" ", f))
log("next: run 03_spatial_signatures.R")
log("  projects PBMC cross-disease DEG signatures onto spatial spots")