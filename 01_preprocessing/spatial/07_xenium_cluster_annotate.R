# sc-triad project
# script: 02_xenium_cluster_annotate.R
# purpose: normalize, cluster, and annotate Xenium lung cells
#          using airway-specific canonical marker genes
#
# xenium normalization note:
#   SCTransform v2 is valid for Xenium (Seurat recommendation for targeted panels)
#   log-normalization also acceptable; we use SCTransform for consistency
#   with PBMC pipeline — enables cross-modality label transfer later
#
# marker gene sources:
#   - Epithelial: Travaglini et al., Nature 2020 (HLCA)
#   - Immune: Sikkema et al., Nature Medicine 2023
#   - Structural: Adams et al., NEJM 2020
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(clustree)
})

set.seed(42)
options(future.globals.maxSize = 16 * 1024^3)

BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
OUT_DIR <- file.path(BASE, "04_spatial", "asthma_lung")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
OBJ_DIR <- file.path(OUT_DIR, "objects")
LOG_FILE <- file.path(OUT_DIR, "02_xenium_cluster.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

save_pdf <- function(plot, fname, width = 12, height = 10) {
  path <- file.path(FIG_DIR, fname)
  pdf(path, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  log(paste("  saved:", fname))
}

log("Xenium cluster + annotate pipeline started")

# ── load QC object ─────────────────────────────────────────────────────────────
obj <- readRDS(file.path(OBJ_DIR, "xenium_merged_qc.rds"))
log(paste("loaded:", ncol(obj), "cells |", nrow(obj), "genes"))

# ── step 1: normalization ──────────────────────────────────────────────────────
# For Xenium targeted panels:
#   - SCTransform v2: recommended when panel is similar in size to HVG subset
#   - vars.to.regress: cell_area accounts for segmentation size differences
#     (larger cells capture more transcripts due to geometry)

DefaultAssay(obj) <- "Xenium"

# join layers if split
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])

log("running SCTransform v2 on Xenium assay...")
obj <- SCTransform(
  obj,
  assay           = "Xenium",
  vst.flavor      = "v2",
  vars.to.regress = if ("cell_area" %in% colnames(obj@meta.data)) "cell_area" else NULL,
  verbose         = FALSE
)
log(paste("SCT done | variable features:", length(VariableFeatures(obj))))

# ── step 2: PCA + UMAP ────────────────────────────────────────────────────────
# Xenium panels are ~300-500 genes — use all available features for PCA
# Do NOT apply 50-PC selection; use ALL genes as features
log("PCA...")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# elbow plot
p_elbow <- ElbowPlot(obj, ndims = 30) + theme_classic(base_size = 11)
save_pdf(p_elbow, "02_pca_elbow.pdf", width = 7, height = 5)

# select PCs: all that explain > 1% variance individually
stdev   <- obj[["pca"]]@stdev
pct_var <- stdev^2 / sum(stdev^2) * 100
n_pcs   <- max(which(pct_var > 1))
n_pcs   <- min(max(n_pcs, 10), 30)
log(paste("PCs selected:", n_pcs))

# Harmony if multiple samples
n_samples <- length(unique(obj$sample_id))
if (n_samples > 1) {
  log("Running Harmony batch correction by sample_id...")

  # FIX 1: Split the SCT assay layers by sample so Harmony knows what the batches are
  obj[["SCT"]] <- split(obj[["SCT"]], f = obj$sample_id)

  obj <- IntegrateLayers(
    object               = obj,
    method               = HarmonyIntegration,
    orig.reduction       = "pca",
    new.reduction        = "harmony",
    assay                = "SCT",
    normalization.method = "SCT",  # FIX 2: Explicitly declare SCT
    verbose              = FALSE
  )

  # FIX 3: Re-join layers so downstream FindClusters and DotPlot run normally
  obj[["SCT"]] <- JoinLayers(obj[["SCT"]])

  red_use <- "harmony"
} else {
  red_use <- "pca"
}

log(paste("UMAP reduction:", red_use))
obj <- RunUMAP(obj, reduction = red_use, dims = 1:n_pcs, verbose = FALSE)

# ── step 3: clustering ────────────────────────────────────────────────────────
obj <- FindNeighbors(obj, reduction = red_use, dims = 1:n_pcs, verbose = FALSE)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (res in resolutions) {
  obj <- FindClusters(obj, resolution = res, algorithm = 1,
                      cluster.name = paste0("snn_res.", res), verbose = FALSE)
}

# FIX: replaced em-dash with hyphen to avoid HPC mbcsToSbcs conversion warnings
p_clustree <- clustree(obj, prefix = "snn_res.") +
  labs(title = "Xenium lung - cluster stability")
save_pdf(p_clustree, "03_clustree.pdf", width = 12, height = 10)

# optimal resolution
counts_per_res <- sapply(paste0("snn_res.", resolutions),
                         function(x) length(unique(obj@meta.data[[x]])))
gains      <- diff(counts_per_res)
stable_idx <- which(gains < 2)[1]
opt_res    <- if (!is.na(stable_idx)) resolutions[stable_idx] else 0.4
Idents(obj) <- paste0("snn_res.", opt_res)
obj$seurat_clusters <- Idents(obj)
log(paste("optimal resolution:", opt_res,
          "| clusters:", length(unique(Idents(obj)))))

# ── step 4: canonical marker gene annotation ──────────────────────────────────
# Curated airway cell type markers
# Sources: HLCA (Sikkema 2023), LungMAP, Adams 2020
# Only genes likely to be in a ~300-500 gene Xenium panel are listed first;
# fallbacks added for each type. Run dotplot to see what's present.

canonical_markers <- list(
  # EPITHELIAL
  "Basal"           = c("KRT5", "KRT17", "TP63", "PDPN"),
  "Club"            = c("SCGB1A1", "SCGB3A2", "CYP2B6", "ALDH3A1"),
  "Goblet"          = c("MUC5AC", "MUC5B", "SPDEF", "FOXA3"),
  "Ciliated"        = c("FOXJ1", "DNAI1", "TUBB4B", "CCDC78"),
  "AT1"             = c("AGER", "PDPN", "HOPX", "CAV1"),
  "AT2"             = c("SFTPB", "SFTPC", "SFTPD", "ABCA3", "NAPSA"),
  "Ionocyte"        = c("FOXI1", "CFTR", "ATP6V1G3"),
  # STROMAL
  "Fibroblast"      = c("COL1A1", "COL1A2", "LUM", "DCN", "PDGFRA"),
  "Smooth_Muscle"   = c("ACTA2", "MYH11", "CNN1", "TAGLN"),
  "Pericyte"        = c("RGS5", "PDGFRB", "MCAM"),
  "Endothelial"     = c("PECAM1", "CDH5", "VWF", "CLDN5"),
  # IMMUNE
  "Macrophage"      = c("CD68", "MARCO", "MRC1", "MSR1", "FCGR3A"),
  "Monocyte"        = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN"),
  "DC"              = c("CD1C", "FCER1A", "CLEC9A", "XCR1"),
  "Mast"            = c("KIT", "TPSAB1", "CPA3", "HDC"),
  "Eosinophil"      = c("SIGLEC8", "CLC", "EPX", "RNASE2"),
  "T_cell"          = c("CD3E", "CD3D", "TRAC", "CD2"),
  "CD4_T"           = c("CD3E", "CD4", "IL7R", "CCR7"),
  "CD8_T"           = c("CD3E", "CD8A", "CD8B", "GZMB"),
  "NK"              = c("NKG7", "GNLY", "NCAM1", "PRF1"),
  "B_cell"          = c("CD19", "MS4A1", "CD79A"),
  "ILC2"            = c("GATA3", "IL1RL1", "KLRG1", "IL13"),
  "Plasma"          = c("IGHG1", "IGHG2", "IGKC", "MZB1"),
  "pDC"             = c("LILRA4", "CLEC4C", "IL3RA", "PLAC8"),
  "Neutrophil"      = c("FCGR3B", "CSF3R", "CXCR2", "S100A9")
)

# Find which markers are present in the panel
panel_genes <- rownames(obj)
log("Marker gene coverage in panel:")
for (ct in names(canonical_markers)) {
  avail <- intersect(canonical_markers[[ct]], panel_genes)
  log(paste0("  ", ct, ": ", length(avail), "/",
             length(canonical_markers[[ct]]), " - ",
             paste(avail, collapse = ", ")))
}

# Dotplot using available markers
all_markers_flat  <- unique(unlist(canonical_markers))
markers_in_panel  <- all_markers_flat[all_markers_flat %in% panel_genes]
markers_missing   <- all_markers_flat[!all_markers_flat %in% panel_genes]
log(paste("Markers in panel:", length(markers_in_panel),
          "| Not in panel:", length(markers_missing)))

if (length(markers_in_panel) >= 5) {
  p_dot <- DotPlot(obj, features = markers_in_panel,
                   group.by = "seurat_clusters", assay = "SCT") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
          axis.text.y = element_text(size = 8)) +
    # FIX: replaced em-dash with hyphen to avoid HPC mbcsToSbcs conversion warnings
    labs(title = "Xenium airway - canonical markers by cluster")
  save_pdf(p_dot, "04_dotplot_canonical_markers.pdf", width = 20, height = 8)
}

# ── step 5: manual annotation from dotplot evidence ──────────────────────────
# !! STUDENT NOTE: Inspect 04_dotplot_canonical_markers.pdf BEFORE running
# this section. Update the mapping below based on your dotplot.
# The mapping below is a biologically informed PRIOR based on expected
# lung cell type frequencies in asthma tissue — not ground truth.
# You MUST validate against your dotplot.

# Placeholder mapping — UPDATE after inspecting dotplot
# Format: "cluster_number" = "CellType"
# Priority annotation logic:
#   High KRT5/TP63 -> Basal
#   High SCGB1A1   -> Club
#   High MUC5AC    -> Goblet
#   High FOXJ1     -> Ciliated
#   High SFTPB/C   -> AT2
#   High COL1A1    -> Fibroblast
#   High ACTA2     -> Smooth Muscle
#   High CD68/MARCO -> Macrophage
#   High TPSAB1/KIT -> Mast
#   High CD3E      -> T cell
#   High NKG7      -> NK
#   High CD19      -> B cell

ANNOTATION_MAP <- c(
  "0"  = "Epithelial_Unknown",
  "1"  = "Epithelial_Unknown",
  "2"  = "Fibroblast",
  "3"  = "Macrophage",
  "4"  = "T_cell",
  "5"  = "Smooth_Muscle",
  "6"  = "Epithelial_Unknown",
  "7"  = "Endothelial",
  "8"  = "Mast",
  "9"  = "B_cell",
  # ADD THESE — inspect dotplot to confirm:
  "10" = "Epithelial_Unknown",
  "11" = "AT1",           # strong PDPN
  "12" = "Epithelial_Unknown",
  "13" = "Epithelial_Unknown",
  "14" = "Epithelial_Unknown",
  "15" = "Pericyte",      # PDGFRB signal
  "16" = "Epithelial_Unknown",
  "17" = "Epithelial_Unknown",
  "18" = "Macrophage",    # CD68-range signal
  "19" = "Epithelial_Unknown",
  "20" = "Epithelial_Unknown",
  "21" = "Epithelial_Unknown",
  "22" = "Endothelial",   # CDH5/VWF
  "23" = "Epithelial_Unknown"
)

# Apply annotation (safe: unmapped clusters -> "Unknown")
clusters_present <- as.character(unique(obj$seurat_clusters))
full_map <- setNames(
  ifelse(clusters_present %in% names(ANNOTATION_MAP),
         ANNOTATION_MAP[clusters_present],
         "Unknown"),
  clusters_present
)
obj$cell_type_xenium <- unname(full_map[as.character(obj$seurat_clusters)])
log("Cell type annotation applied (verify against dotplot)")

ct_counts <- sort(table(obj$cell_type_xenium), decreasing = TRUE)
log("Cell type distribution:")
for (ct in names(ct_counts)) log(paste0("  ", ct, ": ", ct_counts[ct]))

write.csv(as.data.frame(ct_counts),
          file.path(TAB_DIR, "xenium_celltype_counts.csv"))

# ── step 6: UMAP visualizations ───────────────────────────────────────────────
p_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                     label = TRUE, label.size = 3, pt.size = 0.5,
                     raster = FALSE) +
  theme_classic(base_size = 10) + labs(title = "Xenium clusters")

p_sample  <- DimPlot(obj, reduction = "umap", group.by = "sample_id",
                     pt.size = 0.3, raster = FALSE) +
  theme_classic(base_size = 10) + labs(title = "by sample")

p_celltype <- DimPlot(obj, reduction = "umap", group.by = "cell_type_xenium",
                      label = TRUE, label.size = 2.5, repel = TRUE,
                      pt.size = 0.3, raster = FALSE) +
  theme_classic(base_size = 9) + labs(title = "Xenium cell types")

save_pdf(p_cluster + p_sample, "05_umap_cluster_sample.pdf",
         width = 18, height = 8)
save_pdf(p_celltype, "06_umap_celltype.pdf", width = 12, height = 9)

# ── step 7: spatial coordinate plot ───────────────────────────────────────────
# Plot cells in physical tissue space (X, Y coordinates from Xenium)
# This is the key figure unique to spatial data

# FIX: SAMPLES list was defined in 01_xenium_qc_load.R and not available here.
# Replaced with a dynamic loop over unique sample IDs found in obj$sample_id.
for (samp_label in unique(obj$sample_id)) {
  samp_cells <- colnames(obj)[obj$sample_id == samp_label]
  if (length(samp_cells) == 0) next

  obj_samp <- obj[, samp_cells]
  md <- obj_samp@meta.data

  # Get spatial coordinates
  coords <- tryCatch({
    GetTissueCoordinates(obj_samp, which = "centroids")
  }, error = function(e) NULL)

  if (!is.null(coords)) {
    plot_df <- data.frame(
      x         = coords$x,
      y         = coords$y,
      cell_type = obj_samp$cell_type_xenium,
      n_counts  = obj_samp$nCount_Xenium
    )

    p_spatial_ct <- ggplot(plot_df, aes(x = x, y = y, colour = cell_type)) +
      geom_point(size = 0.3, alpha = 0.8) +
      coord_equal() +
      theme_void(base_size = 10) +
      labs(title = paste("Spatial cell types:", samp_label),
           colour = "Cell type") +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1),
                                   ncol = 2)) +
      theme(plot.background = element_rect(fill = "white", colour = NA))

    p_spatial_counts <- ggplot(plot_df, aes(x = x, y = y, colour = log1p(n_counts))) +
      geom_point(size = 0.3, alpha = 0.8) +
      scale_colour_viridis_c(option = "plasma", name = "log(counts+1)") +
      coord_equal() + theme_void(base_size = 10) +
      labs(title = paste("Transcript density:", samp_label)) +
      theme(plot.background = element_rect(fill = "white", colour = NA))

    save_pdf(p_spatial_ct + p_spatial_counts,
             paste0("07_spatial_", samp_label, ".pdf"),
             width = 20, height = 9)
    log(paste("Spatial plot saved:", samp_label))
  }
}

# save checkpoint
saveRDS(obj, file.path(OBJ_DIR, "xenium_annotated.rds"))
log("annotated checkpoint saved: xenium_annotated.rds")
log("next: 03_xenium_pbmc_linkage.R")