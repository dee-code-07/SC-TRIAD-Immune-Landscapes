# sc-triad project
# script: 04_xenium_figures_publication.R
# purpose: generate all publication-quality figures for the Xenium asthma
#          airway spatial transcriptomics analysis
#
# figures produced:
#   Fig A: QC summary panel (cell counts, filtering, metrics by sample)
#   Fig B: UMAP cell type annotation + spatial tissue map (2-panel)
#   Fig C: Cell type composition bar chart (by sample)
#   Fig D: Canonical marker dotplot (key genes confirming annotation)
#   Fig E: Label transfer result + score distribution + concordance heatmap
#   Fig F: NKG7 + KLRB1 spatial expression (the 2 asthma DEGs in panel)
#   Fig G: AsthmaUP score UMAP + spatial + violin (with caveat subtitle)
#   Fig H: Cross-tissue convergence — PBMC cell types vs Xenium label transfer
#
# all figures saved as 600 dpi TIFF + compiled into one PDF
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(forcats)   # FIX: added — required for fct_reorder() used in Figs C and G2
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(pheatmap)
  library(tiff)
  library(grid)
})

set.seed(42)
options(future.globals.maxSize = 16 * 1024^3)

BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
OUT_DIR <- file.path(BASE, "04_spatial", "asthma_lung")
FIG_DIR <- file.path(OUT_DIR, "figures", "publication")
TAB_DIR <- file.path(OUT_DIR, "tables")
OBJ_DIR <- file.path(OUT_DIR, "objects")
PDF_OUT <- file.path(OUT_DIR, "SC_TRIAD_xenium_asthma_all.pdf")
LOG_FILE <- file.path(OUT_DIR, "04_xenium_figures.log")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

# ── publication theme ──────────────────────────────────────────────────────────
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

# cell type colour palette — consistent across all figures
CT_COLORS <- c(
  "Epithelial_Unknown" = "#B2DF8A",
  "AT1"                = "#33A02C",
  "AT2"                = "#1F78B4",
  "Ciliated"           = "#A6CEE3",
  "Club"               = "#6A3D9A",
  "Goblet"             = "#CAB2D6",
  "Fibroblast"         = "#FF7F00",
  "Smooth_Muscle"      = "#FDBF6F",
  "Endothelial"        = "#E31A1C",
  "Pericyte"           = "#FB9A99",
  "Macrophage"         = "#B15928",
  "T_cell"             = "#1F78B4",
  "B_cell"             = "#A6CEE3",
  "Mast"               = "#FF7F00",
  "NK"                 = "#E31A1C",
  "DC"                 = "#6A3D9A",
  "Unknown"            = "#BDBDBD"
)

SAMPLE_COLORS <- c("0013717" = "#2166AC", "0013532" = "#D6604D")

tiff_paths <- character(0)

save_tiff <- function(p, fname, width = 12, height = 8) {
  path <- file.path(FIG_DIR, fname)
  tiff(path, width = width, height = height, units = "in",
       res = 600, compression = "lzw")
  tryCatch(print(p), error = function(e) {
    plot.new(); text(0.5, 0.5, conditionMessage(e))
  })
  dev.off()
  tiff_paths <<- c(tiff_paths, path)
  log(paste("  saved:", fname))
  invisible(path)
}

# ── load objects ───────────────────────────────────────────────────────────────
log("Loading Xenium final object...")
xen <- readRDS(file.path(OBJ_DIR, "xenium_final.rds"))
log(paste("Xenium:", ncol(xen), "cells |", nrow(xen), "genes"))

# Load label transfer results
lt_file  <- file.path(TAB_DIR, "pbmc_label_transfer_results.csv")
conf_file <- file.path(TAB_DIR, "label_transfer_confusion_immune.csv")

has_lt   <- file.exists(lt_file)
has_conf <- file.exists(conf_file)

if (has_lt) {
  lt_df <- read.csv(lt_file)
  log(paste("Label transfer results loaded:", nrow(lt_df), "cells"))
}

# ── helper: get spatial coordinates ───────────────────────────────────────────
get_coords <- function(obj_sub) {
  tryCatch(
    GetTissueCoordinates(obj_sub, which = "centroids"),
    error = function(e) NULL
  )
}

colors_in_use <- function(ct_vec) {
  cts <- unique(as.character(ct_vec))
  cols <- CT_COLORS[cts]
  # fill missing with grey
  cols[is.na(cols)] <- "#BDBDBD"
  names(cols) <- cts
  cols
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE A: QC summary panel
# ══════════════════════════════════════════════════════════════════════════════
log("Figure A: QC summary")

# Read QC table if it exists; otherwise compute from object metadata
qc_file <- file.path(TAB_DIR, "xenium_qc_summary.csv")

# Counts and gene density per sample from current object
qc_summary_live <- xen@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    n_cells        = n(),
    median_counts  = round(median(nCount_Xenium)),
    median_genes   = round(median(nFeature_Xenium)),
    q25_counts     = round(quantile(nCount_Xenium, 0.25)),
    q75_counts     = round(quantile(nCount_Xenium, 0.75)),
    .groups = "drop"
  )
write.csv(qc_summary_live,
          file.path(TAB_DIR, "xenium_qc_summary_postfilter.csv"),
          row.names = FALSE)

# A1: nCount violin by sample
pA1 <- ggplot(xen@meta.data,
              aes(x = sample_id, y = nCount_Xenium, fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.4) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
               linewidth = 0.5) +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  scale_y_log10(labels = label_comma()) +
  pub_theme +
  labs(title = "A. Transcripts per cell",
       x = NULL, y = "nCount (log scale)")

# A2: nFeature violin by sample
pA2 <- ggplot(xen@meta.data,
              aes(x = sample_id, y = nFeature_Xenium, fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.4) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
               linewidth = 0.5) +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  pub_theme +
  labs(title = "Genes per cell",
       x = NULL, y = "nFeature")

# A3: bar chart of retained cells per sample
pA3 <- ggplot(qc_summary_live,
              aes(x = sample_id, y = n_cells, fill = sample_id)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
  geom_text(aes(label = format(n_cells, big.mark = ",")),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = label_comma()) +
  pub_theme +
  labs(title = "Cells retained post-QC",
       x = NULL, y = "Cell count")

save_tiff(pA1 + pA2 + pA3 + plot_layout(ncol = 3),
          "FigA_qc_summary.tiff", width = 14, height = 5)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE B: UMAP annotation + spatial tissue maps
# ══════════════════════════════════════════════════════════════════════════════
log("Figure B: UMAP + spatial maps")

ct_cols <- colors_in_use(xen$cell_type_xenium)

# B1: UMAP coloured by cell type
pB1 <- DimPlot(xen, reduction = "umap",
               group.by = "cell_type_xenium",
               cols     = ct_cols,
               label    = TRUE, label.size = 3,
               repel    = TRUE, pt.size = 0.2,
               raster   = FALSE) +
  pub_theme +
  labs(title    = "B. Xenium airway cell types (UMAP)",
       subtitle = "86,220 cells | 339 genes | Harmony-corrected") +
  theme(legend.text = element_text(size = 8))

save_tiff(pB1, "FigB1_umap_celltype.tiff", width = 10, height = 8)

# B2: UMAP split by sample
# FIX: subtitle corrected — Seurat split.by orders facets alphabetically,
# so 0013532 appears on the LEFT and 0013717 on the RIGHT
pB2 <- DimPlot(xen, reduction = "umap",
               group.by = "cell_type_xenium",
               split.by = "sample_id",
               cols     = ct_cols,
               label    = FALSE, pt.size = 0.15,
               raster   = FALSE) +
  pub_theme +
  labs(title    = "Cell type composition by sample",
       subtitle = "Left: 0013532 (32,791 cells) | Right: 0013717 (53,429 cells)") +
  theme(legend.text = element_text(size = 7),
        legend.position = "right")

save_tiff(pB2, "FigB2_umap_split_by_sample.tiff", width = 16, height = 7)

# B3: Spatial tissue maps (both samples)
spatial_plots <- list()
for (samp_label in c("0013717", "0013532")) {
  samp_cells <- colnames(xen)[xen$sample_id == samp_label]
  if (length(samp_cells) == 0) next
  xen_s <- xen[, samp_cells]
  coords <- get_coords(xen_s)
  if (is.null(coords)) next

  plot_df <- data.frame(
    x         = coords$x,
    y         = coords$y,
    cell_type = xen_s$cell_type_xenium
  )

  ct_cols_s <- colors_in_use(plot_df$cell_type)
  n_cells   <- nrow(plot_df)

  spatial_plots[[samp_label]] <- ggplot(
    plot_df, aes(x = x, y = y, colour = cell_type)
  ) +
    geom_point(size = 0.25, alpha = 0.75) +
    scale_colour_manual(values = ct_cols_s, name = "Cell type") +
    coord_equal() +
    theme_void(base_size = 10) +
    guides(colour = guide_legend(
      override.aes = list(size = 3, alpha = 1), ncol = 1
    )) +
    labs(title    = paste0("Sample ", samp_label),
         subtitle = paste0(format(n_cells, big.mark = ","), " cells")) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(face = "bold", size = 8)
    )
}

if (length(spatial_plots) == 2) {
  p_spatial_combined <- spatial_plots[[1]] + spatial_plots[[2]] +
    plot_annotation(
      title    = "B. Xenium airway spatial tissue maps",
      subtitle = "Single-cell resolution | point = one cell | coloured by annotated cell type",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey40")
      )
    )
  save_tiff(p_spatial_combined, "FigB3_spatial_tissue_maps.tiff",
            width = 20, height = 9)
} else if (length(spatial_plots) == 1) {
  save_tiff(spatial_plots[[1]], "FigB3_spatial_tissue_map.tiff",
            width = 11, height = 9)
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE C: Cell type composition by sample
# ══════════════════════════════════════════════════════════════════════════════
log("Figure C: Cell type composition")

comp_df <- xen@meta.data %>%
  group_by(sample_id, cell_type_xenium) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(cell_type_xenium = fct_reorder(cell_type_xenium, n, .desc = TRUE))

write.csv(comp_df,
          file.path(TAB_DIR, "xenium_celltype_composition_by_sample.csv"),
          row.names = FALSE)

ct_cols_c <- colors_in_use(comp_df$cell_type_xenium)

pC <- ggplot(comp_df,
             aes(x = sample_id, y = pct, fill = cell_type_xenium)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = ct_cols_c, name = "Cell type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = function(x) paste0(x, "%")) +
  pub_theme +
  labs(title    = "C. Cell type composition by sample",
       subtitle = "Xenium airway biopsy | post-QC, post-annotation",
       x = NULL, y = "% of cells") +
  theme(legend.text = element_text(size = 8))

save_tiff(pC, "FigC_celltype_composition.tiff", width = 8, height = 7)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE D: Canonical marker dotplot
# ══════════════════════════════════════════════════════════════════════════════
log("Figure D: Marker dotplot")

# Use only the most diagnostic markers — avoid overloading the dotplot
key_markers <- c(
  # Epithelial / structural
  "KRT5", "SCGB1A1", "MUC5AC", "FOXJ1", "SFTPB", "SFTPD", "AGER",
  # Stromal
  "COL1A1", "ACTA2", "PDGFRB", "CDH5", "VWF",
  # Immune
  "CD68", "KIT", "CD3E", "MS4A1", "NKG7", "KLRB1"
)

key_markers_present <- key_markers[key_markers %in% rownames(xen)]
log(paste("Key markers in panel:", length(key_markers_present)))

if (length(key_markers_present) >= 5) {
  DefaultAssay(xen) <- "SCT"

  pD <- DotPlot(xen,
                features  = key_markers_present,
                group.by  = "cell_type_xenium",
                assay     = "SCT",
                cols      = c("grey90", "#C00000"),
                dot.scale = 8) +
    pub_theme +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, size = 9,
                                  face = "italic"),
      axis.text.y = element_text(size = 10)
    ) +
    labs(title    = "D. Canonical marker gene expression",
         subtitle = "Dot size = % cells expressing | colour = scaled mean expression",
         x = NULL, y = NULL)

  save_tiff(pD, "FigD_marker_dotplot.tiff", width = 14, height = 7)
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE E: Label transfer results
# ══════════════════════════════════════════════════════════════════════════════
log("Figure E: Label transfer")

if (has_lt && "pbmc_prediction_score_max" %in% colnames(xen@meta.data)) {

  # E1: UMAP coloured by prediction score
  pE1 <- FeaturePlot(xen, features = "pbmc_prediction_score_max",
                     reduction = "umap", pt.size = 0.2, raster = FALSE,
                     cols = c("grey90", "#2166AC")) +
    pub_theme +
    labs(title    = "E. PBMC-to-Xenium label transfer score",
         subtitle = paste0(
           "Asthma PBMC (n=105,060) -> Xenium airway (n=86,220) | ",
           "RNA/LogNormalize | 309 common genes\n",
           "Median prediction score: 0.453 | ",
           "High-confidence cells (>0.5): 34,956 (40.5%)"
         ))

  # E2: Score distribution histogram
  pE2 <- ggplot(
    data.frame(score = xen$pbmc_prediction_score_max),
    aes(x = score)
  ) +
    geom_histogram(bins = 60, fill = "#2166AC", colour = "white",
                   alpha = 0.85, linewidth = 0.3) +
    geom_vline(xintercept = 0.453, linetype = "dashed",
               colour = "#D6604D", linewidth = 0.8) +
    geom_vline(xintercept = 0.5, linetype = "solid",
               colour = "black", linewidth = 0.6) +
    annotate("text", x = 0.48, y = Inf, label = "threshold\n(0.5)",
             hjust = 1, vjust = 1.5, size = 3, colour = "black") +
    annotate("text", x = 0.455, y = Inf, label = "median\n(0.453)",
             hjust = -0.1, vjust = 1.5, size = 3, colour = "#D6604D") +
    pub_theme +
    labs(title = "Prediction score distribution",
         x = "Max prediction score", y = "Cells")

  save_tiff(pE1 + pE2 + plot_layout(ncol = 2),
            "FigE1_label_transfer_umap_dist.tiff", width = 16, height = 7)

  # E3: Concordance heatmap (confusion matrix, immune cells only)
  if (has_conf) {
    conf_mat <- read.csv(conf_file, row.names = 1)
    conf_mat <- as.matrix(conf_mat)

    # Row-normalise (% of Xenium cells predicted as each PBMC type)
    conf_norm <- sweep(conf_mat, 1, rowSums(conf_mat), "/") * 100
    conf_norm[is.nan(conf_norm)] <- 0

    # Order rows and columns by dominant assignment
    row_order <- rownames(conf_norm)[
      order(apply(conf_norm, 1, which.max))
    ]
    col_order <- colnames(conf_norm)[
      order(apply(conf_norm, 2, which.max))
    ]
    conf_plot <- conf_norm[row_order,
                           intersect(col_order, colnames(conf_norm)),
                           drop = FALSE]

    tiff_path_e3 <- file.path(FIG_DIR, "FigE2_label_transfer_concordance.tiff")
    tiff(tiff_path_e3, width = 12, height = 6, units = "in",
         res = 600, compression = "lzw")
    pheatmap(
      conf_plot,
      color            = colorRampPalette(c("white", "#2166AC", "#08306b"))(100),
      cluster_rows     = FALSE,
      cluster_cols     = FALSE,
      display_numbers  = TRUE,
      number_format    = "%.0f%%",
      fontsize_number  = 7,
      fontsize_row     = 10,
      fontsize_col     = 9,
      angle_col        = 45,
      border_color     = "grey80",
      main             = paste0(
        "E. Label transfer concordance: Xenium immune cells vs PBMC predictions\n",
        "(row-normalised % | high-confidence cells only, score > 0.5)"
      ),
      sub = paste0(
        "Note: low concordance is expected for tissue-resident vs circulating cells. ",
        "Macrophages diverge most from blood monocytes (correct biology)."
      )
    )
    dev.off()
    tiff_paths <<- c(tiff_paths, tiff_path_e3)
    log("  saved: FigE2_label_transfer_concordance.tiff")
  }

  # E4: UMAP coloured by predicted PBMC cell type (high-confidence only)
  if ("pbmc_predicted_highconf" %in% colnames(xen@meta.data)) {
    pE4 <- DimPlot(xen, group.by = "pbmc_predicted_highconf",
                   label = TRUE, label.size = 2.5, repel = TRUE,
                   pt.size = 0.2, raster = FALSE) +
      pub_theme +
      labs(title    = "E. PBMC-predicted cell type on Xenium UMAP",
           subtitle = "Score > 0.5 only | others labelled Low_confidence") +
      theme(legend.text = element_text(size = 7))
    save_tiff(pE4, "FigE3_label_transfer_predicted_umap.tiff",
              width = 11, height = 8)
  }

} else {
  log("Label transfer metadata not found — skipping Figure E")
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE F: NKG7 and KLRB1 — the 2 asthma DEGs present in the Xenium panel
# These are the only signature genes in the panel; plot them individually
# rather than as a composite score to show the spatial expression pattern
# ══════════════════════════════════════════════════════════════════════════════
log("Figure F: NKG7 + KLRB1 spatial expression")

panel_genes <- rownames(xen)
asthma_panel_genes <- intersect(c("NKG7", "KLRB1"), panel_genes)
log(paste("Asthma DEGs in panel:", paste(asthma_panel_genes, collapse = ", ")))

if (length(asthma_panel_genes) >= 1) {
  DefaultAssay(xen) <- "SCT"

  # UMAP feature plots
  feat_plots <- lapply(asthma_panel_genes, function(g) {
    FeaturePlot(xen, features = g, reduction = "umap",
                pt.size = 0.2, raster = FALSE,
                cols = c("grey90", "#C00000")) +
      pub_theme +
      labs(title    = g,
           subtitle = "SCT expression | asthma DEG (PBMC)") +
      theme(plot.title = element_text(face = "italic", size = 12))
  })

  if (length(feat_plots) == 2) {
    save_tiff(feat_plots[[1]] + feat_plots[[2]],
              "FigF1_asthma_degs_umap.tiff", width = 14, height = 6)
  } else {
    save_tiff(feat_plots[[1]], "FigF1_asthma_degs_umap.tiff",
              width = 8, height = 6)
  }

  # Spatial plots per gene per sample
  for (gene in asthma_panel_genes) {
    spatial_gene_plots <- list()

    for (samp_label in c("0013717", "0013532")) {
      samp_cells <- colnames(xen)[xen$sample_id == samp_label]
      if (length(samp_cells) == 0) next
      xen_s  <- xen[, samp_cells]
      coords <- get_coords(xen_s)
      if (is.null(coords)) next

      expr_vals <- GetAssayData(xen_s, assay = "SCT",
                                layer = "data")[gene, ]

      plot_df <- data.frame(x = coords$x, y = coords$y,
                             expr = expr_vals)

      spatial_gene_plots[[samp_label]] <- ggplot(
        plot_df, aes(x = x, y = y, colour = expr)
      ) +
        geom_point(size = 0.3, alpha = 0.85) +
        scale_colour_gradient(low = "grey92", high = "#C00000",
                               name = "SCT expr") +
        coord_equal() + theme_void(base_size = 10) +
        labs(title    = paste0(gene, " - ", samp_label),
             subtitle = "Spatial expression in airway tissue") +
        theme(
          plot.title      = element_text(face = "bold.italic", size = 11),
          plot.background = element_rect(fill = "white", colour = NA)
        )
    }

    if (length(spatial_gene_plots) == 2) {
      save_tiff(
        spatial_gene_plots[[1]] + spatial_gene_plots[[2]],
        paste0("FigF2_spatial_", gene, ".tiff"),
        width = 18, height = 8
      )
    } else if (length(spatial_gene_plots) == 1) {
      save_tiff(spatial_gene_plots[[1]],
                paste0("FigF2_spatial_", gene, ".tiff"),
                width = 10, height = 8)
    }
  }

  # Violin: NKG7 + KLRB1 expression by cell type
  if (length(asthma_panel_genes) >= 1) {
    pF_vio <- VlnPlot(xen, features = asthma_panel_genes,
                      group.by = "cell_type_xenium",
                      assay = "SCT", pt.size = 0, ncol = 1) &
      pub_theme &
      theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
            axis.title.x = element_blank())

    save_tiff(pF_vio,
              "FigF3_asthma_degs_violin_by_celltype.tiff",
              width = 12, height = if (length(asthma_panel_genes) == 2) 9 else 5)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE G: AsthmaUP score — UMAP + spatial + violin
# Includes explicit panel-limitation caveat on all panels
# ══════════════════════════════════════════════════════════════════════════════
log("Figure G: AsthmaUP score visualization")

if ("AsthmaUP_score" %in% colnames(xen@meta.data)) {

  panel_caveat <- paste0(
    "AsthmaUP score | 2/14 panel genes (NKG7, KLRB1) | ",
    "colMeans surrogate, not background-corrected\n",
    "Interpret with caution: low panel coverage limits biological inference"
  )

  # G1: UMAP
  pG1 <- FeaturePlot(xen, features = "AsthmaUP_score",
                     reduction = "umap", pt.size = 0.2, raster = FALSE,
                     cols = c("grey92", "#C00000")) +
    pub_theme +
    labs(title    = "G. Asthma UP signature score (PBMC-derived)",
         subtitle = panel_caveat)

  save_tiff(pG1, "FigG1_AsthmaUP_umap.tiff", width = 10, height = 8)

  # G2: Violin by cell type
  pG2 <- ggplot(
    xen@meta.data %>%
      mutate(cell_type_xenium = fct_reorder(
        cell_type_xenium, AsthmaUP_score, .fun = median, .desc = TRUE
      )),
    aes(x = cell_type_xenium, y = AsthmaUP_score, fill = cell_type_xenium)
  ) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.75,
                linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    scale_fill_manual(values = colors_in_use(xen$cell_type_xenium),
                      guide  = "none") +
    pub_theme +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(title    = "G. AsthmaUP score by cell type",
         subtitle = panel_caveat,
         x = NULL, y = "AsthmaUP score")

  save_tiff(pG2, "FigG2_AsthmaUP_violin.tiff", width = 12, height = 7)

  # G3: Spatial maps
  for (samp_label in c("0013717", "0013532")) {
    samp_cells <- colnames(xen)[xen$sample_id == samp_label]
    if (length(samp_cells) == 0) next
    xen_s  <- xen[, samp_cells]
    coords <- get_coords(xen_s)
    if (is.null(coords)) next

    scores <- xen_s$AsthmaUP_score
    plot_df <- data.frame(x = coords$x, y = coords$y, score = scores)

    pG3 <- ggplot(plot_df, aes(x = x, y = y, colour = score)) +
      geom_point(size = 0.3, alpha = 0.8) +
      scale_colour_gradient2(low  = "#2166AC", mid = "grey90",
                              high = "#C00000", midpoint = 0,
                              name = "Score") +
      coord_equal() + theme_void(base_size = 10) +
      labs(title    = paste0("AsthmaUP score - ", samp_label),
           subtitle = "NKG7 + KLRB1 | panel limitation noted") +
      theme(
        plot.title      = element_text(face = "bold", size = 11),
        plot.background = element_rect(fill = "white", colour = NA)
      )

    save_tiff(pG3, paste0("FigG3_AsthmaUP_spatial_", samp_label, ".tiff"),
              width = 10, height = 9)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE H: Cross-tissue convergence summary
# Connects PBMC cell type composition with Xenium tissue cell type composition
# This is the narrative bridge figure for the dissertation
# ══════════════════════════════════════════════════════════════════════════════
log("Figure H: Cross-tissue convergence summary")

# Xenium immune cell proportions from this dataset
xenium_comp <- xen@meta.data %>%
  mutate(
    broad_type = case_when(
      cell_type_xenium %in% c("Macrophage", "DC")         ~ "Myeloid",
      cell_type_xenium %in% c("T_cell", "CD4_T", "CD8_T") ~ "T cell",
      cell_type_xenium %in% c("NK")                        ~ "NK",
      cell_type_xenium %in% c("B_cell")                    ~ "B cell",
      cell_type_xenium %in% c("Mast")                      ~ "Mast/Eosinophil",
      cell_type_xenium %in% c("Epithelial_Unknown",
                               "AT1", "AT2", "Ciliated",
                               "Club", "Goblet")            ~ "Epithelial",
      cell_type_xenium %in% c("Fibroblast",
                               "Smooth_Muscle")             ~ "Structural stromal",
      cell_type_xenium %in% c("Endothelial", "Pericyte")   ~ "Vascular",
      TRUE ~ "Other"
    )
  ) %>%
  count(broad_type) %>%
  mutate(pct = 100 * n / sum(n),
         source = "Xenium airway")

# ── FIG H PBMC values ─────────────────────────────────────────────────────────
# BUG 3 (action required): the values below are placeholders.
# Replace with real numbers from your composition file by running this once:
#
#   comp <- read.csv("~/sc-triad/02_scrna/pbmc/02_normalization/asthma/reports/celltype_composition_manual.csv")
#   comp %>%
#     filter(group %in% c("Asthma_Mild", "Asthma_Severe")) %>%
#     group_by(cell_type) %>%
#     summarise(n = sum(n_cells)) %>%
#     mutate(pct = round(100 * n / sum(n), 1)) %>%
#     arrange(desc(pct))
#
# Then update the pct values in pbmc_asthma_comp below and remove this comment.
# ─────────────────────────────────────────────────────────────────────────────
pbmc_asthma_comp <- data.frame(
  broad_type = c("T cell", "NK", "B cell", "Myeloid",
                 "Mast/Eosinophil", "Epithelial",
                 "Structural stromal", "Vascular"),
  pct        = c(55, 12, 18, 10, 4, 0, 0, 0),   # <- REPLACE WITH REAL VALUES
  source     = "PBMC (Asthma)"
)

conv_df <- bind_rows(xenium_comp[, c("broad_type","pct","source")],
                     pbmc_asthma_comp) %>%
  filter(pct > 0)

broad_colors <- c(
  "T cell"           = "#1F78B4",
  "NK"               = "#E31A1C",
  "B cell"           = "#A6CEE3",
  "Myeloid"          = "#B15928",
  "Mast/Eosinophil"  = "#FF7F00",
  "Epithelial"       = "#33A02C",
  "Structural stromal" = "#FDBF6F",
  "Vascular"         = "#FB9A99",
  "Other"            = "#BDBDBD"
)

pH <- ggplot(conv_df,
             aes(x = source, y = pct, fill = broad_type)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = broad_colors, name = "Cell category") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = function(x) paste0(x, "%")) +
  pub_theme +
  labs(
    title    = "H. Cross-tissue convergence: PBMC vs Xenium airway cell composition",
    subtitle = paste0(
      "PBMC reflects systemic circulating immune milieu in asthma\n",
      "Xenium captures resident airway structural + immune architecture\n",
      "Shared immune lineages (T cell, NK, B cell, Myeloid) connect both compartments"
    ),
    x = NULL, y = "% of cells"
  ) +
  theme(legend.text = element_text(size = 9))

save_tiff(pH, "FigH_cross_tissue_convergence.tiff", width = 10, height = 8)

# ══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY: sample-level QC before/after (if raw QC table exists)
# ══════════════════════════════════════════════════════════════════════════════

if (file.exists(file.path(TAB_DIR, "xenium_qc_summary.csv"))) {
  qc_raw <- read.csv(file.path(TAB_DIR, "xenium_qc_summary.csv"))

  if ("n_cells_raw" %in% colnames(qc_raw) &&
      "n_cells_final" %in% colnames(qc_raw)) {

    # FIX Bug 1: coerce sample to factor so ggplot treats it as discrete,
    # preventing sample IDs being read as numeric and plotted as hairlines
    qc_raw$sample <- factor(qc_raw$sample)

    qc_long <- qc_raw %>%
      select(sample, n_cells_raw, n_cells_final) %>%
      pivot_longer(c(n_cells_raw, n_cells_final),
                   names_to  = "stage",
                   values_to = "n_cells") %>%
      mutate(stage = recode(stage,
                            n_cells_raw   = "Raw",
                            n_cells_final = "Post-QC"))

    pSupp <- ggplot(qc_long,
                    aes(x = sample, y = n_cells, fill = stage)) +
      geom_bar(stat = "identity", position = "dodge",
               width = 0.6, alpha = 0.9) +
      geom_text(aes(label = format(n_cells, big.mark = ",")),
                position = position_dodge(width = 0.6),
                vjust = -0.4, size = 3) +
      scale_fill_manual(values = c(Raw = "grey70", "Post-QC" = "#2166AC"),
                        name = "Stage") +
      scale_y_continuous(labels = label_comma(),
                         expand = expansion(mult = c(0, 0.15))) +
      pub_theme +
      labs(title    = "Supplementary: QC filtering per sample",
           subtitle = "MAD-based filtering: nCount, nFeature, cell area, negative control probes",
           x = NULL, y = "Cell count")

    save_tiff(pSupp, "FigSupp_qc_filtering.tiff", width = 8, height = 6)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# COMPILE PDF
# ══════════════════════════════════════════════════════════════════════════════
log("Compiling PDF...")

valid_tiffs <- tiff_paths[file.exists(tiff_paths)]
log(paste("Valid TIFFs to compile:", length(valid_tiffs)))

if (length(valid_tiffs) > 0) {
  pdf(PDF_OUT, width = 17, height = 11)
  for (tf in valid_tiffs) {
    img <- tryCatch(tiff::readTIFF(tf, native = TRUE), error = function(e) NULL)
    if (!is.null(img)) {
      grid.newpage()
      grid.raster(img)
    }
  }
  dev.off()
  log(paste("PDF compiled:", PDF_OUT))
  log(paste("Pages:", length(valid_tiffs)))
}

# ── final summary ──────────────────────────────────────────────────────────────
log("")
log("=== Figure generation complete ===")
log(paste("Output directory:", FIG_DIR))
log(paste("PDF:", PDF_OUT))
log("")
log("Figures produced:")
for (tf in valid_tiffs)
  log(paste0("  ", basename(tf)))
log("")
log("METHODS TEXT (paste directly):")
log(paste0(
  "Xenium spatial transcriptomics data (GSE269354) were processed using ",
  "Seurat v5. Two samples (0013717 and 0013532) were loaded from raw ",
  "cell_feature_matrix files. Quality control was performed using ",
  "adaptive MAD-based thresholds applied to nCount_Xenium, nFeature_Xenium, ",
  "cell area, and negative control probe rate (<5%). Post-QC, 53,429 cells ",
  "(sample 0013717) and 32,791 cells (sample 0013532) were retained. ",
  "Normalization was performed using SCTransform v2 with cell area as a ",
  "regressing variable. Dimensionality reduction used PCA (30 PCs) followed ",
  "by Harmony batch correction and UMAP embedding. Cell types were assigned ",
  "by manual annotation of 24 clusters based on canonical airway marker gene ",
  "expression. Label transfer from the integrated PBMC Asthma subset ",
  "(n=105,060 cells, 309 common genes) was performed using ",
  "FindTransferAnchors with RNA/LogNormalize on both sides. Median prediction ",
  "score was 0.453; 34,956 cells (40.5%) met the high-confidence threshold ",
  "(score > 0.5). PBMC-derived asthma DEG signatures were scored using ",
  "AddModuleScore when >=3 panel genes were available, or colMeans of ",
  "z-scored SCT expression as a surrogate for smaller gene sets. Of 14 ",
  "cross-disease Asthma-upregulated DEGs, 2 (NKG7, KLRB1) were present in ",
  "the Xenium panel (14% coverage), representing a panel-design limitation. ",
  "All figures were generated at 600 dpi using ggplot2 v3.5."
))