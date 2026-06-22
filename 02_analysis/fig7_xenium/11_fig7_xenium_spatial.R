# ==============================================================================
# reproduce_fig7_xenium_pdf.R
# Purpose: Reproduce Fig 7, S12, S13 as native PDFs
# Improvements: Native PDF output, Fixed C/D panel mismatch
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(scales)
  library(viridis)
  library(pheatmap)
})

# Paths
BASE_DIR    <- "/mnt/e/Documents/sctriad"
OBJ_DIR     <- file.path(BASE_DIR, "04_spatial/asthma_lung/objects")
TAB_DIR     <- file.path(BASE_DIR, "04_spatial/asthma_lung/tables")
MAIN_OUT    <- file.path(BASE_DIR, "manuscript_new/main_figures/fig7_xenium")
SUPP_OUT    <- file.path(BASE_DIR, "manuscript_new/supplementary_figures/figs12_s13_xenium")

for (d in c(MAIN_OUT, SUPP_OUT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Constants
CT_COLORS <- c(
  "Epithelial_Unknown" = "#B2DF8A", "AT1" = "#33A02C", "AT2" = "#1F78B4",
  "Ciliated" = "#A6CEE3", "Club" = "#6A3D9A", "Goblet" = "#CAB2D6",
  "Fibroblast" = "#FF7F00", "Smooth_Muscle" = "#FDBF6F", "Endothelial" = "#E31A1C",
  "Pericyte" = "#FB9A99", "Macrophage" = "#B15928", "T_cell" = "#1F78B4",
  "B_cell" = "#A6CEE3", "Mast" = "#FF7F00", "NK" = "#E31A1C", "DC" = "#6A3D9A",
  "Unknown" = "#BDBDBD"
)
SAMPLE_COLORS <- c("0013717" = "#2166AC", "0013532" = "#D6604D")

# Load Data
cat("Loading Xenium final object...\n")
xen <- readRDS(file.path(OBJ_DIR, "xenium_final.rds"))

# --- FIGURE 7: Xenium Airway Atlas ---
cat("Generating Figure 7 Panels...\n")

# A. QC Summary
qc_sum <- xen@meta.data %>%
  group_by(sample_id) %>%
  summarise(n_cells = n(), .groups = "drop")

p7a1 <- ggplot(xen@meta.data, aes(sample_id, nCount_Xenium, fill=sample_id)) +
  geom_violin(scale="width", alpha=0.8) + geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values=SAMPLE_COLORS) + scale_y_log10() + theme_classic() + labs(title="A. QC Metrics", x=NULL)
p7a2 <- ggplot(qc_sum, aes(sample_id, n_cells, fill=sample_id)) +
  geom_bar(stat="identity", width=0.6) + scale_fill_manual(values=SAMPLE_COLORS) + theme_classic() + labs(y="Cell Count")

ggsave(file.path(MAIN_OUT, "fig7a_qc_summary.pdf"), p7a1 + p7a2, width = 10, height = 5)

# B. UMAP + Spatial
p7b1 <- DimPlot(xen, group.by="cell_type_xenium", cols=CT_COLORS, label=TRUE, pt.size=0.1) + theme_classic() + labs(title="B. UMAP")
ggsave(file.path(MAIN_OUT, "fig7b_umap_celltype.pdf"), p7b1, width = 8, height = 7)

# C. Marker Dotplot (WAS Panel D)
cat("Checking markers...\n")
DefaultAssay(xen) <- "Xenium"
all_markers <- rownames(xen)
if (!"SCT" %in% Assays(xen)) {
  cat("SCT assay missing, using Xenium assay.\n")
} else {
  DefaultAssay(xen) <- "SCT"
  all_markers <- unique(c(all_markers, rownames(xen)))
}

key_markers <- c("KRT5", "SCGB1A1", "MUC5AC", "FOXJ1", "SFTPB", "COL1A1", "ACTA2", "CD68", "CD3E", "MS4A1", "NKG7", "KLRB1")
available_markers <- key_markers[key_markers %in% all_markers]

p7c <- DotPlot(xen, features=available_markers, group.by="cell_type_xenium", cols=c("grey90", "red")) + 
       theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(title="C. Marker Gene Expression")
ggsave(file.path(MAIN_OUT, "fig7c_marker_dotplot.pdf"), p7c, width = 12, height = 6)

# D. Composition Bar (WAS Panel C)
comp_df <- xen@meta.data %>% group_by(sample_id, cell_type_xenium) %>% summarise(n=n(), .groups="drop") %>%
  group_by(sample_id) %>% mutate(pct = 100*n/sum(n))
p7d <- ggplot(comp_df, aes(sample_id, pct, fill=cell_type_xenium)) +
  geom_bar(stat="identity", color="white") + scale_fill_manual(values=CT_COLORS) + theme_classic() + labs(title="D. Cell Type Composition")
ggsave(file.path(MAIN_OUT, "fig7d_composition_bar.pdf"), p7d, width = 7, height = 6)

# --- FIGURE S12: Label Transfer ---
cat("Generating Figure S12 Panels...\n")
if ("pbmc_prediction_score_max" %in% colnames(xen@meta.data)) {
  pS12a <- DimPlot(xen, group.by="pbmc_predicted_highconf", label=TRUE, pt.size=0.1) + theme_classic() + labs(title="A. Predicted Cell Type")
  pS12b <- ggplot(xen@meta.data, aes(pbmc_prediction_score_max)) + geom_histogram(bins=50, fill="#2166AC") + theme_classic() + labs(title="B. Prediction Score Dist")
  ggsave(file.path(SUPP_OUT, "figs12_label_transfer.pdf"), wrap_plots(pS12a, pS12b, ncol=2), width = 14, height = 6)
}

# --- FIGURE S13: Spatial DEGs and AsthmaUP Score ---
cat("Generating Figure S13 Panels...\n")
if ("AsthmaUP_score" %in% colnames(xen@meta.data)) {
  pS13a <- FeaturePlot(xen, features="AsthmaUP_score", cols=c("grey90", "red"), pt.size=0.1)
  ggsave(file.path(SUPP_OUT, "figs13a_asthma_up_umap.pdf"), pS13a, width = 10, height = 8)
  
  pS13b <- ggplot(xen@meta.data, aes(x=cell_type_xenium, y=AsthmaUP_score, fill=cell_type_xenium)) +
           geom_violin(scale="width") + geom_boxplot(width=0.1, fill="white") +
           scale_fill_manual(values=CT_COLORS) + theme_classic() +
           theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
           labs(title="B. AsthmaUP by Cell Type", x=NULL)
  ggsave(file.path(SUPP_OUT, "figs13b_asthma_up_violin.pdf"), pS13b, width = 10, height = 6)
}

cat("Xenium PDF reproduction complete.\n")
