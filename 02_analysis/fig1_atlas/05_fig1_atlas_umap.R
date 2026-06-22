# elite reproduction script for figure 1
# goal: cell press (cell genomics) publication quality
# focus: high-vibrancy, bold points, clean aesthetics, vector-editable text

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
  library(dplyr)
  library(ggrepel)
})

# 1. Load Data
obj_path <- "02_scrna/03_integration/triad_integrated_nkt_patched.rds"
if (!file.exists(obj_path)) obj_path <- "02_scrna/03_integration/triad_integrated.rds"
obj <- readRDS(obj_path)

# 2. Elite Color Palettes (Cell Press style)
disease_colors <- c(
  "T2D"    = "#E31A1C", # Red
  "HTN"    = "#1F78B4", # Blue
  "Asthma" = "#33A02C"  # Green
)

celltype_colors <- c(
  "Naive CD4 T"   = "#1F77B4",
  "Memory CD4 T"  = "#AEC7E8",
  "CD8 T"         = "#FF7F0E",
  "NK"            = "#D62728",
  "B cell"        = "#2CA02C",
  "CD14 Monocyte" = "#9467BD",
  "CD16 Monocyte" = "#C5B0D5",
  "DC"            = "#8C564B",
  "Megakaryocyte" = "#E377C2",
  "Basophil"      = "#7F7F7F"
)

# 3. Master Cell Press Theme
theme_cell_press <- function() {
  theme_classic(base_size = 7) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7, color = "black"),
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.25, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# 4. Generate Panels with High Visibility
# pt.size = 1.2 and alpha = 1 for "dark", solid UMAPs
dpi_settings <- c(600, 600)
out_dir <- "manuscript_new/main_figures/fig1_integrated_atlas"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Panel A: Before/After
message("Generating Figure 1A...")
p1a_1 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = "disease", cols = disease_colors, 
                 pt.size = 1.2, raster = TRUE, raster.dpi = dpi_settings) + 
         theme_cell_press() + labs(title = "Before Integration", x = "UMAP 1", y = "UMAP 2")

p1a_2 <- DimPlot(obj, reduction = "umap", group.by = "disease", cols = disease_colors, 
                 pt.size = 1.2, raster = TRUE, raster.dpi = dpi_settings) + 
         theme_cell_press() + labs(title = "After Harmony Integration", x = "UMAP 1", y = "UMAP 2")

pdf(file.path(out_dir, "fig1a_before.pdf"), width = 3, height = 3, useDingbats = FALSE)
print(p1a_1)
dev.off()
pdf(file.path(out_dir, "fig1a_after.pdf"), width = 3, height = 3, useDingbats = FALSE)
print(p1a_2)
dev.off()

# Panel B: Atlas with Boxed Labels
message("Generating Figure 1B...")
cols_b <- celltype_colors[names(celltype_colors) %in% unique(obj$cell_type)]

p_panel_b <- DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cols_b, 
                     pt.size = 1.5, raster = TRUE, raster.dpi = dpi_settings, label = FALSE) + 
             theme_cell_press() + labs(title = "Integrated PBMC Atlas", x = "UMAP 1", y = "UMAP 2")

# Custom Boxed Labels (Professional quality)
label_df <- data.frame(UMAP1 = Embeddings(obj, "umap")[,1], UMAP2 = Embeddings(obj, "umap")[,2], label = obj$cell_type) %>%
            group_by(label) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

p_panel_b <- p_panel_b + 
             geom_label_repel(data = label_df, aes(x = UMAP1, y = UMAP2, label = label, fill = label), 
                              color = "white", fontface = "bold", size = 2, family = "Helvetica",
                              box.padding = 0.2, label.padding = 0.1, segment.size = 0.2, show.legend = FALSE) +
             scale_fill_manual(values = cols_b)

pdf(file.path(out_dir, "fig1b_atlas.pdf"), width = 5, height = 4.5, useDingbats = FALSE)
print(p_panel_b)
dev.off()

# Panel C: Split Condition
message("Generating Figure 1C...")
p_panel_c <- DimPlot(obj, reduction = "umap", group.by = "cell_type", split.by = "disease", 
                     cols = cols_b, pt.size = 1.0, raster = TRUE, raster.dpi = dpi_settings, label = FALSE) + 
             theme_cell_press() + theme(legend.position = "none") +
             labs(title = NULL, x = "UMAP 1", y = "UMAP 2")

pdf(file.path(out_dir, "fig1c_split.pdf"), width = 7, height = 2.5, useDingbats = FALSE)
print(p_panel_c)
dev.off()

# Panel D: Composition
message("Generating Figure 1D...")
comp_df <- obj@meta.data %>% group_by(disease, cell_type) %>% summarise(n = n(), .groups = "drop") %>%
           group_by(disease) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()

p_panel_d <- ggplot(comp_df, aes(x = disease, y = pct, fill = cell_type)) +
             geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
             scale_fill_manual(values = cols_b) +
             scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x, "%")) +
             theme_cell_press() + labs(title = "Cell Composition", x = NULL, y = "% of Cells", fill = "Cell Type")

pdf(file.path(out_dir, "fig1d_composition.pdf"), width = 3, height = 4, useDingbats = FALSE)
print(p_panel_d)
dev.off()

message("Figure 1 Elite Reproduction Complete.")
