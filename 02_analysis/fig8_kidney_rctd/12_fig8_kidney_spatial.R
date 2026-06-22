suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(patchwork)
})

out_dir <- "/mnt/e/Documents/sctriad/manuscript_ft/Main/Figure8_Hypertensive_Kidney_Visium"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

theme_cell_press <- function() {
  theme_classic(base_size = 16) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.size = unit(0.5, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

status_colors <- c("Control" = "#2171B5", "HKD" = "#7A5C99")

# 8B: Stacked bar chart of mean proportions
df_b <- read.csv("/mnt/e/Documents/sctriad/04_spatial/kidney/tables/rctd_summary_by_status.csv")
df_b_long <- df_b %>% select(-n_spots) %>% pivot_longer(cols = -Status, names_to = "CellType", values_to = "Proportion")
# Pick top 15 cell types for legibility, group rest into 'Other'
top_cts <- df_b_long %>% group_by(CellType) %>% summarize(m = max(Proportion)) %>% top_n(15, m) %>% pull(CellType)
df_b_long <- df_b_long %>% mutate(CellType = ifelse(CellType %in% top_cts, CellType, "Other")) %>%
  group_by(Status, CellType) %>% summarize(Proportion = sum(Proportion), .groups="drop")

# Colors
comp_colors <- c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Set2"))[1:length(unique(df_b_long$CellType))]

p8b <- ggplot(df_b_long, aes(x = Status, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width=0.6) + # No border
  scale_fill_manual(values = comp_colors) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(title = "Mean Spot Composition", x = "Disease Status", y = "Proportion") +
  theme_cell_press()
ggsave(file.path(out_dir, "fig8b_stacked_composition.pdf"), p8b, width=7, height=8, useDingbats=FALSE)

# 8C: Volcano plot
df_c <- read.csv("/mnt/e/Documents/sctriad/04_spatial/kidney/tables/rctd_differential_celltypes.csv")
df_c <- df_c %>% mutate(
  sig = case_when(p_adj < 0.05 & log2FC > 0.5 ~ "Up", p_adj < 0.05 & log2FC < -0.5 ~ "Down", TRUE ~ "NS"),
  label = ifelse(sig != "NS", cell_type, NA)
)
p8c <- ggplot(df_c, aes(x = log2FC, y = -log10(p_adj + 1e-10), color = sig)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 6, box.padding = 0.5, max.overlaps = Inf, color="black", fontface="bold") +
  scale_color_manual(values = c("Up" = "#7A5C99", "Down" = "#2171B5", "NS" = "grey80")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "HKD vs Control Proportions", x = "log2(Fold Change)", y = "-log10(padj)") +
  theme_cell_press() + theme(legend.position="none")
ggsave(file.path(out_dir, "fig8c_volcano.pdf"), p8c, width=7, height=7, useDingbats=FALSE)

# 8D: Violin plots of immune cell proportions
df_d <- read.csv("/mnt/e/Documents/sctriad/04_spatial/kidney/tables/rctd_proportions_per_spot.csv")
immune_cells <- c("Mac", "NK", "CD4T", "CD8T", "B_memory", "B_Naive", "cDC", "pDC", "Neutrophil", "Plasma_Cells")
df_d_long <- df_d %>% select(Status, any_of(immune_cells)) %>% pivot_longer(cols = -Status, names_to = "CellType", values_to = "Proportion") %>%
  filter(Proportion > 0.001) # Keep spots with > 0.1% immune content as requested

p8d <- ggplot(df_d_long, aes(x = Status, y = Proportion, fill = Status)) +
  geom_violin(scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  facet_wrap(~CellType, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = status_colors) +
  labs(title = "Immune Enrichment in HKD", x = NULL, y = "Spot Proportion") +
  theme_cell_press() + theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(out_dir, "fig8d_immune_violins.pdf"), p8d, width=16, height=8, useDingbats=FALSE)

cat("Figure 8 Kidney panels perfectly upgraded.\n")
