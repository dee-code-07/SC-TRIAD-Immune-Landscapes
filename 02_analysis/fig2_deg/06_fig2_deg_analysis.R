# elite reproduction script for figure 2
# adherence: cell press / cell genomics standards
# consistency: matches figure 1 palette and theme

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrastr)
  library(ggrepel)
  library(tidyr)
})

# 1. Paths and Theme Setup
BASE_DIR   <- "02_scrna/04_deg"
TABLE_DIR  <- file.path(BASE_DIR, "tables")
OUT_DIR    <- "manuscript_new/main_figures/fig2_differential_gene_expression"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# CONSISTENT PALETTE FROM FIGURE 1
disease_colors <- c(
  "Asthma" = "#4DAF4A", # Green
  "T2D"    = "#E41A1C", # Red
  "HTN"    = "#377EB8"  # Blue
)

# Standard Cell Press Theme (Helvetica, 7pt)
theme_cell_press <- function() {
  theme_classic(base_size = 7) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7, color = "black"),
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 7),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )
}

# 2. Panel A: DEG Lollipop (Consistent with Fig 1 Palette)
message("Generating Figure 2A...")
summary_df <- read.csv(file.path(TABLE_DIR, "DEG_summary.csv"))

lollipop_df <- summary_df |>
  mutate(
    comparison_label = gsub("_vs_Control", "", comparison) |> gsub("Asthma_|T2D_", "", x = _) |> gsub("_", " ", x = _),
    disease = case_when(
      grepl("Asthma", comparison) ~ "Asthma", 
      grepl("T2D", comparison) ~ "T2D", 
      grepl("HTN", comparison) ~ "HTN", 
      TRUE ~ "Other"
    ),
    cell_type = factor(cell_type, levels = rev(c("NK", "CD8 T", "Memory CD4 T", "Naive CD4 T", "B cell", "CD14 Monocyte", "DC"))),
    n_down_neg = -n_down
  )

p_panel_a <- ggplot(lollipop_df) +
  geom_segment(aes(x = cell_type, xend = cell_type, y = 0, yend = n_up, color = disease), linewidth = 0.8) +
  geom_point(aes(x = cell_type, y = n_up, color = disease), size = 2) +
  geom_segment(aes(x = cell_type, xend = cell_type, y = 0, yend = n_down_neg, color = disease), linewidth = 0.8) +
  geom_point(aes(x = cell_type, y = n_down_neg, color = disease), size = 2, shape = 21, fill = "white", stroke = 1) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  scale_color_manual(values = disease_colors) +
  scale_y_continuous(labels = function(x) abs(x)) +
  coord_flip() +
  facet_wrap(~ comparison_label, scales = "free_x", ncol = 4) +
  theme_cell_press() + labs(title = "DEG Distribution", x = NULL, y = "Number of DEGs")

pdf(file.path(OUT_DIR, "fig2a_summary_lollipop.pdf"), width = 7, height = 3.5, useDingbats = FALSE)
print(p_panel_a)
dev.off()

# 3. Panel B: Asthma Correlation (Faceted, High Visibility)
message("Generating Figure 2B...")
cell_types <- c("NK", "CD8 T", "Memory CD4 T", "Naive CD4 T", "B cell", "CD14 Monocyte", "DC")
plot_data <- data.frame()

for (ct in cell_types) {
  mild <- file.path(TABLE_DIR, paste0("Asthma_Mild_vs_Control_", ct, "_DEGs.csv"))
  sev  <- file.path(TABLE_DIR, paste0("Asthma_Severe_vs_Control_", ct, "_DEGs.csv"))
  if (file.exists(mild) && file.exists(sev)) {
    d_m <- read.csv(mild) |> select(gene, log2FoldChange, padj) |> rename(lfc_m = log2FoldChange, p_m = padj)
    d_s <- read.csv(sev)  |> select(gene, log2FoldChange, padj) |> rename(lfc_s = log2FoldChange, p_s = padj)
    d_j <- inner_join(d_m, d_s, by = "gene") |>
      mutate(p_m = if_else(is.na(p_m), 1, p_m), p_s = if_else(is.na(p_s), 1, p_s), cell_type = ct,
             status = case_when(
               (p_m < 0.05 & abs(lfc_m) > 0.5) & (p_s < 0.05 & abs(lfc_s) > 0.5) & (sign(lfc_m) == sign(lfc_s)) ~ "Independent",
               (p_m < 0.05 & abs(lfc_m) > 0.5) ~ "Mild-Only",
               (p_s < 0.05 & abs(lfc_s) > 0.5) ~ "Severe-Only", 
               TRUE ~ "NS"
             ))
    plot_data <- bind_rows(plot_data, d_j)
  }
}

cor_stats <- plot_data |> group_by(cell_type) |> summarize(R = cor(lfc_m, lfc_s, use="complete.obs"), .groups="drop")

# Use Fig 1 derived significance palette
status_cols <- c("NS" = "grey90", "Mild-Only" = "#377EB8", "Severe-Only" = "#E41A1C", "Independent" = "#6A3D9A")

p_panel_b <- ggplot(plot_data, aes(x = lfc_m, y = lfc_s)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point_rast(aes(color = status), alpha = 1, size = 0.6, raster.dpi = 600) +
  scale_color_manual(values = status_cols) +
  geom_text(data = cor_stats, aes(x = -Inf, y = Inf, label = sprintf("R = %.2f", R)), 
            hjust = -0.2, vjust = 1.5, size = 2.5, fontface = "bold", family = "Helvetica") +
  facet_wrap(~ cell_type, ncol = 4) +
  theme_cell_press() + theme(legend.position = "none") +
  labs(title = "Asthma Severity Correlation", x = "Mild (log2FC)", y = "Severe (log2FC)")

pdf(file.path(OUT_DIR, "fig2b_correlation_faceted.pdf"), width = 7, height = 4.5, useDingbats = FALSE)
print(p_panel_b)
dev.off()

# 4. Panel C: BOLD VOLCANO PLOT (Cell Press Standard)
message("Generating Figure 2C...")
# Using Asthma Mild NK as representative
asthma_mild_nk_path <- file.path(TABLE_DIR, "Asthma_Mild_vs_Control_NK_DEGs.csv")
df_nk <- read.csv(asthma_mild_nk_path) |>
         filter(!is.na(padj)) |>
         mutate(
           sig = case_when(
             padj < 0.05 & log2FoldChange > 0.5 ~ "Up", 
             padj < 0.05 & log2FoldChange < -0.5 ~ "Down", 
             TRUE ~ "NS"
           ),
           label = if_else(sig != "NS" & rank(padj) <= 15, gene, NA_character_)
         )

p_panel_c <- ggplot(df_nk, aes(x = log2FoldChange, y = -log10(padj + 1e-300), color = sig)) +
             # alpha = 1 for DARK BOLD points
             geom_point_rast(alpha = 1, size = 1, raster.dpi = 600) +
             # Use Red for Up and Blue for Down (Consistent with Fig 1)
             scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey85")) +
             # Threshold lines
             geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey40", linewidth = 0.3) +
             geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40", linewidth = 0.3) +
             # Professional labels
             geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, 
                              color = "black", fontface = "italic", segment.size = 0.3, family = "Helvetica") +
             theme_cell_press() + theme(legend.position = "none") +
             labs(title = "NK | Asthma Mild vs Control", x = "log2 Fold Change", y = "-log10(adj. p-value)")

pdf(file.path(OUT_DIR, "fig2c_volcano_asthma_mild_nk.pdf"), width = 3.5, height = 3.5, useDingbats = FALSE)
print(p_panel_c)
dev.off()

message("Figure 2 Elite Reproduction Complete.")
