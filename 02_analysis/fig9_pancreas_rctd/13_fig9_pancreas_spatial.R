suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(patchwork)
})

out_dir <- "/mnt/e/Documents/sctriad/manuscript_ft/Main/Figure9_T2D_Pancreas_Visium"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

theme_massive <- function() {
  theme_classic(base_size = 30) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.line = element_line(linewidth = 1, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.title = element_text(size = 32, face = "bold"),
      axis.text = element_text(size = 28, color = "black"),
      legend.title = element_text(size = 32, face = "bold"),
      legend.text = element_text(size = 28),
      legend.key.size = unit(1, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 32, face = "bold"),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(20, 20, 20, 20)
    )
}
status_colors <- c("Control" = "#A9A9A9", "T2D" = "#3D7AB5")

# Use the exact 10-color Xenium palette
CUSTOM_10 <- c("#4E79A7", "#1F4E79", "#76B7B2", "#E15759", "#59A14F", 
               "#F28E2B", "#B07AA1", "#9C755F", "#FF9DA7", "#BAB0AC")

# 9A: Stacked bar chart
df_a <- read.csv("/mnt/e/Documents/sctriad/04_spatial/pancreas/tables/rctd_summary_by_status.csv")
df_a_long <- df_a %>% select(-n_spots, -n_sections) %>% pivot_longer(cols = -status, names_to = "CellType", values_to = "Proportion")

top_cts <- df_a_long %>% group_by(CellType) %>% summarize(m = max(Proportion)) %>% top_n(9, m) %>% pull(CellType)
df_a_long <- df_a_long %>% mutate(CellType = ifelse(CellType %in% top_cts, CellType, "Other")) %>%
  group_by(status, CellType) %>% summarize(Proportion = sum(Proportion), .groups="drop")

p9a <- ggplot(df_a_long, aes(x = status, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width=0.6) +
  scale_fill_manual(values = CUSTOM_10) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "Disease Status", y = "Proportion") +
  theme_massive()
ggsave(file.path(out_dir, "fig9a_stacked_composition.pdf"), p9a, width=12, height=12, useDingbats=FALSE)

# 9B: Effect-size scatter plot
df_b <- read.csv("/mnt/e/Documents/sctriad/04_spatial/pancreas/tables/rctd_differential_celltypes.csv")
compartment_map <- c(alpha="Islet", beta="Islet", delta="Islet", gamma="Islet", epsilon="Islet",
                     acinar="Exocrine", ductal="Exocrine",
                     activated_stellate="Stromal", quiescent_stellate="Stromal", endothelial="Stromal",
                     macrophage="Immune", mast="Immune", t_cell="Immune")
df_b$Compartment <- compartment_map[df_b$cell_type]
df_b$Compartment[is.na(df_b$Compartment)] <- "Other"

p9b <- ggplot(df_b, aes(x = log2FC, y = rank_biserial_r, color = Compartment)) +
  geom_point(size = 6, alpha = 0.8) +
  geom_text_repel(aes(label = cell_type), size = 10, box.padding = 1.5, point.padding = 1, force = 5, min.segment.length = 0, segment.size = 1, max.overlaps = Inf, fontface="bold", show.legend = FALSE) +
  scale_color_brewer(palette = "Set1") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth=1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1.5) +
  labs(x = "log2(Fold Change)", y = "Rank-Biserial r") +
  theme_massive()
ggsave(file.path(out_dir, "fig9b_effect_size_scatter.pdf"), p9b, width=20, height=12, useDingbats=FALSE)

# 9C: Grouped bar chart
df_c <- read.csv("/mnt/e/Documents/sctriad/04_spatial/pancreas/tables/rctd_summary_by_section.csv")
islet_cells <- c("beta", "alpha", "delta", "gamma", "epsilon")
df_c_islet <- df_c %>% select(status, section, any_of(islet_cells)) %>%
  pivot_longer(cols = any_of(islet_cells), names_to = "CellType", values_to = "Proportion") %>%
  group_by(status, CellType) %>%
  summarize(mean_prop = mean(Proportion, na.rm=TRUE), se = sd(Proportion, na.rm=TRUE)/sqrt(n()), .groups="drop")

p9c <- ggplot(df_c_islet, aes(x = CellType, y = mean_prop, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se), width = 0.2, position = position_dodge(0.8), linewidth=1.5) +
  scale_fill_manual(values = status_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Cell Type", y = "Mean Proportion \u00B1 SE") +
  theme_massive()
ggsave(file.path(out_dir, "fig9c_islet_remodeling.pdf"), p9c, width=14, height=12, useDingbats=FALSE)

# 9D: Beta:Alpha ratio
df_d_spots <- read.csv("/mnt/e/Documents/sctriad/04_spatial/pancreas/tables/rctd_proportions_per_spot.csv")
df_d <- df_d_spots %>% select(status, alpha, beta) %>%
  filter(alpha > 0.01 | beta > 0.01) %>%
  mutate(Ratio = beta / (alpha + 1e-4))

p9d <- ggplot(df_d, aes(x = status, y = log2(Ratio + 1), fill = status)) +
  geom_violin(scale = "width", alpha = 1, color = "black", linewidth=1.5) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, linewidth=1.5) +
  scale_fill_manual(values = status_colors) +
  labs(x = NULL, y = "log2(Beta/Alpha + 1)") +
  theme_massive() + theme(legend.position="none")
ggsave(file.path(out_dir, "fig9d_beta_alpha_ratio.pdf"), p9d, width=12, height=12, useDingbats=FALSE)

# 9E: Grouped bar chart of immune cell
immune_cols <- c("macrophage", "t_cell")
df_e <- df_c %>% select(status, section, any_of(immune_cols)) %>%
  pivot_longer(cols = any_of(immune_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(status, CellType) %>%
  summarize(mean_prop = mean(Proportion, na.rm=TRUE), se = sd(Proportion, na.rm=TRUE)/sqrt(n()), .groups="drop")

p9e <- ggplot(df_e, aes(x = CellType, y = mean_prop, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se), width = 0.2, position = position_dodge(0.8), linewidth=1.5) +
  scale_fill_manual(values = status_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Cell Type", y = "Mean Proportion \u00B1 SE") +
  theme_massive()
ggsave(file.path(out_dir, "fig9e_immune_infiltration.pdf"), p9e, width=12, height=12, useDingbats=FALSE)

cat("Figure 9 Pancreas panels flawlessly generated.\n")
