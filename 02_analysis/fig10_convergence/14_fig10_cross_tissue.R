suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

out_dir_10 <- "/mnt/e/Documents/sctriad/manuscript_ft/Main/Figure10_Cross_Tissue_Convergence"
out_dir_s5 <- "/mnt/e/Documents/sctriad/manuscript_ft/Supplementary/FigureS5_Asthmatic_Airway_PBMC_Label_Transfer"
dir.create(out_dir_10, recursive=TRUE, showWarnings=FALSE)
dir.create(out_dir_s5, recursive=TRUE, showWarnings=FALSE)

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

# -----------------------------------------------------------------------------
# FIGURE 10: Cross-tissue immune convergence
# -----------------------------------------------------------------------------
# 1. Kidney
kidney <- read.csv("/mnt/e/Documents/sctriad/04_spatial/kidney/tables/rctd_summary_by_status.csv") %>%
  filter(Status == "Disease") %>% 
  select(any_of(c("Mac", "CD14_Mono", "CD16_Mono", "cDC", "pDC", "NK", "CD4T", "CD8T", "B_memory", "B_Naive"))) %>%
  pivot_longer(cols=everything(), names_to="cell_type", values_to="prop") %>%
  mutate(Broad = case_when(
    grepl("Mac|Mono|DC", cell_type) ~ "Myeloid",
    grepl("NK|CD8", cell_type) ~ "Cytotoxic/NK",
    grepl("CD4|B_", cell_type) ~ "Adaptive",
    TRUE ~ "Other"
  )) %>% group_by(Broad) %>% summarize(prop = sum(prop)) %>% mutate(Tissue = "Kidney (HKD)")

# 2. Pancreas
pancreas <- read.csv("/mnt/e/Documents/sctriad/04_spatial/pancreas/tables/rctd_summary_by_status.csv") %>%
  filter(status == "T2D") %>%
  select(any_of(c("macrophage", "t_cell"))) %>%
  pivot_longer(cols=everything(), names_to="cell_type", values_to="prop") %>%
  mutate(Broad = case_when(
    grepl("mac", cell_type) ~ "Myeloid",
    grepl("t_cell", cell_type) ~ "Adaptive",
    TRUE ~ "Other"
  )) %>% group_by(Broad) %>% summarize(prop = sum(prop)) %>% mutate(Tissue = "Pancreas (T2D)")

# 3. Airway (Xenium)
airway <- read.csv("/mnt/e/Documents/sctriad/04_spatial/asthma_lung/tables/xenium_celltype_composition_by_sample.csv", colClasses = c("sample_id" = "character")) %>%
  filter(grepl("Macrophage|DC|NK|T_cell|B_cell", cell_type_xenium)) %>%
  group_by(cell_type_xenium) %>% summarize(pct = mean(pct)/100) %>% 
  mutate(Broad = case_when(
    grepl("Mac|DC", cell_type_xenium) ~ "Myeloid",
    grepl("NK|CD8", cell_type_xenium) ~ "Cytotoxic/NK", 
    grepl("T_cell|B_cell", cell_type_xenium) ~ "Adaptive", 
    TRUE ~ "Other"
  )) %>% group_by(Broad) %>% summarize(prop = sum(pct)) %>% mutate(Tissue = "Airway (Asthma)")

# 4. PBMC (Asthma) - Missing from previous script, adding it now!
pbmc <- read.csv("/mnt/e/Documents/sctriad/02_scrna/03_integration/reports/celltype_composition_by_disease.csv") %>%
  filter(disease == "Asthma") %>%
  mutate(prop = pct / 100) %>%
  filter(grepl("Monocyte|DC|Macrophage|NK|T|B cell", cell_type)) %>%
  mutate(Broad = case_when(
    grepl("Monocyte|DC|Macrophage", cell_type) ~ "Myeloid",
    grepl("NK|CD8", cell_type) ~ "Cytotoxic/NK",
    grepl("CD4|B cell", cell_type) ~ "Adaptive",
    TRUE ~ "Other"
  )) %>% group_by(Broad) %>% summarize(prop = sum(prop)) %>% mutate(Tissue = "PBMC (Asthma)")

df_10 <- bind_rows(kidney, pancreas, airway, pbmc)
df_10 <- df_10 %>% group_by(Tissue) %>% mutate(prop = prop / sum(prop)) 
df_10$Tissue <- factor(df_10$Tissue, levels = c("PBMC (Asthma)", "Airway (Asthma)", "Kidney (HKD)", "Pancreas (T2D)"))

# Custom beautiful distinct palette for the lineage types
CUSTOM_4 <- c("Myeloid" = "#E4572E", "Cytotoxic/NK" = "#3D7AB5", "Adaptive" = "#7A5C99", "Other" = "grey80")

p10 <- ggplot(df_10, aes(x = Tissue, y = prop, fill = Broad)) +
  geom_bar(stat="identity", position="fill", width=0.6, color="black", linewidth=1.5) +
  scale_fill_manual(values = CUSTOM_4) +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  labs(x = "Tissue Context", y = "Relative Immune Proportion", fill="Lineage") +
  theme_massive() + theme(axis.text.x = element_text(angle=30, hjust=1))

ggsave(file.path(out_dir_10, "fig10_convergence.pdf"), p10, width=16, height=12, useDingbats=FALSE)

# -----------------------------------------------------------------------------
# FIGURE S5: Label Transfer & Signatures
# -----------------------------------------------------------------------------
df_s5a <- read.csv("/mnt/e/Documents/sctriad/04_spatial/asthma_lung/tables/pbmc_label_transfer_results.csv")
if(nrow(df_s5a) > 20000) df_s5a <- df_s5a %>% sample_n(20000)
p_s5a <- ggplot(df_s5a, aes(x = reorder(xenium_cell_type, prediction_score, FUN=median), y = prediction_score, fill=xenium_cell_type)) +
  geom_violin(scale="width", alpha=1, color="black", linewidth=1.5) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=1.5) +
  coord_flip() + theme_massive() + theme(legend.position="none") +
  labs(x = "Xenium Cell Type", y = "Prediction Score")
ggsave(file.path(out_dir_s5, "figS5a_label_transfer.pdf"), p_s5a, width=16, height=14, useDingbats=FALSE)

df_s5b <- read.csv("/mnt/e/Documents/sctriad/04_spatial/asthma_lung/tables/signature_scores_by_celltype.csv")
p_s5b <- ggplot(df_s5b, aes(x = reorder(cell_type_xenium, AsthmaUP_score_mean), y = AsthmaUP_score_mean, fill=cell_type_xenium)) +
  geom_bar(stat="identity", width=0.7, color="black", linewidth=1.5) +
  geom_errorbar(aes(ymin=AsthmaUP_score_mean-AsthmaUP_score_sd, ymax=AsthmaUP_score_mean+AsthmaUP_score_sd), width=0.2, linewidth=1.5) +
  coord_flip() + theme_massive() + theme(legend.position="none") +
  labs(x = "Xenium Cell Type", y = "Mean Module Score")
ggsave(file.path(out_dir_s5, "figS5b_signature_scores.pdf"), p_s5b, width=16, height=14, useDingbats=FALSE)

cat("Figures 10 and S5 successfully recreated with massive typography and PBMC integration.\n")
