# ==============================================================================
# reproduce_fig3_S4_pathway.R
# Purpose: Reproduce Pathway Enrichment figures (Fig 3, S4a) and Table S5
# Input: pathway_enrichment_per_comparison.csv
# Output: Native PDFs in manuscript_new/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(scales)
  library(openxlsx)
})

# Paths
BASE_DIR     <- "/mnt/e/Documents/sctriad"
INPUT_CSV    <- file.path(BASE_DIR, "02_scrna/05_pathway/tables/pathway_enrichment_per_comparison.csv")
OUT_FIG3     <- file.path(BASE_DIR, "manuscript_new/main_figures/fig3_pathway")
OUT_SUPP     <- file.path(BASE_DIR, "manuscript_new/supplementary_figures/figs4_pathway")
OUT_TAB      <- file.path(BASE_DIR, "manuscript_new/supplementary_tables")

for (d in c(OUT_FIG3, OUT_SUPP, OUT_TAB))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Load Data
cat("Loading pathway enrichment data...\n")
df_all <- read.csv(INPUT_CSV)

# Helper: Parse GeneRatio
parse_gene_ratio <- function(r) {
  sapply(r, function(x) {
    p <- as.numeric(strsplit(x, "/")[[1]])
    if (length(p) == 2) p[1]/p[2] else NA
  })
}

# Helper: Wrap labels
wrap_label <- function(x, w = 50) {
  sapply(x, function(s) paste(strwrap(s, width = w), collapse = "\n"))
}

# Helper: Make Dotplot
make_dotplot <- function(df, title, direction = "up") {
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>%
    mutate(GeneRatio_num = parse_gene_ratio(GeneRatio)) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>% # Top 10 for publication clarity
    mutate(Description = fct_reorder(wrap_label(Description), GeneRatio_num))
  
  col_low  <- if(direction == "up") "#fcbf7a" else "#9ecae1"
  col_high <- if(direction == "up") "#a50026" else "#08306b"
  
  ggplot(df, aes(x = GeneRatio_num, y = Description, size = Count, color = -log10(p.adjust))) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = col_low, high = col_high, name = expression(-log[10](p[adj]))) +
    scale_size_continuous(range = c(3, 8)) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(lineheight = 0.9))
}

# --- 1. Figure 3 Panel A (Asthma T/NK Dot Plots) ---
cat("Generating Fig 3 Panel A (Individual Dot Plots)...\n")
# Comparisons requested by user
target_comps <- c(
  "Asthma_Mild_vs_Control_CD14 Monocyte_DEGs_up",
  "Asthma_Mild_vs_Control_CD8 T_DEGs_up",
  "Asthma_Mild_vs_Control_Memory CD4 T_DEGs_up",
  "Asthma_Mild_vs_Control_NK_DEGs_up",
  "Asthma_Mild_vs_Control_Naive CD4 T_DEGs_up",
  "Asthma_Severe_vs_Control_CD8 T_DEGs_up",
  "Asthma_Severe_vs_Control_CD8 T_DEGs_down",
  "Asthma_Severe_vs_Control_Memory CD4 T_DEGs_up",
  "Asthma_Severe_vs_Control_NK_DEGs_up",
  "Asthma_Severe_vs_Control_Naive CD4 T_DEGs_up",
  "Asthma_severity_independent_CD8 T_DEGs_up",
  "Asthma_severity_independent_Memory CD4 T_DEGs_up",
  "Asthma_severity_independent_Memory CD4 T_DEGs_down",
  "Asthma_severity_independent_Naive CD4 T_DEGs_up",
  "HTN_vs_Control_B cell_DEGs_up",
  "HTN_vs_Control_NK_DEGs_up"
)

# Note: The CSV has spaces in cell type names (e.g., "CD14 Monocyte" vs "CD14_Monocyte")
# Standardize based on what's in the CSV
available_comps <- unique(df_all$comparison_direction)
matched_comps <- target_comps[target_comps %in% available_comps]
# Try underscore version if space version not found
if (length(matched_comps) < length(target_comps)) {
  target_comps_underscore <- gsub(" ", "_", target_comps)
  matched_comps <- unique(c(matched_comps, available_comps[available_comps %in% target_comps_underscore]))
}

pdf(file.path(OUT_FIG3, "fig3a_asthma_htn_individual_dotplots.pdf"), width = 10, height = 8)
for (comp in matched_comps) {
  df_sub <- df_all %>% filter(comparison_direction == comp)
  dir <- if(grepl("_down$", comp)) "down" else "up"
  p <- make_dotplot(df_sub, title = gsub("_", " ", comp), direction = dir)
  if (!is.null(p)) print(p)
}
dev.off()

# --- 2. Figure 3 Panel B (Cross-disease Shared) ---
cat("Generating Fig 3 Panel B (Cross-disease Shared)...\n")
df_xd <- df_all %>% filter(grepl("Cross_disease", comparison_direction))
if (nrow(df_xd) > 0) {
  pdf(file.path(OUT_FIG3, "fig3b_cross_disease_shared_enrichment.pdf"), width = 12, height = 10)
  for (comp in unique(df_xd$comparison_direction)) {
    df_sub <- df_xd %>% filter(comparison_direction == comp)
    dir <- if(grepl("_down$", comp)) "down" else "up"
    p <- make_dotplot(df_sub, title = gsub("_", " ", comp), direction = dir)
    if (!is.null(p)) print(p)
  }
  dev.off()
}

# --- 3. Figure S4a (Overlap Summary) ---
cat("Generating Fig S4a (Overlap Summary)...\n")
# Count terms shared across diseases (T2D, Asthma, HTN)
res_overlap <- df_all %>%
  mutate(disease = case_when(
    grepl("^T2D", comparison_direction) ~ "T2D",
    grepl("^Asthma", comparison_direction) ~ "Asthma",
    grepl("^HTN", comparison_direction) ~ "HTN",
    TRUE ~ "Other"
  )) %>%
  filter(disease != "Other") %>%
  group_by(Description) %>%
  summarise(
    diseases = paste(sort(unique(disease)), collapse = "+"),
    n_diseases = n_distinct(disease),
    .groups = "drop"
  ) %>%
  group_by(diseases, n_diseases) %>%
  summarise(term_count = n(), .groups = "drop") %>%
  arrange(desc(n_diseases), desc(term_count))

p_overlap <- ggplot(res_overlap, aes(x = reorder(diseases, -term_count), y = term_count, fill = as.factor(n_diseases))) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = term_count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set1", name = "N Diseases Shared") +
  labs(title = "Pathway Overlap Summary",
       subtitle = "Unique enriched terms (GO BP + KEGG) shared across diseases",
       x = "Disease Overlap", y = "Count of Enriched Terms") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_SUPP, "figs4a_pathway_overlap_summary.pdf"), p_overlap, width = 8, height = 7)

# --- 4. Table S5 Consolidation ---
cat("Saving Table S5...\n")
wb <- createWorkbook()
addWorksheet(wb, "Table_S5_Enrichment")
writeData(wb, "Table_S5_Enrichment", df_all)
saveWorkbook(wb, file.path(OUT_TAB, "Table_S5_Pathway_Enrichment_per_Comparison.xlsx"), overwrite = TRUE)

cat("Pathway Enrichment reproduction complete.\n")
