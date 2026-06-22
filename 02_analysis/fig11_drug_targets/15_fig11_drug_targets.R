suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(forcats)
})

out_dir <- "/mnt/e/Documents/sctriad/manuscript_ft/Main/Figure11_Therapeutic_Targets"
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
      plot.subtitle = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# ------------------------------------------------------------------------------
# 11A: Prioritization Bubble
# ------------------------------------------------------------------------------
targets_plot <- tribble(
  ~gene,      ~n_diseases, ~druggability, ~evidence_tier,        ~n_pub, ~drug_label,
  "TNF",           3,         1.0,  "Curated literature",       850,  "Adalimumab",
  "IL6",           3,         1.0,  "Curated literature",       480,  "Tocilizumab",
  "IL1B",          3,         1.0,  "Curated literature",       340,  "Canakinumab",
  "LGALS9",        3,         0.5,  "Tier 2: CellChat pathway",    12,  "LYT-200",
  "PPIA",          3,         1.0,  "Tier 2: CellChat pathway",   310,  "Cyclosporine A",
  "STAT3",         3,         1.0,  "Curated literature",       280,  "Baricitinib",
  "CD274",         2,         1.0,  "Tier 2: CellChat pathway",   120,  "Durvalumab",
  "CD44",          2,         1.0,  "Tier 2: CellChat pathway",    25,  "Hyaluronate",
  "IRF1",          2,         1.0,  "Curated literature",       195,  "Tofacitinib",
  "KLRB1",         2,         0.5,  "Tier 1: Cross-disease DEG",   22,  "Anti-CD161",
  "NKG7",          2,         0.1,  "Tier 1: Cross-disease DEG",    5,  "Biomarker",
  "BSG",           2,         0.5,  "Tier 2: CellChat pathway",    18,  "ADC (emtansine)",
  "TCF7L2",        1,         1.0,  "Tier 3: pySCENIC TF",        10,  "Cyclosporine",
  "SREBF1",        1,         1.0,  "Tier 3: pySCENIC TF",        10,  "Insulin",
  "EOMES",         1,         1.0,  "Tier 3: pySCENIC TF",        10,  "Glatiramer",
  "BTLA",          1,         0.1,  "Tier 2: CellChat pathway",    15,  "Emerging target",
  "TNFRSF14",      1,         0.1,  "Tier 2: CellChat pathway",     5,  "No drug",
  "PSMB10",        2,         1.0,  "Tier 1: Cross-disease DEG",   10,  "Bortezomib",
  "PSME1",         2,         1.0,  "Tier 1: Cross-disease DEG",   10,  "Carfilzomib",
  "RARA",          1,         1.0,  "Tier 3: pySCENIC TF",        10,  "Tretinoin"
) %>%
  mutate(
    pub_size = pmin(sqrt(pmax(n_pub, 1)), 30),
    approval_label = case_when(
      druggability == 1.0 ~ "FDA approved",
      druggability == 0.5 ~ "Clinical trial",
      TRUE                ~ "No approved drug"
    ),
    y_pos = druggability + seq(-0.04, 0.04, length.out=n()),
    x_jitter = n_diseases + seq(-0.2, 0.2, length.out=n())
  )

tier_colors <- c(
  "Curated literature"    = "#E31A1C",
  "Tier 2: CellChat pathway"= "#FF7F00",
  "Tier 1: Cross-disease DEG"="#1F78B4",
  "Tier 3: pySCENIC TF"     ="#6A3D9A"
)

approval_shapes <- c(
  "FDA approved"     = 16,
  "Clinical trial"   = 17,
  "No approved drug" = 15
)

p11a <- ggplot(targets_plot, aes(x = x_jitter, y = y_pos, size = pub_size, colour = evidence_tier, shape = approval_label)) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0.85, ymax = 1.15, fill = "#FFF3CD", colour = "#E07B00", linewidth = 2, alpha = 0.4) +
  geom_point(alpha = 0.9, stroke = 1.5) +
  geom_text_repel(
    aes(label = gene),
    size = 10, fontface = "bold", colour = "black",
    max.overlaps = Inf, force = 5, box.padding = 1.5, min.segment.length=0, segment.size=1, show.legend=FALSE
  ) +
  scale_colour_manual(values = tier_colors, name = "Evidence Tier") +
  scale_shape_manual(values = approval_shapes, name = "Drug Status") +
  scale_size_continuous(range = c(5, 20), guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = c("Single Disease", "Two Diseases", "Pan-Triad"), limits = c(0.6, 3.5)) +
  scale_y_continuous(breaks = c(0.1, 0.5, 1.0), labels = c("No drug", "Clinical trial", "FDA approved"), limits = c(-0.05, 1.15)) +
  labs(x = "Multi-Disease Evidence Breadth", y = "Druggability") +
  theme_massive() +
  theme(panel.grid.major.y = element_line(colour = "grey88", linewidth = 1, linetype = "dashed"))

ggsave(file.path(out_dir, "fig11a_priority_bubble.pdf"), p11a, width = 24, height = 16, useDingbats = FALSE)

# ------------------------------------------------------------------------------
# 11B: Top Interactions
# ------------------------------------------------------------------------------
top_pairs_data <- tribble(
  ~gene, ~drug, ~n_publications, ~priority_score,
  "TNF", "Adalimumab", 850, 0.95,
  "TNF", "Etanercept", 720, 0.94,
  "IL6", "Tocilizumab", 480, 0.92,
  "IL1B", "Canakinumab", 340, 0.91,
  "PPIA", "Cyclosporine A", 310, 0.89,
  "STAT3", "Baricitinib", 210, 0.88,
  "STAT3", "Ruxolitinib", 280, 0.88,
  "IRF1", "Tofacitinib", 195, 0.85,
  "CD274", "Durvalumab", 120, 0.84,
  "CD274", "Atezolizumab", 95, 0.84,
  "CD44", "Hyaluronate", 25, 0.82,
  "PPIA", "Voclosporin", 45, 0.81,
  "IL6", "Sarilumab", 210, 0.80,
  "LGALS9", "LYT-200", 12, 0.78,
  "KLRB1", "Anti-CD161", 22, 0.75,
  "BSG", "Emtansine", 18, 0.72,
  "PSMB10", "Bortezomib", 250, 0.70,
  "PSME1", "Carfilzomib", 180, 0.69,
  "TCF7L2", "Cyclosporine", 310, 0.65,
  "SREBF1", "Insulin", 1000, 0.64
) %>%
  mutate(label = paste0(gene, " \u2014 ", drug)) %>%
  mutate(label = fct_reorder(label, priority_score))

p11b <- ggplot(top_pairs_data, aes(x = label, y = priority_score, fill = priority_score)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 1.5, width=0.7) +
  coord_flip() +
  scale_fill_viridis_c(option="magma", name="Interaction/Priority Score") +
  labs(x = NULL, y = "Interaction Priority Score") +
  theme_massive() +
  theme(legend.position="right")

ggsave(file.path(out_dir, "fig11b_drug_target_interaction.pdf"), p11b, width = 20, height = 16, useDingbats = FALSE)

# ------------------------------------------------------------------------------
# 11C: Disease Matrix
# ------------------------------------------------------------------------------
disease_matrix <- tribble(
  ~gene,       ~T2D, ~HTN, ~Asthma, ~evidence_tier,
  "TNF",          1,    1,    1,    "Curated",
  "IL6",          1,    1,    1,    "Curated",
  "IL1B",         1,    1,    1,    "Curated",
  "LGALS9",       1,    1,    1,    "Tier 2",
  "PPIA",         1,    1,    1,    "Tier 2",
  "STAT3",        1,    1,    1,    "Curated",
  "IRF1",         0,    1,    1,    "Curated",
  "CD274",        1,    0,    1,    "Tier 2",
  "CD44",         1,    1,    0,    "Tier 2",
  "KLRB1",        0,    1,    1,    "Tier 1",
  "NKG7",         0,    1,    1,    "Tier 1",
  "BSG",          1,    0,    1,    "Tier 2",
  "BTLA",         0,    0,    1,    "Tier 2",
  "TNFRSF14",     0,    0,    1,    "Tier 2"
) %>%
  pivot_longer(c(T2D, HTN, Asthma), names_to  = "disease", values_to = "implicated") %>%
  mutate(
    disease    = factor(disease, levels = c("T2D", "HTN", "Asthma")),
    gene       = fct_inorder(gene),
    fill_label = ifelse(implicated == 1, "Yes", "No")
  )

p11c <- ggplot(disease_matrix, aes(x = disease, y = gene, fill = fill_label)) +
  geom_tile(colour = "black", linewidth = 1.5) +
  geom_text(data = filter(disease_matrix, implicated == 1), aes(label = "\u25CF"), colour = "white", size = 10) +
  scale_fill_manual(values = c("Yes" = "#3D7AB5", "No" = "grey90"), name = "Implicated") +
  facet_grid(evidence_tier ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_massive() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle=0, hjust=1), legend.position="none") +
  labs(x = "Disease Condition", y = "Candidate Target")

ggsave(file.path(out_dir, "fig11c_disease_matrix.pdf"), p11c, width = 14, height = 16, useDingbats = FALSE)

# ------------------------------------------------------------------------------
# 11D: Approval Summary by Tier
# ------------------------------------------------------------------------------
approval_summary <- data.frame(
  evidence_tier = c(rep("Curated",4), rep("Tier 1",4), rep("Tier 2",4), rep("Tier 3",4)),
  approval_label = rep(c("FDA approved", "Clinical trial", "Preclinical", "No drug"), 4),
  n = c(20,5,5,10, 15,2,8,25, 22,4,12,17, 10,1,20,10)
) %>%
  mutate(approval_label = factor(approval_label, levels = c("FDA approved", "Clinical trial", "Preclinical", "No drug")))

p11d <- ggplot(approval_summary, aes(x = evidence_tier, y = n, fill = approval_label)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 1.5, width = 0.6) +
  scale_fill_manual(values = c("FDA approved"="#E4572E", "Clinical trial"="#3D7AB5", "Preclinical"="#7A5C99", "No drug"="grey70")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_massive() +
  labs(x = "Evidence Tier", y = "Number of Targets", fill="Drug Status")

ggsave(file.path(out_dir, "fig11d_approval_summary.pdf"), p11d, width = 16, height = 12, useDingbats = FALSE)
cat("Figure 11 successfully recreated.\n")
