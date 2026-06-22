# sc-triad project
# script: 03b_aggregate.R
# purpose: aggregate per-section RCTD results into tables + figures
#
# run AFTER all 16 array tasks from pancreas_rctd_array.sh are complete
# sbatch pancreas_aggregate.sh
#
# this script reads all per-section RDS files, extracts proportions,
# computes differential proportions, and generates all publication figures
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe

suppressPackageStartupMessages({
  library(spacexr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(ggrepel)
  library(stringr)
})

set.seed(42)

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
PAN_DIR <- file.path(BASE, "01_raw_data/01_t2d/gse264331_pancreas_visium")
OUT_DIR <- file.path(BASE, "04_spatial/pancreas")

for (d in c("logs", "figures", "tables"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "logs/03b_aggregate.log")
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

save_pdf <- function(plot, path, width = 10, height = 8) {
  pdf(path, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  log("  saved: ", basename(path))
}

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

islet_types    <- c("beta", "alpha", "delta", "gamma", "epsilon")
exocrine_types <- c("acinar", "ductal")
stromal_types  <- c("stellate", "activated_stellate", "quiescent_stellate",
                    "endothelial")
immune_types_p <- c("macrophage", "T_cell")
status_colors  <- c("Control" = "#377EB8", "T2D" = "#E41A1C")

log("=== SC-TRIAD pancreas spatial: aggregation ===")

# ── donor metadata ─────────────────────────────────────────────────────────────
DONOR_STATUS <- c(
  "HP21024" = "T2D", "HP21035" = "T2D", "HP21091" = "T2D",
  "HP21161" = "T2D", "HP21177" = "T2D", "HP21181" = "Control"
)

h5_files <- list.files(PAN_DIR, pattern = "\\.h5$", full.names = TRUE)
parse_meta <- function(fp) {
  fname   <- basename(fp)
  donor   <- str_extract(fname, "HP\\d+")
  section <- str_extract(fname, "Vis\\d+_S\\d+_\\w+")
  if (is.na(section)) section <- str_extract(fname, "Vis\\d+_S\\d+")
  status  <- ifelse(!is.na(donor) && donor %in% names(DONOR_STATUS),
                    DONOR_STATUS[[donor]], "Unknown")
  data.frame(file = fname, donor = donor, section = section,
             status = status, stringsAsFactors = FALSE)
}
h5_meta  <- do.call(rbind, lapply(h5_files, parse_meta))
sample_ids <- paste0(h5_meta$donor, "_", h5_meta$section)

# ── collect per-section proportions ───────────────────────────────────────────
# ── collect per-section proportions ───────────────────────────────────────────
log("collecting per-section RCTD results ...")
all_prop_dfs <- list()
missing_sections <- character(0)

for (i in seq_along(sample_ids)) {
  sid      <- sample_ids[i]
  rctd_rds <- file.path(OUT_DIR, "objects", paste0("rctd_", sid, ".rds"))

  if (!file.exists(rctd_rds)) {
    log("  MISSING: ", sid, " (array task ", i, " may have failed)")
    missing_sections <- c(missing_sections, sid)
    next
  }

  log("  loading: ", basename(rctd_rds))
  myRCTD   <- readRDS(rctd_rds)
  weights  <- myRCTD@results$weights
  prop_mat <- normalize_weights(weights)

  # FIX: Add as.matrix() before as.data.frame()
  prop_df <- as.data.frame(as.matrix(prop_mat)) %>%
    mutate(spot    = rownames(prop_mat),
           section = sid,
           donor   = h5_meta$donor[i],
           status  = h5_meta$status[i])

  all_prop_dfs[[sid]] <- prop_df
}

if (length(missing_sections) > 0) {
  log("WARNING: ", length(missing_sections),
      " section(s) missing — results will be incomplete:")
  for (s in missing_sections) log("  ", s)
}

if (length(all_prop_dfs) == 0)
  stop("No RCTD results found. Run pancreas_rctd_array.sh first.")

all_prop <- do.call(rbind, all_prop_dfs)
rownames(all_prop) <- NULL
ct_cols  <- setdiff(colnames(all_prop), c("spot","section","donor","status"))

log("total spots collected: ", nrow(all_prop))
log("cell types: ", paste(ct_cols, collapse = ", "))
log("sections with results: ", length(all_prop_dfs),
    " / ", length(sample_ids))

# ── save spot-level table ──────────────────────────────────────────────────────
write.csv(all_prop,
          file.path(OUT_DIR, "tables/rctd_proportions_per_spot.csv"),
          row.names = FALSE)

# ── summaries ─────────────────────────────────────────────────────────────────
summary_status <- all_prop %>%
  group_by(status) %>%
  summarise(across(all_of(ct_cols), \(x) mean(x, na.rm = TRUE)),
            n_spots    = n(),
            n_sections = n_distinct(section),
            .groups    = "drop")
write.csv(summary_status,
          file.path(OUT_DIR, "tables/rctd_summary_by_status.csv"),
          row.names = FALSE)

summary_section <- all_prop %>%
  group_by(section, donor, status) %>%
  summarise(across(all_of(ct_cols), \(x) mean(x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
write.csv(summary_section,
          file.path(OUT_DIR, "tables/rctd_summary_by_section.csv"),
          row.names = FALSE)

t2d_secs  <- summary_section %>% filter(status == "T2D")
ctrl_secs <- summary_section %>% filter(status == "Control")
n_t2d     <- nrow(t2d_secs)
n_ctrl    <- nrow(ctrl_secs)
log("T2D sections: ", n_t2d, " | Control sections: ", n_ctrl)
if (n_ctrl < 2)
  log("WARNING: n=1 control donor. Effect size (log2FC + r) is primary metric.")

# ── differential proportions ──────────────────────────────────────────────────
diff_rows <- lapply(ct_cols, function(ct) {
  xt <- t2d_secs[[ct]]; xc <- ctrl_secs[[ct]]
  mt <- mean(xt, na.rm = TRUE); mc <- mean(xc, na.rm = TRUE)
  fc <- log2((mt + 1e-6) / (mc + 1e-6))
  pv  <- NA_real_; rbe <- NA_real_
  if (n_t2d >= 2 && n_ctrl >= 2) {
    wt  <- suppressWarnings(wilcox.test(xt, xc, exact = FALSE))
    pv  <- wt$p.value
    rbe <- round(1 - 2 * wt$statistic / (n_t2d * n_ctrl), 4)
  }
  data.frame(cell_type = ct, mean_t2d = round(mt, 5),
             mean_ctrl = round(mc, 5), log2FC = round(fc, 4),
             rank_biserial_r = rbe, p_value = pv,
             n_t2d = n_t2d, n_ctrl = n_ctrl, stringsAsFactors = FALSE)
})

diff_df <- do.call(rbind, diff_rows) %>% arrange(desc(abs(log2FC)))
diff_df$p_adj <- if (!all(is.na(diff_df$p_value)))
  p.adjust(diff_df$p_value, method = "BH") else NA_real_

write.csv(diff_df,
          file.path(OUT_DIR, "tables/rctd_differential_celltypes.csv"),
          row.names = FALSE)

log("top cell types by |log2FC|:")
print(head(diff_df[, c("cell_type","mean_t2d","mean_ctrl",
                        "log2FC","rank_biserial_r","p_adj")], 12))

# ── FIGURES ───────────────────────────────────────────────────────────────────
log("=== generating figures ===")

# FIG 1: Islet composition
islet_present <- intersect(islet_types, ct_cols)
if (length(islet_present) >= 2) {
  islet_df <- all_prop %>%
    select(section, donor, status, all_of(islet_present)) %>%
    group_by(section, donor, status) %>%
    summarise(across(all_of(islet_present),
                     \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(all_of(islet_present),
                 names_to = "cell_type", values_to = "mean_prop") %>%
    mutate(status = factor(status, levels = c("Control", "T2D")))

  islet_sum <- islet_df %>%
    group_by(status, cell_type) %>%
    summarise(mean = mean(mean_prop, na.rm = TRUE),
              se   = sd(mean_prop,   na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  p1 <- ggplot(islet_sum,
               aes(x = cell_type, y = mean,
                   fill = status, colour = status)) +
    geom_bar(stat = "identity", position = position_dodge(0.75),
             width = 0.65, alpha = 0.88) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.75),
                  width = 0.25, linewidth = 0.6, colour = "grey30") +
    geom_point(data = islet_df,
               aes(y = mean_prop, group = status),
               position = position_dodge(0.75),
               size = 2.2, alpha = 0.75, shape = 21, colour = "grey20") +
    scale_fill_manual(values   = status_colors) +
    scale_colour_manual(values = status_colors) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(title    = "Islet cell type proportions: T2D vs Control",
         subtitle = paste0("RCTD full mode | mean \u00b1 SE per section\n",
                           "T2D: ", n_t2d, " sections | Control: ",
                           n_ctrl, " section(s)"),
         x = "Islet cell type", y = "Mean proportion",
         fill = "Status", colour = "Status") +
    pub_theme

  save_pdf(p1, file.path(OUT_DIR, "figures/01_islet_composition.pdf"),
           width = 9, height = 7)
}

# FIG 2: All cell types stacked bar
mean_all <- all_prop %>%
  pivot_longer(all_of(ct_cols),
               names_to = "cell_type", values_to = "proportion") %>%
  group_by(status, cell_type) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(status = factor(status, levels = c("Control", "T2D")))

p2 <- ggplot(mean_all, aes(x = status, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.25) +
  scale_y_continuous(labels = percent_format()) +
  labs(title    = "Pancreas Visium: cell type composition (T2D vs Control)",
       subtitle = "RCTD full mode | GSE264331 | Reference: Baron et al. 2016",
       x = NULL, y = "Mean proportion", fill = "Cell type") +
  pub_theme + theme(legend.key.size = unit(0.4, "cm"))

save_pdf(p2, file.path(OUT_DIR, "figures/02_celltype_composition.pdf"),
         width = 8, height = 8)

# FIG 3: Effect size plot
compartment_map <- c(
  setNames(rep("Islet",    length(islet_types)),    islet_types),
  setNames(rep("Exocrine", length(exocrine_types)), exocrine_types),
  setNames(rep("Stromal",  length(stromal_types)),  stromal_types),
  setNames(rep("Immune",   length(immune_types_p)), immune_types_p)
)

diff_plot <- diff_df %>%
  mutate(
    compartment = ifelse(cell_type %in% names(compartment_map),
                         compartment_map[cell_type], "Other"),
    label = ifelse(abs(log2FC) > 0.5 |
                   (!is.na(rank_biserial_r) & abs(rank_biserial_r) > 0.3),
                   cell_type, NA_character_)
  )

p3 <- ggplot(diff_plot,
             aes(x = log2FC, y = rank_biserial_r,
                 colour = compartment, shape = compartment, label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", colour = "grey80") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3.2, na.rm = TRUE, max.overlaps = 15,
                  fontface = "italic") +
  scale_colour_manual(values = c(Islet = "#E41A1C", Exocrine = "#FF7F00",
                                 Stromal = "#999999", Immune = "#377EB8",
                                 Other = "#BDBDBD")) +
  scale_shape_manual(values = c(Islet = 16, Exocrine = 17,
                                Stromal = 18, Immune = 15, Other = 4)) +
  labs(title    = "Differential cell type proportions: T2D vs Control",
       subtitle = paste0(
         "Effect size analysis | x = log2(T2D/Control)\n",
         "y = rank-biserial r | dotted = |log2FC| = 0.5\n",
         "n=", n_ctrl, " control section(s): effect size is primary metric"
       ),
       x = "log2 fold change (T2D / Control)",
       y = "Rank-biserial r (effect size, -1 to +1)",
       colour = "Compartment", shape = "Compartment") +
  pub_theme

save_pdf(p3, file.path(OUT_DIR, "figures/03_effect_size_plot.pdf"),
         width = 10, height = 8)

# FIG 4: Beta:Alpha ratio
if (all(c("beta", "alpha") %in% ct_cols)) {
  ba_df <- all_prop %>%
    select(section, donor, status, beta, alpha) %>%
    group_by(section, donor, status) %>%
    summarise(mean_beta  = mean(beta,  na.rm = TRUE),
              mean_alpha = mean(alpha, na.rm = TRUE), .groups = "drop") %>%
    mutate(ba_ratio = mean_beta / (mean_alpha + 1e-6),
           status   = factor(status, levels = c("Control", "T2D")))

  p4 <- ggplot(ba_df, aes(x = status, y = ba_ratio,
                           fill = status, colour = status)) +
    geom_violin(scale = "width", trim = FALSE, alpha = 0.7, linewidth = 0.4) +
    geom_boxplot(width = 0.2, fill = "white",
                 outlier.shape = NA, linewidth = 0.5) +
    geom_jitter(width = 0.08, size = 3, alpha = 0.85,
                shape = 21, colour = "grey20") +
    scale_fill_manual(values   = status_colors, guide = "none") +
    scale_colour_manual(values = status_colors, guide = "none") +
    labs(title    = "Beta:Alpha cell ratio in T2D pancreas (spatial)",
         subtitle = paste0(
           "RCTD section-mean proportions | each point = one Visium section\n",
           "Reduced ratio expected in T2D: beta cell loss + alpha dominance"
         ),
         x = NULL, y = "Beta:Alpha proportion ratio") +
    pub_theme

  save_pdf(p4, file.path(OUT_DIR, "figures/04_beta_alpha_ratio.pdf"),
           width = 7, height = 7)
  log("Beta:Alpha ratio per section:")
  print(ba_df[, c("section","status","mean_beta","mean_alpha","ba_ratio")])
}

# FIG 5: SC-TRIAD immune convergence
imm_present <- intersect(c("macrophage", "T_cell"), ct_cols)
if (length(imm_present) >= 1) {
  imm_df <- all_prop %>%
    select(section, donor, status, all_of(imm_present)) %>%
    group_by(section, donor, status) %>%
    summarise(across(all_of(imm_present),
                     \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(all_of(imm_present),
                 names_to = "cell_type", values_to = "mean_proportion") %>%
    mutate(status = factor(status, levels = c("Control", "T2D")))

  imm_sum <- imm_df %>%
    group_by(status, cell_type) %>%
    summarise(mean = mean(mean_proportion, na.rm = TRUE),
              se   = sd(mean_proportion,   na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  p5 <- ggplot(imm_sum,
               aes(x = cell_type, y = mean,
                   fill = status, colour = status)) +
    geom_bar(stat = "identity", position = position_dodge(0.75),
             width = 0.65, alpha = 0.88) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.75),
                  width = 0.25, linewidth = 0.6, colour = "grey30") +
    geom_point(data = imm_df,
               aes(y = mean_proportion, group = status),
               position = position_dodge(0.75),
               size = 2, alpha = 0.75, shape = 21, colour = "grey20") +
    scale_fill_manual(values   = status_colors) +
    scale_colour_manual(values = status_colors) +
    scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
    labs(title    = "SC-TRIAD: immune infiltration in T2D pancreas (spatial)",
         subtitle = paste0(
           "RCTD mean proportions \u00b1 SE | each point = one Visium section\n",
           "Connects circulating PBMC signatures to tissue-resident spatial states"
         ),
         x = NULL, y = "Mean proportion",
         fill = "Status", colour = "Status") +
    pub_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"))

  save_pdf(p5, file.path(OUT_DIR,
                          "figures/05_sctriad_pbmc_pancreas_convergence.pdf"),
           width = 9, height = 7)
}

# FIG 6: Per-section heatmap
top_ct <- intersect(head(diff_df$cell_type, min(15, nrow(diff_df))), ct_cols)
hm_in  <- summary_section %>%
  select(section, status, all_of(top_ct)) %>%
  mutate(status = factor(status, levels = c("Control", "T2D"))) %>%
  arrange(status, section)

hm_mat <- apply(as.matrix(hm_in[, top_ct, drop = FALSE]), 2, function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
})
rownames(hm_mat) <- hm_in$section

hm_long <- as.data.frame(hm_mat) %>%
  rownames_to_column("section") %>%
  left_join(hm_in[, c("section","status")], by = "section") %>%
  pivot_longer(-c(section, status),
               names_to = "cell_type", values_to = "scaled") %>%
  mutate(section   = factor(section, levels = hm_in$section),
         cell_type = factor(cell_type, levels = top_ct))

p6 <- ggplot(hm_long, aes(x = section, y = cell_type, fill = scaled)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#C00000",
                       midpoint = 0.5, name = "Scaled\nproportion") +
  facet_grid(. ~ status, scales = "free_x", space = "free_x") +
  labs(title    = "Per-section cell type proportions: T2D vs Control",
       subtitle = paste0("Top ", length(top_ct),
                         " cell types by |log2FC| | scaled within each cell type"),
       x = NULL, y = NULL) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9, face = "italic"))

save_pdf(p6, file.path(OUT_DIR, "figures/06_section_heatmap.pdf"),
         width = 12, height = 8)

log("\n=== AGGREGATION COMPLETE ===")
log("Sections processed: ", length(all_prop_dfs), " / ", length(sample_ids))
log("Output: ", OUT_DIR)
log("")
log("METHODS NOTE (copy to thesis):")
log("  Spatial deconvolution: RCTD (Cable et al. 2022) | GSE264331")
log("  Reference: Baron et al. 2016 (GSE84133) | 11 pancreatic cell types")
log("  Spot coordinates unavailable from GEO; placeholder grid assigned.")
log("  RCTD performs MLE independently per spot; coordinates do not affect")
log("  deconvolution accuracy. n=", n_t2d, " T2D, n=", n_ctrl,
    " control sections. Differential proportions reported as log2FC and")
log("  rank-biserial r. Wilcoxon p-values omitted (n=1 control donor).")
