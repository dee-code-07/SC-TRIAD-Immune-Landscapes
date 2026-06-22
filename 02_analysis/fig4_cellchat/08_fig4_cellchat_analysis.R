# ==============================================================================
# 06_cellchat_manuscript_robust.R (Complete Version)
# Purpose: CellChat analysis for the T2D–HTN–Asthma triad manuscript.
# Fixes: Seurat v5 Layer Joining, Dynamic DB detection, Isoform cleaning.
# Figures: Main Fig 4, Supplementary Figs S5-S9
# Tables: Supplementary Tables S6-S10
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(grid)
  library(scales)
})

set.seed(42)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 32 * 1024^3)

# Paths
SEURAT_PATH  <- "/home/deekshah/sc-triad/02_scrna/pbmc/03_integration/triad_integrated_nkt_patched.rds"
MAIN_FIG_DIR <- "/home/deekshah/sc-triad/manuscript_new/main_figures/fig4_cellchat"
SUPP_FIG_DIR <- "/home/deekshah/sc-triad/manuscript_new/supplementary_figures/figs5_cellchat"
SUPP_TBL_DIR <- "/home/deekshah/sc-triad/manuscript_new/supplementary_tables"

for (d in c(MAIN_FIG_DIR, SUPP_FIG_DIR, SUPP_TBL_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Constants
MIN_CELLS   <- 10
COND_ORDER  <- c("Control", "T2D_Moderate", "HTN", "Asthma_Mild", "Asthma_Severe")

# Publication Theme
pub_theme <- theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title        = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 5)),
    plot.subtitle     = element_text(size = 10, hjust = 0, colour = "grey40", margin = margin(b = 10), lineheight = 1.2),
    axis.title        = element_text(face = "bold", size = 11),
    axis.text         = element_text(size = 10, colour = "black"),
    axis.line         = element_line(colour = "grey70", linewidth = 0.4),
    axis.ticks        = element_line(colour = "grey70", linewidth = 0.4),
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA)
  )

# Helper: Save PDF
save_pdf <- function(dir, fname, width = 12, height = 10, expr) {
  path <- file.path(dir, paste0(fname, ".pdf"))
  pdf(path, width = width, height = height, useDingbats = FALSE)
  tryCatch(force(expr), error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("error:", conditionMessage(e)), cex = 0.9)
    message("  plot error in ", fname, ": ", conditionMessage(e))
  })
  dev.off()
  cat("  Saved PDF:", path, "\n")
}

# 1. Load and Heal Seurat Object
cat("Loading Seurat object (patched)...\n")
seu <- readRDS(SEURAT_PATH)

if (!"cell_type_clean" %in% colnames(seu@meta.data)) {
  stop("Wrong object loaded. Use triad_integrated_nkt_patched.rds")
}
Idents(seu) <- "cell_type_clean"
seu$cell_type <- seu$cell_type_clean

# Identify condition column
cond_col <- "group" # Verified from diagnostics

# --- SEURAT v5 LAYER HEALING ---
if (length(grep("^data", Layers(seu))) > 1) {
  cat("  Split layers detected. Joining layers for analysis...\n")
  seu <- JoinLayers(seu)
}

# Verify 'data' layer is now populated
raw_data <- LayerData(seu, assay = DefaultAssay(seu), layer = "data")

if (is.null(raw_data) || nrow(raw_data) == 0) {
  cat("  Warning: 'data' layer is still empty. Checking 'counts' layer...\n")
  raw_data <- LayerData(seu, assay = DefaultAssay(seu), layer = "counts")
  
  if (is.null(raw_data) || nrow(raw_data) == 0) {
    stop("FATAL: Both 'data' and 'counts' layers are empty.")
  }
  
  cat("  Normalizing counts since 'data' was missing...\n")
  seu <- NormalizeData(seu)
  raw_data <- LayerData(seu, assay = DefaultAssay(seu), layer = "data")
}

# --- CLEAN GENES ---
cat("Cleaning gene symbols (isoform removal)...\n")
clean_rows <- gsub("\\..*$", "", rownames(raw_data))

if (any(duplicated(clean_rows))) {
  cat("  Handling duplicates...\n")
  total_exp <- rowSums(raw_data)
  keep_idx <- order(total_exp, decreasing = TRUE)
  keep_idx <- keep_idx[!duplicated(clean_rows[keep_idx])]
  clean_matrix <- raw_data[keep_idx, ]
  rownames(clean_matrix) <- clean_rows[keep_idx]
} else {
  clean_matrix <- raw_data
  rownames(clean_matrix) <- clean_rows
}

# 2. ROBUST DATABASE LOADING
cat("\n--- Initializing CellChat Database ---\n")
data("CellChatDB.human", package = "CellChat")
if (!exists("CellChatDB.human")) CellChatDB.human <- get("CellChatDB.human", envir = asNamespace("CellChat"))

int_tab <- CellChatDB.human$interaction
db_genes <- unique(c(
  as.character(int_tab$ligand), as.character(int_tab$receptor),
  as.character(int_tab$ligand.symbol), as.character(int_tab$receptor.symbol)
))
db_genes <- db_genes[!is.na(db_genes) & db_genes != ""]

# Final Overlap Check
overlap <- intersect(rownames(clean_matrix), db_genes)
cat("  FINAL INTERSECTION:", length(overlap), "matching genes found.\n")

db_use <- subsetDB(CellChatDB.human, search = c("Secreted Signaling", "Cell-Cell Contact"))

# 3. Processing Function
make_cellchat <- function(condition) {
  cat("\nProcessing:", condition, "\n")
  cells_idx <- which(as.character(seu@meta.data[[cond_col]]) == condition)
  if (length(cells_idx) == 0) return(NULL)
  
  sub_matrix <- clean_matrix[, cells_idx]
  sub_meta <- seu@meta.data[cells_idx, ]
  
  # Filter cell types
  ct_tab <- table(sub_meta$cell_type_clean)
  keep_ct <- names(ct_tab[ct_tab >= MIN_CELLS])
  if (length(keep_ct) < 2) return(NULL)
  
  keep_cells <- which(sub_meta$cell_type_clean %in% keep_ct)
  
  cc <- createCellChat(object = sub_matrix[, keep_cells], 
                       meta = sub_meta[keep_cells, ], 
                       group.by = "cell_type_clean")
  
  cc@DB <- db_use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", trim = 0.1, population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc
}

# 4. Execution
cat("\n=== Building per-condition CellChat objects ===\n")
cc_list <- list()
for (cond in COND_ORDER) {
  cc <- tryCatch(make_cellchat(cond), error = function(e) {
    cat("  FAILED [", cond, "]:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(cc)) cc_list[[cond]] <- cc
}

available_conds <- names(cc_list)
if (length(available_conds) < 2) stop("Analysis failed: Not enough successful conditions.")

cat("\nMerging Results...\n")
all_ct <- sort(unique(unlist(lapply(cc_list, function(cc) levels(cc@idents)))))
cc_lifted <- lapply(cc_list, function(cc) liftCellChat(cc, group.new = all_ct))
cc_all <- mergeCellChat(cc_lifted, add.names = names(cc_lifted))

# ------------------------------------------------------------------------------
# MAIN FIGURE 4: Landscape Overview
# ------------------------------------------------------------------------------
cat("\nGenerating Main Figure 4...\n")

# (a) Barplots
save_pdf(MAIN_FIG_DIR, "Fig4a_all_interactions", width = 14, height = 6, {
  p1 <- compareInteractions(cc_all, show.legend = FALSE, group = seq_along(available_conds))
  p2 <- compareInteractions(cc_all, show.legend = FALSE, measure = "weight", group = seq_along(available_conds))
  print((p1 + pub_theme) + (p2 + pub_theme) + 
    plot_annotation(title = "Intercellular Communication Landscape", 
                    subtitle = "Left: Interaction count | Right: Aggregate interaction strength",
                    theme = theme(plot.title = element_text(face="bold", size=14),
                                  plot.subtitle = element_text(size=10, colour="grey40"))))
})

# (b) Pathway Presence Heatmap
all_pw <- sort(unique(unlist(lapply(cc_list, function(cc) cc@netP$pathways))))
pw_mat <- do.call(cbind, lapply(cc_list, function(cc) as.integer(all_pw %in% cc@netP$pathways)))
rownames(pw_mat) <- all_pw

ctrl_p <- if ("Control" %in% available_conds) pw_mat[, "Control"] == 1 else rep(FALSE, length(all_pw))
t2d_p  <- if ("T2D_Moderate" %in% available_conds) pw_mat[, "T2D_Moderate"] == 1 else rep(FALSE, length(all_pw))
htn_p  <- if ("HTN" %in% available_conds)          pw_mat[, "HTN"] == 1          else rep(FALSE, length(all_pw))
am_p   <- if ("Asthma_Mild" %in% available_conds)  pw_mat[, "Asthma_Mild"] == 1  else rep(FALSE, length(all_pw))
as_p   <- if ("Asthma_Severe" %in% available_conds) pw_mat[, "Asthma_Severe"] == 1 else rep(FALSE, length(all_pw))
asthma_any <- am_p | as_p

pw_group <- case_when(
  ctrl_p & t2d_p & htn_p & am_p & as_p    ~ "Core",
  !ctrl_p & t2d_p & htn_p & asthma_any    ~ "Triad-gained",
  !ctrl_p & t2d_p & !htn_p & !asthma_any  ~ "T2D-specific",
  !ctrl_p & !t2d_p & !htn_p & asthma_any  ~ "Asthma-specific",
  TRUE                                    ~ "Multi-disease shared"
)
names(pw_group) <- all_pw

hm_df <- as.data.frame(pw_mat) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "condition", values_to = "present") %>%
  mutate(condition = factor(condition, levels = COND_ORDER),
         module = pw_group[pathway])

save_pdf(MAIN_FIG_DIR, "Fig4b_pathway_presence_heatmap", width = 10, height = 14, {
  p <- ggplot(hm_df, aes(x = condition, y = pathway, fill = as.factor(present))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = "#EBEBEB", "1" = "#1A6B6B"), labels = c("Absent", "Present"), name = NULL) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    pub_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Pathway Presence across Conditions", y = NULL, x = NULL)
  print(p)
})

# (c) Info flow
save_pdf(MAIN_FIG_DIR, "Fig4c_pathway_comparison", width = 16, height = 10, {
  p_flow <- rankNet(cc_all, mode = "comparison", stacked = TRUE, do.stat = TRUE)
  print(p_flow + pub_theme)
})

# (d) Priority Chords
priority_pws <- c("GALECTIN", "CypA", "BTLA")
for (pw in priority_pws) {
  save_pdf(MAIN_FIG_DIR, paste0("Fig4d_chord_", pw), width = 14, height = 14, {
    active_conds_pw <- available_conds[sapply(available_conds, function(c) pw %in% cc_list[[c]]@netP$pathways)]
    if (length(active_conds_pw) > 0) {
      n_col <- ceiling(sqrt(length(active_conds_pw)))
      n_row <- ceiling(length(active_conds_pw) / n_col)
      par(mfrow = c(n_row, n_col))
      for (cond in active_conds_pw) {
        netVisual_aggregate(cc_list[[cond]], signaling = pw, layout = "chord", title.name = paste(cond, pw))
      }
    } else {
      plot.new(); text(0.5, 0.5, paste(pw, "not found in any condition"))
    }
  })
}

# ------------------------------------------------------------------------------
# SUPPLEMENTARY FIGURES S5-S9
# ------------------------------------------------------------------------------
cat("\nGenerating Supplementary Figures...\n")

# S8: Sender-receiver scatter plots
save_pdf(SUPP_FIG_DIR, "FigS8_sender_receiver_all", width = 24, height = 8, {
  sr_plots <- lapply(available_conds, function(cond) {
    tryCatch({
      netAnalysis_signalingRole_scatter(cc_list[[cond]], title = cond) + pub_theme
    }, error = function(e) ggplot() + labs(title = paste(cond, "error")) + theme_void())
  })
  print(wrap_plots(sr_plots, nrow = 1) + plot_annotation(title = "Sender-Receiver Centrality across Conditions"))
})

# S9: BTLA chord diagrams
save_pdf(SUPP_FIG_DIR, "FigS9_BTLA_asthma_chords", width = 15, height = 7, {
  asthma_conds <- intersect(available_conds, c("Asthma_Mild", "Asthma_Severe"))
  if (length(asthma_conds) > 0) {
    par(mfrow = c(1, length(asthma_conds)))
    for (cond in asthma_conds) {
      if ("BTLA" %in% cc_list[[cond]]@netP$pathways)
        netVisual_aggregate(cc_list[[cond]], signaling = "BTLA", layout = "chord", title.name = cond)
      else { plot.new(); text(0.5, 0.5, paste(cond, "BTLA absent")) }
    }
  }
})

# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLES (S9-S13)
# ------------------------------------------------------------------------------
cat("\nSaving Supplementary Tables S9-S13...\n")

# S9: Interaction summary
sum_df <- do.call(rbind, lapply(available_conds, function(cond) {
  cc <- cc_list[[cond]]
  data.frame(
    condition      = cond, 
    n_interactions = sum(cc@net$count, na.rm = TRUE), 
    total_strength = sum(cc@net$weight, na.rm = TRUE), 
    n_pathways     = length(cc@netP$pathways), 
    n_cell_types   = length(levels(cc@idents)), 
    exploratory    = (cond == "HTN")
  )
}))
write.xlsx(sum_df, file.path(SUPP_TBL_DIR, "Table_S9_interaction_summary.xlsx"))

# S10: Cell type presence matrix
ct_presence <- as.data.frame(table(seu@meta.data[[cond_col]], seu$cell_type_clean)) %>%
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  rename(cell_type = Var2)
write.xlsx(ct_presence, file.path(SUPP_TBL_DIR, "Table_S10_cell_type_presence.xlsx"))

# S11: Full ligand-receptor interactions
int_df <- do.call(rbind, lapply(available_conds, function(cond) {
  df <- tryCatch(subsetCommunication(cc_list[[cond]]), error = function(e) NULL)
  if (!is.null(df)) df$condition <- cond
  df
}))
write.xlsx(int_df, file.path(SUPP_TBL_DIR, "Table_S11_lr_interactions_per_condition.xlsx"))

# S12: Pathway-level flows
pw_df <- do.call(rbind, lapply(available_conds, function(cond) {
  df <- tryCatch(subsetCommunication(cc_list[[cond]], slot.name = "netP"), error = function(e) NULL)
  if (!is.null(df)) df$condition <- cond
  df
}))
write.xlsx(pw_df, file.path(SUPP_TBL_DIR, "Table_S12_pathway_interactions_per_condition.xlsx"))

# S13: Pathway presence matrix with module annotations
pw_presence <- cbind(data.frame(pathway = all_pw, module = pw_group[all_pw]), as.data.frame(pw_mat))
write.xlsx(pw_presence, file.path(SUPP_TBL_DIR, "Table_S13_pathway_presence_matrix.xlsx"))

cat("\nAnalysis complete. Manuscript PDF figures and Excel tables generated.\n")
