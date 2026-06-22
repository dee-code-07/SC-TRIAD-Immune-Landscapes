# ==============================================================================
# 07_trajectory_manuscript.R (PDF Version)
# Purpose: Trajectory analysis for Fig 5 and Supp Figs S5.1-S5.3
# Output: PDF figures and CSV/XLSX tables
# ==============================================================================

suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(viridis)
  library(openxlsx)
})

set.seed(42)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 32 * 1024^3)

# Paths (Windows format)
BASE_DIR     <- "/home/deekshah/sc-triad"
SEURAT_PATH  <- file.path(BASE_DIR, "02_scrna/03_integration/triad_integrated_nkt_patched.rds")
OUT_DIR      <- file.path(BASE_DIR, "02_scrna/07_trajectory")
FIG_DIR      <- file.path(OUT_DIR, "figures")
TABLE_DIR    <- file.path(OUT_DIR, "tables")
OBJ_DIR      <- file.path(OUT_DIR, "objects")
MAIN_FIG_DIR <- file.path(BASE_DIR, "manuscript_new/main_figures/fig5_trajectory")
SUPP_FIG_DIR <- file.path(BASE_DIR, "manuscript_new/supplementary_figures/figs5_trajectory")

for (d in c(FIG_DIR, TABLE_DIR, OBJ_DIR, MAIN_FIG_DIR, SUPP_FIG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Constants
COND_ORDER  <- c("Control", "T2D_Moderate", "HTN", "Asthma_Mild", "Asthma_Severe")
COND_COLORS <- c(Control       = "#7CAE00",
                 T2D_Moderate  = "#F8766D",
                 HTN           = "#C77CFF",
                 Asthma_Mild   = "#00BFC4",
                 Asthma_Severe = "#00A9FF")

# Trajectory Setup
TRAJECTORIES <- list(
  myeloid = list(name = "Myeloid", cell_types = c("CD14 Monocyte", "CD16 Monocyte", "DC"), start = "CD14 Monocyte", end = "DC", tag = "myeloid"),
  tcell   = list(name = "T Cell", cell_types = c("Naive CD4 T", "Memory CD4 T", "CD8 T"), start = "Naive CD4 T", end = c("Memory CD4 T", "CD8 T"), tag = "tcell"),
  nk      = list(name = "NK", cell_types = c("NK"), start = NULL, tag = "nk")
)

FOCUS_GENES <- list(
  myeloid = c("CD14", "FCGR3A", "FLT3", "ITGAX", "LGALS9", "S100A8", "S100A9", "HLA-DRA", "CCL3", "IL1B"),
  nk      = c("SELL", "XCL1", "XCL2", "CXCR3", "NKG7", "GZMB", "GNLY", "GZMA", "GZMK", "PRF1", "FCGR3A"),
  tcell   = c("CCR7", "SELL", "IL7R", "TCF7", "LEF1", "CD44", "S100A4", "CD74", "HLA-DRA", "GZMB", "GZMA", "PRF1")
)

# Theme
pub_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

# Helper: Save PDF
save_pdf <- function(dir, fname, width = 8, height = 7, expr) {
  path <- file.path(dir, paste0(fname, ".pdf"))
  pdf(path, width = width, height = height, useDingbats = FALSE)
  tryCatch(force(expr), error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("error:", conditionMessage(e)), cex = 0.9)
    message("  plot error in ", fname, ": ", conditionMessage(e))
  })
  dev.off()
  cat("  Saved:", path, "\n")
}

# 1. Load Data
cat("Loading Seurat object (patched)...\n")
seu <- readRDS(SEURAT_PATH)

if (!"cell_type_clean" %in% colnames(seu@meta.data)) {
  stop("Wrong object loaded. Use triad_integrated_nkt_patched.rds")
}
Idents(seu) <- "cell_type_clean"
seu$cell_type <- seu$cell_type_clean

cond_col <- "group"
ct_col   <- "cell_type_clean"

# --- SEURAT v5 LAYER HEALING ---
if (length(grep("^data", Layers(seu))) > 1) {
  cat("  Split layers detected. Joining layers for analysis...\n")
  seu <- JoinLayers(seu)
}
if (!"data" %in% Layers(seu)) {
  cat("  Data layer missing. Normalizing...\n")
  seu <- NormalizeData(seu)
}

# 2. Main Processing Loop
for (traj_id in names(TRAJECTORIES)) {
  traj <- TRAJECTORIES[[traj_id]]
  cat("\n--- Processing:", traj$name, "---\n")
  
  cells <- colnames(seu)[seu@meta.data[[ct_col]] %in% traj$cell_types]
  seu_sub <- subset(seu, cells = cells)
  
  # Clean gene names for correlation later
  rownames(seu_sub) <- gsub("\\..*$", "", rownames(seu_sub))
  
  # Reduced Dim (UMAP)
  umap_mat <- Embeddings(seu_sub, "umap")
  
  # Slingshot Checkpoint
  sds_path <- file.path(OBJ_DIR, paste0("sds_", traj_id, ".rds"))
  if (file.exists(sds_path)) {
    cat("  Loading Slingshot checkpoint...\n")
    sds <- readRDS(sds_path)
  } else {
    cat("  Running Slingshot...\n")
    if (traj_id == "nk") {
      nk_markers <- intersect(c("NKG7", "GNLY", "PRF1"), rownames(seu_sub))
      expr <- colMeans(as.matrix(LayerData(seu_sub, layer = "data")[nk_markers, ]))
      ct_sling <- ifelse(expr <= quantile(expr, 0.25), "NK_immature", "NK_cytotoxic")
      start_clus <- "NK_immature"
    } else {
      ct_sling <- as.character(seu_sub@meta.data[[ct_col]])
      start_clus <- traj$start
    }
    sds <- slingshot(umap_mat, clusterLabels = ct_sling, start.clus = start_clus, end.clus = traj$end)
    saveRDS(sds, sds_path)
  }
  
  # Pseudotime
  pt <- rowMeans(slingPseudotime(sds), na.rm = TRUE)
  seu_sub$pseudotime <- pt
  
  # --- Fig 5 Panels (a, c, e: UMAP + Curves) ---
  save_pdf(MAIN_FIG_DIR, paste0(traj_id, "_pseudotime_umap"), expr = {
    plot_df <- data.frame(dim1=umap_mat[,1], dim2=umap_mat[,2], pt=pt)
    p <- ggplot(plot_df, aes(dim1, dim2, color=pt)) +
      geom_point(size=0.5, alpha=0.6) + scale_color_viridis_c(option="plasma", name="Pseudotime") +
      pub_theme + labs(title=paste(traj$name, "Trajectory"), x="UMAP 1", y="UMAP 2")
    
    curves <- slingCurves(sds)
    for(i in seq_along(curves)) {
      curve_mat <- curves[[i]]$s[curves[[i]]$ord, ]
      curve_df <- data.frame(dim1 = curve_mat[,1], dim2 = curve_mat[,2])
      p <- p + geom_path(data=curve_df, aes(dim1, dim2), color="black", linewidth=1.2, inherit.aes = FALSE)
    }
    print(p)
  })
  
  # --- Fig 5 Panels (b, f: Pseudotime Violin) ---
  if (traj_id != "tcell") {
    save_pdf(MAIN_FIG_DIR, paste0(traj_id, "_pseudotime_violin"), expr = {
      df_vio <- data.frame(pt=pt, cond=factor(seu_sub@meta.data[[cond_col]], levels=COND_ORDER))
      print(ggplot(df_vio, aes(cond, pt, fill=cond)) +
        geom_violin(trim=FALSE, alpha=0.7) + geom_boxplot(width=0.1, fill="white") +
        scale_fill_manual(values=COND_COLORS) + pub_theme +
        labs(title=paste(traj$name, "Pseudotime"), x=NULL, y="Pseudotime") +
        theme(axis.text.x = element_text(angle=45, hjust=1)))
    })
  }
  
  # --- Fig 5 Panel d (T cell Branch Weights) ---
  if (traj_id == "tcell") {
    cat("  Computing branch weights...\n")
    weights <- slingCurveWeights(sds)
    bw_df <- as.data.frame(weights) %>%
      mutate(cond = factor(seu_sub@meta.data[[cond_col]], levels=COND_ORDER)) %>%
      group_by(cond) %>%
      summarise(across(everything(), mean)) %>%
      pivot_longer(-cond, names_to="Lineage", values_to="Weight")
    
    save_pdf(MAIN_FIG_DIR, "tcell_branch_weights", expr = {
      print(ggplot(bw_df, aes(x=cond, y=Weight, fill=Lineage)) +
        geom_bar(stat="identity", position="fill") +
        scale_fill_brewer(palette="Set1", labels=c("Memory CD4", "CD8 T")) + pub_theme +
        labs(title="Lineage Commitment", x=NULL, y="Proportion"))
    })
    write.csv(bw_df, file.path(TABLE_DIR, "tcell_branch_weights.csv"), row.names = FALSE)
  }
  
  # --- Supp Fig S5.1 (Gene Profiles) ---
  fg <- intersect(FOCUS_GENES[[traj_id]], rownames(seu_sub))
  if (length(fg) > 0) {
    save_pdf(SUPP_FIG_DIR, paste0(traj_id, "_focus_genes_pseudotime"), width=12, height=10, expr = {
      df_fg <- FetchData(seu_sub, vars = c(fg, "pseudotime", cond_col)) %>%
        pivot_longer(all_of(fg), names_to="gene", values_to="expr") %>%
        mutate(condition = factor(!!sym(cond_col), levels=COND_ORDER))
      print(ggplot(df_fg, aes(pseudotime, expr, color=condition)) +
        geom_smooth(method="loess", span=0.75) + facet_wrap(~gene, scales="free_y", ncol=3) +
        scale_color_manual(values=COND_COLORS) + pub_theme +
        labs(title=paste(traj$name, "Focus Gene Kinetics")))
    })
  }
  
  # --- Tables ---
  df_sum <- data.frame(pt=pt, cond=seu_sub@meta.data[[cond_col]]) %>%
    group_by(cond) %>%
    summarise(mean_pt = mean(pt, na.rm=TRUE), sd_pt = sd(pt, na.rm=TRUE), n = n(), .groups = "drop") %>%
    mutate(trajectory = traj_id)
  write.csv(df_sum, file.path(TABLE_DIR, paste0(traj_id, "_pseudotime_summary.csv")), row.names = FALSE)

  kw <- kruskal.test(pt ~ seu_sub@meta.data[[cond_col]])
  write.csv(data.frame(test="Kruskal-Wallis", chi2=kw$statistic, p=kw$p.value), 
            file.path(TABLE_DIR, paste0(traj_id, "_kruskal_wallis.csv")))
  
  cat("  Gene correlations...\n")
  expr_mat_sparse <- LayerData(seu_sub, layer="data")
  genes_keep <- rownames(expr_mat_sparse)[!grepl("^RPL|^RPS|^MRPL|^MRPS|^MT-", rownames(expr_mat_sparse))]
  # Sample 2000 genes for speed
  genes_use <- genes_keep[1:min(2000, length(genes_keep))]
  
  cor_res <- sapply(genes_use, function(g) {
    x <- as.numeric(expr_mat_sparse[g, ])
    cor(x, pt, method="spearman", use = "pairwise.complete.obs")
  })
  
  write.csv(data.frame(gene=names(cor_res), rho=cor_res) %>% arrange(desc(abs(rho))), 
            file.path(TABLE_DIR, paste0(traj_id, "_pseudotime_genes.csv")), row.names = FALSE)
}

# 3. Final Table Consolidation (S14-S19)
cat("\nConsolidating Tables S14-S19...\n")
wb <- createWorkbook()
addWorksheet(wb, "S14_Pseudotime_Stats")
all_sum <- do.call(rbind, lapply(list.files(TABLE_DIR, pattern="_pseudotime_summary.csv", full.names = TRUE), read.csv))
writeData(wb, "S14_Pseudotime_Stats", all_sum)
addWorksheet(wb, "S15_Myeloid_Genes")
writeData(wb, "S15_Myeloid_Genes", read.csv(file.path(TABLE_DIR, "myeloid_pseudotime_genes.csv")))
addWorksheet(wb, "S16_NK_Genes")
writeData(wb, "S16_NK_Genes", read.csv(file.path(TABLE_DIR, "nk_pseudotime_genes.csv")))
addWorksheet(wb, "S17_Tcell_Genes")
writeData(wb, "S17_Tcell_Genes", read.csv(file.path(TABLE_DIR, "tcell_pseudotime_genes.csv")))
addWorksheet(wb, "S18_Kruskal_Wallis")
kw_sum <- do.call(rbind, lapply(list.files(TABLE_DIR, pattern="_kruskal_wallis.csv", full.names = TRUE), function(f) { d <- read.csv(f); d$trajectory <- basename(f); d }))
writeData(wb, "S18_Kruskal_Wallis", kw_sum)
addWorksheet(wb, "S19_Branch_Weights")
writeData(wb, "S19_Branch_Weights", read.csv(file.path(TABLE_DIR, "tcell_branch_weights.csv")))
saveWorkbook(wb, file.path(BASE_DIR, "manuscript_new/supplementary_tables/Table_S14_S19_Trajectory_Supplement.xlsx"), overwrite = TRUE)

cat("\nDone. All figures saved as PDF in manuscript_new/.\n")
