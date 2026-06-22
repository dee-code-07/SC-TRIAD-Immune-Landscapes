# sc-triad project
# purpose: reproduce Figure S1 (QC) using EXACTLY the original script logic
# author: agent

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggrastr)
})

set.seed(42)

base <- file.path(Sys.getenv("HOME"), "sc-triad")
base_theme <- theme_classic(base_size = 12) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

datasets <- list(
  t2d = list(
    disease = "T2D",
    path    = file.path(base, "01_raw_data", "01_t2d", "gse255566_pbmc"),
    samples = list(
      GSM8075321 = list(label = "Con1", group = "Control"),
      GSM8075322 = list(label = "Con2", group = "Control"),
      GSM8075323 = list(label = "Con3", group = "Control"),
      GSM8075324 = list(label = "Mod1", group = "T2D_Moderate"),
      GSM8075325 = list(label = "Mod2", group = "T2D_Moderate"),
      GSM8075326 = list(label = "Mod3", group = "T2D_Moderate")
    )
  ),
  htn = list(
    disease = "HTN",
    path    = file.path(base, "01_raw_data", "02_htn", "gse212953_pbmc"),
    samples = list(
      GSM6564434 = list(label = "Control",      group = "Control"),
      GSM6564435 = list(label = "Hypertensive", group = "HTN")
    )
  ),
  asthma = list(
    disease = "Asthma",
    path    = file.path(base, "01_raw_data", "03_asthma", "gse288147_pbmc"),
    samples = list(
      GSM8759883 = list(label = "A01", group = "Asthma_Mild"),
      GSM8759884 = list(label = "A02", group = "Asthma_Mild"),
      GSM8759885 = list(label = "A03", group = "Asthma_Mild"),
      GSM8759886 = list(label = "A04", group = "Asthma_Mild"),
      GSM8759887 = list(label = "A05", group = "Asthma_Mild"),
      GSM8759888 = list(label = "A06", group = "Asthma_Mild"),
      GSM8759889 = list(label = "B02", group = "Asthma_Severe"),
      GSM8759890 = list(label = "B03", group = "Asthma_Severe"),
      GSM8759891 = list(label = "B04", group = "Asthma_Severe"),
      GSM8759892 = list(label = "B05", group = "Asthma_Severe"),
      GSM8759893 = list(label = "C02", group = "Control"),
      GSM8759894 = list(label = "C03", group = "Control"),
      GSM8759895 = list(label = "C04", group = "Control")
    )
  )
)

load_10x_sample <- function(data_path, gsm_id, label, group, disease) {
  all_files  <- list.files(data_path, full.names = TRUE)
  gsm_files  <- all_files[grepl(gsm_id, all_files)]
  tmp <- file.path(tempdir(), gsm_id)
  dir.create(tmp, showWarnings = FALSE)
  file.copy(gsm_files[grepl("barcodes", gsm_files)], file.path(tmp, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(gsm_files[grepl("features", gsm_files)], file.path(tmp, "features.tsv.gz"), overwrite = TRUE)
  file.copy(gsm_files[grepl("matrix",   gsm_files)], file.path(tmp, "matrix.mtx.gz"),   overwrite = TRUE)
  mat <- Read10X(data.dir = tmp)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  obj <- CreateSeuratObject(counts = mat, project = label, min.cells = 3, min.features = 200)
  obj$sample_id    <- gsm_id
  obj$sample_label <- label
  obj$group        <- group
  obj$disease      <- disease
  obj              <- RenameCells(obj, add.cell.id = label)
  return(obj)
}

for (ds_key in names(datasets)) {
  cat("========== Processing", ds_key, "==========\n")
  ds <- datasets[[ds_key]]
  fig_dir <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Load Raw Data (Fast)
  obj_list <- list()
  for (gsm_id in names(ds$samples)) {
    meta <- ds$samples[[gsm_id]]
    obj_list[[gsm_id]] <- load_10x_sample(ds$path, gsm_id, meta$label, meta$group, ds$disease)
  }
  merged <- if (length(obj_list) == 1) obj_list[[1]] else merge(obj_list[[1]], y = obj_list[-1], merge.data = FALSE)
  
  merged[["percent_mt"]]          <- PercentageFeatureSet(merged, pattern = "^MT-")
  merged[["percent_ribo"]]        <- PercentageFeatureSet(merged, pattern = "^RP[SL]")
  merged[["log10_genes_per_umi"]] <- log10(merged$nFeature_RNA) / log10(merged$nCount_RNA)

  # 2. EXACT S1A: Original Violin Plot Logic
  p_violin <- VlnPlot(merged,
    features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"),
    group.by = "sample_label", pt.size = 0, ncol = 4) &
    base_theme &
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  pdf_path_A <- file.path(fig_dir, paste0(ds_key, "_S1A_prefilter_violin.pdf"))
  ggsave(pdf_path_A, p_violin, width = 18, height = 5)
  cat("  Saved S1A:", pdf_path_A, "\n")

  # 3. EXACT S1B: Original Scatter Plot Logic (with ggrastr for performance)
  suppressWarnings({
    p_sc1 <- FeatureScatter(merged, "nCount_RNA", "nFeature_RNA",
      group.by = "sample_label", pt.size = 0.1, raster = FALSE) +
      base_theme + ggtitle("nCount vs nFeature") +
      guides(color=guide_legend(override.aes=list(size=4, alpha=1)))
    p_sc1 <- ggrastr::rasterise(p_sc1, layers="Point", dpi=600)

    p_sc2 <- FeatureScatter(merged, "nCount_RNA", "percent_mt",
      group.by = "sample_label", pt.size = 0.1, raster = FALSE) +
      base_theme +
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.6) +
      ggtitle("nCount vs % Mitochondrial") +
      guides(color=guide_legend(override.aes=list(size=4, alpha=1)))
    p_sc2 <- ggrastr::rasterise(p_sc2, layers="Point", dpi=600)

    p_sc3 <- FeatureScatter(merged, "nFeature_RNA", "percent_mt",
      group.by = "sample_label", pt.size = 0.1, raster = FALSE) +
      base_theme +
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.6) +
      ggtitle("nFeature vs % Mitochondrial") +
      guides(color=guide_legend(override.aes=list(size=4, alpha=1)))
    p_sc3 <- ggrastr::rasterise(p_sc3, layers="Point", dpi=600)
  })

  p_scatter <- p_sc1 + p_sc2 + p_sc3
  pdf_path_B <- file.path(fig_dir, paste0(ds_key, "_S1B_prefilter_scatter.pdf"))
  ggsave(pdf_path_B, p_scatter, width = 18, height = 5)
  cat("  Saved S1B:", pdf_path_B, "\n")
  
  # 4. S1C: Doublets (Fast Reconstruct from CSV)
  csv_path <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key, "reports", "cell_counts_all_stages.csv")
  if (file.exists(csv_path)) {
    counts <- read.csv(csv_path)
    df_prop <- data.frame(
      sample_label = rep(counts$sample_label, 2),
      scDbl_class = rep(c("singlet", "doublet"), each = nrow(counts)),
      count = c(counts$n_cells_post_mad - counts$n_doublets, counts$n_doublets)
    )
    
    p_prop <- ggplot(df_prop, aes(x = sample_label, y = count, fill = scDbl_class)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = c("singlet" = "#4DADE2", "doublet" = "#E84B35")) +
      base_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(ds_key, "— doublet proportions"), x = NULL, y = "proportion", fill = "class")
           
    # Load qc_passed object just to plot the singlet violins
    qc_rds_path <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key, paste0(ds_key, "_qc_passed.rds"))
    if (file.exists(qc_rds_path)) {
      obj_qc <- readRDS(qc_rds_path)
      p_score <- ggplot(obj_qc@meta.data,
        aes(x = sample_label, y = scDbl_score, fill = scDbl_class)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
        scale_fill_manual(values = c("singlet" = "#4DADE2", "doublet" = "#E84B35")) +
        base_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste(ds_key, "— doublet scores (singlets only remaining)"), x = NULL, y = "scDblFinder score", fill = "class")
      rm(obj_qc)
      
      p_doublets <- p_score + p_prop
      pdf_path_C <- file.path(fig_dir, paste0(ds_key, "_S1C_doublets.pdf"))
      ggsave(pdf_path_C, p_doublets, width = 14, height = 5)
      cat("  Saved S1C:", pdf_path_C, "\n")
    }
  }

  rm(merged, obj_list)
  gc()
}

cat("Finished reproducing TRUE Figure S1 pre-filter plots.\n")
