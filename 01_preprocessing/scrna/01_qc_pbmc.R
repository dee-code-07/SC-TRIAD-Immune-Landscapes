# sc-triad project
# script: 01_qc_pbmc.R
# purpose: quality control for pbmc scrnaseq datasets
#          t2d (gse255566), htn (gse212953), asthma (gse288147)
# seurat v5 | r 4.4.3
# author: deeksha h | manipal academy of higher education

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(42)

base <- file.path(Sys.getenv("HOME"), "sc-triad")

# directory layout:
# 02_scrna/pbmc/01_qc/{disease}/              <- rds lives here
# 02_scrna/pbmc/01_qc/{disease}/figures/      <- 600dpi tiffs + overview pdf
# 02_scrna/pbmc/01_qc/{disease}/reports/      <- csv tables + log

make_dirs <- function(disease) {
  root <- file.path(base, "02_scrna", "pbmc", "01_qc", disease)
  dirs <- list(
    root    = root,
    figures = file.path(root, "figures"),
    reports = file.path(root, "reports")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  return(dirs)
}

# global log sits in the pbmc qc root
log_root <- file.path(base, "02_scrna", "pbmc", "01_qc")
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_root, "qc_pipeline.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log("sc-triad pbmc qc pipeline started")
log(paste("seurat version:", packageVersion("Seurat")))


# dataset manifest

datasets <- list(

  t2d = list(
    name    = "T2D_PBMC_GSE255566",
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
    name    = "HTN_PBMC_GSE212953",
    disease = "HTN",
    path    = file.path(base, "01_raw_data", "02_htn", "gse212953_pbmc"),
    samples = list(
      GSM6564434 = list(label = "Control",      group = "Control"),
      GSM6564435 = list(label = "Hypertensive", group = "HTN")
    )
  ),

  asthma = list(
    name    = "Asthma_PBMC_GSE288147",
    disease = "Asthma",
    path    = file.path(base, "01_raw_data", "03_asthma", "gse288147_pbmc"),
    # a = mild asthma | b = severe asthma | c = healthy control
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


# load one 10x sample
# read10x needs a directory with standardized filenames — use temp dir

load_10x_sample <- function(data_path, gsm_id, label, group, disease) {
  all_files  <- list.files(data_path, full.names = TRUE)
  gsm_files  <- all_files[grepl(gsm_id, all_files)]
  barcodes_f <- gsm_files[grepl("barcodes", gsm_files)]
  features_f <- gsm_files[grepl("features", gsm_files)]
  matrix_f   <- gsm_files[grepl("matrix",   gsm_files)]

  if (length(matrix_f) == 0) stop(paste("no matrix file for", gsm_id))

  tmp <- file.path(tempdir(), gsm_id)
  dir.create(tmp, showWarnings = FALSE)
  file.copy(barcodes_f, file.path(tmp, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(features_f, file.path(tmp, "features.tsv.gz"), overwrite = TRUE)
  file.copy(matrix_f,   file.path(tmp, "matrix.mtx.gz"),   overwrite = TRUE)

  mat <- Read10X(data.dir = tmp)

  # cite-seq returns a named list — extract gene expression only
  if (is.list(mat)) {
    log(paste("  cite-seq detected for", gsm_id, "— using gene expression"))
    mat <- mat[["Gene Expression"]]
  }

  obj <- CreateSeuratObject(counts = mat, project = label,
                            min.cells = 3, min.features = 200)
  obj$sample_id    <- gsm_id
  obj$sample_label <- label
  obj$group        <- group
  obj$disease      <- disease
  obj              <- RenameCells(obj, add.cell.id = label)

  log(paste(" ", gsm_id, "(", label, "):", ncol(obj), "cells x", nrow(obj), "genes"))
  return(obj)
}


# compute qc metrics

add_qc_metrics <- function(obj) {
  obj[["percent_mt"]]          <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent_ribo"]]        <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["log10_genes_per_umi"]] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  return(obj)
}


# adaptive mad-based filter per sample
# lun et al. bioconductor approach — more defensible than hardcoded thresholds

mad_filter <- function(x, n_mad = 3, direction = "both") {
  med <- median(x, na.rm = TRUE)
  m   <- mad(x,    na.rm = TRUE)
  if (direction == "both")  return(x >= (med - n_mad * m) & x <= (med + n_mad * m))
  if (direction == "lower") return(x >= (med - n_mad * m))
  if (direction == "upper") return(x <= (med + n_mad * m))
}


# save plot as 600dpi lzw-compressed tiff

save_tiff <- function(plot, path_no_ext, width, height) {
  tiff(paste0(path_no_ext, ".tiff"),
       width = width, height = height, units = "in",
       res = 600, compression = "lzw")
  print(plot)
  dev.off()
}


# generate qc plots for one object/stage, save tiffs, return plot list for pdf

make_qc_plots <- function(obj, label, stage, fig_dir) {

  prefix <- file.path(fig_dir, stage)

  p_violin <- VlnPlot(obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"),
    group.by = "sample_label", pt.size = 0, ncol = 4) &
    theme_classic(base_size = 10) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  save_tiff(p_violin, paste0(prefix, "_violin"), width = 18, height = 5)

  suppressWarnings({
    p_sc1 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA",
      group.by = "sample_label", pt.size = 0.2, raster = FALSE) +
      theme_classic(base_size = 10) + ggtitle("nCount vs nFeature")

    p_sc2 <- FeatureScatter(obj, "nCount_RNA", "percent_mt",
      group.by = "sample_label", pt.size = 0.2, raster = FALSE) +
      theme_classic(base_size = 10) +
      geom_hline(yintercept = 20, linetype = "dashed",
                 color = "red", linewidth = 0.6) +
      ggtitle("nCount vs % Mitochondrial")

    p_sc3 <- FeatureScatter(obj, "nFeature_RNA", "percent_mt",
      group.by = "sample_label", pt.size = 0.2, raster = FALSE) +
      theme_classic(base_size = 10) +
      geom_hline(yintercept = 20, linetype = "dashed",
                 color = "red", linewidth = 0.6) +
      ggtitle("nFeature vs % Mitochondrial")
  })

  p_scatter <- p_sc1 + p_sc2 + p_sc3
  save_tiff(p_scatter, paste0(prefix, "_scatter"), width = 18, height = 5)

  p_mt <- ggplot(obj@meta.data, aes(x = percent_mt, fill = sample_label)) +
    geom_histogram(bins = 60, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
    facet_wrap(~ sample_label, scales = "free_y") +
    theme_classic(base_size = 9) +
    labs(title = paste(label, stage, "— % mitochondrial distribution"),
         x = "% mitochondrial genes", y = "cell count") +
    theme(legend.position = "none")

  save_tiff(p_mt, paste0(prefix, "_mt_histogram"), width = 14, height = 8)

  p_complexity <- ggplot(obj@meta.data,
    aes(x = log10_genes_per_umi, fill = group)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
    theme_classic(base_size = 10) +
    labs(title = paste(label, stage, "— transcriptional complexity"),
         x = "log10(genes) / log10(UMIs)", y = "density")

  save_tiff(p_complexity, paste0(prefix, "_complexity"), width = 8, height = 5)

  log(paste(" ", stage, "tiffs saved to figures/"))
  return(list(violin     = p_violin,
              scatter    = p_scatter,
              mt_hist    = p_mt,
              complexity = p_complexity))
}


# doublet score and proportion plots

make_doublet_plots <- function(obj, label, fig_dir) {

  p_score <- ggplot(obj@meta.data,
    aes(x = sample_label, y = scDbl_score, fill = scDbl_class)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
    scale_fill_manual(values = c("singlet" = "#4DADE2", "doublet" = "#E84B35")) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(label, "— doublet scores"), x = NULL,
         y = "scDblFinder score", fill = "class")

  p_prop <- ggplot(obj@meta.data, aes(x = sample_label, fill = scDbl_class)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = c("singlet" = "#4DADE2", "doublet" = "#E84B35")) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(label, "— doublet proportions"), x = NULL,
         y = "proportion", fill = "class")

  save_tiff(p_score + p_prop,
            file.path(fig_dir, "03_doublets"), width = 14, height = 5)

  return(list(score = p_score, prop = p_prop))
}


# merge all plots into one multi-page pdf for quick review

merge_to_pdf <- function(disease, fig_dir, plots) {
  pdf_path <- file.path(fig_dir, paste0(disease, "_QC_overview.pdf"))
  pdf(pdf_path, width = 18, height = 8, useDingbats = FALSE)
  if (!is.null(plots$pre$violin))      print(plots$pre$violin)
  if (!is.null(plots$pre$scatter))     print(plots$pre$scatter)
  if (!is.null(plots$pre$mt_hist))     print(plots$pre$mt_hist)
  if (!is.null(plots$pre$complexity))  print(plots$pre$complexity)
  if (!is.null(plots$post$violin))     print(plots$post$violin)
  if (!is.null(plots$post$scatter))    print(plots$post$scatter)
  if (!is.null(plots$post$mt_hist))    print(plots$post$mt_hist)
  if (!is.null(plots$post$complexity)) print(plots$post$complexity)
  if (!is.null(plots$dbl$score))       print(plots$dbl$score + plots$dbl$prop)
  dev.off()
  log(paste("  overview pdf:", pdf_path))
}


# main qc function

run_qc <- function(ds_key) {

  ds   <- datasets[[ds_key]]
  dirs <- make_dirs(ds_key)

  log(paste("processing:", ds$name))

  # load and merge all samples
  obj_list <- list()
  for (gsm_id in names(ds$samples)) {
    meta <- ds$samples[[gsm_id]]
    tryCatch({
      obj_list[[gsm_id]] <- load_10x_sample(
        ds$path, gsm_id, meta$label, meta$group, ds$disease)
    }, error = function(e) log(paste("  error:", gsm_id, e$message)))
  }

  if (length(obj_list) == 0) stop("no samples loaded")

  merged <- if (length(obj_list) == 1) obj_list[[1]] else
              merge(obj_list[[1]], y = obj_list[-1], merge.data = FALSE)

  n_raw <- ncol(merged)
  log(paste("raw cells:", n_raw, "| genes:", nrow(merged)))

  merged            <- add_qc_metrics(merged)
  pre_median_genes  <- median(merged$nFeature_RNA)
  pre_median_counts <- median(merged$nCount_RNA)
  pre_median_mt     <- round(median(merged$percent_mt), 2)

  log(paste("pre-filter: median genes =", pre_median_genes,
            "| median %mt =", pre_median_mt))

  pre_plots <- make_qc_plots(merged, ds_key, "01_prefilter", dirs$figures)

  # per-sample adaptive mad filtering
  md           <- merged@meta.data
  pass_filters <- rep(TRUE, ncol(merged))
  thresh_rows  <- list()

  for (samp in unique(md$sample_label)) {
    idx             <- which(md$sample_label == samp)
    feat_pass       <- mad_filter(md$nFeature_RNA[idx],   n_mad = 3, direction = "both")
    count_pass      <- mad_filter(md$nCount_RNA[idx],     n_mad = 3, direction = "both")
    mt_hard         <- md$percent_mt[idx] < 20
    mt_mad          <- mad_filter(md$percent_mt[idx],     n_mad = 3, direction = "upper")
    complexity_pass <- md$log10_genes_per_umi[idx] > 0.8
    sample_pass     <- feat_pass & count_pass & mt_hard & mt_mad & complexity_pass
    pass_filters[idx] <- sample_pass

    thresh_rows[[samp]] <- data.frame(
      dataset              = ds$name,
      disease              = ds$disease,
      sample               = samp,
      group                = unique(md$group[idx]),
      n_cells_raw          = length(idx),
      n_cells_post_mad     = sum(sample_pass),
      pct_kept_after_mad   = round(100 * sum(sample_pass) / length(idx), 1),
      feat_min_kept        = round(min(md$nFeature_RNA[idx][feat_pass])),
      feat_max_kept        = round(max(md$nFeature_RNA[idx][feat_pass])),
      median_nfeature      = round(median(md$nFeature_RNA[idx]), 1),
      median_ncount        = round(median(md$nCount_RNA[idx]),   1),
      median_mt_pct        = round(median(md$percent_mt[idx]),   2),
      n_removed_low_feat   = sum(!feat_pass),
      n_removed_high_count = sum(!count_pass),
      n_removed_high_mt    = sum(!mt_hard | !mt_mad),
      n_removed_low_compl  = sum(!complexity_pass)
    )
  }

  thresh_df <- do.call(rbind, thresh_rows)
  rownames(thresh_df) <- NULL
  log("per-sample mad filter results:")
  print(thresh_df[, c("sample", "group", "n_cells_raw",
                      "n_cells_post_mad", "pct_kept_after_mad", "median_mt_pct")])

  write.csv(thresh_df,
            file.path(dirs$reports, "per_sample_filter_thresholds.csv"),
            row.names = FALSE)

  filtered  <- merged[, pass_filters]
  n_postmad <- ncol(filtered)
  log(paste("post-mad:", n_postmad, "cells (",
            round(100 * n_postmad / n_raw, 1), "% of raw)"))

  post_plots <- make_qc_plots(filtered, ds_key, "02_postfilter", dirs$figures)

  # doublet detection — join layers first (seurat v5 stores samples as layers)
  log("running scdblfinder...")
  filtered_joined <- JoinLayers(filtered)
  sce             <- as.SingleCellExperiment(filtered_joined)
  sce             <- scDblFinder(sce, samples = "sample_label",
                                 BPPARAM = BiocParallel::SerialParam())

  filtered$scDbl_class <- sce$scDblFinder.class
  filtered$scDbl_score <- sce$scDblFinder.score

  n_doublets <- sum(filtered$scDbl_class == "doublet")
  dbl_rate   <- round(100 * n_doublets / ncol(filtered), 2)
  log(paste("doublets detected:", n_doublets, "(", dbl_rate, "%)"))

  dbl_plots <- make_doublet_plots(filtered, ds_key, dirs$figures)

  # snapshot doublet counts per sample before removal
  dbl_counts <- filtered@meta.data %>%
    group_by(sample_label) %>%
    summarise(
      n_doublets       = sum(scDbl_class == "doublet"),
      doublet_rate_pct = round(100 * mean(scDbl_class == "doublet"), 2),
      .groups = "drop"
    )

  # remove doublets
  filtered  <- filtered[, filtered$scDbl_class == "singlet"]
  n_final   <- ncol(filtered)
  pct_final <- round(100 * n_final / n_raw, 1)
  log(paste("final cells:", n_final, "(", pct_final, "% of raw)"))

  merge_to_pdf(ds_key, dirs$figures,
               list(pre = pre_plots, post = post_plots, dbl = dbl_plots))

  # detailed per-sample cell count table tracking all stages
  pre_counts <- thresh_df[, c("sample", "group", "disease",
                               "n_cells_raw", "n_cells_post_mad")]
  colnames(pre_counts)[1] <- "sample_label"

  final_counts <- filtered@meta.data %>%
    group_by(sample_label, group) %>%
    summarise(n_cells_final = n(), .groups = "drop")

  cell_counts <- pre_counts %>%
    left_join(dbl_counts,   by = "sample_label") %>%
    left_join(final_counts, by = c("sample_label", "group")) %>%
    mutate(
      n_removed_mad        = n_cells_raw - n_cells_post_mad,
      n_removed_doublets   = n_doublets,
      pct_final_of_raw     = round(100 * n_cells_final / n_cells_raw, 1)
    ) %>%
    select(sample_label, group, disease,
           n_cells_raw,
           n_removed_mad,
           n_cells_post_mad,
           n_doublets,
           doublet_rate_pct,
           n_cells_final,
           pct_final_of_raw)

  write.csv(cell_counts,
            file.path(dirs$reports, "cell_counts_all_stages.csv"),
            row.names = FALSE)

  log("cell_counts_all_stages.csv saved")

  # dataset-level summary row for master report
  summary_row <- data.frame(
    dataset                 = ds$name,
    disease                 = ds$disease,
    n_samples               = length(ds$samples),
    n_cells_raw             = n_raw,
    n_cells_post_mad        = n_postmad,
    pct_retained_after_mad  = round(100 * n_postmad / n_raw, 1),
    n_doublets_detected     = n_doublets,
    doublet_rate_pct        = dbl_rate,
    n_cells_final           = n_final,
    pct_cells_final_of_raw  = pct_final,
    median_genes_raw        = pre_median_genes,
    median_genes_final      = median(filtered$nFeature_RNA),
    median_counts_raw       = pre_median_counts,
    median_counts_final     = median(filtered$nCount_RNA),
    median_mt_raw           = pre_median_mt,
    median_mt_final         = round(median(filtered$percent_mt), 2)
  )

  # save rds in disease root
  rds_path <- file.path(dirs$root, paste0(ds_key, "_qc_passed.rds"))
  saveRDS(filtered, rds_path)
  log(paste("rds saved:", rds_path))

  return(list(summary = summary_row, object = filtered))
}


# run pipeline — skip datasets with existing rds checkpoint

all_results <- list()

for (ds_key in names(datasets)) {
  rds_path <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key,
                        paste0(ds_key, "_qc_passed.rds"))
  if (file.exists(rds_path)) {
    log(paste("checkpoint found — loading:", ds_key))
    obj <- readRDS(rds_path)
    all_results[[ds_key]] <- list(
      summary = data.frame(
        dataset                = datasets[[ds_key]]$name,
        disease                = datasets[[ds_key]]$disease,
        n_samples              = length(datasets[[ds_key]]$samples),
        n_cells_raw            = NA_integer_,
        n_cells_post_mad       = NA_integer_,
        pct_retained_after_mad = NA_real_,
        n_doublets_detected    = NA_integer_,
        doublet_rate_pct       = NA_real_,
        n_cells_final          = ncol(obj),
        pct_cells_final_of_raw = NA_real_,
        median_genes_raw       = NA_real_,
        median_genes_final     = median(obj$nFeature_RNA),
        median_counts_raw      = NA_real_,
        median_counts_final    = median(obj$nCount_RNA),
        median_mt_raw          = NA_real_,
        median_mt_final        = round(median(obj$percent_mt), 2)
      ),
      object = obj
    )
  } else {
    tryCatch({
      all_results[[ds_key]] <- run_qc(ds_key)
    }, error = function(e) {
      log(paste("failed:", ds_key, "-", e$message))
    })
  }
}


# master summary across all three datasets

log("generating master qc summary")

summary_rows <- lapply(all_results, function(r) r$summary)
summary_rows <- summary_rows[!sapply(summary_rows, is.null)]

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  print(summary_df)
  write.csv(summary_df,
            file.path(log_root, "qc_summary_all_datasets.csv"),
            row.names = FALSE)
  log("qc_summary_all_datasets.csv saved")
}

# print final output structure to log
log("final output structure:")
for (ds_key in names(datasets)) {
  root <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key)
  if (dir.exists(root)) {
    log(paste0("02_scrna/pbmc/01_qc/", ds_key, "/"))
    for (f in list.files(root, recursive = TRUE))
      log(paste0("  ", f))
  }
}

log("qc pipeline complete")