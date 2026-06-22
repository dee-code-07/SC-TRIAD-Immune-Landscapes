# sc-triad project
# script: 02_norm_cluster_annotate.R
# purpose: normalization, clustering, and cell type annotation for pbmc datasets
# workflow: sctransform v2 -> pca -> harmony (t2d/asthma) -> umap ->
#           louvain clustering -> singler annotation -> marker validation
# references:
#   sctransform v2: choudhary & satija, genome biology 2022
#   harmony: korsunsky et al., nature methods 2019
#   singler: aran et al., nature immunology 2019
#   seurat v5 vignettes: satijalab.org/seurat/articles/
#   pc selection method: github.com/satijalab/seurat/issues/8288
# seurat v5 | r 4.4.3
# author: deeksha h | manipal academy of higher education

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(SingleR)
  library(celldex)
  library(clustree)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(42)

# increase memory limit for large matrix operations
options(future.globals.maxSize = 3e+09)

base     <- file.path(Sys.getenv("HOME"), "sc-triad")
qc_root  <- file.path(base, "02_scrna", "pbmc", "01_qc")
out_root <- file.path(base, "02_scrna", "pbmc", "02_normalization")

make_dirs <- function(disease) {
  root <- file.path(out_root, disease)
  dirs <- list(
    root    = root,
    figures = file.path(root, "figures"),
    reports = file.path(root, "reports")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  return(dirs)
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_root, "norm_pipeline.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log("sc-triad normalization pipeline started")
log(paste("seurat:", packageVersion("Seurat")))


# dataset config
# harmony_batch: metadata column for batch correction | NULL = skip harmony
# rationale for htn: n=1 per group, correction would destroy biological signal

datasets <- list(

  t2d = list(
    name          = "T2D_PBMC_GSE255566",
    disease       = "T2D",
    harmony_batch = "sample_label"
  ),

  htn = list(
    name          = "HTN_PBMC_GSE212953",
    disease       = "HTN",
    harmony_batch = NULL
  ),

  asthma = list(
    name          = "Asthma_PBMC_GSE288147",
    disease       = "Asthma",
    harmony_batch = "sample_label"
  )
)


# canonical pbmc marker genes for manual validation post-singler
# source: seurat pbmc vignette + monaco immune reference labels

canonical_markers <- list(
  "CD4 T"         = c("CD3E", "CD4", "IL7R", "CCR7"),
  "CD8 T"         = c("CD3E", "CD8A", "CD8B"),
  "NK"            = c("NCAM1", "NKG7", "GNLY", "FCGR3A"),
  "B cell"        = c("CD19", "MS4A1", "CD79A"),
  "CD14 Mono"     = c("CD14", "LYZ", "S100A8", "S100A9"),
  "CD16 Mono"     = c("FCGR3A", "MS4A7"),
  "DC"            = c("FCER1A", "CST3", "LILRA4"),
  "Platelet"      = c("PPBP", "PF4")
)


# save 600 dpi lzw-compressed tiff

save_tiff <- function(plot, path_no_ext, width, height) {
  tiff(paste0(path_no_ext, ".tiff"),
       width = width, height = height, units = "in",
       res = 600, compression = "lzw")
  print(plot)
  dev.off()
}


# pc selection using cumulative variance method
# source: github.com/satijalab/seurat/issues/8288
# co1: first pc where cumulative variance > 90% and individual variance < 5%
# co2: first pc where incremental drop in variance > 0.1%
# n_pcs: min of co1 and co2, clamped between 15 and 40

select_pcs <- function(pca_obj) {
  stdev   <- pca_obj@stdev
  pct     <- stdev / sum(stdev) * 100
  cumu    <- cumsum(pct)

  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
              decreasing = TRUE)[1] + 1

  n_pcs <- min(co1, co2, na.rm = TRUE)
  n_pcs <- max(n_pcs, 25)
  n_pcs <- min(n_pcs, 50)

  return(list(
    n_pcs   = n_pcs,
    cumvar  = round(cumu[n_pcs], 1),
    pct_tbl = data.frame(pc = seq_along(pct),
                         stdev = round(stdev, 4),
                         var_pct = round(pct, 3),
                         cum_var_pct = round(cumu, 3))
  ))
}


# main processing function

process_dataset <- function(ds_key) {

  ds   <- datasets[[ds_key]]
  dirs <- make_dirs(ds_key)

  log(paste("========== processing:", ds$name, "=========="))

  rds_in <- file.path(qc_root, ds_key, paste0(ds_key, "_qc_passed.rds"))
  if (!file.exists(rds_in)) stop(paste("qc rds not found:", rds_in))

  obj <- readRDS(rds_in)
  log(paste("loaded:", ncol(obj), "cells |", nrow(obj), "genes"))


  # qc rds files come out of merge() with layers already split per sample
  # must join first, then re-split so sctransform gets clean per-sample layers
  # reference: satijalab.org/seurat/articles/integration_introduction

  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_label)
  log(paste("rna layers split:", length(Layers(obj, assay = "RNA")), "layers"))


  # sctransform v2 normalization
  # vst.flavor = "v2" is default in seurat v5 (uses glmgampoi backend)
  # replaces normalizedata + findvariablefeatures + scaledata
  # vars.to.regress = "percent_mt" removes apoptotic/damaged cell signal
  # returns 3000 variable features by default (vs 2000 in lognormalize)
  # reference: choudhary & satija, genome biology 2022

  log("running sctransform v2...")
  obj <- SCTransform(obj,
                     vst.flavor      = "v2",
                     vars.to.regress = "percent_mt",
                     verbose         = FALSE)

  log(paste("sct complete | variable features:", length(VariableFeatures(obj))))


  # pca on sct-normalized data
  # sctransform benefits from higher npcs than lognormalize workflow
  # because technical variation is more effectively removed, so higher
  # pcs capture real biological signal rather than technical noise
  # reference: satijalab.org/seurat/articles/sctransform_vignette

  log("running pca (50 pcs)...")
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

  # select optimal number of pcs
  pc_res <- select_pcs(obj[["pca"]])
  n_pcs  <- pc_res$n_pcs
  log(paste("pcs selected:", n_pcs,
            "| cumulative variance:", pc_res$cumvar, "%"))

  write.csv(pc_res$pct_tbl,
            file.path(dirs$reports, "pca_variance_explained.csv"),
            row.names = FALSE)

  # elbow plot
  p_elbow <- ElbowPlot(obj, ndims = 50) +
    geom_vline(xintercept = n_pcs, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    theme_classic(base_size = 11) +
    labs(title = paste(ds_key, "- PCA elbow (selected:", n_pcs, "PCs)"))

  save_tiff(p_elbow, file.path(dirs$figures, "01_pca_elbow"), width = 8, height = 5)


  # harmony batch correction (t2d and asthma only)
  # integrates pca embeddings across samples to remove donor-level batch effects
  # correct syntax for sct+harmony in seurat v5:
  #   IntegrateLayers(..., assay = "SCT") - NOT normalization.method = "SCT"
  #   normalization.method = "SCT" is only for cca/rpca integration
  # reference: satijalab.org/seurat/reference/harmonyintegration
  # reference: github.com/satijalab/seurat/issues/8288

  if (!is.null(ds$harmony_batch)) {
    log(paste("running harmony batch correction by:", ds$harmony_batch))
    obj <- IntegrateLayers(
      object         = obj,
      method         = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction  = "harmony",
      assay          = "SCT",
      verbose        = FALSE
    )
    reduction_use <- "harmony"
    log("harmony complete")
  } else {
    log("harmony skipped - htn has n=1 per group, correction = signal loss")
    reduction_use <- "pca"
  }


  # umap using selected reduction and pcs

  log(paste("running umap | reduction:", reduction_use, "| dims: 1 to", n_pcs))
  obj <- RunUMAP(obj,
                 reduction = reduction_use,
                 dims      = 1:n_pcs,
                 verbose   = FALSE)


  # knn graph for clustering

  obj <- FindNeighbors(obj,
                       reduction = reduction_use,
                       dims      = 1:n_pcs,
                       verbose   = FALSE)


  # louvain clustering at multiple resolutions for clustree analysis
  # louvain (algorithm=1) is the seurat default and always available
  # leiden (algorithm=4) requires the leiden package - not used here

  log("clustering at resolutions: 0.2 0.4 0.6 0.8 1.0")
  resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)

  for (res in resolutions) {
    obj <- FindClusters(obj,
                        algorithm        = 1,
                        resolution       = res,
                        cluster.name     = paste0("SCT_snn_res.", res),
                        verbose          = FALSE)
  }

  # clustree: visualize cluster stability across resolutions
  # helps select resolution where clusters stop splitting meaningfully
  # reference: zappia & oshlack, gigascience 2018

  p_clustree <- clustree(obj, prefix = "SCT_snn_res.") +
    theme(legend.position = "right") +
    labs(title = paste(ds_key, "- cluster tree"))

  save_tiff(p_clustree, file.path(dirs$figures, "02_clustree"),
            width = 12, height = 10)

  # record cluster counts per resolution
  cluster_counts <- sapply(paste0("SCT_snn_res.", resolutions),
                           function(col) length(unique(obj[[col]][, 1])))

  cluster_df <- data.frame(resolution = resolutions,
                           n_clusters = cluster_counts)
  log("clusters per resolution:")
  print(cluster_df)
  write.csv(cluster_df,
            file.path(dirs$reports, "cluster_counts_by_resolution.csv"),
            row.names = FALSE)

  # select optimal resolution: first where gain drops to fewer than 2 new clusters
  gains      <- diff(cluster_counts)
  stable_idx <- which(gains < 2)[1]
  opt_res    <- if (!is.na(stable_idx)) resolutions[stable_idx] else 0.4
  n_clusters <- cluster_counts[resolutions == opt_res]

  log(paste("optimal resolution:", opt_res, "| clusters:", n_clusters))

  Idents(obj) <- paste0("SCT_snn_res.", opt_res)
  obj$seurat_clusters <- Idents(obj)


  # umap plots

  p_cluster <- DimPlot(obj, reduction = "umap",
                       group.by = "seurat_clusters",
                       label = TRUE, label.size = 3.5,
                       repel = TRUE, pt.size = 0.3,
                       raster = FALSE) +
    theme_classic(base_size = 11) +
    labs(title = paste(ds_key, "clusters | res", opt_res))

  p_group <- DimPlot(obj, reduction = "umap",
                     group.by = "group",
                     pt.size = 0.3, raster = FALSE) +
    theme_classic(base_size = 11) +
    labs(title = paste(ds_key, "by group"))

  p_sample <- DimPlot(obj, reduction = "umap",
                      group.by = "sample_label",
                      pt.size = 0.3, raster = FALSE) +
    theme_classic(base_size = 11) +
    labs(title = paste(ds_key, "by sample"))

  save_tiff(p_cluster,
            file.path(dirs$figures, "03_umap_clusters"), width = 10, height = 8)
  save_tiff(p_group,
            file.path(dirs$figures, "04_umap_group"),    width = 10, height = 8)
  save_tiff(p_sample,
            file.path(dirs$figures, "05_umap_sample"),   width = 10, height = 8)
  save_tiff(p_cluster + p_group,
            file.path(dirs$figures, "06_umap_panel"),    width = 18, height = 8)


  # singler cell type annotation
  # reference: monaco immune data - best resolution for pbmc subtypes
  # correct seurat v5 approach: extract sct data layer after joining layers
  # reference: aran et al., nature immunology 2019

  log("running singler annotation (monaco immune reference)...")

  ref <- celldex::MonacoImmuneData()

  # extract log-normalized rna counts for singler input
  # using rna assay lognorm rather than sct data layer - more robust in v5
  # sct data layer can be degenerate after split/join/harmony operations
  # rna lognorm is the standard singler input used in the singler vignette
  obj_j    <- JoinLayers(obj, assay = "RNA")
  obj_j    <- NormalizeData(obj_j, assay = "RNA", verbose = FALSE)
  norm_mat <- GetAssayData(obj_j, assay = "RNA", layer = "data")

  # run singler twice:
  # pass 1: label.main (8 broad types) - robust, used for majority vote
  # pass 2: label.fine (29 subtypes) - detailed, saved for reference
  # rationale: fine labels require high confidence to distinguish subtypes,
  # often collapsing to one dominant type. main labels are more reliable
  # for initial annotation and manuscript figures.
  # reference: aran et al., nature immunology 2019

  singler_main <- SingleR(
    test      = norm_mat,
    ref       = ref,
    labels    = ref$label.main,
    de.method = "wilcox"
  )

  singler_fine <- SingleR(
    test      = norm_mat,
    ref       = ref,
    labels    = ref$label.fine,
    de.method = "wilcox"
  )

  obj$singler_label_main   <- singler_main$labels
  obj$singler_label_fine   <- singler_fine$labels
  obj$singler_pruned_main  <- singler_main$pruned.labels
  obj$singler_pruned_fine  <- singler_fine$pruned.labels

  # use main labels as primary annotation
  obj$singler_label  <- obj$singler_label_main
  obj$singler_pruned <- obj$singler_pruned_main

  n_types <- length(unique(na.omit(obj$singler_label_main)))
  n_fine  <- length(unique(na.omit(obj$singler_label_fine)))
  log(paste("singler complete | main:", n_types, "types | fine:", n_fine, "subtypes"))

  # per-cluster annotation summary for main labels
  annot_main <- as.data.frame(table(
    cluster       = obj$seurat_clusters,
    singler_label = obj$singler_label_main
  )) %>%
    filter(Freq > 0) %>%
    arrange(cluster, desc(Freq))

  write.csv(annot_main,
            file.path(dirs$reports, "singler_main_labels_per_cluster.csv"),
            row.names = FALSE)

  # per-cluster annotation summary for fine labels (for reference)
  annot_fine <- as.data.frame(table(
    cluster       = obj$seurat_clusters,
    singler_label = obj$singler_label_fine
  )) %>%
    filter(Freq > 0) %>%
    arrange(cluster, desc(Freq))

  write.csv(annot_fine,
            file.path(dirs$reports, "singler_fine_labels_per_cluster.csv"),
            row.names = FALSE)

  # majority vote on main labels
  cluster_celltypes <- annot_main %>%
    group_by(cluster) %>%
    slice_max(Freq, n = 1, with_ties = FALSE) %>%
    rename(cell_type = singler_label) %>%
    select(cluster, cell_type, n_cells = Freq)

  write.csv(cluster_celltypes,
            file.path(dirs$reports, "cluster_celltype_final.csv"),
            row.names = FALSE)

  log("cluster -> cell type assignments:")
  print(as.data.frame(cluster_celltypes))

  # add majority-vote label to metadata
  cmap        <- setNames(cluster_celltypes$cell_type,
                          as.character(cluster_celltypes$cluster))
  # unname() required - cmap has cluster ids as names, not cell barcodes
  # without unname() seurat tries to match by name and throws "no cell overlap"
  obj$cell_type <- unname(cmap[as.character(obj$seurat_clusters)])


  # umap by cell type
  p_celltype <- DimPlot(obj, reduction = "umap",
                        group.by = "cell_type",
                        label = TRUE, label.size = 2.8,
                        repel = TRUE, pt.size = 0.3,
                        raster = FALSE) +
    theme_classic(base_size = 9) +
    labs(title = paste(ds_key, "- SingleR cell types"))

  save_tiff(p_celltype,
            file.path(dirs$figures, "07_umap_celltype"), width = 12, height = 8)


  # canonical marker validation - dotplot + featureplot
  # confirms singler annotations using known lineage markers

  all_markers     <- unique(unlist(canonical_markers))
  markers_present <- all_markers[all_markers %in% rownames(obj)]
  missing         <- all_markers[!all_markers %in% rownames(obj)]

  if (length(missing) > 0)
    log(paste("markers absent from dataset:", paste(missing, collapse = ", ")))

  if (length(markers_present) > 0) {

    p_dot <- DotPlot(obj,
                     features = markers_present,
                     group.by = "seurat_clusters",
                     assay    = "SCT") +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
      labs(title = paste(ds_key, "- canonical markers by cluster"))

    save_tiff(p_dot,
              file.path(dirs$figures, "08_dotplot_markers"), width = 16, height = 8)

    # feature plots for primary lineage markers
    lineage_markers <- c("CD3E", "CD4", "CD8A", "CD19",
                         "CD14", "FCGR3A", "NCAM1", "NKG7")
    lineage_present <- lineage_markers[lineage_markers %in% rownames(obj)]

    if (length(lineage_present) > 0) {
      p_feat <- FeaturePlot(obj,
                            features  = lineage_present,
                            reduction = "umap",
                            pt.size   = 0.2,
                            ncol      = 4,
                            raster    = FALSE) &
        theme_classic(base_size = 8)

      save_tiff(p_feat,
                file.path(dirs$figures, "09_featureplot_lineage"),
                width = 18, height = 9)
    }
  }


  # cell type composition by group
  composition <- obj@meta.data %>%
    group_by(group, cell_type) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
    ungroup() %>%
    arrange(group, desc(n_cells))

  write.csv(composition,
            file.path(dirs$reports, "celltype_composition_by_group.csv"),
            row.names = FALSE)

  p_comp <- ggplot(composition, aes(x = group, y = pct, fill = cell_type)) +
    geom_bar(stat = "identity") +
    theme_classic(base_size = 10) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.text  = element_text(size = 7),
          legend.key.size = unit(0.4, "cm")) +
    labs(title = paste(ds_key, "- cell type composition by group"),
         x = NULL, y = "% of cells", fill = "cell type")

  save_tiff(p_comp,
            file.path(dirs$figures, "10_celltype_composition"),
            width = 10, height = 7)


  # overview pdf - all key plots in one file for quick review

  pdf(file.path(dirs$figures, paste0(ds_key, "_processing_overview.pdf")),
      width = 14, height = 8, useDingbats = FALSE)
  print(p_elbow)
  print(p_clustree)
  print(p_cluster + p_group)
  print(p_celltype)
  print(p_dot)
  print(p_comp)
  dev.off()

  log("overview pdf saved")


  # dataset-level summary row

  summary_row <- data.frame(
    dataset           = ds$name,
    disease           = ds$disease,
    n_cells           = ncol(obj),
    n_variable_genes  = length(VariableFeatures(obj)),
    n_pcs_used        = n_pcs,
    pca_cumvar_pct    = pc_res$cumvar,
    harmony_applied   = !is.null(ds$harmony_batch),
    reduction_used    = reduction_use,
    optimal_res       = opt_res,
    n_clusters        = n_clusters,
    n_singler_main_types = n_types,
    n_singler_fine_types = n_fine
  )


  # save annotated rds checkpoint
  # note: we call joinlayers before saving so the object is ready for
  # downstream deg analysis without needing to rejoin layers

  obj <- JoinLayers(obj, assay = "RNA")
  rds_out <- file.path(dirs$root, paste0(ds_key, "_annotated.rds"))
  saveRDS(obj, rds_out)
  log(paste("annotated rds saved:", rds_out))

  return(summary_row)
}


# run for all datasets - skip if annotated rds exists

all_summaries <- list()

for (ds_key in names(datasets)) {
  rds_out <- file.path(out_root, ds_key, paste0(ds_key, "_annotated.rds"))
  if (file.exists(rds_out)) {
    log(paste("checkpoint found, skipping:", ds_key))
    next
  }
  tryCatch({
    all_summaries[[ds_key]] <- process_dataset(ds_key)
  }, error = function(e) {
    log(paste("failed:", ds_key, "-", e$message))
  })
}


# master summary

if (length(all_summaries) > 0) {
  summary_df <- do.call(rbind, all_summaries)
  rownames(summary_df) <- NULL
  print(summary_df)
  write.csv(summary_df,
            file.path(out_root, "normalization_summary.csv"),
            row.names = FALSE)
  log("normalization_summary.csv saved")
}

# print output structure
log("output structure:")
for (ds_key in names(datasets)) {
  d <- file.path(out_root, ds_key)
  if (dir.exists(d)) {
    log(paste0("02_scrna/pbmc/02_normalization/", ds_key, "/"))
    for (f in list.files(d, recursive = TRUE)) log(paste0("  ", f))
  }
}

log("normalization pipeline complete")