# sc-triad project
# script: 03_integration.R
# purpose: integrate t2d, htn, asthma pbmc datasets into a single object
#          for cross-disease analysis of shared cell states
#
# workflow (seurat v5 + harmony):
#   merge 3 annotated rds objects
#   join rna layers -> split by sample_label
#   sctransform v2 on merged object (new joint model)
#   pca -> harmony (correct by disease + sample_label)
#   umap -> clustering -> figures
#
# rationale for re-running sctransform:
#   per-disease sct models are not transferable to merged object
#   joint sct learns a single regularization model across all cells
#   this is the correct approach for cross-dataset integration
#   reference: choudhary & satija, genome biology 2022
#
# rationale for harmony group.by.vars = c("disease", "sample_label"):
#   disease: corrects dataset-level batch (different geo studies, protocols)
#   sample_label: corrects donor-level technical variation
#   harmony aligns shared cell types across batches while preserving
#   transcriptional differences between disease states within cell types
#   reference: korsunsky et al., nature methods 2019
#
# output: 02_scrna/pbmc/03_integration/
# author: deeksha h | manipal academy of higher education

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(clustree)
})

options(future.globals.maxSize = 20 * 1024^3)
set.seed(42)

base    <- file.path(Sys.getenv("HOME"), "sc-triad")
norm    <- file.path(base, "02_scrna", "pbmc", "02_normalization")
out_dir <- file.path(base, "02_scrna", "pbmc", "03_integration")
fig_dir <- file.path(out_dir, "figures")
rep_dir <- file.path(out_dir, "reports")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "integration.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

save_tiff <- function(plot, path_no_ext, width, height) {
  tiff(paste0(path_no_ext, ".tiff"),
       width  = width, height = height, units = "in",
       res    = 600, compression = "lzw")
  print(plot)
  dev.off()
}

log("sc-triad integration pipeline started")
log(paste("seurat:", packageVersion("Seurat")))
log(paste("harmony:", packageVersion("harmony")))


# load individual annotated rds objects

log("loading per-disease annotated objects...")

t2d    <- readRDS(file.path(norm, "t2d",    "t2d_annotated.rds"))
htn    <- readRDS(file.path(norm, "htn",    "htn_annotated.rds"))
asthma <- readRDS(file.path(norm, "asthma", "asthma_annotated.rds"))

log(paste("t2d cells:",    ncol(t2d)))
log(paste("htn cells:",    ncol(htn)))
log(paste("asthma cells:", ncol(asthma)))

# add disease label to metadata (will be used for harmony batch correction)
t2d$disease    <- "T2D"
htn$disease    <- "HTN"
asthma$disease <- "Asthma"

# keep only essential metadata to reduce object size after merge
keep_cols <- c("sample_label", "disease", "group",
               "cell_type", "cell_type_manual", "cell_type_singler",
               "singler_label_main", "singler_label_fine",
               "nCount_RNA", "nFeature_RNA", "percent_mt",
               "seurat_clusters")

for (col in keep_cols) {
  if (!col %in% colnames(t2d@meta.data))    t2d[[col]]    <- NA
  if (!col %in% colnames(htn@meta.data))    htn[[col]]    <- NA
  if (!col %in% colnames(asthma@meta.data)) asthma[[col]] <- NA
}


# drop sct assay from each object before merging
# per-disease sct models have different gene sets (different variable features
# selected per dataset), so merge.SCTAssay fails with "number of columns must match"
# we re-run sctransform on the merged object anyway, so per-disease sct is not needed
# only the RNA assay (raw counts) is needed for the joint sct model

log("removing per-disease SCT assays before merge...")
DefaultAssay(t2d)    <- "RNA"
DefaultAssay(htn)    <- "RNA"
DefaultAssay(asthma) <- "RNA"
if ("SCT" %in% Assays(t2d))    t2d[["SCT"]]    <- NULL
if ("SCT" %in% Assays(htn))    htn[["SCT"]]    <- NULL
if ("SCT" %in% Assays(asthma)) asthma[["SCT"]] <- NULL
log("SCT assays removed - RNA assay retained for joint SCTransform")

# merge all three objects
# seurat v5: merge keeps rna layers split per sample automatically
# reference: satijalab.org/seurat/articles/essential_commands

log("merging t2d + htn + asthma...")

merged <- merge(
  x            = t2d,
  y            = list(htn, asthma),
  add.cell.ids = c("T2D", "HTN", "Asthma"),
  merge.data   = FALSE
)

log(paste("merged object: cells =", ncol(merged),
          "| genes =", nrow(merged)))
log(paste("rna layers after merge:", length(Layers(merged, assay = "RNA"))))

# remove per-disease objects to free memory
rm(t2d, htn, asthma)
gc()


# join then split rna layers by sample_label for sctransform
# sctransform learns a separate regularization model per sample
# this is required for correct variance stabilization in multi-sample objects

log("joining then splitting rna layers by sample_label...")
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample_label)
log(paste("rna layers after split:", length(Layers(merged, assay = "RNA"))))

# sct assay already removed before merge - nothing to clean up here


# sctransform v2 on merged object
# joint model learns regularization across all cells and samples
# vars.to.regress: mitochondrial percentage only
# n_cells = 5000 per sample for faster glmgampoi fitting on large object

log("running sctransform v2 on merged object...")
log(paste("total cells:", ncol(merged), "- this will take ~20-30 min"))

merged <- SCTransform(
  merged,
  vst.flavor       = "v2",
  vars.to.regress  = "percent_mt",
  variable.features.n = 3000,
  verbose          = FALSE
)

log(paste("sct complete | variable features:", length(VariableFeatures(merged))))


# pca on sct variable features

log("running pca (50 pcs)...")
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

# variance explained
pct  <- merged[["pca"]]@stdev^2 / sum(merged[["pca"]]@stdev^2) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]
co2  <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1),
             decreasing = TRUE)[1] + 1
n_pcs <- min(co1, co2, na.rm = TRUE)
n_pcs <- max(n_pcs, 25)
n_pcs <- min(n_pcs, 50)

log(paste("pcs selected:", n_pcs, "| cumulative variance:", round(cumu[n_pcs], 1), "%"))

write.csv(
  data.frame(pc = seq_along(pct), stdev = merged[["pca"]]@stdev,
             var_pct = round(pct, 3), cumvar = round(cumu, 3)),
  file.path(rep_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

p_elbow <- ElbowPlot(merged, ndims = 50) +
  geom_vline(xintercept = n_pcs, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme_classic(base_size = 10) +
  labs(title = paste("integrated - PCA elbow (selected:", n_pcs, "PCs)"))

save_tiff(p_elbow, file.path(fig_dir, "01_pca_elbow"), width = 8, height = 5)


# umap before integration (for comparison figure)

log("running unintegrated umap for comparison...")
merged <- RunUMAP(merged,
                  reduction  = "pca",
                  dims       = 1:n_pcs,
                  reduction.name = "umap.unintegrated",
                  verbose    = FALSE)


# harmony integration
# group.by.vars:
#   disease: removes cross-study technical batch (different geo datasets)
#   sample_label: removes donor-level technical variation
# theta: default (2) - moderate correction strength
# reference: IntegrateLayers with HarmonyIntegration, assay = "SCT"
# from satijalab docs: for sct data, pass assay = "SCT" to IntegrateLayers

log("running harmony integration (batch: disease + sample_label)...")

merged <- IntegrateLayers(
  object        = merged,
  method        = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  assay         = "SCT",
  group.by.vars = c("disease", "sample_label"),
  verbose       = FALSE
)

log("harmony integration complete")
log(paste("harmony dims available:", ncol(merged[["harmony"]]@cell.embeddings)))


# umap on harmony embedding

log(paste("running integrated umap | dims: 1 to", n_pcs))
merged <- RunUMAP(merged,
                  reduction  = "harmony",
                  dims       = 1:n_pcs,
                  reduction.name = "umap",
                  verbose    = FALSE)


# clustering on harmony embedding
# test multiple resolutions for clustree stability analysis

log("clustering at resolutions: 0.2 0.4 0.6 0.8 1.0...")
merged <- FindNeighbors(merged,
                        reduction = "harmony",
                        dims      = 1:n_pcs,
                        verbose   = FALSE)

for (res in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
  merged <- FindClusters(merged, resolution = res,
                         algorithm = 1, verbose = FALSE)
}

# cluster counts per resolution
res_cols <- grep("SCT_snn_res", colnames(merged@meta.data), value = TRUE)
cluster_counts <- sapply(res_cols, function(x) length(unique(merged@meta.data[[x]])))
res_df <- data.frame(
  resolution = as.numeric(gsub("SCT_snn_res.", "", res_cols)),
  n_clusters = cluster_counts
)
log("clusters per resolution:")
print(res_df)
write.csv(res_df, file.path(rep_dir, "cluster_counts_by_resolution.csv"),
          row.names = FALSE)

# optimal resolution: first where gain < 2 new clusters vs previous
n_vec   <- res_df$n_clusters
gains   <- c(n_vec[1], diff(n_vec))
opt_idx <- which(gains < 2)[1]
if (is.na(opt_idx)) opt_idx <- which.min(gains)
opt_res     <- res_df$resolution[opt_idx]
n_clusters  <- res_df$n_clusters[opt_idx]

log(paste("optimal resolution:", opt_res, "| clusters:", n_clusters))
Idents(merged) <- paste0("SCT_snn_res.", opt_res)
merged$seurat_clusters_integrated <- Idents(merged)


# figures

log("generating figures...")

celltype_colors <- c(
  "Naive CD4 T"    = "#4E9DC4",
  "Memory CD4 T"   = "#2171B5",
  "CD4 T"          = "#08519C",
  "CD8 T"          = "#6BAED6",
  "NK"             = "#FB6A4A",
  "CD16 NK"        = "#EF3B2C",
  "B cell"         = "#74C476",
  "CD14 Monocyte"  = "#FD8D3C",
  "CD16 Monocyte"  = "#D94801",
  "DC"             = "#9E9AC8",
  "Megakaryocyte"  = "#E7298A",
  "Basophil"       = "#A1D99B",
  "Neutrophil"     = "#FDAE6B",
  "Other"          = "#BDBDBD"
)

disease_colors <- c(
  "T2D"    = "#E41A1C",
  "HTN"    = "#377EB8",
  "Asthma" = "#4DAF4A"
)

# before vs after integration comparison
p_before <- DimPlot(merged, reduction = "umap.unintegrated",
                    group.by = "disease",
                    cols = disease_colors,
                    pt.size = 0.1, raster = FALSE) +
  theme_classic(base_size = 9) +
  labs(title = "Before integration") +
  theme(legend.text = element_text(size = 8))

p_after <- DimPlot(merged, reduction = "umap",
                   group.by = "disease",
                   cols = disease_colors,
                   pt.size = 0.1, raster = FALSE) +
  theme_classic(base_size = 9) +
  labs(title = "After Harmony integration") +
  theme(legend.text = element_text(size = 8))

save_tiff(p_before + p_after,
          file.path(fig_dir, "02_before_after_integration"),
          width = 20, height = 8)

# umap by disease
p_disease <- DimPlot(merged, reduction = "umap",
                     group.by = "disease",
                     cols = disease_colors,
                     pt.size = 0.1, raster = FALSE) +
  theme_classic(base_size = 10) +
  labs(title = "Integrated PBMC - by disease")

save_tiff(p_disease, file.path(fig_dir, "03_umap_disease"),
          width = 10, height = 8)

# umap by cell type
colors_use <- celltype_colors[names(celltype_colors) %in% unique(merged$cell_type)]

p_celltype <- DimPlot(merged, reduction = "umap",
                      group.by = "cell_type",
                      cols = colors_use,
                      label = TRUE, label.size = 3,
                      repel = TRUE, pt.size = 0.1,
                      raster = FALSE) +
  theme_classic(base_size = 10) +
  labs(title = "Integrated PBMC - cell types")

save_tiff(p_celltype, file.path(fig_dir, "04_umap_celltype"),
          width = 12, height = 8)

# umap by sample
p_sample <- DimPlot(merged, reduction = "umap",
                    group.by = "sample_label",
                    pt.size = 0.1, raster = FALSE) +
  theme_classic(base_size = 9) +
  labs(title = "Integrated PBMC - by sample") +
  theme(legend.text = element_text(size = 6))

save_tiff(p_sample, file.path(fig_dir, "05_umap_sample"),
          width = 12, height = 8)

# umap by integrated cluster
p_cluster <- DimPlot(merged, reduction = "umap",
                     group.by = "seurat_clusters_integrated",
                     label = TRUE, label.size = 3,
                     pt.size = 0.1, raster = FALSE) +
  theme_classic(base_size = 10) +
  labs(title = paste("Integrated clusters | res", opt_res))

save_tiff(p_cluster, file.path(fig_dir, "06_umap_clusters"),
          width = 10, height = 8)

# split umap by disease (key figure for midsem presentation)
p_split <- DimPlot(merged, reduction = "umap",
                   group.by = "cell_type",
                   split.by = "disease",
                   cols = colors_use,
                   pt.size = 0.1, raster = FALSE,
                   label = TRUE, label.size = 2.5,
                   repel = TRUE) +
  theme_classic(base_size = 9) +
  labs(title = "Cell type composition by disease") +
  theme(legend.text = element_text(size = 7))

save_tiff(p_split, file.path(fig_dir, "07_umap_split_by_disease"),
          width = 24, height = 8)

# clustree
p_clustree <- clustree(merged, prefix = "SCT_snn_res.") +
  labs(title = "Integration - cluster tree") +
  theme(legend.text = element_text(size = 7))

save_tiff(p_clustree, file.path(fig_dir, "08_clustree"),
          width = 18, height = 12)

# cell type composition across diseases
composition <- merged@meta.data %>%
  group_by(disease, cell_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(disease) %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  ungroup()

write.csv(composition,
          file.path(rep_dir, "celltype_composition_by_disease.csv"),
          row.names = FALSE)

p_comp <- ggplot(composition,
                 aes(x = disease, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors, na.value = "#BDBDBD") +
  theme_classic(base_size = 11) +
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(title = "Cell type composition across diseases",
       x = NULL, y = "% of cells", fill = "cell type")

save_tiff(p_comp, file.path(fig_dir, "09_celltype_composition_by_disease"),
          width = 10, height = 7)

# overview pdf
log("saving overview pdf...")
pdf(file.path(fig_dir, "integration_overview.pdf"),
    width = 16, height = 9, useDingbats = FALSE)
print(p_before + p_after)
print(p_disease)
print(p_celltype)
print(p_split)
print(p_comp)
print(p_clustree)
dev.off()
log("overview pdf saved")


# integration quality metrics
# lisi score (local inverse simpson index) requires lisi package
# kbet requires kbet package
# computing basic mixing metrics instead: per-cluster disease entropy

log("computing integration quality metrics...")

# disease mixing entropy per integrated cluster
entropy_df <- merged@meta.data %>%
  group_by(seurat_clusters_integrated, disease) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters_integrated) %>%
  mutate(
    total    = sum(n),
    prop     = n / total,
    entropy  = -sum(prop * log2(prop + 1e-10))
  ) %>%
  slice(1) %>%
  select(cluster = seurat_clusters_integrated, total_cells = total, entropy) %>%
  arrange(desc(entropy))

write.csv(entropy_df,
          file.path(rep_dir, "integration_disease_mixing_entropy.csv"),
          row.names = FALSE)

log(paste("mean disease mixing entropy per cluster:",
          round(mean(entropy_df$entropy), 3),
          "(max possible:", round(log2(3), 3), ")"))


# save integrated rds
rds_out <- file.path(out_dir, "triad_integrated.rds")
log(paste("saving integrated rds:", rds_out))
saveRDS(merged, rds_out)
log(paste("integrated rds saved | cells:", ncol(merged),
          "| clusters:", n_clusters))


# output summary
log("output structure:")
all_files <- c(
  list.files(fig_dir, full.names = FALSE),
  list.files(rep_dir, full.names = FALSE),
  "triad_integrated.rds"
)
for (f in sort(all_files)) log(paste(" ", f))

log("integration pipeline complete")