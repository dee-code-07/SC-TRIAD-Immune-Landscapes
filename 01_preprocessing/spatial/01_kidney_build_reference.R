# sc-triad project
# script: 01_build_reference_v2.R
# purpose: build seurat kidney reference — v2 fixes harmony layer-split error
#
# fix log vs v1:
#   v1 error: "attempt to set an attribute on NULL" in IntegrateLayers
#   root cause: seurat v5 IntegrateLayers requires RNA assay to be split into
#               per-batch layers BEFORE normalization. v1 called NormalizeData
#               on a joined object, destroying the layer structure harmony needs.
#   fix: split RNA by tech → normalize (per-layer) → PCA → IntegrateLayers
#        → JoinLayers → UMAP/clustering
#
#   v1 warning: "data layer not found, using counts"
#   root cause: same issue — NormalizeData writes to "data" layer per split,
#               but if assay is not split, it writes to a single "data" layer
#               that DimPlot and other functions couldn't locate correctly.
#
# reference: satijalab.org/seurat/articles/seurat5_integration
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

options(future.globals.maxSize = 100 * 1024^3)
set.seed(42)

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
SCRNA   <- file.path(BASE, "01_raw_data/02_htn/gse211785_kidney/scrna")
OUT_DIR <- file.path(BASE, "04_spatial/kidney")

for (d in c("logs", "objects", "figures", "reports", "tables"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "logs/01_build_reference_v2.log")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

save_tiff <- function(plot, path, width, height) {
  tiff(path, width = width, height = height, units = "in",
       res = 300, compression = "lzw")
  print(plot)
  dev.off()
}

log("=== kidney reference builder v2 ===")
log(paste("Seurat:", packageVersion("Seurat")))
log(paste("harmony:", packageVersion("harmony")))

# ── step 1: load data ──────────────────────────────────────────────────────────
log("loading count matrix...")
counts <- readRDS(file.path(SCRNA, "GSE211785_counts.rds"))
log(paste("counts dim:", nrow(counts), "x", ncol(counts)))

log("loading metadata...")
meta <- read.csv(
  gzfile(file.path(SCRNA,
    "GSE211785_scRNA-seq_snRNA-seq_snATAC-seq_metadata.txt.gz")),
  row.names = 1
)
log(paste("metadata rows:", nrow(meta)))

# ── step 2: fix barcode mismatch ───────────────────────────────────────────────
log("fixing barcode mismatch (integer indices → full barcodes)...")

if (ncol(counts) != nrow(meta))
  stop(paste("dimension mismatch:", ncol(counts), "vs", nrow(meta)))

col_as_int   <- suppressWarnings(as.integer(colnames(counts)))
is_sequential <- all(!is.na(col_as_int)) &&
                  all(col_as_int == seq(0, ncol(counts) - 1))

if (is_sequential) {
  colnames(counts) <- rownames(meta)
  log("positional mapping applied successfully")
} else {
  common <- intersect(colnames(counts), rownames(meta))
  if (length(common) < 10000)
    stop(paste("only", length(common), "matching barcodes — check format"))
  counts <- counts[, common]
  meta   <- meta[common, ]
  log(paste("direct name match:", length(common), "cells"))
}

# verify spot check
log("spot checks after barcode fix:")
for (i in c(1, 1000, 50000, 200000)) {
  if (i <= ncol(counts))
    log(paste0("  pos ", i, ": ", colnames(counts)[i],
               " | group=", meta$group[i],
               " | cell_type=", meta$Cluster_Idents[i]))
}

# ── step 3: filter RNA-only (exclude SN_ATAC) ─────────────────────────────────
log("filtering to SC_RNA + SN_RNA only (removing SN_ATAC)...")
log("tech counts before filter:")
print(table(meta$tech))

rna_mask   <- meta$tech %in% c("SC_RNA", "SN_RNA")
counts_rna <- counts[, rna_mask]
meta_rna   <- meta[rna_mask, ]

log(paste("cells after RNA filter:", ncol(counts_rna)))
log("group distribution:")
print(table(meta_rna$group))

rm(counts); gc()

# ── step 4: create seurat object ───────────────────────────────────────────────
log("creating seurat object...")
obj <- CreateSeuratObject(
  counts       = counts_rna,
  meta.data    = meta_rna,
  min.cells    = 3,
  min.features = 200
)

obj$cell_type <- obj$Cluster_Idents
Idents(obj)   <- "cell_type"

log(paste("object created:", ncol(obj), "cells |", nrow(obj), "genes"))
rm(counts_rna); gc()

# ── step 5: QC metrics ────────────────────────────────────────────────────────
log("computing QC metrics...")
obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# soft filter: remove extreme outliers only
# we are deliberately conservative here — this is a published dataset
# the authors already QC'd; aggressive refiltering would distort cell type
# proportions we need for deconvolution

before <- ncol(obj)
obj <- subset(obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 10000 &
           percent_mt   < 25)
log(paste("cells after soft QC filter:", ncol(obj),
          "(removed:", before - ncol(obj), ")"))

# ── step 6: CRITICAL — split layers BEFORE normalization ─────────────────────
# this is the fix for the v1 harmony error
# seurat v5 IntegrateLayers requires the RNA assay to have one layer per batch
# splitting here means NormalizeData will write normalized values per-layer
# harmony then reads these separate layers to compute batch-corrected embeddings
# reference: https://satijalab.org/seurat/articles/seurat5_integration

log("splitting RNA layers by tech (required for Harmony in Seurat v5)...")
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$tech)
log(paste("RNA layers after split:",
          length(Layers(obj, assay = "RNA"))))

# ── step 7: normalize PER LAYER ───────────────────────────────────────────────
# NormalizeData on a split assay writes a "data" layer per batch
# this is what DimPlot and downstream functions need

log("running NormalizeData (per-layer, LogNormalize)...")
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor         = 10000,
                     verbose              = FALSE)

log("finding variable features (3000)...")
obj <- FindVariableFeatures(obj,
                            selection.method = "vst",
                            nfeatures        = 3000,
                            verbose          = FALSE)

log("scaling data...")
# ScaleData on split layers: only scale variable features
# do NOT regress percent_mt here — it can cause issues with split layers
# instead, pass it to Harmony as a covariate if needed
obj <- ScaleData(obj,
                 features = VariableFeatures(obj),
                 verbose  = FALSE)

# ── step 8: PCA ───────────────────────────────────────────────────────────────
log("running PCA (50 PCs)...")
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# PC selection
pct  <- obj[["pca"]]@stdev^2 / sum(obj[["pca"]]@stdev^2) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]
co2  <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1),
             decreasing = TRUE)[1] + 1
n_pcs <- min(max(min(co1, co2, na.rm = TRUE), 25), 50)
log(paste("PCs selected:", n_pcs, "| cumvar:", round(cumu[n_pcs], 1), "%"))

p_elbow <- ElbowPlot(obj, ndims = 50) +
  geom_vline(xintercept = n_pcs, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  theme_classic(base_size = 10) +
  labs(title = paste("Kidney reference: PCA elbow (selected:", n_pcs, "PCs)"))

save_tiff(p_elbow,
          file.path(OUT_DIR, "figures/01_pca_elbow.tiff"),
          width = 8, height = 5)

# ── step 9: Harmony integration ───────────────────────────────────────────────
# corrects SC_RNA vs SN_RNA batch effect
# SC_RNA: cytoplasmic + nuclear RNA → higher mitochondrial, more total counts
# SN_RNA: nuclear RNA only → lower counts, different gene detection pattern
# not correcting this would cause SC_RNA and SN_RNA cells to form separate
# clusters by modality rather than by cell type — ruining the reference

log("running Harmony integration (batch: tech)...")
log("  this corrects SC_RNA vs SN_RNA technical batch")

obj <- IntegrateLayers(
  object          = obj,
  method          = HarmonyIntegration,
  orig.reduction  = "pca",
  new.reduction   = "harmony",
  group.by.vars   = "tech",   # batch variable = tech modality
  verbose         = FALSE
)

log("Harmony complete")
log(paste("harmony dims:", ncol(obj[["harmony"]]@cell.embeddings)))

# ── step 10: UMAP and clustering ──────────────────────────────────────────────
log(paste("running UMAP on harmony (dims 1 to", n_pcs, ")..."))
obj <- RunUMAP(obj,
               reduction      = "harmony",
               dims           = 1:n_pcs,
               reduction.name = "umap",
               verbose        = FALSE)

obj <- FindNeighbors(obj,
                     reduction = "harmony",
                     dims      = 1:n_pcs,
                     verbose   = FALSE)

obj <- FindClusters(obj,
                    resolution = 0.5,
                    algorithm  = 1,
                    verbose    = FALSE)

log(paste("clusters found:", length(unique(obj$seurat_clusters))))

# ── step 11: join layers before saving ────────────────────────────────────────
# CRITICAL: join layers before saving and before RCTD
# RCTD reads the normalized count matrix directly from the "data" layer
# if layers are split, GetAssayData will only return the first batch's data
# joining merges all per-batch "data" layers into one complete matrix

log("joining RNA layers (required before RCTD and saving)...")
obj <- JoinLayers(obj, assay = "RNA")
log("layers joined")

# verify data layer exists now
test_layer <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data")[1:3, 1:3],
  error = function(e) NULL
)
if (is.null(test_layer)) {
  stop("data layer still missing after JoinLayers — investigate")
} else {
  log("data layer verified: exists and non-empty")
}

# ── step 12: UMAP figures ─────────────────────────────────────────────────────
log("generating UMAP figures...")

disease_colors <- c(
  Control = "#7CAE00",
  HKD     = "#C77CFF",
  DKD     = "#F8766D"
)

tech_colors <- c(
  SC_RNA = "#E41A1C",
  SN_RNA = "#377EB8"
)

# before vs after harmony: check tech mixing
p_before <- DimPlot(obj,
  reduction = "umap",
  group.by  = "tech",
  cols      = tech_colors,
  pt.size   = 0.05) +
  theme_classic(base_size = 9) +
  labs(title = "Post-Harmony: tech mixing (should be mixed)")

p_group <- DimPlot(obj,
  reduction = "umap",
  group.by  = "group",
  cols      = disease_colors,
  pt.size   = 0.05) +
  theme_classic(base_size = 9) +
  labs(title = "Disease group distribution")

save_tiff(p_before + p_group,
          file.path(OUT_DIR, "figures/02_umap_tech_group.tiff"),
          width = 18, height = 7)

# cell type UMAP — the key QC figure
# label concordance with published Cluster_Idents validates
# that our normalization/integration hasn't scrambled the biology

p_ct <- DimPlot(obj,
  reduction  = "umap",
  group.by   = "cell_type",
  label      = TRUE,
  label.size = 2.5,
  repel      = TRUE,
  pt.size    = 0.05) +
  theme_classic(base_size = 9) +
  labs(title = "Kidney reference: published cell type labels",
       subtitle = paste0(ncol(obj), " cells | 36 cell types | post-Harmony")) +
  theme(legend.text   = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

save_tiff(p_ct,
          file.path(OUT_DIR, "figures/03_umap_celltypes.tiff"),
          width = 16, height = 10)

# split by disease — shows cell type loss in DKD (PT loss) and
# immune infiltration in HKD (CD8T, NK gain)
p_split <- DimPlot(obj,
  reduction  = "umap",
  group.by   = "cell_type",
  split.by   = "group",
  pt.size    = 0.05,
  label      = TRUE,
  label.size = 2,
  repel      = TRUE) +
  theme_classic(base_size = 8) +
  labs(title = "Cell type distribution by disease group") +
  theme(legend.position = "none")

save_tiff(p_split,
          file.path(OUT_DIR, "figures/04_umap_split_disease.tiff"),
          width = 24, height = 8)

# ── step 13: validate published cell type labels with canonical markers ────────
# key kidney markers — if UMAP separates these correctly, reference is valid

kidney_markers <- list(
  "Proximal Tubule" = c("SLC34A1", "CUBN", "LRP2"),
  "Loop of Henle"   = c("UMOD", "SLC12A1"),
  "Distal Tubule"   = c("SLC12A3", "CALB1"),
  "Podocyte"        = c("NPHS1", "NPHS2", "PODXL"),
  "Endothelial"     = c("PECAM1", "CDH5", "ENG"),
  "Immune"          = c("PTPRC", "CD3E", "CD14", "NKG7"),
  "Fibroblast"      = c("COL1A1", "PDGFRA", "DCN")
)

all_markers <- unique(unlist(kidney_markers))
present     <- all_markers[all_markers %in% rownames(obj)]
missing     <- all_markers[!all_markers %in% rownames(obj)]

if (length(missing) > 0)
  log(paste("canonical markers absent:", paste(missing, collapse = ", ")))

if (length(present) >= 5) {
  p_dot <- DotPlot(obj,
    features = present,
    group.by = "cell_type",
    assay    = "RNA") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) +
    labs(title = "Canonical marker validation",
         subtitle = "confirms published cell type label accuracy")

  save_tiff(p_dot,
            file.path(OUT_DIR, "figures/05_canonical_marker_dotplot.tiff"),
            width = 16, height = 14)
}

# ── step 14: save reference objects ───────────────────────────────────────────
log("saving reference objects...")

# full reference
saveRDS(obj,
        file.path(OUT_DIR, "objects/kidney_ref_annotated.rds"))
log("full reference saved")

# HKD + Control — primary HTN spatial deconvolution reference
ref_hkd <- subset(obj, subset = group %in% c("HKD", "Control"))
saveRDS(ref_hkd,
        file.path(OUT_DIR, "objects/kidney_ref_HKD_Control.rds"))
log(paste("HKD+Control reference:", ncol(ref_hkd), "cells"))

# DKD + Control — T2D-kidney reference
ref_dkd <- subset(obj, subset = group %in% c("DKD", "Control"))
saveRDS(ref_dkd,
        file.path(OUT_DIR, "objects/kidney_ref_DKD_Control.rds"))
log(paste("DKD+Control reference:", ncol(ref_dkd), "cells"))

# ── step 15: cell composition report ──────────────────────────────────────────
composition <- obj@meta.data %>%
  group_by(group, cell_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 3)) %>%
  ungroup()

write.csv(composition,
          file.path(OUT_DIR, "reports/celltype_composition_by_group.csv"),
          row.names = FALSE)

# summary
summary_df <- data.frame(
  metric = c("Total cells", "Genes", "Cell types",
             "Control", "HKD", "DKD",
             "SC_RNA", "SN_RNA", "PCs used"),
  value  = c(ncol(obj), nrow(obj), length(unique(obj$cell_type)),
             sum(obj$group == "Control"),
             sum(obj$group == "HKD"),
             sum(obj$group == "DKD"),
             sum(obj$tech == "SC_RNA"),
             sum(obj$tech == "SN_RNA"),
             n_pcs)
)

write.csv(summary_df,
          file.path(OUT_DIR, "reports/reference_summary.csv"),
          row.names = FALSE)

log("reference build summary:")
print(summary_df)
log("=== reference builder v2 complete ===")
log("next: submit 02_rctd_deconvolution.sh")