# sc-triad project
# purpose: reproduce supplementary figures for section 2.1 in pdf format
#          rasterize umaps using ggrastr for publication quality
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
datasets <- c("t2d", "htn", "asthma")

# Exact theme from 01_umap_celltype.R
base_theme <- theme_classic(base_size = 12) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

celltype_colors <- c(
  "Naive CD4 T"   = "#4E79A7",
  "Memory CD4 T"  = "#1F4E79",
  "CD8 T"         = "#76B7B2",
  "NK"            = "#E15759",
  "B cell"        = "#59A14F",
  "CD14 Monocyte" = "#F28E2B",
  "CD16 Monocyte" = "#B07AA1",
  "DC"            = "#9C755F",
  "Megakaryocyte" = "#FF9DA7",
  "Basophil"      = "#BAB0AC"
)

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

for (ds_key in datasets) {
  cat("========== Processing", ds_key, "==========\n")
  
  # ---------------------------------------------------------
  # 1. QC Figures (S1A, S1B, S1C)
  # ---------------------------------------------------------
  qc_rds_path <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key, paste0(ds_key, "_qc_passed.rds"))
  qc_fig_dir <- file.path(base, "02_scrna", "pbmc", "01_qc", ds_key, "figures")
  dir.create(qc_fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(qc_rds_path)) {
    cat("Loading", qc_rds_path, "\n")
    obj_qc <- readRDS(qc_rds_path)
    
    # S1A: QC violins (using theme from 01_umap_celltype.R)
    p_violin <- VlnPlot(obj_qc,
      features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"),
      group.by = "sample_label", pt.size = 0, ncol = 4) &
      base_theme &
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    
    pdf_path <- file.path(qc_fig_dir, paste0(ds_key, "_01_prefilter_violin.pdf"))
    ggsave(pdf_path, p_violin, width = 18, height = 5)
    cat("  Saved:", pdf_path, "\n")
    
    # S1B: QC scatter (matching the vector style, but rasterized for performance)
    suppressWarnings({
      p_sc1 <- FeatureScatter(obj_qc, "nCount_RNA", "nFeature_RNA",
        group.by = "sample_label", pt.size = 0.1, raster = FALSE) +
        base_theme + ggtitle("nCount vs nFeature") +
        guides(color=guide_legend(override.aes=list(size=4, alpha=1)))
      p_sc1 <- ggrastr::rasterise(p_sc1, layers="Point", dpi=600)

      p_sc2 <- FeatureScatter(obj_qc, "nCount_RNA", "percent_mt",
        group.by = "sample_label", pt.size = 0.1, raster = FALSE) +
        base_theme +
        geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.6) +
        ggtitle("nCount vs % Mitochondrial") +
        guides(color=guide_legend(override.aes=list(size=4, alpha=1)))
      p_sc2 <- ggrastr::rasterise(p_sc2, layers="Point", dpi=600)

      p_sc3 <- FeatureScatter(obj_qc, "nFeature_RNA", "percent_mt",
        group.by = "sample_label", pt.size = 0.1, raster = FALSE) +
        base_theme +
        geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.6) +
        ggtitle("nFeature vs % Mitochondrial") +
        guides(color=guide_legend(override.aes=list(size=4, alpha=1)))
      p_sc3 <- ggrastr::rasterise(p_sc3, layers="Point", dpi=600)
    })
    
    p_scatter <- p_sc1 + p_sc2 + p_sc3
    pdf_path <- file.path(qc_fig_dir, paste0(ds_key, "_01_prefilter_scatter.pdf"))
    ggsave(pdf_path, p_scatter, width = 18, height = 5)
    cat("  Saved:", pdf_path, "\n")
    
    # S1C: Doublet scores
    # Since qc_passed.rds has removed doublets, the meta.data only has singlets.
    # We reconstruct the proportions from the saved csv report to accurately show doublets vs singlets.
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
        labs(title = paste(ds_key, "— doublet proportions"), x = NULL,
             y = "proportion", fill = "class")
             
      p_score <- ggplot(obj_qc@meta.data,
        aes(x = sample_label, y = scDbl_score, fill = scDbl_class)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
        scale_fill_manual(values = c("singlet" = "#4DADE2", "doublet" = "#E84B35")) +
        base_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste(ds_key, "— doublet scores (singlets only remaining)"), x = NULL,
             y = "scDblFinder score", fill = "class")
             
      p_doublets <- p_score + p_prop
      pdf_path <- file.path(qc_fig_dir, paste0(ds_key, "_03_doublets.pdf"))
      ggsave(pdf_path, p_doublets, width = 14, height = 5)
      cat("  Saved:", pdf_path, "\n")
    }
    
    rm(obj_qc)
    gc()
  } else {
    cat("Warning:", qc_rds_path, "not found\n")
  }
  
  # ---------------------------------------------------------
  # 2. Normalization Figures (S2A, S2B, S2C)
  # ---------------------------------------------------------
  norm_rds_path <- file.path(base, "02_scrna", "pbmc", "02_normalization", ds_key, paste0(ds_key, "_annotated.rds"))
  norm_fig_dir <- file.path(base, "02_scrna", "pbmc", "02_normalization", ds_key, "figures")
  dir.create(norm_fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(norm_rds_path)) {
    cat("Loading", norm_rds_path, "\n")
    obj_norm <- readRDS(norm_rds_path)
    
    cols <- celltype_colors[names(celltype_colors) %in% unique(obj_norm$cell_type)]
    if (length(cols) == 0) cols <- NULL
    
    # S2A: Per-disease UMAPs
    # Using exact logic from 01_umap_celltype.R: pt.size=0.08, but with ggrastr for performance
    p_celltype <- DimPlot(obj_norm, reduction = "umap",
                          group.by = "cell_type",
                          cols = cols,
                          label = TRUE, label.size = 3.5,
                          label.box = TRUE, label.color = "white",
                          repel = TRUE, pt.size = 0.08,
                          raster = FALSE) +
      base_theme +
      labs(title = paste(ds_key, "- SingleR cell types"), x="UMAP 1", y="UMAP 2", color="Cell type") +
      guides(color=guide_legend(override.aes=list(size=4, alpha=1), ncol=1))
      
    p_celltype <- ggrastr::rasterise(p_celltype, layers="Point", dpi=600)
      
    pdf_path <- file.path(norm_fig_dir, paste0(ds_key, "_07_umap_celltype.pdf"))
    ggsave(pdf_path, p_celltype, width = 12, height = 8)
    cat("  Saved:", pdf_path, "\n")
    
    # S2B: Marker dot plots
    all_markers <- unique(unlist(canonical_markers))
    markers_present <- all_markers[all_markers %in% rownames(obj_norm)]
    if (length(markers_present) > 0) {
      p_dot <- DotPlot(obj_norm,
                       features = markers_present,
                       group.by = "seurat_clusters",
                       assay    = "SCT") +
        base_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
        labs(title = paste(ds_key, "- canonical markers by cluster"))
        
      pdf_path <- file.path(norm_fig_dir, paste0(ds_key, "_08_dotplot_markers.pdf"))
      ggsave(pdf_path, p_dot, width = 16, height = 8)
      cat("  Saved:", pdf_path, "\n")
    }
    
    # S2C: Per-sample composition
    composition <- obj_norm@meta.data %>%
      group_by(group, cell_type) %>%
      summarise(n_cells = n(), .groups = "drop") %>%
      group_by(group) %>%
      mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
      ungroup() %>%
      arrange(group, desc(n_cells))

    # Apply composition styling
    p_comp <- ggplot(composition, aes(x = group, y = pct, fill = cell_type)) +
      geom_bar(stat = "identity", width = 0.6, colour = "white", linewidth = 0.3) +
      (if(!is.null(cols)) scale_fill_manual(values = cols) else NULL) +
      scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x,"%")) +
      base_theme +
      theme(axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
            legend.text  = element_text(size = 7),
            legend.key.size = unit(0.4, "cm")) +
      labs(title = paste(ds_key, "- cell type composition by group"),
           x = NULL, y = "% of cells", fill = "cell type") +
      guides(fill = guide_legend(ncol = 1, reverse = TRUE))
      
    pdf_path <- file.path(norm_fig_dir, paste0(ds_key, "_10_celltype_composition.pdf"))
    ggsave(pdf_path, p_comp, width = 10, height = 7)
    cat("  Saved:", pdf_path, "\n")
    
    rm(obj_norm)
    gc()
  } else {
    cat("Warning:", norm_rds_path, "not found\n")
  }
}

cat("Finished reproducing supplementary figures for Section 2.1 in PDF format.\n")
