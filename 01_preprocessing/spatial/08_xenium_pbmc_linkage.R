# sc-triad project
# script: 03_xenium_pbmc_linkage.R  (v5 — label transfer class-check fix)
#
# Changes from v4:
#   FIX 4 (revised): class check now covers BOTH the PBMC SCT assay AND the
#   Xenium SCT assay before choosing the transfer path. In v4 the check only
#   covered Xenium; today's run showed PBMC SCT is the legacy SCTAssay class
#   (not Assay5), so the SCT-to-SCT branch was entered, then crashed because
#   FindTransferAnchors detected the mismatch internally. The corrected guard:
#     if (xen_sct_class == "SCTAssay" || pbmc_sct_class == "SCTAssay")
#   ensures RNA/LogNormalize is used whenever EITHER side is legacy SCTAssay.
#
# All other fixes from v4 are unchanged:
#   FIX 1: explicit features/assay/normalization.method in FindTransferAnchors
#   FIX 2: colMeans fallback when <3 panel genes
#   FIX 3: safe_join() skips JoinLayers on SCTAssay objects
#   FIX 5: winsorise colMeans scores at 99th percentile
#   FIX 6: deduplicate byte-for-byte identical score vectors
#   FIX 7: suppress Inf/-Inf log warnings on NA-only scores
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

set.seed(42)
options(future.globals.maxSize = 32 * 1024^3)

BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
OUT_DIR <- file.path(BASE, "04_spatial", "asthma_lung")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
OBJ_DIR <- file.path(OUT_DIR, "objects")
LOG_FILE <- file.path(OUT_DIR, "03_xenium_pbmc_linkage_v5.log")

PBMC_RDS   <- file.path(BASE, "02_scrna", "pbmc", "03_integration",
                         "triad_integrated.rds")
XDEG_CSV   <- file.path(BASE, "02_scrna", "pbmc", "04_deg", "tables",
                         "cross_disease_shared_DEGs.csv")
ASTHMA_CSV <- file.path(BASE, "02_scrna", "pbmc", "04_deg", "tables",
                         "Asthma_severity_independent_DEGs.csv")

# ── helpers ────────────────────────────────────────────────────────────────────
log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

save_tiff <- function(plot, fname, width = 12, height = 10) {
  path <- file.path(FIG_DIR, fname)
  tiff(path, width = width, height = height, units = "in",
       res = 300, compression = "lzw")
  print(plot)
  dev.off()
  log(paste("  saved:", fname))
}

# FIX 3: safe JoinLayers — skip for legacy SCTAssay, apply only to Assay5
safe_join <- function(obj, assay_name = NULL) {
  assay_name <- if (is.null(assay_name)) DefaultAssay(obj) else assay_name
  assay_obj  <- obj[[assay_name]]
  if (inherits(assay_obj, "Assay5")) {
    obj[[assay_name]] <- JoinLayers(assay_obj)
    log(paste("  JoinLayers applied to", assay_name, "(Assay5)"))
  } else {
    log(paste("  JoinLayers skipped for", assay_name,
              "(class:", class(assay_obj)[1], "- layers already unified)"))
  }
  return(obj)
}

# FIX 2 + FIX 5: scoring with colMeans fallback and outlier winsorisation
score_signature <- function(obj, genes, score_name, min_ams = 3,
                            winsor_pct = 0.99) {
  genes_in_panel <- intersect(genes, rownames(obj))
  n_in           <- length(genes_in_panel)
  n_total        <- length(genes)

  log(paste0("  ", score_name, ": ", n_in, "/", n_total, " genes in panel"))

  if (n_in == 0) {
    log("    -> 0 genes: score set to NA")
    obj[[score_name]]                     <- NA_real_
    obj[[paste0(score_name, "_n_genes")]] <- 0L
    return(obj)
  }

  if (n_in >= min_ams) {
    obj <- AddModuleScore(
      obj,
      features = list(genes_in_panel),
      name     = score_name,
      ctrl     = min(n_in, 20),
      assay    = "SCT",
      seed     = 42
    )
    obj[[score_name]]  <- obj@meta.data[[paste0(score_name, "1")]]
    obj@meta.data[[paste0(score_name, "1")]] <- NULL
    log(paste0("    -> AddModuleScore | range: ",
               round(min(obj[[score_name]], na.rm = TRUE), 3), " to ",
               round(max(obj[[score_name]], na.rm = TRUE), 3)))

  } else {
    # colMeans of per-gene z-scored SCT data — not background corrected
    expr_mat    <- GetAssayData(obj, assay = "SCT", layer = "data")[
      genes_in_panel, , drop = FALSE
    ]
    expr_scaled <- t(scale(t(as.matrix(expr_mat))))
    expr_scaled[is.nan(expr_scaled)] <- 0
    raw_scores  <- colMeans(expr_scaled, na.rm = TRUE)

    # FIX 5: winsorise
    cap       <- quantile(raw_scores, winsor_pct,     na.rm = TRUE)
    floor_val <- quantile(raw_scores, 1 - winsor_pct, na.rm = TRUE)
    n_capped  <- sum(raw_scores > cap | raw_scores < floor_val, na.rm = TRUE)
    raw_scores <- pmin(pmax(raw_scores, floor_val), cap)

    obj[[score_name]] <- raw_scores
    log(paste0("    -> colMeans fallback (", n_in, " genes) | ",
               "range after winsorise: ",
               round(min(raw_scores, na.rm = TRUE), 3), " to ",
               round(max(raw_scores, na.rm = TRUE), 3),
               " | ", n_capped, " outlier cells capped"))
    log("    WARNING: colMeans score, not background-corrected. Note in methods.")
  }

  obj[[paste0(score_name, "_n_genes")]] <- n_in
  return(obj)
}

# ── start ──────────────────────────────────────────────────────────────────────
log("Xenium-PBMC linkage v5 started")

xen <- readRDS(file.path(OBJ_DIR, "xenium_annotated.rds"))
log(paste("Xenium:", ncol(xen), "cells |", nrow(xen), "genes"))
panel_genes <- rownames(xen)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Label transfer (PBMC Asthma -> Xenium airway)
#
# FIX 4 (REVISED v5): Check SCT class on BOTH objects before choosing path.
#
# Seurat version history produces two possible SCT assay classes:
#   SCTAssay  — produced by SCTransform in Seurat v4 / early v5.
#               Internal layers already unified. JoinLayers crashes on it.
#               FindTransferAnchors does NOT support mixing SCTAssay with Assay5.
#   Assay5    — produced by SCTransform in Seurat v5.1+.
#               Layers may be split per sample; JoinLayers is needed.
#               FindTransferAnchors supports Assay5 on both sides.
#
# The safe strategy:
#   If EITHER object carries SCTAssay -> force RNA/LogNormalize on both sides.
#   If BOTH carry Assay5 -> attempt SCT-to-SCT transfer.
#
# This is the only change from v4. The RNA fallback worked in all prior test
# runs where it was reached; the issue was that v4 only checked Xenium.
# ══════════════════════════════════════════════════════════════════════════════

log("Loading PBMC object for label transfer...")
pbmc <- readRDS(PBMC_RDS)
pbmc_asthma <- subset(pbmc, subset = disease == "Asthma")
log(paste("PBMC asthma subset:", ncol(pbmc_asthma), "cells"))

# ── FIX 4 (v5): check both SCT assay classes ──────────────────────────────────
xen_sct_class  <- class(xen[["SCT"]])[1]
pbmc_sct_class <- class(pbmc_asthma[["SCT"]])[1]

log(paste("Xenium SCT class:", xen_sct_class))
log(paste("PBMC   SCT class:", pbmc_sct_class))

if (xen_sct_class == "SCTAssay" || pbmc_sct_class == "SCTAssay") {
  # At least one side is legacy SCTAssay.
  # Mixing SCTAssay with Assay5 in FindTransferAnchors is not supported.
  # Force RNA/LogNormalize on both sides.
  log("FIX 4 (v5): SCTAssay detected on at least one side.")
  log("  Forcing RNA/LogNormalize for label transfer on both objects.")

  # PBMC: switch to RNA, normalize
  DefaultAssay(pbmc_asthma) <- "RNA"
  pbmc_asthma <- safe_join(pbmc_asthma, "RNA")
  pbmc_asthma <- NormalizeData(pbmc_asthma, verbose = FALSE)

  # Xenium: build RNA assay from SCT counts, normalize
  xen_raw_counts <- GetAssayData(xen, assay = "SCT", layer = "counts")
  xen[["RNA"]]   <- CreateAssayObject(counts = xen_raw_counts)
  DefaultAssay(xen) <- "RNA"
  xen <- NormalizeData(xen, verbose = FALSE)

  common_genes <- intersect(rownames(xen), rownames(pbmc_asthma))
  log(paste("Common genes (RNA x RNA):", length(common_genes)))
  ref_assay   <- "RNA"
  query_assay <- "RNA"
  norm_method <- "LogNormalize"

} else {
  # Both are Assay5 — attempt SCT-to-SCT transfer
  log("Both objects have Assay5 SCT — attempting SCT-to-SCT label transfer")

  DefaultAssay(pbmc_asthma) <- "SCT"
  pbmc_asthma <- safe_join(pbmc_asthma, "SCT")

  DefaultAssay(xen) <- "SCT"
  xen <- safe_join(xen, "SCT")

  common_genes <- intersect(rownames(xen), rownames(pbmc_asthma))
  log(paste("Common genes (SCT x SCT):", length(common_genes)))

  if (length(common_genes) < 20) {
    log("WARNING: <20 SCT common genes. Falling back to RNA...")

    DefaultAssay(pbmc_asthma) <- "RNA"
    pbmc_asthma <- safe_join(pbmc_asthma, "RNA")
    pbmc_asthma <- NormalizeData(pbmc_asthma, verbose = FALSE)

    xen_raw_counts <- GetAssayData(xen, assay = "SCT", layer = "counts")
    xen[["RNA"]]   <- CreateAssayObject(counts = xen_raw_counts)
    DefaultAssay(xen) <- "RNA"
    xen <- NormalizeData(xen, verbose = FALSE)

    common_genes <- intersect(rownames(xen), rownames(pbmc_asthma))
    log(paste("Common genes (RNA fallback):", length(common_genes)))
    ref_assay   <- "RNA"
    query_assay <- "RNA"
    norm_method <- "LogNormalize"
  } else {
    ref_assay   <- "SCT"
    query_assay <- "SCT"
    norm_method <- "SCT"
  }
}

do_label_transfer <- length(common_genes) >= 20

if (do_label_transfer) {
  log(paste("Running FindTransferAnchors:",
            ref_assay, "->", query_assay,
            "| features:", length(common_genes)))

  tryCatch({
    anchors <- FindTransferAnchors(
      reference            = pbmc_asthma,
      query                = xen,
      features             = common_genes,    # FIX 1: explicit feature list
      reference.assay      = ref_assay,       # FIX 1: explicit assay
      query.assay          = query_assay,     # FIX 1: explicit assay
      normalization.method = norm_method,     # FIX 1: explicit norm method
      dims                 = 1:30,
      verbose              = FALSE
    )

    log(paste("Anchors found:", nrow(anchors@anchors)))

    predictions <- TransferData(
      anchorset = anchors,
      refdata   = pbmc_asthma$cell_type,
      dims      = 1:30,
      verbose   = FALSE
    )

    xen$pbmc_predicted_celltype   <- predictions$predicted.id
    xen$pbmc_prediction_score_max <- predictions$prediction.score.max
    xen$pbmc_predicted_highconf   <- ifelse(
      xen$pbmc_prediction_score_max > 0.5,
      xen$pbmc_predicted_celltype,
      "Low_confidence"
    )

    log(paste("Label transfer done | median score:",
              round(median(xen$pbmc_prediction_score_max), 3),
              "| high-conf cells:", sum(xen$pbmc_prediction_score_max > 0.5)))

    # Save per-cell predictions
    write.csv(
      data.frame(
        cell_id          = colnames(xen),
        xenium_cell_type = xen$cell_type_xenium,
        pbmc_predicted   = xen$pbmc_predicted_celltype,
        prediction_score = xen$pbmc_prediction_score_max
      ),
      file.path(TAB_DIR, "pbmc_label_transfer_results.csv"),
      row.names = FALSE
    )

    # UMAP coloured by label transfer result
    p_lt_umap <- DimPlot(xen, group.by = "pbmc_predicted_highconf",
                         label = TRUE, label.size = 2.5, repel = TRUE,
                         pt.size = 0.3, raster = FALSE) +
      theme_classic(base_size = 9) +
      labs(title    = "PBMC label transfer onto Xenium",
           subtitle = "Asthma PBMC -> airway tissue | score > 0.5 shown")
    save_tiff(p_lt_umap, "08_pbmc_label_transfer_umap.tiff",
              width = 12, height = 9)

    # Prediction score distribution
    p_score <- ggplot(
      data.frame(score = xen$pbmc_prediction_score_max),
      aes(x = score)
    ) +
      geom_histogram(bins = 50, fill = "#00BFC4", colour = "white",
                     linewidth = 0.3) +
      geom_vline(xintercept = 0.5, linetype = "dashed",
                 colour = "red", linewidth = 0.8) +
      theme_classic(base_size = 11) +
      labs(title    = "Label transfer prediction scores",
           subtitle = "Dashed = 0.5 high-confidence threshold",
           x = "Max prediction score", y = "Cells")
    save_tiff(p_score, "08b_label_transfer_score_dist.tiff",
              width = 8, height = 5)

    # Concordance table: Xenium annotation vs PBMC prediction (immune only)
    immune_types <- c("T_cell", "CD4_T", "CD8_T", "NK", "B_cell",
                      "Macrophage", "Monocyte", "DC", "Mast")
    immune_mask  <- xen$cell_type_xenium %in% immune_types &
                    xen$pbmc_prediction_score_max > 0.5

    if (sum(immune_mask) > 50) {
      conf_mat <- table(
        Xenium    = xen$cell_type_xenium[immune_mask],
        PBMC_pred = xen$pbmc_predicted_celltype[immune_mask]
      )
      write.csv(as.data.frame.matrix(conf_mat),
                file.path(TAB_DIR, "label_transfer_confusion_immune.csv"))
      log("Concordance table saved")
    }

    # Spatial plot: prediction score by tissue location (per sample)
    for (samp_label in unique(xen$sample_id)) {
      samp_cells <- colnames(xen)[xen$sample_id == samp_label]
      if (length(samp_cells) == 0) next

      xen_samp <- xen[, samp_cells]
      coords   <- tryCatch(
        GetTissueCoordinates(xen_samp, which = "centroids"),
        error = function(e) NULL
      )
      if (is.null(coords)) next

      plot_df <- data.frame(
        x     = coords$x,
        y     = coords$y,
        score = xen_samp$pbmc_prediction_score_max
      )

      p_spatial_lt <- ggplot(plot_df, aes(x = x, y = y, colour = score)) +
        geom_point(size = 0.4, alpha = 0.8) +
        scale_colour_viridis_c(option = "plasma", name = "Pred. score") +
        coord_equal() + theme_void(base_size = 10) +
        labs(title    = paste("Label transfer score -", samp_label),
             subtitle = "PBMC -> Xenium prediction confidence") +
        theme(plot.background = element_rect(fill = "white", colour = NA))

      save_tiff(p_spatial_lt,
                paste0("08c_spatial_label_transfer_", samp_label, ".tiff"),
                width = 10, height = 9)
    }

  }, error = function(e) {
    log(paste("Label transfer failed:", conditionMessage(e)))
    log("Proceeding with signature scoring only")
  })

} else {
  log("Skipping label transfer: insufficient common genes")
}

# Reset to SCT for all signature scoring
DefaultAssay(xen) <- "SCT"

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Signature scoring (unchanged from v4)
# ══════════════════════════════════════════════════════════════════════════════

log("Loading PBMC DEG signatures...")

triad_up <- triad_dn <- asthma_up <- asthma_dn <- character(0)

if (file.exists(XDEG_CSV)) {
  xdeg <- read.csv(XDEG_CSV)
  log(paste("Cross-disease DEGs loaded:", nrow(xdeg), "genes"))

  triad_up  <- xdeg$gene[xdeg$n_diseases == 3 & xdeg$direction == "Up"]
  triad_dn  <- xdeg$gene[xdeg$n_diseases == 3 & xdeg$direction == "Down"]
  asthma_up <- xdeg$gene[grepl("Asthma", xdeg$diseases) & xdeg$direction == "Up"]
  asthma_dn <- xdeg$gene[grepl("Asthma", xdeg$diseases) & xdeg$direction == "Down"]

  log(paste("Triad UP:", length(triad_up), "| DOWN:", length(triad_dn)))
  log(paste("Asthma UP:", length(asthma_up), "| DOWN:", length(asthma_dn)))

  log("=== PANEL COVERAGE REPORT (paste into methods section) ===")
  for (sig in list(
    list(name = "TriadUP",  genes = triad_up),
    list(name = "TriadDN",  genes = triad_dn),
    list(name = "AsthmaUP", genes = asthma_up),
    list(name = "AsthmaDN", genes = asthma_dn)
  )) {
    n_in   <- length(intersect(sig$genes, panel_genes))
    n_tot  <- length(sig$genes)
    pct    <- if (n_tot > 0) round(100 * n_in / n_tot) else 0
    log(paste0("  ", sig$name, ": ", n_in, "/", n_tot, " (", pct, "%)"))
    if (n_in > 0 && n_tot > 0)
      log(paste0("    genes in panel: ",
                 paste(intersect(sig$genes, panel_genes), collapse = ", ")))
  }
  log("==========================================================")
}

asthma_sev_up_nk <- asthma_sev_up_mono <- character(0)
if (file.exists(ASTHMA_CSV)) {
  asthma_sev <- read.csv(ASTHMA_CSV)
  asthma_sev_up_nk   <- asthma_sev$gene[
    asthma_sev$cell_type == "NK" & asthma_sev$direction == "Up"
  ]
  asthma_sev_up_mono <- asthma_sev$gene[
    asthma_sev$cell_type %in% c("CD14 Monocyte", "CD16 Monocyte") &
    asthma_sev$direction == "Up"
  ]
  log(paste("NK sev-independent UP:", length(asthma_sev_up_nk)))
  log(paste("Monocyte sev-independent UP:", length(asthma_sev_up_mono)))
}

nk_signature <- c("NKG7", "TRBC1", "KLRB1", "HCST", "IFITM2", "PSMB10",
                  "PFN1", "SH3BGRL3", "LY6E", "PSME2", "GZMM", "PSME1",
                  "PSMB9", "MYL6")

log("Scoring signatures...")

if (length(triad_up)           > 0) xen <- score_signature(xen, triad_up,           "TriadUP_score")
if (length(triad_dn)           > 0) xen <- score_signature(xen, triad_dn,           "TriadDN_score")
if (length(asthma_up)          > 0) xen <- score_signature(xen, asthma_up,          "AsthmaUP_score")
if (length(asthma_dn)          > 0) xen <- score_signature(xen, asthma_dn,          "AsthmaDN_score")
                                    xen <- score_signature(xen, nk_signature,        "NK_CrossDisease_score")
if (length(asthma_sev_up_nk)   > 0) xen <- score_signature(xen, asthma_sev_up_nk,   "NK_AsthmaUP_score")
if (length(asthma_sev_up_mono) > 0) xen <- score_signature(xen, asthma_sev_up_mono, "Mono_AsthmaUP_score")

# ── verify + deduplicate ───────────────────────────────────────────────────────
score_cols <- grep("_score$", colnames(xen@meta.data), value = TRUE)

log("Score range check:")
for (sc in score_cols) {
  vals    <- xen@meta.data[[sc]]
  n_genes <- xen@meta.data[[paste0(sc, "_n_genes")]][1]
  if (all(is.na(vals))) {
    log(paste0("  ", sc, " (", n_genes, " genes): all NA - skipped"))
  } else {
    lo <- round(suppressWarnings(min(vals, na.rm = TRUE)), 4)
    hi <- round(suppressWarnings(max(vals, na.rm = TRUE)), 4)
    log(paste0("  ", sc, " (", n_genes, " genes): ", lo, " to ", hi))
  }
}

# Active: non-NA and non-flat
active_scores <- score_cols[sapply(score_cols, function(sc) {
  vals <- xen@meta.data[[sc]]
  !all(is.na(vals)) &&
    suppressWarnings(diff(range(vals, na.rm = TRUE))) > 1e-6
})]

# FIX 6: drop byte-for-byte duplicates
if (length(active_scores) > 1) {
  rounded_vecs <- lapply(active_scores,
                         function(sc) round(xen@meta.data[[sc]], 6))
  is_dup <- duplicated(rounded_vecs)
  if (any(is_dup)) {
    log("FIX 6: Duplicate score vectors excluded from plots:")
    for (sc in active_scores[is_dup])
      log(paste0("  '", sc, "' is identical to a previous score"))
  }
  active_scores <- active_scores[!is_dup]
}

log(paste("Active unique scores:", paste(active_scores, collapse = ", ")))

# ── visualizations ─────────────────────────────────────────────────────────────
for (score in active_scores) {
  n_genes      <- xen@meta.data[[paste0(score, "_n_genes")]][1]
  subtitle_txt <- paste0(
    n_genes, " panel genes | ",
    if (n_genes < 3)
      "colMeans surrogate (winsorised 99th pct) - interpret cautiously"
    else
      "AddModuleScore background-corrected"
  )

  p_umap <- FeaturePlot(xen, features = score, reduction = "umap",
                        pt.size = 0.3, raster = FALSE,
                        cols = c("grey92", "#C00000")) +
    theme_classic(base_size = 10) +
    labs(title = gsub("_score", "", score), subtitle = subtitle_txt)
  save_tiff(p_umap, paste0("09_umap_", score, ".tiff"), width = 9, height = 8)

  p_vio <- VlnPlot(xen, features = score, group.by = "cell_type_xenium",
                   pt.size = 0, sort = "increasing") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(title = paste(gsub("_score", "", score), "- by cell type"),
         subtitle = subtitle_txt)
  save_tiff(p_vio, paste0("10_violin_", score, ".tiff"),
            width = 12, height = 6)
}

# Spatial score plots per sample (top 2 unique scores)
for (samp_label in unique(xen$sample_id)) {
  samp_cells <- colnames(xen)[xen$sample_id == samp_label]
  if (length(samp_cells) == 0) next

  xen_samp <- xen[, samp_cells]
  coords   <- tryCatch(
    GetTissueCoordinates(xen_samp, which = "centroids"),
    error = function(e) NULL
  )
  if (is.null(coords)) next

  for (score in head(active_scores, 2)) {
    scores_vec <- xen_samp@meta.data[[score]]
    n_genes    <- xen_samp@meta.data[[paste0(score, "_n_genes")]][1]

    plot_df <- data.frame(x = coords$x, y = coords$y, score = scores_vec)

    p_spatial <- ggplot(plot_df, aes(x = x, y = y, colour = score)) +
      geom_point(size = 0.4, alpha = 0.8) +
      scale_colour_gradient2(
        low = "#2166AC", mid = "grey90", high = "#C00000", midpoint = 0
      ) +
      coord_equal() + theme_void(base_size = 10) +
      labs(title    = paste(gsub("_score", "", score), "-", samp_label),
           subtitle = paste(n_genes, "panel genes"),
           colour   = "Score") +
      theme(plot.background = element_rect(fill = "white", colour = NA))

    save_tiff(p_spatial,
              paste0("11_spatial_", samp_label, "_", score, ".tiff"),
              width = 10, height = 9)
  }
}

# ── summary table ──────────────────────────────────────────────────────────────
score_summary <- xen@meta.data %>%
  group_by(cell_type_xenium) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(active_scores),
      list(
        mean = \(x) round(mean(x, na.rm = TRUE), 4),
        sd   = \(x) round(sd(x,   na.rm = TRUE), 4)
      )
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.csv(score_summary,
          file.path(TAB_DIR, "signature_scores_by_celltype.csv"),
          row.names = FALSE)
log("Score summary table saved")
print(as.data.frame(score_summary), row.names = FALSE)

# Methods-ready panel coverage table
panel_coverage <- data.frame(
  signature         = score_cols,
  n_genes_panel     = sapply(score_cols, function(sc)
    xen@meta.data[[paste0(sc, "_n_genes")]][1]),
  method_used       = sapply(score_cols, function(sc) {
    n <- xen@meta.data[[paste0(sc, "_n_genes")]][1]
    if (is.na(n) || n == 0) "none"
    else if (n >= 3)         "AddModuleScore"
    else                     "colMeans_scaled_winsorised"
  }),
  included_in_plots = score_cols %in% active_scores
)
write.csv(panel_coverage,
          file.path(TAB_DIR, "signature_panel_coverage.csv"),
          row.names = FALSE)
log("Panel coverage table saved")
print(panel_coverage)

# ── save ───────────────────────────────────────────────────────────────────────
saveRDS(xen, file.path(OBJ_DIR, "xenium_final.rds"))
log("xenium_final.rds saved")
log("Next: 04_xenium_figures_publication.R")