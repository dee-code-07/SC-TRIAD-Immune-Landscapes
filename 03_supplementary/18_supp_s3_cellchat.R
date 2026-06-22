suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
})

obj_dir <- "~/sc-triad/02_scrna/pbmc/06_cellchat/objects"
out_dir <- "~/sc-triad/02_scrna/pbmc/06_cellchat/figures"

conds <- c("Control", "T2D_Moderate", "HTN", "Asthma_Mild", "Asthma_Severe")
cc_list <- list()
for (cond in conds) {
  fname <- file.path(obj_dir, paste0("cellchat_", cond, ".rds"))
  if(file.exists(fname)) cc_list[[cond]] <- readRDS(fname)
}

cat("Generating Figure S3A (Heatmaps)...\n")
for(cond in names(cc_list)) {
  fname <- file.path(out_dir, paste0("FigS3a_heatmap_", cond, ".pdf"))
  pdf(fname, width = 14, height = 9)
  ht1 <- netAnalysis_signalingRole_heatmap(cc_list[[cond]], pattern = "outgoing", title = "Outgoing", font.size=12, font.size.title=16)
  ht2 <- netAnalysis_signalingRole_heatmap(cc_list[[cond]], pattern = "incoming", title = "Incoming", font.size=12, font.size.title=16)
  draw(ht1 + ht2, column_title = paste0("Fig S3A: ", cond, " - Sender & Receiver Roles"), column_title_gp = gpar(fontsize = 22, fontface = "bold"))
  dev.off()
}

cat("Generating Figure S3B (BTLA Chord)...\n")
for(cond in c("Asthma_Mild", "Asthma_Severe")) {
  if("BTLA" %in% cc_list[[cond]]@netP$pathways) {
    fname <- file.path(out_dir, paste0("FigS3b_chord_BTLA_", cond, ".pdf"))
    pdf(fname, width=9, height=9)
    tryCatch({
      circos.clear()
      netVisual_aggregate(cc_list[[cond]], signaling = "BTLA", layout = "chord", title.name = paste0("Fig S3B: ", cond, " - BTLA"))
    }, error=function(e){})
    dev.off()
  }
}

cat("Figure S3 generation complete!\n")
