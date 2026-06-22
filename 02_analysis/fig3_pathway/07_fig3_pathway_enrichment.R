
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(stringr)
  library(scales)
  library(patchwork)
})

# Path Setup
base_dir <- "/mnt/e/Documents/sctriad"
deg_dir  <- file.path(base_dir, "02_scrna/04_deg/tables")
out_dir  <- file.path(base_dir, "manuscript_new/main_figures/fig3_pathway_enrichment")
supp_dir <- file.path(base_dir, "manuscript_new/supplementary_figures/figs4_pathway_supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

# Constants from original script
PADJ_CUTOFF   <- 0.05
QVAL_CUTOFF   <- 0.20
MIN_GSSIZE    <- 5L
MAX_GSSIZE    <- 300L
SIMPLIFY_CUT  <- 0.5
TOP_N_DISPLAY <- 15L # Reduced for cleaner plots
MIN_COUNT     <- 3L
MIN_TERMS     <- 3L # Slightly more permissive to catch shared ones

# Blacklists
KEGG_BL <- c(
  "Pathways of neurodegeneration - multiple diseases", "Huntington disease", "Parkinson disease",
  "Prion disease", "Alzheimer disease", "Amyotrophic lateral sclerosis", "Thermogenesis",
  "Non-alcoholic fatty liver disease", "Cardiac muscle contraction", "Metabolic pathways",
  "Ribosome", "Herpes simplex virus 1 infection", "Human cytomegalovirus infection",
  "Influenza A", "Viral myocarditis", "Toxoplasmosis", "Human T-cell leukemia virus 1 infection",
  "Spinocerebellar ataxia", "Diabetic cardiomyopathy", "Chemical carcinogenesis - reactive oxygen species",
  "Retrograde endocannabinoid signaling", "RNA polymerase", "Ribosome biogenesis in eukaryotes",
  "Nucleotide excision repair", "Spliceosome", "Protein export", "Leishmaniasis",
  "Salmonella infection", "Yersinia infection", "Epstein-Barr virus infection",
  "Shigellosis", "Bacterial invasion of epithelial cells", "Staphylococcus aureus infection",
  "Tuberculosis", "Tight junction", "Fc epsilon RI signaling pathway",
  "Complement and coagulation cascades", "Fc gamma R-mediated phagocytosis",
  "T cell differentiation in thymus"
)
GO_BL <- c(
  "cytoplasmic translation", "formation of cytoplasmic translation initiation complex",
  "rRNA metabolic process", "ribonucleoprotein complex biogenesis", "maturation of SSU-rRNA",
  "protein-RNA complex organization", "protein folding", "protein targeting",
  "protein localization to endoplasmic reticulum", "establishment of protein localization to membrane",
  "RNA 5'-end processing", "mRNA splicing, via spliceosome", "RNA splicing, via transesterification reactions",
  "axon development", "glial cell differentiation", "positive regulation of neuron projection development",
  "negative regulation of amyloid fibril formation", "connective tissue development",
  "regulation of chondrocyte differentiation", "erythrocyte differentiation", "erythrocyte homeostasis",
  "myeloid cell homeostasis", "response to BMP", "cellular response to BMP stimulus",
  "response to ketone", "regulation of platelet activation", "platelet aggregation",
  "T cell differentiation", "regulation of viral life cycle", "negative regulation of viral life cycle",
  "negative regulation of viral process", "regulation of viral entry into host cell",
  "negative regulation of viral entry into host cell", "modulation by symbiont of entry into host",
  "regulation of biological process involved in symbiotic interaction", "disruption of anatomical structure in another organism",
  "disruption of cell in another organism", "killing of cells of another organism",
  "spliceosomal snRNP assembly", "neuron projection development", "oligodendrocyte differentiation",
  "intracellular chemical homeostasis", "anatomical structure homeostasis", "tissue homeostasis"
)
ALL_BL <- c(KEGG_BL, GO_BL)

# Theme
pub_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle     = element_text(size = 9, hjust = 0, colour = "grey40"),
    axis.title        = element_text(face = "bold", size = 10),
    axis.text         = element_text(size = 9, colour = "black"),
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.2)
  )

# Helper functions
symbols_to_entrez <- function(symbols) {
  if (length(symbols) == 0L) return(character(0))
  tbl <- suppressMessages(suppressWarnings(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ))
  as.character(unique(tbl$ENTREZID))
}

run_enrichment <- function(entrez_genes, universe = NULL) {
  if (length(entrez_genes) < 3) return(NULL)
  
  go_res <- tryCatch(
    enrichGO(gene = entrez_genes, universe = universe, OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             minGSSize = 3, maxGSSize = 500, readable = TRUE),
    error = function(e) NULL
  )
  
  kegg_res <- tryCatch(
    enrichKEGG(gene = entrez_genes, universe = universe, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               minGSSize = 3, maxGSSize = 500),
    error = function(e) NULL
  )
  
  res_list <- list()
  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    go_res <- clusterProfiler::simplify(go_res, cutoff = 0.7, by = "p.adjust", select_fun = min)
    res_list$go <- as.data.frame(go_res) %>% mutate(source = "GO_BP")
  }
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    res_list$kegg <- as.data.frame(kegg_res) %>% mutate(source = "KEGG")
  }
  
  if (length(res_list) == 0) return(NULL)
  
  bind_rows(res_list) %>%
    filter(!Description %in% ALL_BL) %>%
    mutate(GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text=x))))
}

make_dotplot <- function(df, title, color_scheme = "up") {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = TOP_N_DISPLAY) %>%
    mutate(Description = fct_reorder(Description, GeneRatio_num))
  
  col_low <- if(color_scheme == "up") "#fcbf7a" else "#9ecae1"
  col_high <- if(color_scheme == "up") "#a50026" else "#08306b"
  
  ggplot(df, aes(x = GeneRatio_num, y = Description, size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_gradient(low = col_low, high = col_high) +
    scale_size_continuous(range = c(3, 8)) +
    labs(title = title, x = "Gene Ratio", y = NULL, color = "-log10(p.adj)", size = "Count") +
    pub_theme
}

# 1. Asthma T cells (Panel 3a candidate)
message("Processing Asthma T cells...")
asthma_files <- list.files(deg_dir, pattern = "Asthma_.*_T_DEGs.csv", full.names = TRUE)
asthma_results <- list()
for (f in asthma_files) {
  df <- read.csv(f)
  up_genes <- df %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)
  if (length(up_genes) >= 5) {
    entrez <- symbols_to_entrez(up_genes)
    enr <- run_enrichment(entrez)
    if (!is.null(enr)) {
      stem <- sub("_DEGs.csv", "", basename(f))
      enr$Comparison <- stem
      asthma_results[[stem]] <- enr
    }
  }
}

if (length(asthma_results) > 0) {
  # Focus on a representative one: Asthma_Severe_vs_Control_CD8_T
  rep_key <- names(asthma_results)[grep("Severe_vs_Control_CD8_T", names(asthma_results))[1]]
  p3a <- make_dotplot(asthma_results[[rep_key]], paste("Pathway Enrichment: Asthma CD8+ T cells (Up)"))
  ggsave(file.path(out_dir, "fig3a_asthma_cd8t_enrichment.pdf"), p3a, width = 8, height = 6)
}

# 2. Cross-disease Shared (Panel 3b candidate)
message("Processing Cross-disease Shared...")
xd_file <- file.path(deg_dir, "cross_disease_shared_DEGs.csv")
if (file.exists(xd_file)) {
  xd_df <- read.csv(xd_file)
  # Group by cell type and direction
  xd_summary <- xd_df %>%
    group_by(cell_type, direction) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  
  for (i in 1:nrow(xd_summary)) {
    genes <- unlist(xd_summary$genes[i])
    ct <- xd_summary$cell_type[i]
    dir <- xd_summary$direction[i]
    message("  ", ct, " ", dir, ": ", length(genes), " genes")
    entrez <- symbols_to_entrez(genes)
    enr <- run_enrichment(entrez)
    if (!is.null(enr)) {
      p <- make_dotplot(enr, paste("Shared Pathways:", ct, dir), color_scheme = tolower(dir))
      fname <- paste0("fig3b_shared_", tolower(ct), "_", tolower(dir), "_enrichment.pdf")
      ggsave(file.path(out_dir, fname), p, width = 8, height = 6)
    }
  }
}

# 3. NK Signature (Panel 3c candidate)
message("Processing NK Signature...")
NK_SIGNATURE <- c("NKG7","TRBC1","KLRB1","HCST","IFITM2","PSMB10", "PFN1","SH3BGRL3","LY6E","PSME2","GZMM","PSME1", "PSMB9","MYL6")
nk_entrez <- symbols_to_entrez(NK_SIGNATURE)
nk_enr <- run_enrichment(nk_entrez)
if (!is.null(nk_enr)) {
  p3c <- make_dotplot(nk_enr, "Enrichment of Cross-disease NK Signature")
  ggsave(file.path(out_dir, "fig3c_nk_signature_enrichment.pdf"), p3c, width = 8, height = 5)
}

# 4. Supplementary: All comparisons
message("Processing all comparisons for supplementary...")
all_files <- list.files(deg_dir, pattern = "_DEGs.csv", full.names = TRUE)
all_files <- all_files[!basename(all_files) %in% c("all_significant_DEGs.csv", "cross_disease_shared_DEGs.csv", "Asthma_severity_independent_DEGs.csv")]

pdf(file.path(supp_dir, "figS5_all_pathway_enrichments.pdf"), width = 10, height = 8)
for (f in all_files) {
  stem <- sub("_DEGs.csv", "", basename(f))
  df <- read.csv(f)
  # Use generic column names
  lfc_col <- if("avg_log2FC" %in% colnames(df)) "avg_log2FC" else "log2FoldChange"
  padj_col <- if("p_val_adj" %in% colnames(df)) "p_val_adj" else "padj"
  
  for (dir in c("Up", "Down")) {
    genes <- if(dir == "Up") {
      df %>% filter(.data[[padj_col]] < 0.05, .data[[lfc_col]] > 0.25) %>% pull(gene)
    } else {
      df %>% filter(.data[[padj_col]] < 0.05, .data[[lfc_col]] < -0.25) %>% pull(gene)
    }
    
    if (length(genes) >= 5) {
      entrez <- symbols_to_entrez(genes)
      enr <- run_enrichment(entrez)
      if (!is.null(enr)) {
        p <- make_dotplot(enr, paste(stem, dir), color_scheme = tolower(dir))
        print(p)
      }
    }
  }
}
dev.off()

message("Done!")
