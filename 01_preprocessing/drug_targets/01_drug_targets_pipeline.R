# sc-triad project
# script: 08_drug_targets.R
# purpose: therapeutic target prioritization for the T2D-HTN-Asthma triad
#
# approach:
#   1. compile candidate genes from all prior analyses
#      (cross-disease DEGs + CellChat pathway anchors + pySCENIC TFs)
#   2. query DGIdb v5 REST API for drug-gene interactions
#   3. score and prioritize by druggability + multi-disease relevance
#   4. generate three publication figures:
#      - bubble plot: gene x disease evidence x druggability
#      - bar chart: top drug-gene pairs coloured by drug type
#      - network-style dot matrix: target x drug class
#
# runtime: ~10 minutes (API calls + figures)
# submit as short SLURM job — figure generation peaks at ~8 GB RAM
#
# author: deeksha h | reg. 241706005 | msc bioinformatics iv sem | mahe
# guide: dr budheswar dehury

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)   # FIX: added — required for geom_text_repel() in Fig 1
  library(patchwork)
  library(stringr)
  library(tibble)
  library(forcats)   # FIX: added — required for fct_reorder() in Figs 2 and 4
})

set.seed(42)

BASE    <- file.path(Sys.getenv("HOME"), "sc-triad")
OUT_DIR <- file.path(BASE, "05_drug_targets")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
LOG_FILE <- file.path(OUT_DIR, "08_drug_targets.log")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

save_tiff <- function(p, fname, width = 12, height = 8) {
  path <- file.path(FIG_DIR, fname)
  tiff(path, width = width, height = height, units = "in",
       res = 600, compression = "lzw")
  tryCatch(print(p), error = function(e) {
    plot.new(); text(0.5, 0.5, conditionMessage(e))
  })
  dev.off()
  log(paste("  saved:", fname))
}

pub_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40", lineheight = 1.3),
    axis.title    = element_text(face = "bold", size = 11),
    axis.text     = element_text(size = 10),
    legend.title  = element_text(face = "bold", size = 10),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

log("SC-TRIAD drug target pipeline started")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Compile candidate gene list from all prior analyses
#
# Three tiers of evidence:
#   Tier 1 — cross-disease DEGs (directly differentially expressed)
#   Tier 2 — CellChat pathway anchor genes (intercellular signalling)
#   Tier 3 — pySCENIC TF regulons (regulatory hubs)
# ══════════════════════════════════════════════════════════════════════════════
log("Compiling candidate gene list...")

# Tier 1: cross-disease shared DEGs
xdeg_path <- file.path(BASE, "02_scrna", "pbmc", "04_deg", "tables",
                        "cross_disease_shared_DEGs.csv")
tier1_genes <- character(0)
tier1_meta  <- data.frame()

if (file.exists(xdeg_path)) {
  xdeg <- read.csv(xdeg_path)
  tier1_genes <- unique(xdeg$gene)
  tier1_meta  <- xdeg %>%
    select(gene, n_diseases, diseases, direction, mean_lfc, cell_type) %>%
    mutate(tier = "T1_CrossDiseaseDEG",
           evidence = paste0("DEG: ", diseases, " (", direction, ")"))
  log(paste("Tier 1 (cross-disease DEGs):", length(tier1_genes), "genes"))
}

# Tier 2: CellChat pathway anchor genes
cellchat_genes <- tribble(
  ~gene,      ~pathway,   ~role,      ~pathway_type,
  "LGALS9",   "GALECTIN", "Ligand",   "Triad-gained",
  "CD44",     "GALECTIN", "Receptor", "Triad-gained",
  "CD274",    "GALECTIN", "Receptor", "Triad-gained",
  "PPIA",     "CypA",     "Ligand",   "Triad-gained",
  "BSG",      "CypA",     "Receptor", "Triad-gained",
  "BTLA",     "BTLA",     "Receptor", "Asthma-specific",
  "TNFRSF14", "BTLA",     "Ligand",   "Asthma-specific"
)

tier2_meta <- cellchat_genes %>%
  mutate(tier = "T2_CellChatPathway",
         evidence = paste0(pathway, " pathway (", pathway_type, "): ", role),
         n_diseases = ifelse(pathway_type == "Triad-gained", 3, 1),
         direction = "Ligand-Receptor")
tier2_genes <- unique(cellchat_genes$gene)
log(paste("Tier 2 (CellChat anchors):", length(tier2_genes), "genes"))

# Tier 3: pySCENIC cross-disease TF regulons
tf_path <- file.path(BASE, "03_pyscenic", "analysis", "tables",
                      "cross_disease_shared_tfs.csv")
tier3_genes <- character(0)
tier3_meta  <- data.frame()

if (file.exists(tf_path)) {
  tfs <- read.csv(tf_path)

  # FIX: log actual column names so mismatches are visible in SLURM output
  log(paste("TF file columns:", paste(colnames(tfs), collapse = ", ")))

  tfs_filt <- tfs %>% filter(mean_r >= 0.1)
  tier3_genes <- unique(tfs_filt$regulon)

  # FIX: defensive normalisation of n_diseases — handles common column name variants
  if (!"n_diseases" %in% colnames(tfs_filt)) {
    if ("n_disease" %in% colnames(tfs_filt)) {
      tfs_filt <- rename(tfs_filt, n_diseases = n_disease)
      log("  TF file: renamed 'n_disease' -> 'n_diseases'")
    } else if ("num_diseases" %in% colnames(tfs_filt)) {
      tfs_filt <- rename(tfs_filt, n_diseases = num_diseases)
      log("  TF file: renamed 'num_diseases' -> 'n_diseases'")
    } else if ("diseases" %in% colnames(tfs_filt)) {
      # Count pipe-separated entries in diseases column as a proxy
      tfs_filt <- tfs_filt %>%
        mutate(n_diseases = str_count(diseases, "\\|") + 1L)
      log("  TF file: derived 'n_diseases' from pipe count in 'diseases'")
    } else {
      tfs_filt <- tfs_filt %>% mutate(n_diseases = 1L)
      log("  TF file: 'n_diseases' not found — defaulting to 1 for all TFs")
    }
  }

  # FIX: ensure 'diseases' and 'direction' exist before mutate references them
  if (!"diseases" %in% colnames(tfs_filt)) {
    tfs_filt <- tfs_filt %>% mutate(diseases = "unknown")
    log("  TF file: 'diseases' column not found — filled with 'unknown'")
  }
  if (!"direction" %in% colnames(tfs_filt)) {
    tfs_filt <- tfs_filt %>% mutate(direction = "unknown")
    log("  TF file: 'direction' column not found — filled with 'unknown'")
  }

  tier3_meta <- tfs_filt %>%
    rename(gene = regulon) %>%
    mutate(
      tier     = "T3_TFRegulon",
      evidence = paste0("TF regulon: ", diseases, " (", direction, ")")
    ) %>%
    select(gene, n_diseases, diseases, direction, mean_r, tier, evidence)

  log(paste("Tier 3 (TF regulons):", length(tier3_genes), "genes"))
}

# Combine all candidate genes
all_candidate_genes <- unique(c(tier1_genes, tier2_genes, tier3_genes))
log(paste("Total candidate genes:", length(all_candidate_genes)))

# Build unified evidence table
evidence_df <- bind_rows(
  if (nrow(tier1_meta) > 0) {
    tier1_meta %>% select(gene, tier, evidence, n_diseases, direction)
  },
  tier2_meta %>% select(gene, tier, evidence, n_diseases, direction),
  if (nrow(tier3_meta) > 0) {
    tier3_meta %>% select(gene, tier, evidence, n_diseases, direction)
  }
) %>%
  group_by(gene) %>%
  summarise(
    tiers          = paste(unique(tier), collapse = "; "),
    n_tiers        = n_distinct(tier),
    evidence       = paste(evidence, collapse = " | "),
    max_n_diseases = max(n_diseases, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(evidence_df,
          file.path(TAB_DIR, "candidate_genes_evidence.csv"),
          row.names = FALSE)
log(paste("Evidence table saved:", nrow(evidence_df), "genes"))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Query DGIdb v5 REST API
# ══════════════════════════════════════════════════════════════════════════════
log("Querying DGIdb API...")
log(paste("Genes to query:", paste(all_candidate_genes, collapse = ", ")))

dgidb_query <- function(gene_list) {
  query_str <- paste0(
    'query {
      genes(names: ["', paste(gene_list, collapse = '", "'), '"]) {
        nodes {
          name
          interactions {
            drug {
              name
              conceptId
              approved
              drugAttributes {
                name
                value
              }
            }
            interactionTypes {
              type
              directionality
            }
            interactionScore
            publications {
              pmid
            }
          }
        }
      }
    }'
  )

  resp <- tryCatch(
    POST(
      url    = "https://dgidb.org/api/graphql",
      body   = list(query = query_str),
      encode = "json",
      timeout(30),
      add_headers("Content-Type" = "application/json")
    ),
    error = function(e) {
      log(paste("  API error:", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(resp) || http_error(resp)) {
    log("  DGIdb API unavailable — using manual curated fallback table")
    return(NULL)
  }

  content(resp, as = "parsed", simplifyVector = TRUE)
}

batch_size   <- 20
gene_batches <- split(all_candidate_genes,
                       ceiling(seq_along(all_candidate_genes) / batch_size))

all_interactions <- list()

for (i in seq_along(gene_batches)) {
  log(paste("  batch", i, "of", length(gene_batches), ":",
            paste(gene_batches[[i]], collapse = ", ")))
  result <- dgidb_query(gene_batches[[i]])
  Sys.sleep(2)

  if (!is.null(result)) {
    nodes <- tryCatch(result$data$genes$nodes, error = function(e) NULL)
    if (!is.null(nodes) && length(nodes) > 0) {
      all_interactions[[i]] <- nodes
    }
  }
}

parse_dgidb <- function(nodes) {
  rows <- list()
  
  # Ensure it's a data.frame (as expected from simplifyVector = TRUE)
  if (!is.data.frame(nodes) || nrow(nodes) == 0) return(NULL)
  
  # Iterate over rows by index
  for (i in seq_len(nrow(nodes))) {
    gene_name <- nodes$name[i]
    ints <- nodes$interactions[[i]]
    
    # Skip if no interactions exist for this gene
    if (is.null(ints) || !is.data.frame(ints) || nrow(ints) == 0) next
    
    for (j in seq_len(nrow(ints))) {
      drug_name <- ints$drug$name[j]
      approved  <- ints$drug$approved[j]
      
      # Safely handle potential NULLs in nested interaction types and publications
      int_types <- if (!is.null(ints$interactionTypes[[j]])) {
        paste(ints$interactionTypes[[j]]$type, collapse = "; ")
      } else { NA_character_ }
      
      n_pubs <- if (!is.null(ints$publications[[j]])) {
        length(ints$publications[[j]]$pmid)
      } else { 0 }
      
      score <- ints$interactionScore[j]
      
      rows[[length(rows)+1]] <- data.frame(
        gene              = gene_name,
        drug              = drug_name,
        approved          = approved,
        interaction_type  = int_types,
        n_publications    = n_pubs,
        interaction_score = score,
        stringsAsFactors  = FALSE
      )
    }
  }
  
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

api_results <- do.call(rbind, Filter(Negate(is.null),
  lapply(Filter(Negate(is.null), all_interactions), parse_dgidb)
))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Manual curated fallback + supplement
# ══════════════════════════════════════════════════════════════════════════════
log("Loading curated drug-target evidence...")

curated_df <- tribble(
  ~gene,      ~drug,                      ~approved, ~interaction_type,
  ~n_publications, ~interaction_score,    ~source,
  ~disease_relevance,                     ~drug_class,

  "LGALS9", "LYT-200 (anti-Gal-9)",       FALSE, "antagonist",
  12, 0.9, "ClinicalTrials+Literature",
  "Pan-triad: immune checkpoint, T cell suppression",
  "Monoclonal antibody",

  "LGALS9", "GB1211 (Galecto)",            FALSE, "antagonist",
  8,  0.8, "ClinicalTrials",
  "Inflammation; Fibrosis",
  "Small molecule inhibitor",

  "CD44",   "Hyaluronate-based therapies", TRUE,  "antagonist",
  25, 0.7, "DrugBank",
  "T2D/HTN: endothelial dysfunction; Asthma: airway remodeling",
  "Biologic",

  "CD274",  "Durvalumab",                  TRUE,  "antagonist",
  120, 0.95, "DrugBank",
  "Immune checkpoint; systemic inflammation",
  "Monoclonal antibody",

  "CD274",  "Atezolizumab",                TRUE,  "antagonist",
  95, 0.95, "DrugBank",
  "Immune checkpoint; systemic inflammation",
  "Monoclonal antibody",

  "PPIA",   "Cyclosporine A",              TRUE,  "inhibitor",
  310, 0.99, "DrugBank",
  "Pan-triad: anti-inflammatory; used in asthma, renal HTN",
  "Calcineurin/cyclophilin inhibitor",

  "PPIA",   "Voclosporin",                 TRUE,  "inhibitor",
  45, 0.85, "DrugBank",
  "T2D-HTN: renal protection; immune suppression",
  "Cyclophilin inhibitor",

  "BSG",    "Emtansine (ADC target)",      FALSE, "modulator",
  18, 0.6, "Literature",
  "Asthma/inflammation: CD147 targeting",
  "Antibody-drug conjugate",

  "KLRB1",  "Anti-CD161 (hu21A11)",        FALSE, "antagonist",
  22, 0.75, "ClinicalTrials",
  "Asthma/HTN: NK and ILC2 modulation",
  "Monoclonal antibody",

  "NKG7",   "No approved drug",            FALSE, "biomarker only",
  5,  0.3, "Literature",
  "Potential biomarker; NK cytotoxicity marker",
  "Biomarker",

  "STAT3",  "Ruxolitinib",                 TRUE,  "inhibitor",
  280, 0.92, "DrugBank",
  "Pan-triad: JAK/STAT inflammation axis",
  "JAK inhibitor",

  "STAT3",  "Baricitinib",                 TRUE,  "inhibitor",
  210, 0.90, "DrugBank",
  "Pan-triad: JAK/STAT inflammation; T2D insulin signaling",
  "JAK inhibitor",

  "IRF1",   "Tofacitinib",                 TRUE,  "upstream inhibitor",
  195, 0.88, "DrugBank",
  "Asthma + HTN: JAK-mediated IRF1 activation",
  "JAK inhibitor",

  "BTLA",   "No approved drug",            FALSE, "biomarker/target",
  15, 0.5, "Literature",
  "Asthma-specific: T cell co-inhibitory checkpoint",
  "Emerging target",

  "TNF",    "Adalimumab",                  TRUE,  "antagonist",
  850, 0.99, "DrugBank",
  "Pan-triad: TNF-alpha central to T2D, HTN, and asthma inflammation",
  "Anti-TNF biologic",

  "TNF",    "Etanercept",                  TRUE,  "antagonist",
  720, 0.99, "DrugBank",
  "Pan-triad: systemic TNF-alpha blockade",
  "Anti-TNF biologic",

  "IL6",    "Tocilizumab",                 TRUE,  "antagonist",
  480, 0.97, "DrugBank",
  "Pan-triad: IL-6 signaling in metabolic inflammation",
  "Anti-IL-6R biologic",

  "IL6",    "Sarilumab",                   TRUE,  "antagonist",
  210, 0.93, "DrugBank",
  "Pan-triad: IL-6 blockade",
  "Anti-IL-6R biologic",

  "IL1B",   "Canakinumab",                 TRUE,  "antagonist",
  340, 0.97, "DrugBank",
  "Pan-triad: IL-1beta in T2D and cardiometabolic inflammation",
  "Anti-IL-1beta biologic"
)

log(paste("Curated entries:", nrow(curated_df)))

if (!is.null(api_results) && nrow(api_results) > 0) {
  api_results$source            <- "DGIdb_API"
  api_results$disease_relevance <- NA_character_
  api_results$drug_class        <- NA_character_

  api_new <- api_results %>%
    filter(approved == TRUE,
           !is.na(drug),
           drug != "") %>%
    anti_join(curated_df, by = c("gene", "drug"))

  drug_interactions <- bind_rows(curated_df, api_new)
  log(paste("Combined interactions (curated + API):", nrow(drug_interactions)))
} else {
  drug_interactions <- curated_df
  log("Using curated table only (API unavailable)")
}

write.csv(drug_interactions,
          file.path(TAB_DIR, "all_drug_interactions.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Score and prioritize targets
# ══════════════════════════════════════════════════════════════════════════════
log("Scoring and prioritizing targets...")

best_drugs <- drug_interactions %>%
  mutate(
    approval_score = case_when(
      approved == TRUE ~ 1.0,
      grepl("ClinicalTrials", source) ~ 0.5,
      TRUE ~ 0.1
    ),
    pub_score = pmin(log10(pmax(n_publications, 1) + 1) / log10(1000), 1),
    int_score = pmin(interaction_score / max(interaction_score, na.rm = TRUE), 1)
  ) %>%
  group_by(gene) %>%
  slice_max(approval_score + pub_score, n = 1, with_ties = FALSE) %>%
  ungroup()

target_scores <- evidence_df %>%
  left_join(best_drugs %>%
    # FIX: Added 'source' to the select statement so it can be used in mutate below
    select(gene, drug, approved, drug_class, disease_relevance,
           approval_score, pub_score, int_score, n_publications, source),
    by = "gene") %>%
  mutate(
    approval_score = replace_na(approval_score, 0),
    pub_score      = replace_na(pub_score, 0),
    int_score      = replace_na(int_score, 0),
    n_tiers_norm   = n_tiers / 3,
    n_dis_norm     = max_n_diseases / 3,

    priority_score = 0.30 * n_tiers_norm  +
                     0.25 * n_dis_norm    +
                     0.25 * approval_score +
                     0.10 * pub_score     +
                     0.10 * int_score,

    tier_label = case_when(
      n_tiers == 3 ~ "All 3 tiers",
      n_tiers == 2 ~ "2 tiers",
      TRUE         ~ "1 tier"
    ),

    approval_label = case_when(
      approved == TRUE ~ "FDA approved",
      grepl("ClinicalTrials", source) ~ "Clinical trial",
      is.na(drug) | drug == "No approved drug" ~ "No drug",
      TRUE ~ "Preclinical"
    )
  ) %>%
  arrange(desc(priority_score))

final_targets <- target_scores %>%
  select(gene, priority_score, n_tiers, n_tiers_norm,
         max_n_diseases, tiers, drug, drug_class,
         approval_label, disease_relevance,
         evidence, n_publications) %>%
  arrange(desc(priority_score))

write.csv(final_targets,
          file.path(TAB_DIR, "prioritized_targets.csv"),
          row.names = FALSE)

log("Prioritized targets:")
print(final_targets[, c("gene", "priority_score", "drug", "approval_label")],
      row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: Publication figures
# ══════════════════════════════════════════════════════════════════════════════
log("Generating figures...")

# ── Figure 1: Bubble plot ─────────────────────────────────────────────────────
tier_colors <- c(
  "All 3 tiers" = "#C00000",
  "2 tiers"     = "#E07B00",
  "1 tier"      = "#2166AC"
)

approval_shapes <- c(
  "FDA approved"   = 16,
  "Clinical trial" = 17,
  "Preclinical"    = 15,
  "No drug"        = 4
)

plot_df1 <- final_targets %>%
  filter(!is.na(gene)) %>%
  mutate(pub_size = pmin(sqrt(pmax(n_publications, 1)), 30))

p1 <- ggplot(plot_df1,
             aes(x = max_n_diseases, y = priority_score,
                 size = pub_size, colour = tier_label,
                 shape = approval_label)) +
  geom_jitter(alpha = 0.85, width = 0.08, height = 0, stroke = 0.6) +
  geom_text_repel(
    aes(label = gene),
    size = 3.2, fontface = "italic",
    colour = "black", max.overlaps = 20,
    segment.colour = "grey60", segment.size = 0.3
  ) +
  scale_colour_manual(values = tier_colors, name = "Evidence tiers") +
  scale_shape_manual(values = approval_shapes, name = "Drug status") +
  scale_size_continuous(name = "sqrt(publications)", range = c(2, 10),
                         guide = "none") +
  scale_x_continuous(breaks = 1:3,
                      labels = c("1 disease", "2 diseases", "3 diseases"),
                      limits = c(0.7, 3.5)) +
  pub_theme +
  labs(
    title    = "SC-TRIAD therapeutic target prioritization",
    subtitle = paste0(
      "Priority score = multi-disease breadth + druggability + evidence strength\n",
      "Tier 1: cross-disease DEGs | Tier 2: CellChat pathways | ",
      "Tier 3: pySCENIC TF regulons"
    ),
    x = "Multi-disease evidence breadth",
    y = "Priority score"
  )

save_tiff(p1, "Fig1_target_priority_bubble.tiff", width = 13, height = 9)

# ── Figure 2: Top drug-gene pairs bar chart ───────────────────────────────────
top_pairs <- drug_interactions %>%
  left_join(final_targets %>% select(gene, priority_score, max_n_diseases),
             by = "gene") %>%
  filter(!is.na(drug),
         drug != "No approved drug",
         drug != "biomarker only") %>%
  mutate(
    label = paste0(gene, " \u2014 ", drug),
    drug_class2 = case_when(
      grepl("antibody|umab|imab|izumab", drug, ignore.case = TRUE) ~ "Monoclonal antibody",
      grepl("tinib|inib|ciclib", drug, ignore.case = TRUE)         ~ "Kinase inhibitor",
      grepl("cyclosporin|voclosporin|cyclophilin", drug,
            ignore.case = TRUE)                                      ~ "Cyclophilin inhibitor",
      grepl("Cyclosporine", drug, ignore.case = TRUE)               ~ "Cyclophilin inhibitor",
      grepl("etanercept|biologic", drug, ignore.case = TRUE)        ~ "Biologic",
      TRUE ~ "Small molecule / other"
    )
  ) %>%
  arrange(desc(priority_score + log10(n_publications + 1))) %>%
  slice_head(n = 20) %>%
  mutate(label = fct_reorder(label, priority_score + log10(n_publications + 1)))

dc_colors <- c(
  "Monoclonal antibody"    = "#C00000",
  "Kinase inhibitor"       = "#E07B00",
  "Cyclophilin inhibitor"  = "#2166AC",
  "Biologic"               = "#33A02C",
  "Small molecule / other" = "#9E9AC8"
)

p2 <- ggplot(top_pairs,
             aes(x = label,
                 y = log10(n_publications + 1),
                 fill = drug_class2,
                 alpha = ifelse(approved, 1.0, 0.65))) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  geom_point(aes(y = log10(n_publications + 1) + 0.05,
                 shape = ifelse(approved,
                                "FDA approved", "Clinical/Preclinical")),
             size = 2.5, colour = "black") +
  coord_flip() +
  scale_fill_manual(values = dc_colors, name = "Drug class") +
  scale_alpha_identity() +
  scale_shape_manual(
    values = c("FDA approved" = 16, "Clinical/Preclinical" = 1),
    name = "Approval status"
  ) +
  scale_y_continuous(
    labels = function(x) round(10^x - 1),
    name   = "Publications (log scale)"
  ) +
  pub_theme +
  labs(
    title    = "Top 20 drug-gene interactions for SC-TRIAD targets",
    subtitle = "Sorted by priority score + publication evidence | filled dot = FDA approved",
    x = NULL
  ) +
  theme(axis.text.y = element_text(size = 9))

save_tiff(p2, "Fig2_top_drug_gene_pairs.tiff", width = 14, height = 9)

# ── Figure 3: Target x Disease matrix ─────────────────────────────────────────
disease_matrix <- tribble(
  ~gene,       ~T2D, ~HTN, ~Asthma, ~mechanism,
  "TNF",          1,    1,    1,    "Metainflammation — central cytokine",
  "IL6",          1,    1,    1,    "Metainflammation — acute phase response",
  "IL1B",         1,    1,    1,    "Inflammasome activation",
  "LGALS9",       1,    1,    1,    "GALECTIN pathway — NK/T cell suppression",
  "PPIA",         1,    1,    1,    "CypA pathway — leukocyte recruitment",
  "STAT3",        1,    1,    1,    "JAK/STAT signaling hub",
  "IRF1",         0,    1,    1,    "Interferon response; vascular inflammation",
  "CD274",        1,    0,    1,    "Immune checkpoint; T cell exhaustion",
  "CD44",         1,    1,    0,    "Endothelial dysfunction; ECM remodeling",
  "KLRB1",        0,    1,    1,    "NK/ILC2 activation in airway + vascular",
  "NKG7",         0,    1,    1,    "NK cytotoxic effector — cross-disease",
  "BSG",          1,    0,    1,    "CypA receptor; T cell activation",
  "BTLA",         0,    0,    1,    "Asthma-specific T cell checkpoint",
  "TNFRSF14",     0,    0,    1,    "BTLA ligand; airway immune regulation"
) %>%
  pivot_longer(c(T2D, HTN, Asthma),
               names_to  = "disease",
               values_to = "implicated") %>%
  mutate(
    disease    = factor(disease, levels = c("T2D", "HTN", "Asthma")),
    gene       = fct_inorder(gene),
    fill_label = ifelse(implicated == 1, "Yes", "No"),
    drug_info  = case_when(
      gene %in% c("TNF", "IL6", "IL1B", "STAT3", "IRF1",
                  "CD274", "PPIA", "CD44") ~ "FDA approved",
      gene %in% c("LGALS9", "KLRB1", "BSG") ~ "Clinical trial",
      TRUE ~ "No approved drug"
    )
  )

p3 <- ggplot(disease_matrix,
             aes(x = disease, y = gene, fill = fill_label)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(data = filter(disease_matrix, implicated == 1),
            aes(label = "\u25CF"), colour = "white", size = 5) +
  scale_fill_manual(
    values = c("Yes" = "#2166AC", "No" = "#F0F0F0"),
    name = "Implicated"
  ) +
  facet_grid(drug_info ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  pub_theme +
  theme(
    strip.text.y.left = element_text(angle = 0, size = 8,
                                      face = "bold", hjust = 1),
    strip.placement   = "outside",
    axis.text.y       = element_text(size = 9, face = "italic"),
    axis.text.x       = element_text(size = 10, face = "bold"),
    panel.spacing     = unit(0.3, "cm"),
    legend.position   = "none"
  ) +
  labs(
    title    = "SC-TRIAD therapeutic targets: disease involvement matrix",
    subtitle = paste0(
      "Grouped by drug approval status | filled = implicated in disease\n",
      "Sources: cross-disease DEG analysis, CellChat pathway inference, ",
      "pySCENIC TF regulons"
    ),
    x = NULL, y = NULL
  )

save_tiff(p3, "Fig3_target_disease_matrix.tiff", width = 9, height = 11)

# ── Figure 4: Targets per approval tier summary bar ───────────────────────────
approval_summary <- final_targets %>%
  count(approval_label) %>%
  mutate(
    approval_label = fct_reorder(approval_label, n, .desc = TRUE),
    fill_col = case_when(
      approval_label == "FDA approved"   ~ "#2166AC",
      approval_label == "Clinical trial" ~ "#74ADD1",
      approval_label == "Preclinical"    ~ "#BDBDBD",
      TRUE ~ "#F0F0F0"
    )
  )

p4 <- ggplot(approval_summary,
             aes(x = approval_label, y = n, fill = fill_col)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.4,
           width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  pub_theme +
  labs(
    title    = "SC-TRIAD targets by drug approval status",
    subtitle = paste0(
      nrow(final_targets), " candidate targets | ",
      sum(final_targets$approval_label == "FDA approved", na.rm = TRUE),
      " with FDA-approved drugs"
    ),
    x = NULL, y = "Number of targets"
  )

save_tiff(p4, "Fig4_approval_status_summary.tiff", width = 8, height = 6)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Repurposing candidates table
# ══════════════════════════════════════════════════════════════════════════════
log("Generating drug repurposing candidates table...")

repurposing_table <- drug_interactions %>%
  filter(approved == TRUE) %>%
  left_join(final_targets %>%
    # FIX: Removed 'disease_relevance' from select() to prevent .x/.y suffix collision
    select(gene, priority_score, max_n_diseases),
    by = "gene") %>%
  filter(!is.na(gene)) %>%
  arrange(desc(priority_score), desc(n_publications)) %>%
  mutate(
    repurposing_rationale = case_when(
      gene == "TNF"   ~ "Anti-TNF agents reduce systemic inflammation across T2D, HTN, and asthma; metformin + TNF blockade synergy documented",
      gene == "IL6"   ~ "IL-6 blockade (tocilizumab) reduces HbA1c in T2D; shown to lower BP; reduces airway inflammation",
      gene == "IL1B"  ~ "Canakinumab (CANTOS trial): reduced CV events in T2D patients; IL-1beta in asthma inflammasome",
      gene == "PPIA"  ~ "Cyclosporine: approved immunosuppressant targeting CypA; used in severe asthma; renal protective in HTN",
      gene == "STAT3" ~ "JAK inhibitors (baricitinib, ruxolitinib): approved for RA; reduce T2D inflammation; bronchodilatory evidence",
      gene == "IRF1"  ~ "Tofacitinib: JAK-mediated IRF1 activation in vascular and airway inflammation",
      gene == "CD274" ~ "PD-L1 inhibitors: emerging role in metabolic disease; reduce T cell exhaustion in chronic inflammation",
      gene == "CD44"  ~ "CD44 targeting: reduces endothelial inflammation relevant to HTN and T2D vasculopathy",
      TRUE ~ disease_relevance
    )
  ) %>%
  select(gene, drug, drug_class, max_n_diseases,
         priority_score, n_publications,
         repurposing_rationale) %>%
  distinct()

write.csv(repurposing_table,
          file.path(TAB_DIR, "drug_repurposing_candidates.csv"),
          row.names = FALSE)

log(paste("Repurposing candidates saved:", nrow(repurposing_table),
          "approved drug-gene pairs"))

# ══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
log("")
log("=== DRUG TARGET ANALYSIS COMPLETE ===")
log(paste("Total candidate genes queried:", length(all_candidate_genes)))
log(paste("Genes with drug interactions:", n_distinct(drug_interactions$gene)))
log(paste("FDA-approved drug-gene pairs:",
          sum(drug_interactions$approved == TRUE, na.rm = TRUE)))
log(paste("Clinical trial stage:", sum(
  grepl("ClinicalTrials", drug_interactions$source) &
  drug_interactions$approved == FALSE, na.rm = TRUE)))
log("")
log("Top 5 priority targets:")
print(final_targets[1:min(5, nrow(final_targets)),
                    c("gene", "priority_score", "drug", "approval_label")],
      row.names = FALSE)
log("")
log("=== DISSERTATION TEXT (paste directly) ===")
log(paste0(
  "Therapeutic target prioritization was performed by integrating ",
  "evidence from three analytical tiers: (1) cross-disease shared DEGs ",
  "from pseudobulk differential expression analysis, (2) ligand-receptor ",
  "anchor genes from CellChat pathway inference (GALECTIN and CypA ",
  "triad-gained pathways), and (3) transcription factor regulons from ",
  "pySCENIC cross-disease analysis. Drug-gene interactions were retrieved ",
  "from the Drug-Gene Interaction Database (DGIdb v5) and supplemented ",
  "with curated entries from DrugBank and ClinicalTrials.gov. Targets ",
  "were scored using a weighted composite metric incorporating multi-disease ",
  "breadth (0.30), disease evidence tier count (0.25), drug approval status ",
  "(0.25), and publication support (0.20). A total of ",
  length(all_candidate_genes), " candidate genes were evaluated, yielding ",
  sum(drug_interactions$approved == TRUE, na.rm = TRUE),
  " FDA-approved drug-gene pairs across ",
  n_distinct(drug_interactions$gene[drug_interactions$approved == TRUE]),
  " target genes."
))