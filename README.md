# SC-TRIAD-Immune-Landscapes

**Construction and cell-type annotation of the SC-TRIAD cross-disease PBMC atlas**

This repository contains the official code for the SC-TRIAD (Single-Cell integration of T2D, HTN, and Asthma) project, which establishes a cross-disease PBMC atlas and validates systemic immune reprogramming across local tissue microenvironments (hypertensive kidney, T2D pancreas, and asthmatic airway).

## Overview

The T2D–HTN–Asthma triad represents a highly prevalent yet mechanistically fragmented disease cluster. In this project, we integrated single-cell transcriptomics across 154,797 PBMCs to resolve shared and disease-specific immune trajectories. 

Key computational components include:
* **Single-Cell PBMC Integration**: Joint integration (Harmony, SCTransform) and manual annotation of major immune lineages across disease conditions.
* **Transcriptional Reprogramming**: Pseudobulk differential expression profiling uncovering a unique severity-independent asthma trajectory and an HTN-asthma convergent NK cytotoxic signature.
* **Intercellular Communication**: `CellChat` mapping of pan-triad signaling nodes (e.g., CypA and GALECTIN).
* **Pseudotime Trajectory & TF Regulon Inference**: Monocyte-to-DC maturation dynamics, T-cell effector polarization, and `pySCENIC` regulatory network mapping defining a shared myeloid inflammatory axis.
* **Spatial Transcriptomics Deconvolution**: Robust spacial validation using Xenium and Visium (RCTD) across three target tissues (asthmatic airway, hypertensive kidney, T2D pancreas).
* **Drug Repurposing**: In silico prioritization of FDA-approved therapeutics (e.g., Cyclosporine A, checkpoint modulators, and proteasome inhibitors) targeting cross-disease hubs.

## Repository Structure

The analysis pipeline is categorized into three major components:

```
├── 01_data_processing
│   ├── 01_qc_pbmc.R
│   ├── 02_norm_cluster_annotate.R
│   ├── 03_integration.R
│   └── 04_patch_celltype.R
├── 02_analysis
│   ├── fig1_atlas
│   │   └── 05_fig1_atlas_umap.R
│   ├── fig2_deg
│   │   └── 06_fig2_deg_analysis.R
│   ├── fig3_pathway
│   │   └── 07_fig3_pathway_enrichment.R
│   ├── fig4_cellchat
│   │   └── 08_fig4_cellchat_analysis.R
│   ├── fig5_pseudotime
│   │   └── 09_fig5_trajectory_analysis.R
│   ├── fig6_pyscenic
│   │   └── 10_fig6_pyscenic_regulons.R
│   ├── fig7_xenium
│   │   └── 11_fig7_xenium_spatial.R
│   ├── fig8_kidney_rctd
│   │   └── 12_fig8_kidney_spatial.R
│   ├── fig9_pancreas_rctd
│   │   └── 13_fig9_pancreas_spatial.R
│   ├── fig10_convergence
│   │   └── 14_fig10_cross_tissue.R
│   └── fig11_drug_targets
│       └── 15_fig11_drug_targets.R
└── 03_supplementary
    ├── 16_supp_s1_qc.R
    ├── 17_supp_s2_annotation.R
    ├── 18_supp_s3_cellchat.R
    └── 19_supp_s5_label_transfer.R
```

## Data Availability
Raw sequencing data for the PBMC cohorts are available under GEO accessions GSE255566, GSE212953, and GSE288147. Spatial datasets were obtained from GSE269354 (Xenium), GSE211785 (Visium), and GSE264331 (Visium). All processed matrices and Seurat objects required to reproduce the figures will be deposited in Zenodo.

## Dependencies
- R (>= 4.2.0)
- Seurat v5.0+
- Harmony
- SCTransform v2
- pySCENIC
- CellChat v2
- spacexr (RCTD)

## Citation
*Manuscript in preparation.*
