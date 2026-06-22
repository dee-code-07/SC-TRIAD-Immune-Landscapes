# SC-TRIAD: Single-cell and spatial integration of Type 2 Diabetes, Hypertension, and Asthma

**Authors:** Deeksha H and Budheswar Dehury*  
**Affiliation:** Department of Bioinformatics, Manipal School of Life Sciences, Manipal Academy of Higher Education, Manipal-576104, India  
**Correspondence:** budheswar.dehury@manipal.edu

---

## Overview
This repository contains the official codebase for the SC-TRIAD (Single-Cell integration of T2D, HTN, and Asthma) project. By integrating large-scale PBMC single-cell RNA sequencing with matched spatial transcriptomics (Visium and Xenium), this study establishes a cross-disease immune atlas to map shared and disease-specific immune reprogramming across cardiometabolic and allergic respiratory conditions.

### Key Findings
- **Convergent Inflammatory States:** Identified shared myeloid inflammatory programs across T2D, HTN, and Asthma.
- **Disease-Specific Trajectories:** Resolved a unique severity-independent asthma trajectory and an HTN-asthma convergent NK cytotoxic signature.
- **Spatial Validation:** Validated single-cell observations within target local tissue microenvironments (hypertensive kidney, T2D pancreas, and asthmatic airway).
- **Drug Repurposing Hubs:** Prioritized FDA-approved therapeutics targeting cross-disease hubs, emphasizing modulators of the CypA and GALECTIN signaling pathways.

---

## Data Availability
Raw sequencing data and spatial datasets are publicly available from the Gene Expression Omnibus (GEO):
- **T2D PBMCs:** GSE255566
- **HTN PBMCs:** GSE212953
- **Asthma PBMCs:** GSE288147
- **Asthmatic Airway (Xenium):** GSE269354
- **Hypertensive Kidney (Visium):** GSE211785
- **T2D Pancreas (Visium):** GSE264331

All processed matrices and Seurat objects required to reproduce the figures will be deposited in Zenodo.

---

## Repository Structure

The scripts are organized sequentially to reflect the workflow presented in the manuscript:

### 📂 `01_preprocessing/`
Clean, step-by-step pipelines for initial data processing.
- **`scrna/`**: Quality control, normalization, Harmony integration, and cluster annotation.
- **`pyscenic/`**: Export to loom, Gene Regulatory Network (GRN) inference, motif enrichment, and AUCell scoring.
- **`spatial/`**: Visium and Xenium processing, label transfer, and RCTD spatial deconvolution.
- **`drug_targets/`**: Pipeline for the identification and scoring of repurposable therapeutic candidates.

### 📂 `02_analysis/`
The core analysis scripts used to generate primary figures:
- **`fig1_atlas/`**: UMAP visualizations and cell-type composition analysis.
- **`fig2_deg/`**: Differential gene expression across cohorts.
- **`fig3_pathway/`**: Pathway and functional enrichment analysis.
- **`fig4_cellchat/`**: Intercellular communication mapping.
- **`fig5_pseudotime/`**: Trajectory inference and lineage analysis.
- **`fig6_pyscenic/`**: Regulatory network and TF regulon mapping.
- **`fig7_xenium/`**: Xenium airway spatial transcriptomics analysis.
- **`fig8_kidney_rctd/`**: Kidney RCTD spatial mapping.
- **`fig9_pancreas_rctd/`**: Pancreas RCTD spatial mapping.
- **`fig10_convergence/`**: Cross-tissue convergence evaluation.
- **`fig11_drug_targets/`**: Network-based drug target prioritization.

### 📂 `03_supplementary/`
Additional validation scripts covering QC metrics, label transfer accuracy, and secondary analyses.

---

## Dependencies
To run these scripts, you will need **R (>= 4.2.0)** and **Python 3.9+** with the following key packages:
- **Core**: `Seurat (v5)`, `tidyverse`, `patchwork`, `ggplot2`, `Harmony`, `SCTransform v2`
- **Network & Signaling**: `CellChat v2`, `pySCENIC`
- **Spatial**: `spacexr` (RCTD)

---

## Usage
Scripts are numbered sequentially. For best results, follow the order:
1. Run `01_preprocessing` modules to generate the integrated objects and baseline models.
2. Execute `02_analysis` folders to reproduce manuscript figures.
3. Consult `03_supplementary` for deeper metric validation.
