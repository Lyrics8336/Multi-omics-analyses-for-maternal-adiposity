# Multi-omics profiling of maternal adiposity: molecular signatures and risk stratification for maternal and offspring outcomes

This repository contains the R source code for the computational analyses and visualizations of our study on gut microbial and serum metabolomic signatures of maternal adiposity and their relevance to maternal and offspring outcomes.

## Project Overview

This study integrates longitudinal serum metabolomics, gut microbiome profiling (16S rRNA gene sequencing, fungal ITS sequencing and shotgun metagenomics), host genetic data and clinical information from the THSBC discovery cohort (*n* = 5,212) and the independent HBC cohort (*n* = 393). Omics measurements span early (T1), middle (T2) and late (T3) pregnancy; availability varies by platform and timepoint.

The analyses characterize signatures associated with pre-pregnancy body mass index (pBMI) and overweight/obesity (OWO), examine their associations with maternal and offspring outcomes, and evaluate incremental outcome discrimination beyond clinical variables and observed pBMI. Additional analyses develop multi-omics-inferred pBMI, which identified women with divergent maternal and offspring risk profiles despite similar observed pBMI and captured risk information not fully reflected by conventional anthropometric measures. 

## Repository Structure

Scripts are grouped into computational analyses, visualizations and shared utilities. The numbering follows the main analytical workflow.

### 1. Core Computational Analysis (`1_analysis/`)

* **`1_multiomics_profiles.R`**: Evaluates factors associated with overall omics variation, gut microbial diversity and the Firmicutes/Bacteroidetes ratio, and multi-omics distances across pBMI groups and gestational timepoints (Figure 2).
* **`2_metabolomic_analysis.R`**: Characterizes pre-pregnancy OWO-associated serum metabolites using regression and mixed-effects models; prioritizes metabolites using repeated cross-validated OPLS-DA; examines longitudinal association patterns, clinical correlates, cross-cohort consistency and determinants of metabolomic variation (Figure 3).
* **`3_microbial_analysis.R`**: Characterizes pre-pregnancy OWO-associated bacterial and fungal features and their clinical correlates in the THSBC and HBC (Figure 4).
* **`4_outcome_associations_mediation.R`**: Evaluates pBMI-outcome associations and performs exploratory score-level and member-level mediation analyses using previously derived T1 DIABLO bacterial-metabolite modules for gestational diabetes mellitus (GDM) and hypertensive disorders of pregnancy (HDP) (Figure 5).

### 2. Predictive Modeling and Adiposity Phenotyping (`1_analysis/`)

* **`5_multiomics_outcome_prediction.R`**: Uses repeated nested cross-validation to compare clinical, clinical-plus-pBMI and omics-augmented outcome models. Evaluates discrimination, paired AUROC differences and feature-selection stability across gestational sampling windows (Figure 5).
* **`6_multiomics_pbmi_modeling.R`**: Develops and internally evaluates LASSO models for pBMI using clinical covariates, the BMI polygenic risk score (PRS), metabolites and microbial genera. Exports out-of-fold predictions and model-performance summaries (Figure 6).
* **`7_pgm_outcome_associations.R`**: Constructs PGM scores from observed or multi-omics-inferred pBMI and BMI-PRS, and evaluates outcome incidence and adjusted associations across mismatch groups, including subgroup and offspring-sex interaction analyses (Figure 6).

### 3. Visualization Scripts (`2_visualization/`)

* **`plot_fig2.R` to `plot_fig6.R`**: Generate plots corresponding to Figures 2-6 from processed data and analytical outputs. Some statistical calculations are included in the visualization scripts; for example, the PCA and PERMANOVA analyses for Figure 2b are implemented in `plot_fig2.R`.

### 4. Utilities

* **`01_help_func.R`**: Shared functions for statistical analysis, data transformations and other operations used by the analysis scripts.
* **`14_module_mediation_helpers.R`**: Implements DIABLO score-level mediation analyses.
* **`14_pair_mediation_helpers.R`**: Implements bacterium-metabolite pair and OWO-single-feature mediation analyses.
* **`14_omics_cohort_helpers.R`**: Loads the full T1 16S cohort and projects previously derived DIABLO bacterial scores for mediation analyses.

## Software Requirements

The scripts use R. Main package dependencies include:

| Purpose | Packages |
| --- | --- |
| Data handling and export | `tidyverse`, `data.table`, `readxl`, `writexl` |
| Statistical and multivariate analysis | `vegan`, `permute`, `pairwiseAdonis`, `lme4`, `lmerTest`, `survey`, `RNOmni`, `ropls`, `latrend`, `mediation`, `performance` |
| Predictive modeling and parallel computing | `caret`, `glmnet`, `pROC`, `PRROC`, `doParallel`, `foreach` |
| Visualization | `ggplot2`, `patchwork`, `cowplot`, `ggpubr`, `ggh4x`, `ggforce`, `ggrepel`, `ggridges`, `ggnewscale`, `circlize`, `ComplexHeatmap`, `ggsankey`, `ggradar`, `colorspace`, `cols4all`, `MetBrewer`, `scales` |

Install the packages required by each script before execution. Package calls within the scripts and helper functions provide the detailed dependencies. Metabolite pathway enrichment was performed using the MetaboAnalyst 6.0 web platform; `plot_fig3.R` reads its exported enrichment results.

## Data Availability

The raw fungal ITS, bacterial 16S rRNA and shotgun metagenomic sequencing data have been deposited in the [Genome Sequence Archive (GSA), National Genomics Data Center](https://ngdc.cncb.ac.cn/gsa/) under accession numbers **CRA014764**, **CRA014766** and **CRA014529**, respectively. Other raw genotype and phenotype data are available from the corresponding authors (X.-F.P. or A.P.) upon reasonable request, as described in the manuscript's Data Availability statement.

## License

The code is distributed under the MIT License. See [LICENSE](LICENSE) for details.
