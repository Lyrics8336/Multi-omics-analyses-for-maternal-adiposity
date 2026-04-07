# Multi-omics profiling of maternal adiposity: mechanistic drivers and enhanced risk stratification for perinatal and offspring outcomes

This repository contains the R source code and plot data for the visualization presented in our study. The project leverages longitudinal multi-omics profiling in the Tongji-Huaxi-Shuangliu Birth Cohort (THSBC) and Huaxi Birth Cohort (HBC) to systematically delineate the molecular signatures, mechanistic pathways, and clinical heterogeneity underlying maternal pre-pregnancy adiposity.

## Software Requirements:
Data visualization was performed in R (v4.4.1) using core packages including circlize (v0.4.17), ComplexHeatmap (v2.20.0), cowplot (v1.2.0), ggbreak (v0.1.6), ggExtra (v0.11.0), ggforce (v0.5.0), ggh4x (v0.3.1), ggplot2 (v4.0.1), ggpubr (v0.6.2), ggradar (v0.2), ggrepel (v0.9.6), ggridges (v0.5.7), ggsankey (v0.0.99999), ggsankeyfier (v0.1.8), networkD3 (v0.4.1), patchwork (v1.3.2), and pheatmap (v1.0.13).

## Installation Guide:


Package Installation
From an `R` session, type:
```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")

remotes::install_github("ricardo-bion/ggradar")
remotes::install_github("davidsjoberg/ggsankey")

cran_list <- c("cowplot", "ggbreak", "ggExtra", "ggforce", "ggh4x", 
               "ggplot2", "ggpubr", "ggrepel", "ggridges", "ggsankeyfier", 
               "networkD3", "patchwork", "pheatmap")
install.packages(cran_list)
```

## Repository Structure
The repository is organized by figure order as presented in the manuscript:

- Figure2.R: Beta-diversity, PCA, Alpha-diversity (Richness/Shannon), and F/B ratio analysis.

- Figure3.R: Circular heatmaps for serum metabolites, Sankey diagrams of metabolic pathways, and cluster trajectories.

- Figure4.R: Taxonomic heatmaps (Genus/Species) and association of multi-omics with clinical outcomes.

- Figure5.R: Logistic/Linear regression models, Polar plots for outcome proportion, and Chord diagrams for mediation analysis.

- Figure6.R: Machine learning performance (R²), Radar plots for phenotypic vs. genetic mismatch, and health outcome frequency analysis.

/data/: Contains the raw data sheets required for all visualization.


## Instructions for use
The scripts are organized and numbered according to the logical flow of the visualization. For full reproduction, please execute them in numerical order. 

Running the scripts will generate high-resolution visualization files or plot them directly in the R IDE.

Total run time for all figures (Fig 2-6): < 5 minutes on a standard desktop computer.


