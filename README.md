# An agent-interactive human immune cell noise landscape reveals an unrecognized functional dimension of health and diseases
This repository contains the analysis pipeline and code for the study "An agent-interactive human immune cell noise landscape reveals an unrecognized functional dimension of health and diseases"
# Data
Healthy PBMC scRNA-seq data were composed of two population-scale cohorts, AIDA Phase 1 Data Freeze v1 and OneK1K.
PBMC scRNA-seq data from 7 diseases including COVID-19, influenza, sepsis, post-COVID-19, SLE, RA, and NSCLC were collected from 6 previously published studies.
# Pipeline
# Healthy data processing
To ensure cross-cohort comparability, cell-type annotations were harmonized by consolidating the original fine-grained subtypes into unified major immune cell categories, including monocytes, CD4+ T cells, CD8+ T cells, other T cells, B cells, NK cells, dendritic cells, and other cells. For each donor and each harmonized cell type, QC procedures were applied to ensure data reliability. At the cellular level, cells were retained only if they satisfied the following criteria: (i) the number of detected genes ranged from 300 to 6,000, thereby excluding low-quality cells and potential doublets. (ii) the proportion of mitochondrial gene counts was ≤ 20%, with cells exceeding this threshold removed as low-viability or apoptotic cells. At the gene level, only genes detected in at least 5 cells were retained to ensure robust expression estimation. UMI count matrices were normalized using the R package SCTransform (v0.4.1).
# Expression-decoupled noise quantification framework
To distinguish biological noise from technical noise, the squared coefficient of variation (CV²) was modeled as a function of the mean expression level (μ). The observed CV² (CVobs²) was decomposed into technical (CVtech²) and biological (CVbio²) components
# Multi-metric calibration factor (α_total) for correction and validation of correction effectiveness

# Healthy baseline noise atlas analyses
# Disease data processing
# Disease data analysis
# Interactive multi-agent system
We developed a multi-agent noise analysis system based on the model context protocol (MCP), comprising four functional layers.
