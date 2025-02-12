# SeuratWorkflow
ATTENTION: The user needs to source the custom FindMinimumPCs function before running the SeuratWorkflow, as this function is required to compute the minimum PCs. Please see the documentation for FindMinimumPCs. 

## Description
The SeuratWorkflow function is a custom analysis pipeline designed to preprocess, integrate, and visualize single-cell RNA sequencing (scRNA-seq) data using Seurat. It provides a flexible and comprehensive approach to analyzing scRNA-seq datasets, using normalization, feature selection, dimensionality reduction, and integration methods such as Harmony and Canonical Correlation Analysis (CCA). It also includes clustering, dimensionality reduction visualizations (t-SNE and UMAP), and allows for the handling of batch effects during integration. This function is for researchers who wish to integrate multiple datasets from different experimental batches or need a streamlined tool to help them perform common analysis steps with ease.

## Input Parameters
- **seurat**: The Seurat object containing your single-cell RNA sequencing data.
- **batch_column**: The metadata column name in the Seurat object representing the batch or condition variable (default is "orig.ident").
- **nfeatures**: The number of variable features to select for downstream analysis (default is 2000).
- **seed**: A random seed for reproducibility (default is 123).
- **integration_method**: The integration method(s) to use. Options are "harmony" and "cca.integration" (default is "harmony").
- **min.pc_cca**: Minimum number of PCs to use for CCA integration. This parameter should be defined by the user, as the CCA integration reduction does not provide a standard deviation for automatic selection.
- **run_clustering**: A boolean parameter to control whether clustering should be run. By default, clustering is performed unless both Harmony and CCA integration methods are selected.
- **resolutions**: A vector of resolutions to be used in clustering.
