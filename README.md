# SeuratWorkflow
ATTENTION: The user needs to source the custom FindMinimumPCs function before running the SeuratWorkflow, as this function is required to compute the minimum PCs. Please see the documentation for FindMinimumPCs (https://github.com/vickinformatics/FindMinimumPCs). 

## Description
The SeuratWorkflow function is a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline designed to preprocess, integrate, and visualize data using the Seurat workflow. It includes steps for normalization, feature selection, dimensionality reduction, and data integration using methods such as Harmony, Canonical Correlation Analysis (CCA), Reciprocal PCA (RPCA), and Joint PCA (JPCA). It also supports clustering and dimensionality reduction visualizations like t-SNE and UMAP.

## Arguments
- **seurat**: The Seurat object containing your single-cell RNA sequencing data.
- **batch_column**: The metadata column name in the Seurat object representing the batch or condition variable (default is "orig.ident").
- **nfeatures**: The number of variable features to select for downstream analysis (default is 2000).
- **seed**: A random seed for reproducibility (default is 123).
- **PCs_pca**: The number of principal components (PCs) to use for PCA. If not specified, it will be calculated using the FindMinimumPCs() function.
- **PCs_harmony**: The number of principal components (PCs) to use for Harmony integration. If not specified, it will be calculated using the FindMinimumPCs() function.
- **PCs_cca**: The number of principal components (PCs) to use for Canonical Correlation Analysis (CCA) integration (default is 30).
- **PCs_rpca**: The number of principal components (PCs) to use for Reciprocal PCA integration (default is 30).
- **PCs_jpca**: The number of principal components (PCs) to use for Joint PCA integration (default is 30).
- **run_integration**: Specifies whether integration methods such as Harmony, CCA, RPCA, and Joint PCA should be applied. When set to TRUE, integration will be performed using the specified integration methods (default is FALSE).
- **integration_method**: The integration method(s) to use. Options include "harmony", "CCAIntegration", "RPCAIntegration", and "JointPCAIntegration". Multiple methods can be specified (default is "harmony").
- **run_clustering**: Controls whether clustering is performed. If set to TRUE, clustering will be performed on the integrated data unless more than one integration method is specified. If integration is not performed, clustering will default to using PCA.
- **resolutions**: A numeric vector of resolutions to be used in clustering.

## Updates
