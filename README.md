# SeuratWorkflow
Note: If the user chooses not to define the number of PCs for PCA or Harmony (i.e., set `PCs_pca` or `PCs_harmony` as NULL), they must source the custom `FindMinimumPCs` function before running `SeuratWorkflow`, as this will automatically determine the minimum number of PCs. Please refer to the documentation for [FindMinimumPCs](https://github.com/vickinformatics/FindMinimumPCs).

## Description
The `SeuratWorkflow` function is a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline designed to preprocess, integrate, and visualize data using the Seurat workflow. It includes steps for normalization, feature selection, dimensionality reduction, and data integration using methods such as Harmony, Canonical Correlation Analysis (CCA), Reciprocal PCA (RPCA), and Joint PCA (JPCA). It also supports clustering and dimensionality reduction visualizations like t-SNE and UMAP.

## Arguments
- `seurat`: The Seurat object containing the single-cell RNA sequencing data.
- `batch_column`: The metadata column name in the Seurat object representing the batch or condition variable (default is "orig.ident").
- `nfeatures`: The number of variable features to select for downstream analysis (default is 2000).
- `seed`: A random seed for reproducibility (default is 123).
- `PCs_pca`: The number of principal components (PCs) to use for PCA. If not specified, it will be calculated using the FindMinimumPCs() function.
- `PCs_harmony`: The number of principal components (PCs) to use for Harmony integration. If not specified, it will be calculated using the `FindMinimumPCs` function.
- `PCs_cca`: The number of principal components (PCs) to use for CCA integration (default is 30).
- `PCs_rpca`: The number of principal components (PCs) to use for RPCA integration (default is 30).
- `PCs_jpca`: The number of principal components (PCs) to use for JPCA integration (default is 30).
- `run_integration`: Specifies whether integration methods such as Harmony, CCA, RPCA, and JPCA should be applied. When set to TRUE, integration will be performed using the specified integration methods (default is FALSE).
- `integration_method`: The integration method(s) to use. Options include "harmony", "CCAIntegration", "RPCAIntegration", and "JointPCAIntegration". Multiple methods can be specified (default is "harmony").
- `run_clustering`: Controls whether clustering is performed. If set to TRUE, clustering will be performed on the integrated data unless more than one integration method is specified. If integration is not performed, clustering will default to using PCA.
- `resolutions`: A numeric vector of resolutions to be used in clustering.

## Updates
<ins>March 5<sup>th</sup>, 2025</ins>
- The user can specify the number of principal components (`PCs_pca`, `PCs_harmony`, etc.) to use for each reduction method.
- Added `run_integration` option to allow users to choose whether to perform data integration.
- Introduced more integration methods, including Reciprocal PCA and Joint PCA (RPCA and JPCA), alongside the existing Harmony and CCA methods.
