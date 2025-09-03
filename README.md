# SeuratWorkflow
Note: If the user chooses not to define the number of PCs for PCA or Harmony (i.e., setting `PCs_pca` or `PCs_harmony` as NULL), they must source the custom `FindMinimumPCs` function before running `SeuratWorkflow`, as this will automatically determine the minimum number of PCs. Please refer to the documentation for [FindMinimumPCs](https://github.com/vickinformatics/FindMinimumPCs).

## Description
The `SeuratWorkflow` function is a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline designed to preprocess, integrate, and visualize data using the Seurat workflow. It includes steps for normalization, feature selection, dimensionality reduction, and data integration using methods such as Harmony, Canonical Correlation Analysis (CCA), Reciprocal PCA (RPCA), and Joint PCA (JPCA). It also supports clustering and dimensionality reduction visualizations like t-SNE and UMAP.

## Arguments
- `seurat`: A Seurat object containing the single-cell RNA sequencing data.
- `batch_column`: The metadata column name in the Seurat object representing the batch or condition variable (default is "orig.ident").
- `nfeatures`: The number of variable features to select for downstream analysis (default is 2000).
- `variables_to_regress`: A character vector of metadata variables to regress out during scaling (e.g., c("percent.mt", "nCount_RNA")). If NULL, no regression is performed (default is NULL).
- `seed`: A random seed for reproducibility (default is 123).
- `PCs_pca`: The number of principal components (PCs) to use for PCA. If not specified, it will be calculated using the `FindMinimumPCs` function.
- `PCs_harmony`: The number of principal components (PCs) to use for Harmony integration. If not specified, it will be calculated using the `FindMinimumPCs` function.
- `PCs_cca`: The number of principal components (PCs) to use for CCA integration (default is 30).
- `PCs_rpca`: The number of principal components (PCs) to use for RPCA integration (default is 30).
- `PCs_jpca`: The number of principal components (PCs) to use for JPCA integration (default is 30).
- `run_integration`: Specifies whether integration methods such as Harmony, CCA, RPCA, and JPCA should be applied. When set to TRUE, integration will be performed using the specified integration methods (default is FALSE).
- `integration_method`: The integration method(s) to use. Options include "harmony", "CCAIntegration", "RPCAIntegration", and "JointPCAIntegration". Multiple methods can be specified (default is "harmony").
- `run_tSNE_UMAP`: Specifies whether tSNE and UMAP dimensionality reductions should be performed (default is FALSE).
- `check_duplicates`: Checks for duplicated cell names in `RunTSNE()` (default is FALSE).
- `run_clustering`: Controls whether clustering is performed. If set to TRUE, clustering will be performed on the integrated data unless more than one integration method is specified. If integration is not performed, clustering will default to using PCA.
- `resolutions`: A numeric vector of resolutions to be used in clustering.

## Updates
<ins>September 3<sup>rd</sup>, 2025</ins>
- Added `check_duplicates` argument to allow control over the `RunTSNE()` parameter `check_duplicates`.

<ins>April 3<sup>rd</sup>, 2025</ins>
- Added `variables_to_regress` argument to allow users to regress out unwanted sources of variation (e.g., mitochondrial content, RNA count) during scaling via `ScaleData(vars.to.regress = ...)`.

<ins>March 5<sup>th</sup>, 2025</ins>
- The user can specify the number of principal components (`PCs_pca`, `PCs_harmony`, etc.) to use for each reduction method.
- Added `run_integration` option to allow users to choose whether to perform data integration.
- Introduced more integration methods, including Reciprocal PCA and Joint PCA (RPCA and JPCA), alongside the existing Harmony and CCA methods.
- Added `run_tSNE_UMAP` option to allow users to choose whether to perform tSNE and UMAP dimensionality reductions.

## References
<ins>Seurat</ins>
- Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, Srivastava A, Molla G, Madad S, Fernandez-Granda C, Satija R (2023). “Dictionary learning for integrative, multimodal and scalable single-cell analysis.” _Nature Biotechnology_. doi:10.1038/s41587-023-01767-y, https://doi.org/10.1038/s41587-023-01767-y.
- Hao Y, Hao S, Andersen-Nissen E, III WMM, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zagar M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LB, Yeung B, Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R (2021). “Integrated analysis of multimodal single-cell data.” _Cell_. doi:10.1016/j.cell.2021.04.048, https://doi.org/10.1016/j.cell.2021.04.048.
- Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, III WMM, Hao Y, Stoeckius M, Smibert P, Satija R (2019). “Comprehensive Integration of Single-Cell Data.” _Cell_, **177**, 1888-1902. doi:10.1016/j.cell.2019.05.031, https://doi.org/10.1016/j.cell.2019.05.031.
- Butler A, Hoffman P, Smibert P, Papalexi E, Satija R (2018). “Integrating single-cell transcriptomic data across different conditions, technologies, and species.” _Nature Biotechnology_, **36**, 411-420. doi:10.1038/nbt.4096, https://doi.org/10.1038/nbt.4096.
- Satija R, Farrell JA, Gennert D, Schier AF, Regev A (2015). “Spatial reconstruction of single-cell gene expression data.” _Nature Biotechnology_, **33**, 495-502. doi:10.1038/nbt.3192, https://doi.org/10.1038/nbt.3192.

<ins>Harmony</ins>
- Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh P, Raychaudhuri S (2019). “Fast, sensitive and accurate integration of single-cell data with Harmony.” Nature Methods, 16(12), 1289–1296. https://doi.org/10.1038/s41592-019-0619-0.
