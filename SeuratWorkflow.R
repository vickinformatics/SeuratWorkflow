#' SeuratWorkflow
#' 
#' The `SeuratWorkflow` function is a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline designed to preprocess, integrate, and visualize data using the Seurat workflow.
#' It includes steps for normalization, feature selection, dimensionality reduction, and data integration using methods such as Harmony, Canonical Correlation Analysis (CCA), Reciprocal PCA (RPCA), and Joint PCA (JPCA).
#' It also supports clustering and dimensionality reduction visualizations like t-SNE and UMAP.
#'
#' @param seurat A Seurat object containing the single-cell RNA sequencing data.
#' @param batch_column The metadata column name in the Seurat object representing the batch or condition variable (default is "orig.ident").
#' @param nfeatures The number of variable features to select for downstream analysis (default is 2000).
#' @param seed A random seed for reproducibility (default is 123).
#' @param PCs_pca The number of principal components (PCs) to use for PCA. If not specified, it will be calculated using the `FindMinimumPCs` function.
#' @param PCs_harmony The number of principal components (PCs) to use for Harmony integration. If not specified, it will be calculated using the `FindMinimumPCs` function.
#' @param PCs_cca The number of principal components (PCs) to use for CCA integration (default is 30).
#' @param PCs_rpca The number of principal components (PCs) to use for RPCA integration (default is 30).
#' @param PCs_jpca The number of principal components (PCs) to use for JPCA integration (default is 30).
#' @param run_integration Specifies whether integration methods such as Harmony, CCA, RPCA, and JPCA should be applied. When set to TRUE, integration will be performed using the specified integration methods (default is FALSE).
#' @param integration_method The integration method(s) to use. Options include "harmony", "CCAIntegration", "RPCAIntegration", and "JointPCAIntegration". Multiple methods can be specified (default is "harmony").
#' @param run_tSNE_UMAP Specifies whether tSNE and UMAP dimensionality reductions should be performed (default is FALSE).
#' @param run_clustering Controls whether clustering is performed. If set to TRUE, clustering will be performed on the integrated data unless more than one integration method is specified. If integration is not performed, clustering will default to using PCA.
#' @param resolutions A numeric vector of resolutions to be used in clustering.
#' 
#' @return A Seurat object with the results of preprocessing, integration, and clustering.
#'
#' @note If the user chooses not to define the number of PCs for PCA or Harmony (i.e., setting `PCs_pca` or `PCs_harmony` to NULL), they must source the custom `FindMinimumPCs` function
#' before running `SeuratWorkflow`, as this will automatically determine the minimum number of PCs.
#' 
#' @citations
#' Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, Srivastava A, Molla G, Madad S, Fernandez-Granda C, Satija R (2023). "Dictionary learning for integrative, multimodal and scalable single-cell analysis." *Nature Biotechnology*. doi:10.1038/s41587-023-01767-y, https://doi.org/10.1038/s41587-023-01767-y.
#' 
#' Hao Y, Hao S, Andersen-Nissen E, III WMM, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zagar M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LB, Yeung B, Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R (2021). "Integrated analysis of multimodal single-cell data." *Cell*. doi:10.1016/j.cell.2021.04.048, https://doi.org/10.1016/j.cell.2021.04.048.
#' 
#' Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, III WMM, Hao Y, Stoeckius M, Smibert P, Satija R (2019). "Comprehensive Integration of Single-Cell Data." *Cell*, 177, 1888-1902. doi:10.1016/j.cell.2019.05.031, https://doi.org/10.1016/j.cell.2019.05.031.
#' 
#' Butler A, Hoffman P, Smibert P, Papalexi E, Satija R (2018). "Integrating single-cell transcriptomic data across different conditions, technologies, and species." *Nature Biotechnology*, 36, 411-420. doi:10.1038/nbt.4096, https://doi.org/10.1038/nbt.4096.
#' 
#' Satija R, Farrell JA, Gennert D, Schier AF, Regev A (2015). "Spatial reconstruction of single-cell gene expression data." *Nature Biotechnology*, 33, 495-502. doi:10.1038/nbt.3192, https://doi.org/10.1038/nbt.3192.
#' 
#' Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh P, Raychaudhuri S (2019). "Fast, sensitive and accurate integration of single-cell data with Harmony." *Nature Methods*, 16(12), 1289–1296. https://doi.org/10.1038/s41592-019-0619-0.
#'
#' @author Vicki Do
#' @lastUpdated 2025-3-5
#'
#' @examples
#' seurat_obj <- SeuratWorkflow(seurat = seurat_obj, run_integration = TRUE, integration_method = c("harmony", "CCAIntegration), run_clustering = FALSE)
#' seurat_obj <- SeuratWorkflow(seurat = seurat_obj, PCs_pca = 30, run_integration = FALSE, run_clustering = TRUE, resolutions = c(0.4, 0.6, 0.8))

SeuratWorkflow <- function(seurat,
                           batch_column = "orig.ident",
                           nfeatures = 2000,
                           seed = 123,
                           PCs_pca = NULL, 
                           PCs_harmony = NULL,
                           PCs_cca = 30,
                           PCs_rpca = 30,
                           PCs_jpca = 30,
                           run_integration = FALSE,
                           integration_method = "harmony",
                           run_tSNE_UMAP = FALSE,
                           run_clustering = FALSE,
                           resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0)) {
  
  # Normalization, scaling, and PCA
  seurat <- seurat %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
    ScaleData() %>%
    RunPCA()
  print(ElbowPlot(seurat))
  
  # Use provided PCs_pca or calculate using FindMinimumPCs (see documentation)
  if (is.null(PCs_pca)) {
    PCs_pca <- FindMinimumPCs(seurat, reduction_type = "pca")
  }
  
  # Run Harmony integration if specified and run_integration is TRUE
  if (run_integration && "harmony" %in% integration_method) {
    set.seed(seed)
    seurat <- seurat %>% 
      RunHarmony(batch_column, plot_convergence = TRUE, reduction = "pca", reduction.save = "harmony")
  }
  
  # Use provided PCs_harmony or calculate using FindMinimumPCs (see documentation)
  if (is.null(PCs_harmony) && run_integration && "harmony" %in% integration_method) {
    PCs_harmony <- FindMinimumPCs(seurat, reduction_type = "harmony")
  }
  
  # Save PCs_pca and PCs_harmony to the global environment
  assign("PCs_pca", PCs_pca, envir = .GlobalEnv)
  assign("PCs_harmony", PCs_harmony, envir = .GlobalEnv)
  
  # Run other integration methods if specified and run_integration is TRUE
  if (run_integration && any(c("CCAIntegration", "RPCAIntegration", "JointPCAIntegration") %in% integration_method)) {
    set.seed(seed)
    
    if ("CCAIntegration" %in% integration_method) {
      seurat <- IntegrateLayers(
        object = seurat, method = "CCAIntegration",
        orig.reduction = "pca", new.reduction = "cca",
        verbose = TRUE)
    }
    
    if ("RPCAIntegration" %in% integration_method) {
      seurat <- IntegrateLayers(
        object = seurat, method = "RPCAIntegration",
        orig.reduction = "pca", new.reduction = "rpca",
        verbose = TRUE)
    }
    
    if ("JointPCAIntegration" %in% integration_method) {
      seurat <- IntegrateLayers(
        object = seurat, method = "JointPCAIntegration",
        new.reduction = "jpca",
        verbose = TRUE)
    }
  }
  
  # Run tSNE and UMAP using PCA
  if (run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    seurat <- seurat %>%
      RunTSNE(reduction = "pca", dims = 1:PCs_pca, reduction.name = "tsne.unintegrated") %>%
      RunUMAP(reduction = "pca", dims = 1:PCs_pca, reduction.name = "umap.unintegrated")
    print(DimPlot(seurat, reduction = "tsne.unintegrated", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.unintegrated", raster = FALSE))
  }
  
  # Run tSNE and UMAP using Harmony if integration is performed
  if (run_integration && "harmony" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    seurat <- seurat %>%
      RunTSNE(reduction = "harmony", dims = 1:PCs_harmony, reduction.name = "tsne.harmony") %>%
      RunUMAP(reduction = "harmony", dims = 1:PCs_harmony, reduction.name = "umap.harmony")
    print(DimPlot(seurat, reduction = "tsne.harmony", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.harmony", raster = FALSE))
  }
  
  # Run tSNE and UMAP using CCA integration if integration is performed
  if (run_integration && "CCAIntegration" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    seurat <- seurat %>%
      RunTSNE(reduction = "cca", dims = 1:PCs_cca, reduction.name = "tsne.cca") %>%
      RunUMAP(reduction = "cca", dims = 1:PCs_cca, reduction.name = "umap.cca")
    print(DimPlot(seurat, reduction = "tsne.cca", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.cca", raster = FALSE))
  }
  
  # Run tSNE and UMAP using RPCA integration if integration is performed
  if (run_integration && "RPCAIntegration" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    seurat <- seurat %>%
      RunTSNE(reduction = "rpca", dims = 1:PCs_rpca, reduction.name = "tsne.rpca") %>%
      RunUMAP(reduction = "rpca", dims = 1:PCs_rpca, reduction.name = "umap.rpca")
    print(DimPlot(seurat, reduction = "tsne.rpca", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.rpca", raster = FALSE))
  }
  
  # Run tSNE and UMAP using JointPCA integration if integration is performed
  if (run_integration && "jpca" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    seurat <- seurat %>%
      RunTSNE(reduction = "jpca", dims = 1:PCs_jpca, reduction.name = "tsne.jpca") %>%
      RunUMAP(reduction = "jpca", dims = 1:PCs_jpca, reduction.name = "umap.jpca")
    print(DimPlot(seurat, reduction = "tsne.jpca", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.jpca", raster = FALSE))
  }
  
  # Prioritize the "run_clustering" check, even if more than one integration method is defined
  if (!run_clustering) {
    message("Clustering skipped because run_clustering is FALSE.")
    
  } else if (run_integration && length(integration_method) > 1) {
    message("Clustering skipped because more than one integration method was defined.")
    
  } else if (run_clustering) {
    
    # If run_integration is FALSE, always run clustering on PCA
    reduction_for_clustering <- ifelse(run_integration, "harmony", "pca")
    
    # Avoid looking for 'cca', 'rpca', or 'jpca' if no integration is performed
    if (run_integration) {
      if ("CCAIntegration" %in% integration_method) {
        reduction_for_clustering <- "cca"
      } else if ("RPCAIntegration" %in% integration_method) {
        reduction_for_clustering <- "rpca"
      } else if ("JointPCAIntegration" %in% integration_method) {
        reduction_for_clustering <- "jpca"
      }
    }
    
    # Run neighbors and clustering based on selected reduction method (default to PCA if no integration method is defined)
    seurat <- FindNeighbors(seurat, reduction = reduction_for_clustering, dims = 1:get(paste0("PCs_", reduction_for_clustering)))
    seurat <- FindClusters(seurat, resolution = resolutions)  
  }
  
  # Return object
  return(seurat)
}