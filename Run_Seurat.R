#' Run_Seurat
#' 
#' The `Run_Seurat()` function is a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline designed to preprocess, integrate, and visualize data using the Seurat workflow.
#' It includes steps for normalization, feature selection, dimensionality reduction, and data integration using methods such as Harmony, Canonical Correlation Analysis (CCA), Reciprocal PCA (RPCA), and Joint PCA (JPCA).
#' It also supports clustering and dimensionality reduction visualizations like t-SNE and UMAP.
#'
#' @param seurat A Seurat object containing the single-cell RNA sequencing data.
#' @param batch_column The metadata column name in the Seurat object representing the batch or condition variable (default is "orig.ident").
#' @param nfeatures The number of variable features to select for downstream analysis (default is 2000).
#' @param variables_to_regress A character vector of metadata variables to regress out during scaling (e.g., c("percent.mt", "nCount_RNA")). If NULL, no regression is performed (default is NULL).
#' @param seed A random seed for reproducibility (default is 123).
#' @param PCs_pca The number of principal components (PCs) to use for PCA. If not specified, it will be calculated using the `FindMinimumPCs()` function.
#' @param PCs_harmony The number of principal components (PCs) to use for Harmony integration. If not specified, it will be calculated using the `FindMinimumPCs()` function.
#' @param PCs_cca The number of principal components (PCs) to use for CCA integration (default is 30).
#' @param PCs_rpca The number of principal components (PCs) to use for RPCA integration (default is 30).
#' @param PCs_jpca The number of principal components (PCs) to use for JPCA integration (default is 30).
#' @param run_integration Specifies whether integration methods such as Harmony, CCA, RPCA, and JPCA should be applied. When set to TRUE, integration will be performed using the specified integration methods (default is FALSE).
#' @param integration_method The integration method(s) to use. Options include "harmony", "CCAIntegration", "RPCAIntegration", and "JointPCAIntegration". Multiple methods can be specified (default is "harmony").
#' @param run_tSNE_UMAP Specifies whether t-SNE and UMAP dimensionality reductions should be performed (default is FALSE).
#' @param check_duplicates Checks for duplicated cell names in `RunTSNE()` (default is FALSE).
#' @param run_clustering Controls whether clustering is performed. If set to TRUE, clustering will be performed on the integrated data unless more than one integration method is specified. If integration is not performed, clustering will default to using PCA.
#' @param resolutions A numeric vector of resolutions to be used in clustering.
#' @param normalize.args A named list of additional arguments to pass to NormalizeData() (default is NULL).
#' @param variable.features.args A named list of additional arguments to pass to FindVariableFeatures() (default is NULL).
#' @param scale.args A named list of additional arguments to pass to ScaleData() (default is NULL).
#' @param pca.args A named list of additional arguments to pass to RunPCA() (default is NULL).
#' @param harmony.args A named list of additional arguments to pass to RunHarmony() (default is NULL).
#' @param cca.args A named list of additional arguments to pass to IntegrateLayers() for CCA (default is NULL).
#' @param rpca.args A named list of additional arguments to pass to IntegrateLayers() for RPCA (default is NULL).
#' @param jpca.args A named list of additional arguments to pass to IntegrateLayers() for JPCA (default is NULL).
#' @param tsne.args A named list of additional arguments to pass to RunTSNE() (default is NULL).
#' @param umap.args A named list of additional arguments to pass to RunUMAP() (default is NULL).
#' @param neighbors.args A named list of additional arguments to pass to FindNeighbors() (default is NULL).
#' @param clusters.args A named list of additional arguments to pass to FindClusters() (default is NULL).
#' 
#' @return A Seurat object with the results of preprocessing, integration, and clustering.
#'
#' @note If the user chooses not to define the number of PCs for PCA or Harmony (i.e., setting `PCs_pca` or `PCs_harmony` to NULL), they must source the custom `FindMinimumPCs()` function
#' before running `Run_Seurat()`, as this will automatically determine the minimum number of PCs.
#' 
#' @references
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
#' @lastUpdated 2025-11-7
#'
#' @examples
#' seurat_obj <- Run_Seurat(seurat = seurat_obj, run_integration = TRUE, integration_method = c("harmony", "CCAIntegration"), run_clustering = FALSE)
#' seurat_obj <- Run_Seurat(seurat = seurat_obj, PCs_pca = 30, run_integration = FALSE, run_clustering = TRUE, resolutions = c(0.4, 0.6, 0.8))
#' seurat_obj <- Run_Seurat(seurat = seurat_obj, run_integration = TRUE, integration_method = "RPCAIntegration", rpca.args = list(k.weight = 50))
#' seurat_obj <- Run_Seurat(seurat = seurat_obj, run_tSNE_UMAP = TRUE, umap.args = list(min.dist = 0.1, n.neighbors = 50))
#' seurat_obj <- Run_Seurat(seurat = seurat_obj, run_integration = TRUE, integration_method = c("RPCAIntegration", "CCAIntegration"), rpca.args = list(k.weight = 25), cca.args = list(k.weight = 150))

Run_Seurat <- function(seurat,
                       batch_column = "orig.ident",
                       nfeatures = 2000,
                       variables_to_regress = NULL,
                       seed = 123,
                       PCs_pca = NULL, 
                       PCs_harmony = NULL,
                       PCs_cca = 30,
                       PCs_rpca = 30,
                       PCs_jpca = 30,
                       run_integration = FALSE,
                       integration_method = "harmony",
                       run_tSNE_UMAP = FALSE,
                       check_duplicates = FALSE,
                       run_clustering = FALSE,
                       resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0),
                       normalize.args = NULL,
                       variable.features.args = NULL,
                       scale.args = NULL,
                       pca.args = NULL,
                       harmony.args = NULL,
                       cca.args = NULL,
                       rpca.args = NULL,
                       jpca.args = NULL,
                       tsne.args = NULL,
                       umap.args = NULL,
                       neighbors.args = NULL,
                       clusters.args = NULL) {
  
  # Normalization
  if (is.null(normalize.args)) {
    seurat <- NormalizeData(seurat)
  } else {
    norm_args <- modifyList(list(object = seurat), normalize.args)
    seurat <- do.call(NormalizeData, norm_args)
  }
  
  # Find Variable Features
  if (is.null(variable.features.args)) {
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = nfeatures)
  } else {
    fvf_args <- modifyList(list(object = seurat, selection.method = "vst", nfeatures = nfeatures), variable.features.args)
    seurat <- do.call(FindVariableFeatures, fvf_args)
  }
  
  # Scale Data
  if (is.null(scale.args)) {
    seurat <- ScaleData(seurat, vars.to.regress = variables_to_regress)
  } else {
    scale_args <- modifyList(list(object = seurat, vars.to.regress = variables_to_regress), scale.args)
    seurat <- do.call(ScaleData, scale_args)
  }
  
  # Run PCA
  if (is.null(pca.args)) {
    seurat <- RunPCA(seurat)
  } else {
    pca_args <- modifyList(list(object = seurat), pca.args)
    seurat <- do.call(RunPCA, pca_args)
  }
  print(ElbowPlot(seurat))
  
  # Use provided PCs_pca or calculate using FindMinimumPCs function (see documentation)
  if (is.null(PCs_pca)) {
    PCs_pca <- FindMinimumPCs(seurat, reduction_type = "pca")
  }
  
  # Run Harmony integration if specified and run_integration is TRUE
  if (run_integration && "harmony" %in% integration_method) {
    set.seed(seed)
    if (is.null(harmony.args)) {
      seurat <- RunHarmony(seurat, batch_column, plot_convergence = TRUE, reduction = "pca", reduction.save = "harmony")
    } else {
      harm_args <- modifyList(list(object = seurat, group.by.vars = batch_column, plot_convergence = TRUE, 
                                   reduction = "pca", reduction.save = "harmony"), harmony.args)
      seurat <- do.call(RunHarmony, harm_args)
    }
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
      if (is.null(cca.args)) {
        seurat <- IntegrateLayers(
          object = seurat, method = "CCAIntegration",
          orig.reduction = "pca", new.reduction = "cca",
          verbose = TRUE)
      } else {
        # Build the call string and evaluate it to avoid do.call issues
        args_str <- paste(names(cca.args), "=", sapply(cca.args, deparse), collapse = ", ")
        eval_str <- sprintf("IntegrateLayers(object = seurat, method = 'CCAIntegration', orig.reduction = 'pca', new.reduction = 'cca', verbose = TRUE, %s)", args_str)
        seurat <- eval(parse(text = eval_str))
      }
    }
    
    if ("RPCAIntegration" %in% integration_method) {
      if (is.null(rpca.args)) {
        seurat <- IntegrateLayers(
          object = seurat, method = "RPCAIntegration",
          orig.reduction = "pca", new.reduction = "rpca",
          verbose = TRUE)
      } else {
        # Build the call string and evaluate it to avoid do.call issues
        args_str <- paste(names(rpca.args), "=", sapply(rpca.args, deparse), collapse = ", ")
        eval_str <- sprintf("IntegrateLayers(object = seurat, method = 'RPCAIntegration', orig.reduction = 'pca', new.reduction = 'rpca', verbose = TRUE, %s)", args_str)
        seurat <- eval(parse(text = eval_str))
      }
    }
    
    if ("JointPCAIntegration" %in% integration_method) {
      if (is.null(jpca.args)) {
        seurat <- IntegrateLayers(
          object = seurat, method = "JointPCAIntegration",
          new.reduction = "jpca",
          verbose = TRUE)
      } else {
        # Build the call string and evaluate it to avoid do.call issues
        args_str <- paste(names(jpca.args), "=", sapply(jpca.args, deparse), collapse = ", ")
        eval_str <- sprintf("IntegrateLayers(object = seurat, method = 'JointPCAIntegration', new.reduction = 'jpca', verbose = TRUE, %s)", args_str)
        seurat <- eval(parse(text = eval_str))
      }
    }
  }
  
  # Run tSNE and UMAP using PCA
  if (run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    
    # tSNE on PCA
    if (is.null(tsne.args)) {
      seurat <- RunTSNE(seurat, reduction = "pca", dims = 1:PCs_pca, 
                        reduction.name = "tsne.unintegrated", check_duplicates = check_duplicates)
    } else {
      tsne_pca_args <- modifyList(list(object = seurat, reduction = "pca", dims = 1:PCs_pca, 
                                       reduction.name = "tsne.unintegrated", check_duplicates = check_duplicates), tsne.args)
      seurat <- do.call(RunTSNE, tsne_pca_args)
    }
    
    # UMAP on PCA
    if (is.null(umap.args)) {
      seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCs_pca, 
                        reduction.name = "umap.unintegrated")
    } else {
      umap_pca_args <- modifyList(list(object = seurat, reduction = "pca", dims = 1:PCs_pca, 
                                       reduction.name = "umap.unintegrated"), umap.args)
      seurat <- do.call(RunUMAP, umap_pca_args)
    }
    
    print(DimPlot(seurat, reduction = "tsne.unintegrated", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.unintegrated", raster = FALSE))
  }
  
  # Run tSNE and UMAP using Harmony if integration is performed
  if (run_integration && "harmony" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    
    # tSNE on Harmony
    if (is.null(tsne.args)) {
      seurat <- RunTSNE(seurat, reduction = "harmony", dims = 1:PCs_harmony, 
                        reduction.name = "tsne.harmony", check_duplicates = check_duplicates)
    } else {
      tsne_harmony_args <- modifyList(list(object = seurat, reduction = "harmony", dims = 1:PCs_harmony, 
                                           reduction.name = "tsne.harmony", check_duplicates = check_duplicates), tsne.args)
      seurat <- do.call(RunTSNE, tsne_harmony_args)
    }
    
    # UMAP on Harmony
    if (is.null(umap.args)) {
      seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:PCs_harmony, 
                        reduction.name = "umap.harmony")
    } else {
      umap_harmony_args <- modifyList(list(object = seurat, reduction = "harmony", dims = 1:PCs_harmony, 
                                           reduction.name = "umap.harmony"), umap.args)
      seurat <- do.call(RunUMAP, umap_harmony_args)
    }
    
    print(DimPlot(seurat, reduction = "tsne.harmony", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.harmony", raster = FALSE))
  }
  
  # Run tSNE and UMAP using CCA integration if integration is performed
  if (run_integration && "CCAIntegration" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    
    # tSNE on CCA
    if (is.null(tsne.args)) {
      seurat <- RunTSNE(seurat, reduction = "cca", dims = 1:PCs_cca, 
                        reduction.name = "tsne.cca", check_duplicates = check_duplicates)
    } else {
      tsne_cca_args <- modifyList(list(object = seurat, reduction = "cca", dims = 1:PCs_cca, 
                                       reduction.name = "tsne.cca", check_duplicates = check_duplicates), tsne.args)
      seurat <- do.call(RunTSNE, tsne_cca_args)
    }
    
    # UMAP on CCA
    if (is.null(umap.args)) {
      seurat <- RunUMAP(seurat, reduction = "cca", dims = 1:PCs_cca, 
                        reduction.name = "umap.cca")
    } else {
      umap_cca_args <- modifyList(list(object = seurat, reduction = "cca", dims = 1:PCs_cca, 
                                       reduction.name = "umap.cca"), umap.args)
      seurat <- do.call(RunUMAP, umap_cca_args)
    }
    
    print(DimPlot(seurat, reduction = "tsne.cca", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.cca", raster = FALSE))
  }
  
  # Run tSNE and UMAP using RPCA integration if integration is performed
  if (run_integration && "RPCAIntegration" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    
    # tSNE on RPCA
    if (is.null(tsne.args)) {
      seurat <- RunTSNE(seurat, reduction = "rpca", dims = 1:PCs_rpca, 
                        reduction.name = "tsne.rpca", check_duplicates = check_duplicates)
    } else {
      tsne_rpca_args <- modifyList(list(object = seurat, reduction = "rpca", dims = 1:PCs_rpca, 
                                        reduction.name = "tsne.rpca", check_duplicates = check_duplicates), tsne.args)
      seurat <- do.call(RunTSNE, tsne_rpca_args)
    }
    
    # UMAP on RPCA
    if (is.null(umap.args)) {
      seurat <- RunUMAP(seurat, reduction = "rpca", dims = 1:PCs_rpca, 
                        reduction.name = "umap.rpca")
    } else {
      umap_rpca_args <- modifyList(list(object = seurat, reduction = "rpca", dims = 1:PCs_rpca, 
                                        reduction.name = "umap.rpca"), umap.args)
      seurat <- do.call(RunUMAP, umap_rpca_args)
    }
    
    print(DimPlot(seurat, reduction = "tsne.rpca", raster = FALSE))
    print(DimPlot(seurat, reduction = "umap.rpca", raster = FALSE))
  }
  
  # Run tSNE and UMAP using JointPCA integration if integration is performed
  if (run_integration && "JointPCAIntegration" %in% integration_method && run_tSNE_UMAP) {
    Idents(seurat) <- batch_column
    
    # tSNE on JPCA
    if (is.null(tsne.args)) {
      seurat <- RunTSNE(seurat, reduction = "jpca", dims = 1:PCs_jpca, 
                        reduction.name = "tsne.jpca", check_duplicates = check_duplicates)
    } else {
      tsne_jpca_args <- modifyList(list(object = seurat, reduction = "jpca", dims = 1:PCs_jpca, 
                                        reduction.name = "tsne.jpca", check_duplicates = check_duplicates), tsne.args)
      seurat <- do.call(RunTSNE, tsne_jpca_args)
    }
    
    # UMAP on JPCA
    if (is.null(umap.args)) {
      seurat <- RunUMAP(seurat, reduction = "jpca", dims = 1:PCs_jpca, 
                        reduction.name = "umap.jpca")
    } else {
      umap_jpca_args <- modifyList(list(object = seurat, reduction = "jpca", dims = 1:PCs_jpca, 
                                        reduction.name = "umap.jpca"), umap.args)
      seurat <- do.call(RunUMAP, umap_jpca_args)
    }
    
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
    
    # Run neighbors
    if (is.null(neighbors.args)) {
      seurat <- FindNeighbors(seurat, reduction = reduction_for_clustering, 
                              dims = 1:get(paste0("PCs_", reduction_for_clustering)))
    } else {
      fn_args <- modifyList(list(object = seurat, reduction = reduction_for_clustering, 
                                 dims = 1:get(paste0("PCs_", reduction_for_clustering))), neighbors.args)
      seurat <- do.call(FindNeighbors, fn_args)
    }
    
    # Run clustering
    if (is.null(clusters.args)) {
      seurat <- FindClusters(seurat, resolution = resolutions)
    } else {
      fc_args <- modifyList(list(object = seurat, resolution = resolutions), clusters.args)
      seurat <- do.call(FindClusters, fc_args)
    }
  }
  
  # Return object
  return(seurat)
}
