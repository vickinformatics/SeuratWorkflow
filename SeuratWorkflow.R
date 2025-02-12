### Description: The SeuratWorkflow function is a custom analysis pipeline designed to preprocess, integrate, and visualize single-cell RNA sequencing (scRNA-seq) data using Seurat. It provides a flexible and comprehensive approach to analyzing scRNA-seq datasets, using normalization, feature selection, dimensionality reduction, and integration methods such as Harmony and Canonical Correlation Analysis (CCA). It also includes clustering, dimensionality reduction visualizations (t-SNE and UMAP), and allows for the handling of batch effects during integration. This function is for researchers who wish to integrate multiple datasets from different experimental batches or need a streamlined tool to help them perform common analysis steps with ease.

### Last updated: 2025-2-7

### Author: Vicki Do

SeuratWorkflow <- function(seurat,
                           batch_column = "orig.ident",
                           nfeatures = 2000,
                           seed = 123,
                           integration_method = "harmony",
                           min.pc_cca = 30, # User needs to define this because cca.integration reduction does not provide a stdev
                           run_clustering = TRUE,
                           resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2)) {
  
  # Normalization, scaling, and PCA
  seurat <- seurat %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
    ScaleData() %>%
    RunPCA()
  print(ElbowPlot(seurat))
  
  # Run Harmony integration if specified
  if ("harmony" %in% integration_method) {
    set.seed(seed)
    seurat_integrated <- seurat %>% 
      RunHarmony(batch_column, plot_convergence = TRUE, reduction = "pca", reduction.save = "harmony")
  }
  
  # Run CCAIntegration if specified
  if ("cca.integration" %in% integration_method) {
    set.seed(seed)
    seurat_integrated <- IntegrateLayers(seurat_integrated, method = "CCAIntegration", orig.reduction = "pca", new.reduction = "cca.integration", verbose = FALSE)
  }
  
  # Find minimum PCs (see FindMinimumPCs documentation for details) and save to global environment
  min.pc_pca <- FindMinimumPCs(seurat, reduction_type = "pca")
  min.pc_harmony <- FindMinimumPCs(seurat_integrated, reduction_type = "harmony")
  
  assign("min.pc_pca", min.pc_pca, envir = .GlobalEnv)
  assign("min.pc_harmony", min.pc_harmony, envir = .GlobalEnv)
  assign("min.pc_cca", min.pc_cca, envir = .GlobalEnv)
  
  # Run tSNE and UMAP using PCA
  Idents(seurat_integrated) <- batch_column
  seurat_integrated <- seurat_integrated %>%
    RunTSNE(reduction = "pca", dims = 1:min.pc_pca, reduction.name = "tsne.unintegrated") %>%
    RunUMAP(reduction = "pca", dims = 1:min.pc_pca, reduction.name = "umap.unintegrated")
  print(DimPlot(seurat_integrated, reduction = "tsne.unintegrated", raster = FALSE))
  print(DimPlot(seurat_integrated, reduction = "umap.unintegrated", raster = FALSE))
  
  # Run tSNE and UMAP using Harmony
  if ("harmony" %in% integration_method) {
    Idents(seurat_integrated) <- batch_column
    seurat_integrated <- seurat_integrated %>%
      RunTSNE(reduction = "harmony", dims = 1:min.pc_harmony, reduction.name = "tsne.harmony") %>%
      RunUMAP(reduction = "harmony", dims = 1:min.pc_harmony, reduction.name = "umap.harmony")
    print(DimPlot(seurat_integrated, reduction = "tsne.harmony", raster = FALSE))
    print(DimPlot(seurat_integrated, reduction = "umap.harmony", raster = FALSE))
  }
  
  # Run tSNE and UMAP using CCAIntegration
  if ("cca.integration" %in% integration_method) {
    Idents(seurat_integrated) <- batch_column
    seurat_integrated <- seurat_integrated %>%
      RunTSNE(reduction = "cca.integration", dims = 1:min.pc_cca, reduction.name = "tsne.cca") %>%
      RunUMAP(reduction = "cca.integration", dims = 1:min.pc_cca, reduction.name = "umap.cca")
    print(DimPlot(seurat_integrated, reduction = "tsne.cca", raster = FALSE))
    print(DimPlot(seurat_integrated, reduction = "umap.cca", raster = FALSE))
  }
  
  # Run neighbors and clustering
  if (run_clustering && !(all(c("harmony", "cca.integration") %in% integration_method))) { # Skip clustering if both integration methods are selected
    seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:min.pc_harmony)
    seurat_integrated <- FindClusters(seurat_integrated, resolution = resolutions)
  }
  
  # Return integrated object
  return(seurat_integrated)
}