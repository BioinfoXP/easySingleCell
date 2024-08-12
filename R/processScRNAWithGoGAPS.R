#' @title Run CoGAPS Analysis on Single-cell Data
#' @description This function runs CoGAPS analysis on single-cell data for a specified pattern number and saves the results.
#' @param scRNA A Seurat object containing single-cell data.
#' @param nPatterns A numeric value specifying the number of patterns to be used in CoGAPS analysis.
#' @param out_dir A character string specifying the output directory where results will be saved.
#' @param nIterations A numeric value specifying the number of iterations for CoGAPS (default is 1000).
#' @param seed A numeric value specifying the random seed for reproducibility (default is 42).
#' @param nThreads A numeric value specifying the number of threads to be used (default is 4).
#' @return The CoGAPS result object.
#' @export
#' @import Seurat
#' @import CoGAPS
#' @import parallel
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(CoGAPS)
#' library(parallel)
#'
#' # Example Seurat object
#' seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' seurat_obj <- NormalizeData(seurat_obj)
#' seurat_obj <- FindVariableFeatures(seurat_obj)
#' seurat_obj <- ScaleData(seurat_obj)
#' seurat_obj <- RunPCA(seurat_obj)
#' seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#'
#' # Run CoGAPS analysis
#' runCogapsAnalysis(
#'   scRNA = seurat_obj,
#'   nPatterns = 4,
#'   out_dir = './output_data/Figure7/cogaps',
#'   nIterations = 1000,
#'   seed = 42,
#'   nThreads = 4
#' )
#' }

runCogapsAnalysis <- function(scRNA, nPatterns, out_dir, nIterations = 1000, seed = 42, nThreads = 4) {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(CoGAPS)
  library(parallel)
  library(magrittr)

  # Extract count matrix and normalize
  counts <- GetAssayData(object = scRNA, slot = "counts") %>% as.matrix()
  norm_counts <- log1p(counts)

  # Create output directory if it does not exist
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Set CoGAPS parameters
  params <- CogapsParams(
    nIterations = nIterations,
    seed = seed,
    nPatterns = nPatterns,
    sparseOptimization = TRUE,
    distributed = "genome-wide"
  )

  # Run CoGAPS analysis
  cogaps_result <- CoGAPS(norm_counts, params, nThreads = nThreads)

  # Save the result
  saveRDS(cogaps_result, file = paste0(out_dir, "/cogaps_result_k_", nPatterns, ".rds"))

  # Return the result object
  return(cogaps_result)
}

# Example usage
# library(Seurat)
# library(CoGAPS)
# library(parallel)
#
# # Example Seurat object
# seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#
# # Run CoGAPS analysis
# cogaps_result <- runCogapsAnalysis(
#   scRNA = seurat_obj,
#   nPatterns = 4,
#   out_dir = './output_data/cogaps',
#   nIterations = 1000,
#   seed = 42,
#   nThreads = 4
# )
