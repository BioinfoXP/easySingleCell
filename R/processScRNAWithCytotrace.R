#' Run CytoTRACE Analysis
#'
#' This function runs CytoTRACE analysis on a Seurat object and saves the results and plots.
#'
#' @param scRNA A Seurat object.
#' @param celltype A character string specifying the column name in `scRNA@meta.data` that contains cell type information. Default is 'celltype'.
#' @param output_data_dir Directory to save the output data. Default is "./output_data/".
#' @param output_figure_dir Directory to save the output figures. Default is "./output_data/".
#' @param ncores Number of cores to use for CytoTRACE analysis. Default is 5.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' runCytoTRACEAnalysis(scRNA = your_seurat_object, celltype = 'celltype', output_data_dir = "./output_data/", output_figure_dir = "./output_data/", ncores = 5)
#' }
runCytoTRACEAnalysis <- function(scRNA, celltype = 'celltype',
                                 output_data_dir = "./output_data/", output_figure_dir = "./output_data/", ncores = 5) {
  # Check if necessary packages are installed
  required_packages <- c("Seurat", "CytoTRACE")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Create output directories if they do not exist
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }
  if (!dir.exists(output_figure_dir)) {
    dir.create(output_figure_dir, recursive = TRUE)
  }

  # Load required libraries
  library(Seurat)
  library(CytoTRACE)

  # Validate inputs
  if (!inherits(scRNA, "Seurat")) {
    stop("The 'scRNA' parameter must be a Seurat object.")
  }
  if (!is.character(celltype) || length(celltype) != 1) {
    stop("The 'celltype' parameter must be a single character string.")
  }
  if (!celltype %in% colnames(scRNA@meta.data)) {
    stop(paste("The specified celltype column '", celltype, "' does not exist in scRNA@meta.data.", sep = ""))
  }
  if (!is.character(output_data_dir) || length(output_data_dir) != 1) {
    stop("The 'output_data_dir' parameter must be a single character string.")
  }
  if (!is.character(output_figure_dir) || length(output_figure_dir) != 1) {
    stop("The 'output_figure_dir' parameter must be a single character string.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0) {
    stop("The 'ncores' parameter must be a single positive number.")
  }

  # Run CytoTRACE analysis
  phe <- scRNA@meta.data[, celltype]
  phe <- as.character(phe)
  names(phe) <- rownames(scRNA@meta.data)

  mat_3k <- as.matrix(scRNA@assays$RNA@counts)

  results <- CytoTRACE(mat = mat_3k, ncores = ncores)
  save(results, file = file.path(output_data_dir, "Cytotrace_results.Rdata"))

  plotCytoTRACE(results, outputDir = output_figure_dir, emb = scRNA@reductions$umap@cell.embeddings)
}
