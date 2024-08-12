#' Run CellChat Analysis
#'
#' This function runs CellChat analysis on a Seurat object, normalizes and clusters the data, and saves the results.
#'
#' @param sce A Seurat object.
#' @param celltype A character string specifying the column name in `sce@meta.data` that contains cell type information. Default is 'cell_type'.
#' @param groupby A character string specifying the column name in `sce@meta.data` to split the object by. If NULL, the object will not be split. Default is NULL.
#' @param species A character string specifying the species. Must be either 'human' or 'mouse'.
#' @param output_data_dir Directory to save the output data. Default is "./output_data/".
#' @param ncores Number of cores to use for parallel processing. Default is 1.
#' @param workers Number of workers to use for future processing. Default is 5.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' runCellChatAnalysis(sce = your_seurat_object, celltype = 'cell_type', groupby = 'group', species = 'human', output_data_dir = "./output_data/", ncores = 4, workers = 5)
#' }
runCellChatAnalysis <- function(sce, celltype = 'cell_type', groupby = NULL, species,
                                output_data_dir = "./output_data/", ncores = 1, workers = 5) {
  # Check if necessary packages are installed
  required_packages <- c("Seurat", "CellChat", "future", "dplyr", "qs")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Create output directory if it does not exist
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }

  # Load required libraries
  library(Seurat)
  library(CellChat)
  library(future)
  library(dplyr)
  library(qs)

  # Validate inputs
  if (!inherits(sce, "Seurat")) {
    stop("The 'sce' parameter must be a Seurat object.")
  }
  if (!is.character(celltype) || length(celltype) != 1) {
    stop("The 'celltype' parameter must be a single character string.")
  }
  if (!celltype %in% colnames(sce@meta.data)) {
    stop(paste("The specified celltype column '", celltype, "' does not exist in sce@meta.data.", sep = ""))
  }
  if (!is.null(groupby) && (!is.character(groupby) || length(groupby) != 1)) {
    stop("The 'groupby' parameter must be a single character string or NULL.")
  }
  if (!is.null(groupby) && !groupby %in% colnames(sce@meta.data)) {
    stop(paste("The specified groupby column '", groupby, "' does not exist in sce@meta.data.", sep = ""))
  }
  if (!is.character(species) || length(species) != 1 || !(species %in% c("human", "mouse"))) {
    stop("The 'species' parameter must be a single character string and must be either 'human' or 'mouse'.")
  }
  if (!is.character(output_data_dir) || length(output_data_dir) != 1) {
    stop("The 'output_data_dir' parameter must be a single character string.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0) {
    stop("The 'ncores' parameter must be a single positive number.")
  }
  if (!is.numeric(workers) || length(workers) != 1 || workers <= 0) {
    stop("The 'workers' parameter must be a single positive number.")
  }

  # Select appropriate CellChat database and PPI data
  if (species == "human") {
    CellChatDB <- CellChatDB.human
    PPI <- PPI.human
  } else if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
    PPI <- PPI.mouse
  }

  # Split object by group if specified
  if (!is.null(groupby)) {
    sce_list <- SplitObject(sce, split.by = groupby)
    names_list <- names(sce_list)
  } else {
    sce_list <- list(sce)
    names_list <- "all"
  }

  # Normalize and cluster
  sce_list <- mclapply(1:length(sce_list), function(i) {
    sce_list[[i]] %>%
      NormalizeData() %>%
      FindClusters()
  }, mc.cores = ncores)

  # Run CellChat analysis
  cellchat_res <- mclapply(1:length(sce_list), function(i) {
    data.input <- sce_list[[i]]@assays$RNA@data
    identity <- sce_list[[i]]@meta.data
    Idents(sce_list[[i]]) <- celltype

    cellchat <- createCellChat(data.input, meta = identity, group.by = celltype)
    cellchat <- addMeta(cellchat, meta = identity)
    cellchat <- setIdent(cellchat, ident.use = celltype)
    groupSize <- as.numeric(table(cellchat@idents))

    # Use the selected CellChat database
    cellchat@DB <- CellChatDB

    # Subset data and identify overexpressed genes and interactions
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = workers)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    return(cellchat)
  }, mc.cores = ncores)

  # Save results
  names(cellchat_res) <- names_list
  save(sce_list, cellchat_res, file = file.path(output_data_dir, "cellchat_res.Rdata"))
}
