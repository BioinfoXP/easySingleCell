# =========== 1. PrepareCpdb =============

#' @title Preprocess and Save Single-cell Data for CellPhoneDB Analysis
#' @description This function preprocesses single-cell data, normalizes it, and saves the results for CellPhoneDB analysis.
#' @param scRNA A Seurat object containing the single-cell RNA data.
#' @param output_dir A character string specifying the output directory where results will be saved. Default is './output_data/'.
#' @param celltype_column A character string specifying the column name in the metadata that contains cell type information. Default is 'celltype'.
#' @param normalization_threshold A numeric value specifying the threshold for filtering genes based on mean expression. Default is 0.001.
#' @param prefix A character string specifying the prefix to be added to the output file names. Default is "".
#' @return None. The function saves the normalized counts and metadata files in the specified output directory.
#' @export
#' @import Seurat
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' PrepareCpdb(
#'   scRNA = sce,
#'   output_dir = './output_data/',
#'   celltype_column = 'celltype',
#'   normalization_threshold = 0.001,
#'   prefix = ''
#' )
#' }

PrepareCpdb <- function(scRNA, output_dir = './output_data/', celltype_column = 'celltype', normalization_threshold = 0.001, prefix = "") {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(magrittr)

  # Create the output directory for CellPhoneDB if it doesn't exist
  cpdb_dir <- file.path(output_dir, 'cpdb')
  dir.create(cpdb_dir, showWarnings = FALSE, recursive = TRUE)

  # Normalize the counts and filter genes with low expression
  Normalized_counts <- GetAssayData(scRNA, slot = 'data') %>% as.data.frame()
  Normalized_counts$mean_exp <- rowMeans(Normalized_counts)
  Normalized_counts <- Normalized_counts[Normalized_counts$mean_exp > normalization_threshold, ]
  Normalized_counts <- Normalized_counts[, -ncol(Normalized_counts)]
  Normalized_counts <- cbind(Gene = rownames(Normalized_counts), Normalized_counts)

  # Create metadata file
  metadata <- data.frame(Cell = rownames(scRNA@meta.data), celltype = scRNA@meta.data[, celltype_column])
  metadata$celltype <- gsub(' ', '_', metadata$celltype)

  # Save the normalized counts and metadata files
  write.table(Normalized_counts, file.path(cpdb_dir, paste0(prefix, 'Normalized_counts.txt')), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(metadata, file.path(cpdb_dir, paste0(prefix, 'cellphonedb_meta.txt')), row.names = FALSE, sep = '\t', quote = FALSE)
}

# Example usage
# PrepareCpdb(
#   scRNA = sce,
#'   output_dir = './output_data/',
#'   celltype_column = 'celltype',
#'   normalization_threshold = 0.001,
#'   prefix = ''
#' )


# ============ 2. PrepareCell2loc ==========

#' @title Prepare Data for Cell2loc
#' @description This function prepares data for Cell2loc analysis by exporting counts matrix, UMAP embeddings, and metadata from a Seurat object.
#' @param input_file A character string specifying the path to the Seurat object file.
#' @param output_dir A character string specifying the directory to save the output files. Default is "output_data".
#' @param prefix A character string specifying the prefix to be added to the output file names. Default is "".
#' @return NULL
#' @export
#' @import Seurat
#' @import qs
#' @examples
#' \dontrun{
#' PrepareCell2loc(input_file = './output_data/Figure11/sce.tumor.所有细胞注释.qs', output_dir = 'output_data', prefix = 'Figure12_')
#' }

PrepareCell2loc <- function(input_file, output_dir = "output_data", prefix = "") {
  # Load required libraries
  library(Seurat)
  library(qs)

  # Ensure the output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load the Seurat object
  scRNA <- qread(input_file)

  # Export counts matrix
  counts_matrix <- t(as.matrix(scRNA@assays$RNA@counts))
  write.table(counts_matrix, file = file.path(output_dir, paste0(prefix, "scRNA_counts.csv")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Export UMAP embeddings
  embedding <- Embeddings(scRNA, "umap")
  write.table(embedding, file = file.path(output_dir, paste0(prefix, "scRNA_embedding.csv")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Export metadata
  metadata <- scRNA@meta.data
  write.table(metadata, file = file.path(output_dir, paste0(prefix, "metadata.csv")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  return(NULL)
}

# Example usage
# PrepareCell2loc(input_file = './output_data/Figure11/sce.tumor.所有细胞注释.qs', output_dir = 'output_data', prefix = 'Figure12_')



# =========== 3. PreparePyscenic =============

#' @title Prepare Data for pySCENIC Analysis
#' @description This function prepares the data for pySCENIC analysis by splitting the Seurat object by cell type, creating meta cells, and generating the necessary files for pySCENIC.
#' @param scRNA A Seurat object containing single-cell RNA data.
#' @param output_dir A character string specifying the output directory where results will be saved. Default is './output_data/'.
#' @param min_cells An integer specifying the minimum number of cells for meta cell creation. Default is 0.
#' @param reduction A character string specifying the reduction method to use. Default is 'umap'.
#' @param dims A numeric vector specifying the dimensions to use for the reduction. Default is 1:2.
#' @param k_param An integer specifying the k parameter for meta cell creation. Default is 10.
#' @param cores An integer specifying the number of cores to use. Default is 10.
#' @return None. The function saves the meta cell matrix and other necessary files for pySCENIC analysis in the specified output directory.
#' @export
#' @import Seurat
#' @importFrom magrittr %>%
#' @import data.table
#' @import arrow
#' @import SCopeLoomR
#' @examples
#' \dontrun{
#' PreparePyscenic(
#'   scRNA = scRNA,
#'   output_dir = './output_data/',
#'   min_cells = 0,
#'   reduction = 'umap',
#'   dims = 1:2,
#'   k_param = 10,
#'   cores = 10
#' )
#' }

PreparePyscenic <- function(scRNA, output_dir = './output_data/', min_cells = 0, reduction = 'umap', dims = 1:2, k_param = 10, cores = 10) {
  # Load necessary libraries
  library(Seurat)
  library(magrittr)
  library(data.table)
  library(arrow)
  library(SCopeLoomR)

  # Load pre-prepared data
  load(system.file("data", "pyscenic_database.Rdata", package = "easySingleCell"))

  # Load pre-definced function
  source(system.file("data", "pyscenic.R", package = "easySingleCell"))

  # Split the Seurat object by cell type
  seu.list <- SplitObject(scRNA, split.by = "celltype")
  for (i in 1:length(names(seu.list))) {
    seu.list[[i]]@project.name <- names(seu.list)[i]
  }

  # Create meta cells for each cell type
  metacells.list <- lapply(seq_along(seu.list), function(ii) {
    makeMetaCells(
      seu       = seu.list[[ii]],
      min.cells = min_cells,
      reduction = reduction,
      dims      = dims,
      k.param   = k_param,
      cores     = cores
    )
  })

  # Combine meta cell matrices and metadata
  mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
  mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)

  # Save the meta cell matrix
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(mc.mat, file.path(output_dir, "00-1.mc.mat.rds"))

  # Read the saved meta cell matrix
  mc.mat <- readRDS(file.path(output_dir, "00-1.mc.mat.rds"))

  # Filter low expression genes
  expr.in.cells <- rowSums(mc.mat > 0)
  mc.mat <- mc.mat[expr.in.cells >= 5, ]

  # Intersect genes with cisDB
  genes.use <- intersect(cisdb.genes, rownames(mc.mat))
  mc.mat <- mc.mat[genes.use, ]

  # Create a loom file for pySCENIC
  loom <- SCopeLoomR::build_loom(
    file.name         = file.path(output_dir, "00-2.mc_mat_for_step1.loom"),
    dgem              = mc.mat,
    default.embedding = NULL
  )
  loom$close()
}

# Example usage
# PreparePyscenic(
#   scRNA = scRNA,
#   output_dir = './output_data/',
#   min_cells = 0,
#   reduction = 'umap',
#   dims = 1:2,
#   k_param = 10,
#   cores = 10
# )
