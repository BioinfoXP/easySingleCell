#' @title Convert Gene Symbols to Ensembl IDs
#' @description Converts Seurat object row names from Symbols to Ensembl IDs.
#' Optimized for Geneformer and other Ensembl-based downstream tools.
#'
#' @param sce A Seurat object.
#' @param assay Name of the assay to use. Defaults to the current DefaultAssay.
#' @param db Annotation database object. Defaults to org.Hs.eg.db.
#' @param verbose Logical. Whether to print progress. Default TRUE.
#'
#' @return A new Seurat object with Ensembl ID row names.
#'
#' @examples
#' \dontrun{
#'   # Example with Seurat built-in dataset
#'   library(Seurat)
#'   data("pbmc_small")
#'   # Convert symbols to Ensembl
#'   pbmc_ens <- symbol2ensembl(sce = pbmc_small)
#' }
#'
#' @export
symbol2ensembl <- function(sce,
                           assay = NULL,
                           db = NULL,
                           verbose = TRUE) {

  # Check for necessary namespaces
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package not found.")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("AnnotationDbi package not found.")

  # Default to Human database if none provided
  if (is.null(db)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("org.Hs.eg.db not found.")
    db <- org.Hs.eg.db::org.Hs.eg.db
  }

  if (is.null(assay)) assay <- Seurat::DefaultAssay(sce)

  # 1. Extract raw counts (Compatible with V4 slots and V5 layers)
  counts_matrix <- Seurat::GetAssayData(sce, assay = assay, slot = "counts")
  symbols <- rownames(counts_matrix)

  if (verbose) message(paste("Mapping IDs for assay:", assay))

  # 2. Map Symbols to Ensembl using AnnotationDbi
  mapping <- AnnotationDbi::mapIds(
    x = db,
    keys = symbols,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  # 3. Filter out unmapped genes
  valid_indices <- !is.na(mapping)
  clean_mapping <- mapping[valid_indices]

  # 4. Handle duplicated Ensembl IDs (keep first occurrence)
  duplicate_check <- !duplicated(clean_mapping)
  final_mapping <- clean_mapping[duplicate_check]

  # 5. Subset matrix and assign new names
  final_counts <- counts_matrix[names(final_mapping), ]
  rownames(final_counts) <- final_mapping

  # 6. Reconstruct Seurat Object and inherit metadata
  new_obj <- Seurat::CreateSeuratObject(
    counts = final_counts,
    project = "Ensembl_Converted",
    meta.data = sce@meta.data[colnames(final_counts), , drop = FALSE]
  )

  if (verbose) message(paste("Successfully mapped", length(final_mapping), "genes."))
  return(new_obj)
}
