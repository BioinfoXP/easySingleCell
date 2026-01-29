#' @title Convert Gene Symbols to Ensembl IDs (Multi-species)
#' @description Converts Seurat object row names from Symbols to Ensembl IDs.
#' Optimized for Geneformer. Supports Human and Mouse, including cross-species conversion.
#'
#' @param sce A Seurat object.
#' @param assay Name of the assay to use. Defaults to `Seurat::DefaultAssay(sce)`.
#' @param species Input species. Options: "human" (default), "mouse".
#' @param target_species Output species for Ensembl IDs. Options: "same" (default), "human", "mouse".
#' @param db_human Human annotation database. Default `org.Hs.eg.db`.
#' @param db_mouse Mouse annotation database. Default `org.Mm.eg.db`.
#' @param verbose Logical. Whether to print progress messages. Default TRUE.
#'
#' @return A new Seurat object with Ensembl ID row names.
#'
#' @examples
#' \dontrun{
#'   # Human Symbol to Human Ensembl
#'   sce_ens <- symbol2ensembl(sce = pbmc_small)
#'
#'   # Human Symbol to Mouse Ensembl (Cross-species)
#'   sce_mouse_ens <- symbol2ensembl(sce = pbmc_small, target_species = "mouse")
#' }
#'
#' @export
symbol2ensembl <- function(sce,
                           assay = NULL,
                           species = "human",
                           target_species = "same",
                           db_human = NULL,
                           db_mouse = NULL,
                           verbose = TRUE) {

  # --- 1. Dependencies Check ---
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' required.")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Package 'AnnotationDbi' required.")

  # --- 2. Data Extraction ---
  if (is.null(assay)) assay <- Seurat::DefaultAssay(sce)
  counts_matrix <- Seurat::GetAssayData(sce, assay = assay, slot = "counts")
  current_symbols <- rownames(counts_matrix)

  if (length(current_symbols) == 0) stop("No genes found in the specified assay.")

  # --- 3. Resolve Species Parameters ---
  species <- tolower(species)
  if (target_species == "same") target_species <- species
  target_species <- tolower(target_species)

  if (verbose) {
    message(paste0(">>> Input Species: ", species))
    if (species != target_species) {
      message(paste0(">>> Cross-species Mode: ", species, " -> ", target_species))
    }
  }

  # --- 4. Load Target Database ---
  load_db <- function(db_arg, pkg_name) {
    if (!is.null(db_arg)) return(db_arg)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      stop(paste0("Package '", pkg_name, "' required for mapping."))
    }
    return(getExportedValue(pkg_name, pkg_name))
  }

  current_db <- if (target_species == "human") {
    load_db(db_human, "org.Hs.eg.db")
  } else if (target_species == "mouse") {
    load_db(db_mouse, "org.Mm.eg.db")
  } else {
    stop("Unsupported target species. Use 'human' or 'mouse'.")
  }

  # --- 5. Symbol Transformation ---
  query_keys <- current_symbols
  if (species == "human" && target_species == "mouse") {
    # TP53 -> Tp53
    to_title <- function(x) {
      s <- tolower(x); substr(s, 1, 1) <- toupper(substr(s, 1, 1)); s
    }
    query_keys <- sapply(current_symbols, to_title)
  } else if (species == "mouse" && target_species == "human") {
    # Tp53 -> TP53
    query_keys <- toupper(current_symbols)
  }

  # --- 6. Perform Mapping ---
  if (verbose) message(">>> Querying database...")
  mapping <- AnnotationDbi::mapIds(
    x = current_db,
    keys = query_keys,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  # --- 7. Reconstruction ---
  map_df <- data.frame(Original = current_symbols, Ensembl = as.character(mapping), stringsAsFactors = FALSE)
  map_df <- map_df[!is.na(map_df$Ensembl), ]
  map_df <- map_df[!duplicated(map_df$Ensembl), ]

  if (nrow(map_df) == 0) stop("No genes mapped to Ensembl.")

  final_counts <- counts_matrix[map_df$Original, ]
  rownames(final_counts) <- map_df$Ensembl

  new_obj <- Seurat::CreateSeuratObject(
    counts = final_counts,
    project = "Ensembl_Converted",
    meta.data = sce@meta.data[colnames(final_counts), , drop = FALSE]
  )

  if (verbose) message(paste0(">>> Success: ", nrow(map_df), " genes retained."))
  return(new_obj)
}
