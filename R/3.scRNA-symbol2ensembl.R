# =========== convert_sce_id ===========
# ========== 1.convert_sce_id ==========

#' @title Convert Feature IDs (Symbol <-> Ensembl)
#' @description A tool to convert Seurat row names between Gene Symbols and Ensembl IDs.
#' Compatible with Geneformer requirements (Raw counts + Ensembl IDs).
#'
#' @param sce A Seurat object.
#' @param assay Name of the assay to use. Defaults to Seurat::DefaultAssay(sce).
#' @param species Input species. Options: "human", "mouse".
#' @param target_species Output species. Options: "same" (default), "human", "mouse".
#' @param from_type The current format of rownames. Options: "symbol", "ensembl".
#' @param to_type The desired output format. Options: "symbol", "ensembl".
#' @param db_human Human annotation database. Default org.Hs.eg.db.
#' @param db_mouse Mouse annotation database. Default org.Mm.eg.db.
#' @param verbose Logical. Whether to print progress. Default TRUE.
#'
#' @return A new Seurat object with converted row names.
#'
#' @export
convert_sce_id <- function(sce,
                               assay = NULL,
                               species = "human",
                               target_species = "same",
                               from_type = "symbol",
                               to_type = "ensembl",
                               db_human = NULL,
                               db_mouse = NULL,
                               verbose = TRUE) {

  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Package 'AnnotationDbi' is required.")

  if (is.null(assay)) assay <- Seurat::DefaultAssay(sce)
  counts_matrix <- Seurat::GetAssayData(sce, assay = assay, slot = "counts")
  current_features <- rownames(counts_matrix)

  if (length(current_features) == 0) stop("No genes found in the specified assay.")

  species <- tolower(species)
  if (target_species == "same") target_species <- species
  target_species <- tolower(target_species)
  from_type <- toupper(from_type)
  to_type <- toupper(to_type)

  if (verbose) {
    message(paste0(">>> Mode: ", species, " (", from_type, ") -> ", target_species, " (", to_type, ")"))
  }

  load_db <- function(db_arg, pkg_name) {
    if (!is.null(db_arg)) return(db_arg)
    if (!requireNamespace(pkg_name, quietly = TRUE)) stop(paste("Package", pkg_name, "is required."))
    return(getExportedValue(pkg_name, pkg_name))
  }

  db_in <- if (species == "human") load_db(db_human, "org.Hs.eg.db") else load_db(db_mouse, "org.Mm.eg.db")
  db_out <- if (target_species == "human") load_db(db_human, "org.Hs.eg.db") else load_db(db_mouse, "org.Mm.eg.db")

  # Step A: Normalize to Symbol
  if (from_type == "ENSEMBL") {
    if (verbose) message(">>> Step 1: Mapping Input Ensembl -> Source Symbol...")
    src_symbols <- tryCatch({
      AnnotationDbi::mapIds(db_in, keys = current_features,
                            column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    }, error = function(e) stop("Mapping failed. Check if from_type matches your data."))
  } else {
    src_symbols <- current_features
    names(src_symbols) <- current_features
  }

  # Step B: Cross-species Transformation
  query_symbols <- src_symbols
  if (species != target_species) {
    if (verbose) message(">>> Step 2: Transforming Symbols for target species...")
    if (species == "human" && target_species == "mouse") {
      to_title <- function(x) {
        if (is.na(x)) return(NA)
        s <- tolower(x); substr(s, 1, 1) <- toupper(substr(s, 1, 1)); s
      }
      query_symbols <- sapply(src_symbols, to_title)
    } else if (species == "mouse" && target_species == "human") {
      query_symbols <- toupper(src_symbols)
    }
  }

  # Step C: Final Mapping
  if (to_type == "ENSEMBL") {
    if (verbose) message(">>> Step 3: Mapping Target Symbol -> Target Ensembl...")
    valid_mask <- !is.na(query_symbols)
    final_mapping <- rep(NA, length(query_symbols))
    if (any(valid_mask)) {
      mapped_vals <- AnnotationDbi::mapIds(db_out, keys = as.character(query_symbols[valid_mask]),
                                           column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
      final_mapping[valid_mask] <- as.character(mapped_vals)
    }
  } else {
    final_mapping <- query_symbols
  }

  res_df <- data.frame(Old = current_features, New = as.character(final_mapping), stringsAsFactors = FALSE)
  res_df <- res_df[!is.na(res_df$New), ]
  res_df <- res_df[!duplicated(res_df$New), ]

  if (nrow(res_df) == 0) stop("Conversion resulted in 0 genes.")

  final_counts <- counts_matrix[res_df$Old, ]
  rownames(final_counts) <- res_df$New

  new_obj <- Seurat::CreateSeuratObject(
    counts = final_counts,
    project = "Geneformer_Ready",
    meta.data = sce@meta.data[colnames(final_counts), , drop = FALSE]
  )

  if (verbose) message(paste0(">>> Success: ", nrow(res_df), " features converted."))
  return(new_obj)
}

# =========== convert_id ===========
# ========== 2.convert_id ==========

#' @title Convert Gene IDs (Vector -> DataFrame)
#' @description Converts a vector of gene IDs to a target format/species.
#'
#' @param features Character vector of gene IDs.
#' @param species Input species ("human", "mouse").
#' @param target_species Output species ("same", "human", "mouse").
#' @param from_type Input format ("symbol", "ensembl").
#' @param to_type Output format ("symbol", "ensembl").
#' @param db_human Default org.Hs.eg.db.
#' @param db_mouse Default org.Mm.eg.db.
#' @param verbose Print progress. Default TRUE.
#'
#' @return A data frame with columns: 'Original_ID', 'Converted_ID'.
#'
#' @examples
#' \dontrun{
#'   genes <- c("TP53", "EGFR")
#'   # Human Symbol -> Human Ensembl
#'   df <- convert_id(genes, from_type = "symbol", to_type = "ensembl")
#'
#'   # Human Symbol -> Mouse Symbol
#'   df_ms <- convert_id(genes, target_species = "mouse", to_type = "symbol")
#' }
#'
#' @export
convert_id <- function(features,
                       species = "human",
                       target_species = "same",
                       from_type = "symbol",
                       to_type = "ensembl",
                       db_human = NULL,
                       db_mouse = NULL,
                       verbose = TRUE) {

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Package 'AnnotationDbi' is required.")

  features <- as.character(unique(features))
  if (length(features) == 0) stop("Input vector is empty.")

  species <- tolower(species)
  if (target_species == "same") target_species <- species
  target_species <- tolower(target_species)
  from_type <- toupper(from_type)
  to_type <- toupper(to_type)

  if (verbose) {
    message(paste0(">>> Mode: ", species, " (", from_type, ") -> ", target_species, " (", to_type, ")"))
  }

  load_db <- function(db_arg, pkg_name) {
    if (!is.null(db_arg)) return(db_arg)
    if (!requireNamespace(pkg_name, quietly = TRUE)) stop(paste("Package", pkg_name, "is required."))
    return(getExportedValue(pkg_name, pkg_name))
  }

  db_in <- if (species == "human") load_db(db_human, "org.Hs.eg.db") else load_db(db_mouse, "org.Mm.eg.db")
  db_out <- if (target_species == "human") load_db(db_human, "org.Hs.eg.db") else load_db(db_mouse, "org.Mm.eg.db")

  # Step A: Normalize to Symbol
  if (from_type == "ENSEMBL") {
    src_symbols <- tryCatch({
      AnnotationDbi::mapIds(db_in, keys = features,
                            column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    }, error = function(e) stop("Mapping failed. Check if from_type matches data."))
  } else {
    src_symbols <- features
    names(src_symbols) <- features
  }

  # Step B: Cross-species
  query_symbols <- src_symbols
  if (species != target_species) {
    if (species == "human" && target_species == "mouse") {
      to_title <- function(x) {
        if (is.na(x)) return(NA)
        s <- tolower(x); substr(s, 1, 1) <- toupper(substr(s, 1, 1)); s
      }
      query_symbols <- sapply(src_symbols, to_title)
    } else if (species == "mouse" && target_species == "human") {
      query_symbols <- toupper(src_symbols)
    }
  }

  # Step C: Final Mapping
  if (to_type == "ENSEMBL") {
    valid_mask <- !is.na(query_symbols)
    final_mapping <- rep(NA, length(query_symbols))
    if (any(valid_mask)) {
      mapped_vals <- AnnotationDbi::mapIds(db_out, keys = as.character(query_symbols[valid_mask]),
                                           column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
      final_mapping[valid_mask] <- as.character(mapped_vals)
    }
  } else {
    final_mapping <- query_symbols
  }

  # Output DataFrame
  res_df <- data.frame(
    Original_ID = features,
    Converted_ID = as.character(final_mapping),
    stringsAsFactors = FALSE
  )

  res_df <- res_df[!is.na(res_df$Converted_ID), ]

  if (verbose) message(paste0(">>> Success: ", nrow(res_df), "/", length(features), " converted."))

  return(res_df)
}
