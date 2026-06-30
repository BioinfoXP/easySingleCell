# =========== Gene ID conversion ===========

.esc_normalize_keytype <- function(x) {
  if (is.null(x)) return(NULL)
  x <- toupper(as.character(x))
  aliases <- c(
    "GENEID" = "ENTREZID",
    "ENTREZ" = "ENTREZID",
    "ENTREZGENE" = "ENTREZID",
    "ENSEMBLID" = "ENSEMBL",
    "SYMBOLS" = "SYMBOL",
    "GENESYMBOL" = "SYMBOL",
    "REFSEQID" = "REFSEQ",
    "UNIPROTID" = "UNIPROT"
  )
  if (x %in% names(aliases)) aliases[[x]] else x
}

.esc_load_orgdb <- function(species,
                            db = NULL,
                            db_human = NULL,
                            db_mouse = NULL,
                            db_rat = NULL) {
  if (!is.null(db)) return(db)

  species <- tolower(species)
  pkg_name <- switch(
    species,
    "human" = "org.Hs.eg.db",
    "homo sapiens" = "org.Hs.eg.db",
    "mouse" = "org.Mm.eg.db",
    "mus musculus" = "org.Mm.eg.db",
    "rat" = "org.Rn.eg.db",
    "rattus norvegicus" = "org.Rn.eg.db",
    stop("Unsupported species: ", species, ". Use 'human', 'mouse', 'rat', or pass an OrgDb object via 'db'.")
  )

  if (pkg_name == "org.Hs.eg.db" && !is.null(db_human)) return(db_human)
  if (pkg_name == "org.Mm.eg.db" && !is.null(db_mouse)) return(db_mouse)
  if (pkg_name == "org.Rn.eg.db" && !is.null(db_rat)) return(db_rat)

  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("Package '", pkg_name, "' is required. Install it with BiocManager::install('", pkg_name, "').")
  }
  getExportedValue(pkg_name, pkg_name)
}

.esc_validate_keytype <- function(db, keytype, arg_name = "keytype") {
  available <- AnnotationDbi::keytypes(db)
  if (!keytype %in% available) {
    stop(
      "'", arg_name, "' = '", keytype, "' is not available in this OrgDb. Available keytypes: ",
      paste(available, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.esc_symbol_bridge <- function(symbols, species, target_species) {
  species <- tolower(species)
  target_species <- tolower(target_species)
  if (identical(species, target_species)) return(symbols)

  warning(
    "Cross-species conversion uses a simple symbol-name bridge, not an ortholog database. ",
    "For publication-grade ortholog mapping, use a dedicated ortholog resource.",
    call. = FALSE
  )

  if (species %in% c("human", "homo sapiens") && target_species %in% c("mouse", "mus musculus", "rat", "rattus norvegicus")) {
    return(vapply(symbols, function(x) {
      if (is.na(x)) return(NA_character_)
      x <- tolower(x)
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }, character(1)))
  }

  if (species %in% c("mouse", "mus musculus", "rat", "rattus norvegicus") && target_species %in% c("human", "homo sapiens")) {
    return(toupper(symbols))
  }

  symbols
}

.esc_map_ids <- function(features,
                         species = "human",
                         target_species = "same",
                         from_type = "SYMBOL",
                         to_type = "ENSEMBL",
                         db = NULL,
                         target_db = NULL,
                         db_human = NULL,
                         db_mouse = NULL,
                         db_rat = NULL,
                         multiVals = "first") {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required.")
  }

  features <- as.character(features)
  species <- tolower(species)
  if (identical(target_species, "same")) target_species <- species
  target_species <- tolower(target_species)

  from_type <- .esc_normalize_keytype(from_type)
  to_type <- .esc_normalize_keytype(to_type)

  db_in <- .esc_load_orgdb(species, db = db, db_human = db_human, db_mouse = db_mouse, db_rat = db_rat)
  db_out <- .esc_load_orgdb(target_species, db = target_db, db_human = db_human, db_mouse = db_mouse, db_rat = db_rat)

  .esc_validate_keytype(db_in, from_type, "from_type")
  .esc_validate_keytype(db_out, to_type, "to_type")

  if (identical(species, target_species)) {
    mapped <- AnnotationDbi::mapIds(
      db_in,
      keys = features,
      column = to_type,
      keytype = from_type,
      multiVals = multiVals
    )
    return(as.character(mapped))
  }

  .esc_validate_keytype(db_in, "SYMBOL", "bridge SYMBOL")
  .esc_validate_keytype(db_out, "SYMBOL", "bridge SYMBOL")

  source_symbols <- if (from_type == "SYMBOL") {
    stats::setNames(features, features)
  } else {
    AnnotationDbi::mapIds(
      db_in,
      keys = features,
      column = "SYMBOL",
      keytype = from_type,
      multiVals = multiVals
    )
  }

  target_symbols <- .esc_symbol_bridge(as.character(source_symbols), species, target_species)
  valid <- !is.na(target_symbols)
  mapped <- rep(NA_character_, length(features))

  if (any(valid)) {
    if (to_type == "SYMBOL") {
      mapped[valid] <- target_symbols[valid]
    } else {
      mapped_vals <- AnnotationDbi::mapIds(
        db_out,
        keys = target_symbols[valid],
        column = to_type,
        keytype = "SYMBOL",
        multiVals = multiVals
      )
      mapped[valid] <- as.character(mapped_vals)
    }
  }

  mapped
}

.esc_get_assay_data <- function(object, assay, layer) {
  if ("LayerData" %in% getNamespaceExports("Seurat")) {
    mat <- tryCatch(
      getExportedValue("Seurat", "LayerData")(object, assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (!is.null(mat)) return(mat)
  }

  Seurat::GetAssayData(object, assay = assay, slot = layer)
}

#' @title Convert Feature IDs in a Seurat Object
#' @description Convert Seurat feature row names with AnnotationDbi keytypes.
#' Supports common IDs such as `SYMBOL`, `ENSEMBL`, `ENTREZID`, `ALIAS`,
#' `UNIPROT`, and `REFSEQ` when the selected OrgDb provides them.
#'
#' @param object A Seurat object. Old argument name `sce` is still supported.
#' @param sce Deprecated alias for `object`.
#' @param assay Name of the assay. Defaults to `Seurat::DefaultAssay(object)`.
#' @param layer Matrix layer/slot to convert. Default `"counts"`.
#' @param species Input species. Options include `"human"`, `"mouse"`, and `"rat"`.
#' @param target_species Output species. Use `"same"` to keep the same species.
#' @param from_type Input AnnotationDbi keytype. Lowercase aliases such as
#'   `"symbol"` and `"ensembl"` are accepted.
#' @param to_type Output AnnotationDbi keytype.
#' @param keytype Alias for `from_type`.
#' @param column Alias for `to_type`.
#' @param db Optional source OrgDb object.
#' @param target_db Optional target OrgDb object.
#' @param db_human Human annotation database. Default org.Hs.eg.db.
#' @param db_mouse Mouse annotation database. Default org.Mm.eg.db.
#' @param db_rat Rat annotation database. Default org.Rn.eg.db.
#' @param multiVals Passed to `AnnotationDbi::mapIds`. Default `"first"`.
#' @param project Project name for the returned Seurat object.
#' @param keep_unmapped Keep unmapped features with their original IDs. Default `FALSE`.
#' @param verbose Logical. Whether to print progress.
#'
#' @return A new Seurat object with converted row names and conversion metadata in `@misc$id_conversion`.
#' @export
convert_sce_id <- function(object = NULL,
                           sce = object,
                           assay = NULL,
                           layer = "counts",
                           species = "human",
                           target_species = "same",
                           from_type = "symbol",
                           to_type = "ensembl",
                           keytype = NULL,
                           column = NULL,
                           db = NULL,
                           target_db = NULL,
                           db_human = NULL,
                           db_mouse = NULL,
                           db_rat = NULL,
                           multiVals = "first",
                           project = NULL,
                           keep_unmapped = FALSE,
                           verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Package 'AnnotationDbi' is required.")

  object <- sce
  if (is.null(object)) stop("A Seurat object must be provided via 'object' or 'sce'.")

  if (!is.null(keytype)) from_type <- keytype
  if (!is.null(column)) to_type <- column

  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  counts_matrix <- .esc_get_assay_data(object, assay = assay, layer = layer)
  current_features <- rownames(counts_matrix)
  if (length(current_features) == 0) stop("No features found in the specified assay.")

  from_type_norm <- .esc_normalize_keytype(from_type)
  to_type_norm <- .esc_normalize_keytype(to_type)

  if (verbose) {
    message(">>> Mode: ", species, " (", from_type_norm, ") -> ", target_species, " (", to_type_norm, ")")
  }

  mapped <- .esc_map_ids(
    features = current_features,
    species = species,
    target_species = target_species,
    from_type = from_type_norm,
    to_type = to_type_norm,
    db = db,
    target_db = target_db,
    db_human = db_human,
    db_mouse = db_mouse,
    db_rat = db_rat,
    multiVals = multiVals
  )

  res_df <- data.frame(
    Old = current_features,
    New = as.character(mapped),
    stringsAsFactors = FALSE
  )

  if (keep_unmapped) {
    res_df$New[is.na(res_df$New) | res_df$New == ""] <- res_df$Old[is.na(res_df$New) | res_df$New == ""]
  } else {
    res_df <- res_df[!is.na(res_df$New) & res_df$New != "", , drop = FALSE]
  }

  res_df <- res_df[!duplicated(res_df$New), , drop = FALSE]
  if (nrow(res_df) == 0) stop("Conversion resulted in 0 features.")

  final_counts <- counts_matrix[res_df$Old, , drop = FALSE]
  rownames(final_counts) <- res_df$New

  if (is.null(project)) {
    project <- tryCatch(object@project.name, error = function(e) "ID_Converted")
    if (is.null(project) || project == "") project <- "ID_Converted"
  }

  new_obj <- Seurat::CreateSeuratObject(
    counts = final_counts,
    project = project,
    assay = assay,
    meta.data = object@meta.data[colnames(final_counts), , drop = FALSE]
  )

  new_obj@misc$id_conversion <- list(
    species = species,
    target_species = target_species,
    from_type = from_type_norm,
    to_type = to_type_norm,
    mapping = res_df
  )

  if (verbose) message(">>> Success: ", nrow(res_df), "/", length(current_features), " features converted.")
  new_obj
}

#' @title Convert Gene IDs
#' @description Convert a vector of IDs with AnnotationDbi keytypes. This is not
#' limited to symbol/Ensembl conversion; any keytype supported by the chosen
#' OrgDb can be used.
#'
#' @param features Character vector of gene IDs.
#' @param species Input species (`"human"`, `"mouse"`, or `"rat"`).
#' @param target_species Output species. Use `"same"` to keep the same species.
#' @param from_type Input AnnotationDbi keytype.
#' @param to_type Output AnnotationDbi keytype.
#' @param keytype Alias for `from_type`.
#' @param column Alias for `to_type`.
#' @param db Optional source OrgDb object.
#' @param target_db Optional target OrgDb object.
#' @param db_human Human annotation database. Default org.Hs.eg.db.
#' @param db_mouse Mouse annotation database. Default org.Mm.eg.db.
#' @param db_rat Rat annotation database. Default org.Rn.eg.db.
#' @param multiVals Passed to `AnnotationDbi::mapIds`. Default `"first"`.
#' @param keep_unmapped Keep unmapped features with original IDs. Default `FALSE`.
#' @param verbose Print progress. Default TRUE.
#'
#' @return A data frame with original IDs, converted IDs, keytypes, and species.
#'
#' @examples
#' \dontrun{
#' genes <- c("TP53", "EGFR")
#' convert_id(genes, from_type = "symbol", to_type = "ensembl")
#' convert_id(genes, from_type = "symbol", to_type = "entrezid")
#' convert_id(c("7157", "1956"), from_type = "entrezid", to_type = "symbol")
#' }
#' @export
convert_id <- function(features,
                       species = "human",
                       target_species = "same",
                       from_type = "symbol",
                       to_type = "ensembl",
                       keytype = NULL,
                       column = NULL,
                       db = NULL,
                       target_db = NULL,
                       db_human = NULL,
                       db_mouse = NULL,
                       db_rat = NULL,
                       multiVals = "first",
                       keep_unmapped = FALSE,
                       verbose = TRUE) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Package 'AnnotationDbi' is required.")

  features <- as.character(features)
  features <- features[!is.na(features) & features != ""]
  features <- unique(features)
  if (length(features) == 0) stop("Input vector is empty.")

  if (!is.null(keytype)) from_type <- keytype
  if (!is.null(column)) to_type <- column

  from_type_norm <- .esc_normalize_keytype(from_type)
  to_type_norm <- .esc_normalize_keytype(to_type)
  target_species_norm <- if (identical(target_species, "same")) tolower(species) else tolower(target_species)

  if (verbose) {
    message(">>> Mode: ", species, " (", from_type_norm, ") -> ", target_species, " (", to_type_norm, ")")
  }

  mapped <- .esc_map_ids(
    features = features,
    species = species,
    target_species = target_species,
    from_type = from_type_norm,
    to_type = to_type_norm,
    db = db,
    target_db = target_db,
    db_human = db_human,
    db_mouse = db_mouse,
    db_rat = db_rat,
    multiVals = multiVals
  )

  if (keep_unmapped) {
    mapped[is.na(mapped) | mapped == ""] <- features[is.na(mapped) | mapped == ""]
  }

  res_df <- data.frame(
    Original_ID = features,
    Converted_ID = as.character(mapped),
    FromType = from_type_norm,
    ToType = to_type_norm,
    SourceSpecies = tolower(species),
    TargetSpecies = target_species_norm,
    stringsAsFactors = FALSE
  )

  if (!keep_unmapped) {
    res_df <- res_df[!is.na(res_df$Converted_ID) & res_df$Converted_ID != "", , drop = FALSE]
  }

  if (verbose) message(">>> Success: ", nrow(res_df), "/", length(features), " converted.")
  res_df
}
