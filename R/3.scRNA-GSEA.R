# =============== Single-cell GSEA ================

#' Run GO GSEA from single-cell differential expression
#'
#' @description
#' A Seurat-oriented wrapper that runs `Seurat::FindMarkers()`, converts the
#' marker table into a pre-ranked gene vector, and then runs
#' `clusterProfiler::gseGO()`. If a marker table is already available, pass it
#' through `markers =` to reuse the ranking and GSEA step without rerunning
#' differential expression.
#'
#' @param object A Seurat object. Required when `markers = NULL`.
#' @param ident.1 Identity class, group name, or cell vector passed to
#'   `Seurat::FindMarkers()`.
#' @param ident.2 Optional second identity class, group name, or cell vector.
#' @param group.by Optional metadata column used to regroup cells before
#'   differential expression.
#' @param subset.ident Optional identity class to subset before regrouping.
#' @param assay Assay used by `Seurat::FindMarkers()`. Default uses Seurat's
#'   active assay.
#' @param slot Expression slot used by `Seurat::FindMarkers()`. Default
#'   `"data"`.
#' @param min.pct Minimum detection fraction passed to `FindMarkers()`.
#' @param logfc.threshold LogFC threshold passed to `FindMarkers()`. Default
#'   `0` keeps genes available for pre-ranked GSEA.
#' @param test.use Differential expression test passed to `FindMarkers()`.
#' @param latent.vars Optional latent variables passed to `FindMarkers()`.
#' @param markers Optional marker table from `FindMarkers()` or similar.
#' @param gene_col Optional column containing gene symbols or IDs. If `NULL`,
#'   row names are used when available, otherwise common columns such as
#'   `gene` or `symbol` are detected.
#' @param rank.by Ranking column. Use `"auto"` to prefer Seurat logFC columns
#'   such as `avg_log2FC`, then `avg_logFC`, then `stat`.
#' @param species One of `"human"`, `"mouse"`, or `"rat"` when `OrgDb` is not
#'   supplied.
#' @param org_db Optional organism annotation database.
#' @param OrgDb ClusterProfiler-style alias for `org_db`.
#' @param ont GO ontology: `"BP"`, `"MF"`, `"CC"`, or `"ALL"`.
#' @param key_type Input gene ID key type. Default `"SYMBOL"`.
#' @param keyType ClusterProfiler-style alias for `key_type`.
#' @param p_cutoff GSEA p-value cutoff.
#' @param pvalueCutoff ClusterProfiler-style alias for `p_cutoff`.
#' @param min_size Minimum gene set size.
#' @param minGSSize ClusterProfiler-style alias for `min_size`.
#' @param max_size Maximum gene set size.
#' @param maxGSSize ClusterProfiler-style alias for `max_size`.
#' @param seed Random seed used by `gseGO()`.
#' @param verbose Logical. Passed to Seurat and clusterProfiler.
#' @param ... Additional arguments passed to `Seurat::FindMarkers()` when
#'   `markers = NULL`.
#'
#' @return An `easy_scGSEA` list with `gsea`, `markers`, `gene_rank`, and
#'   `params`.
#' @export
#' @importFrom Seurat FindMarkers
#' @importFrom clusterProfiler gseGO
#'
#' @examples
#' \dontrun{
#' gsea_res <- scGSEA(
#'   object = sce,
#'   ident.1 = "Tumor",
#'   ident.2 = "Normal",
#'   group.by = "group",
#'   species = "human",
#'   ont = "BP"
#' )
#'
#' head(gsea_res$gsea)
#' head(gsea_res$markers)
#' }
scGSEA <- function(object = NULL,
                   ident.1 = NULL,
                   ident.2 = NULL,
                   group.by = NULL,
                   subset.ident = NULL,
                   assay = NULL,
                   slot = "data",
                   min.pct = 0.1,
                   logfc.threshold = 0,
                   test.use = "wilcox",
                   latent.vars = NULL,
                   markers = NULL,
                   gene_col = NULL,
                   rank.by = "auto",
                   species = c("human", "mouse", "rat"),
                   org_db = NULL,
                   OrgDb = NULL,
                   ont = "BP",
                   key_type = "SYMBOL",
                   keyType = NULL,
                   p_cutoff = 0.05,
                   pvalueCutoff = NULL,
                   min_size = 10,
                   minGSSize = NULL,
                   max_size = 500,
                   maxGSSize = NULL,
                   seed = 123,
                   verbose = TRUE,
                   ...) {
  if (!is.null(OrgDb)) org_db <- OrgDb
  if (is.null(org_db)) {
    species <- match.arg(species)
  } else {
    species <- as.character(species)[1]
    if (is.na(species) || !nzchar(species)) species <- "custom"
  }
  if (!is.null(keyType)) key_type <- keyType
  if (!is.null(pvalueCutoff)) p_cutoff <- pvalueCutoff
  if (!is.null(minGSSize)) min_size <- minGSSize
  if (!is.null(maxGSSize)) max_size <- maxGSSize

  if (is.null(markers)) {
    if (is.null(object)) {
      stop("A Seurat object must be provided via 'object' when 'markers' is NULL.", call. = FALSE)
    }
    if (is.null(ident.1)) {
      stop("'ident.1' is required when running scGSEA from a Seurat object.", call. = FALSE)
    }

    markers <- Seurat::FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group.by,
      subset.ident = subset.ident,
      assay = assay,
      slot = slot,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    )
  }

  .scGSEA_from_markers(
    markers = markers,
    gene_col = gene_col,
    rank.by = rank.by,
    species = species,
    org_db = org_db,
    ont = ont,
    key_type = key_type,
    p_cutoff = p_cutoff,
    min_size = min_size,
    max_size = max_size,
    seed = seed,
    verbose = verbose,
    params = list(
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group.by,
      subset.ident = subset.ident,
      assay = assay,
      slot = slot,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      latent.vars = latent.vars,
      species = species,
      ont = ont,
      key_type = key_type,
      p_cutoff = p_cutoff,
      min_size = min_size,
      max_size = max_size,
      seed = seed
    )
  )
}

.scGSEA_from_markers <- function(markers,
                                 gene_col = NULL,
                                 rank.by = "auto",
                                 species = "human",
                                 org_db = NULL,
                                 ont = "BP",
                                 key_type = "SYMBOL",
                                 p_cutoff = 0.05,
                                 min_size = 10,
                                 max_size = 500,
                                 seed = 123,
                                 verbose = TRUE,
                                 params = list(),
                                 gsea_fun = clusterProfiler::gseGO) {
  gene_rank <- .scGSEA_rank_markers(markers, rank.by = rank.by, gene_col = gene_col)
  if (length(gene_rank) < min_size) {
    stop(
      "Not enough ranked genes for GSEA: ", length(gene_rank),
      " genes after filtering, but min_size is ", min_size, ".",
      call. = FALSE
    )
  }

  org_db <- .scGSEA_orgdb(species = species, OrgDb = org_db)

  if (isTRUE(verbose)) {
    message("Running single-cell GO GSEA with ", length(gene_rank), " ranked genes...")
  }
  set.seed(seed)
  gsea_res <- gsea_fun(
    geneList = gene_rank,
    OrgDb = org_db,
    ont = ont,
    keyType = key_type,
    minGSSize = min_size,
    maxGSSize = max_size,
    pvalueCutoff = p_cutoff,
    verbose = verbose,
    seed = TRUE
  )

  gsea_n <- if (is.null(gsea_res)) 0 else nrow(as.data.frame(gsea_res))
  if (gsea_n == 0) {
    warning("No significant GSEA terms found.", call. = FALSE)
  } else if (isTRUE(verbose)) {
    message("single-cell GSEA done. Significant terms: ", gsea_n)
  }

  structure(
    list(
      gsea = gsea_res,
      markers = markers,
      gene_rank = gene_rank,
      params = params
    ),
    class = "easy_scGSEA"
  )
}

.scGSEA_rank_markers <- function(markers, rank.by = "auto", gene_col = NULL) {
  if (!is.data.frame(markers)) {
    stop("'markers' must be a data frame returned by Seurat::FindMarkers() or similar.", call. = FALSE)
  }
  if (nrow(markers) == 0) {
    stop("'markers' is empty.", call. = FALSE)
  }

  genes <- .scGSEA_marker_genes(markers, gene_col = gene_col)
  rank_col <- .scGSEA_rank_column(markers, rank.by = rank.by)
  scores <- markers[[rank_col]]
  if (!is.numeric(scores)) {
    stop("Ranking column '", rank_col, "' must be numeric.", call. = FALSE)
  }

  keep <- !is.na(scores) & is.finite(scores) & !is.na(genes) & nzchar(genes)
  genes <- as.character(genes[keep])
  scores <- scores[keep]
  if (length(scores) == 0) {
    stop("No valid genes remain after removing missing or non-finite ranking scores.", call. = FALSE)
  }

  ord <- order(scores, decreasing = TRUE, na.last = NA)
  genes <- genes[ord]
  scores <- scores[ord]
  keep_unique <- !duplicated(genes)
  stats::setNames(scores[keep_unique], genes[keep_unique])
}

.scGSEA_marker_genes <- function(markers, gene_col = NULL) {
  if (!is.null(gene_col)) {
    if (!gene_col %in% colnames(markers)) {
      stop("Gene column '", gene_col, "' was not found in markers.", call. = FALSE)
    }
    return(markers[[gene_col]])
  }

  rn <- rownames(markers)
  default_rn <- identical(rn, as.character(seq_len(nrow(markers))))
  if (!is.null(rn) && !default_rn && all(nzchar(rn))) {
    return(rn)
  }

  candidates <- c("gene", "symbol", "features", "feature", "Gene", "SYMBOL")
  hit <- candidates[candidates %in% colnames(markers)]
  if (length(hit) > 0) {
    return(markers[[hit[[1]]]])
  }

  stop(
    "Could not detect gene names. Provide 'gene_col' or keep genes as marker row names.",
    call. = FALSE
  )
}

.scGSEA_rank_column <- function(markers, rank.by = "auto") {
  if (length(rank.by) != 1 || is.na(rank.by) || !nzchar(rank.by)) {
    stop("'rank.by' must be one column name or 'auto'.", call. = FALSE)
  }
  if (!identical(rank.by, "auto")) {
    if (!rank.by %in% colnames(markers)) {
      stop("Ranking column '", rank.by, "' was not found in markers.", call. = FALSE)
    }
    return(rank.by)
  }

  candidates <- c("avg_log2FC", "avg_logFC", "avg_log2fc", "log2FC", "logFC", "stat")
  hit <- candidates[candidates %in% colnames(markers)]
  if (length(hit) == 0) {
    stop(
      "Could not detect a ranking column. Use 'rank.by' with a numeric column such as avg_log2FC.",
      call. = FALSE
    )
  }
  hit[[1]]
}

.scGSEA_orgdb <- function(species = "human", OrgDb = NULL) {
  if (!is.null(OrgDb)) {
    return(OrgDb)
  }

  species <- tolower(as.character(species)[1])
  if (species == "human") {
    return(org.Hs.eg.db::org.Hs.eg.db)
  }
  if (species == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package 'org.Mm.eg.db' is required for species = 'mouse'.", call. = FALSE)
    }
    return(getExportedValue("org.Mm.eg.db", "org.Mm.eg.db"))
  }
  if (species == "rat") {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
      stop("Package 'org.Rn.eg.db' is required for species = 'rat'.", call. = FALSE)
    }
    return(getExportedValue("org.Rn.eg.db", "org.Rn.eg.db"))
  }

  stop("Unsupported species: ", species, ". Use 'human', 'mouse', 'rat', or provide OrgDb.", call. = FALSE)
}

#' @export
print.easy_scGSEA <- function(x, ...) {
  gsea_n <- if (is.null(x$gsea)) 0 else nrow(as.data.frame(x$gsea))
  cat("<easySingleCell scGSEA>\n")
  cat("  ranked genes: ", length(x$gene_rank), "\n", sep = "")
  cat("  GSEA terms:   ", gsea_n, "\n", sep = "")
  invisible(x)
}
