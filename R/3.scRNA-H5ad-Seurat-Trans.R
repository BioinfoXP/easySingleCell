# =============== readH5AD ================
#' @title Read H5AD file to Seurat Object (Fixed Version)
#' @description Imports an AnnData object (.h5ad) and converts it to a Seurat object,
#' preserving counts, metadata, and dimensionality reductions.
#'
#' @param file Path to the .h5ad file.
#' @param data_type Character. Specifies which matrix to use:
#'   - "X": uses `adata.X` (default).
#'   - "raw": uses `adata.raw.X`.
#'   - "counts", "logcounts", etc.: uses `adata.layers['name']`.
#'   If the requested type is not found, it automatically falls back to "X".
#' @param assay Name of the assay to create in Seurat object. Default is "RNA".
#' @param project Project name for the returned Seurat object. If NULL, uses the input file name.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#' @param import_obsm Logical. Whether to import AnnData `obsm` entries.
#' @param import_uns Logical. Whether to import AnnData `uns` entries into `@misc$uns`.
#' @param store_non_reduction_obsm_as_meta Logical. Store non-reduction `obsm` matrices in `meta.data`.
#' @param skip_obsm_patterns Character vector of `obsm` name patterns to skip as reductions.
#' @param reduction_patterns Character vector of `obsm` name patterns to import as reductions.
#' @param make_unique_features Logical. Make duplicated feature names unique before creating the Seurat object.
#'
#' @return A Seurat object.
#' @export
readH5AD <- function(
    file,
    data_type = "X",
    assay = "RNA",
    project = NULL,
    verbose = TRUE,
    import_obsm = TRUE,
    import_uns = FALSE,
    store_non_reduction_obsm_as_meta = TRUE,
    skip_obsm_patterns = c("^spatial$", "fate_prob", "trend", "gene", "impute", "expression"),
    reduction_patterns = c("pca", "umap", "tsne", "diffmap", "diffusion", "dm_", "lsi", "ica"),
    make_unique_features = TRUE
) {

  # ============================================================================
  # 0. \u68C0\u67E5\u4F9D\u8D56
  # ============================================================================
  required_pkgs <- c("reticulate", "Matrix", "Seurat")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # ============================================================================
  # 1. \u5DE5\u5177\u51FD\u6570
  # ============================================================================
  .msg <- function(...) {
    if (isTRUE(verbose)) message(...)
  }

  .warn <- function(...) {
    warning(..., call. = FALSE)
  }

  .safe_py_to_r <- function(x, default = NULL) {
    tryCatch(
      reticulate::py_to_r(x),
      error = function(e) default
    )
  }

  .check_raw_exists <- function(a) {
    tryCatch({
      !is.null(a$raw) && !is.null(a$raw$X)
    }, error = function(e) FALSE)
  }

  .to_sparse_feature_by_cell <- function(x) {
    # x \u9884\u671F\u662F cells x genes\uFF0C\u9700\u8981\u8F6C\u6210 genes x cells
    xr <- reticulate::py_to_r(x)

    if (inherits(xr, "dgRMatrix")) {
      xr <- as(xr, "CsparseMatrix")
    } else if (inherits(xr, "dgCMatrix")) {
      # do nothing
    } else if (inherits(xr, "sparseMatrix")) {
      xr <- as(xr, "CsparseMatrix")
    } else {
      xr <- as(as.matrix(xr), "CsparseMatrix")
    }

    Matrix::t(xr)
  }

  .to_numeric_matrix <- function(x) {
    if (is.null(x)) return(NULL)

    if (is.data.frame(x)) {
      bad_cols <- vapply(x, is.list, logical(1))
      if (any(bad_cols)) return(NULL)
      x <- data.matrix(x)
    } else if (is.vector(x) && !is.list(x)) {
      x <- matrix(x, ncol = 1)
    } else if (!is.matrix(x)) {
      x <- tryCatch(as.matrix(x), error = function(e) NULL)
    }

    if (is.null(x)) return(NULL)
    if (length(dim(x)) != 2) return(NULL)

    suppressWarnings(storage.mode(x) <- "numeric")
    x <- as.matrix(x)

    if (!is.numeric(x)) return(NULL)
    if (nrow(x) == 0 || ncol(x) == 0) return(NULL)

    x
  }

  .make_unique_case_insensitive <- function(name, existing_names) {
    out <- name
    i <- 1
    while (tolower(out) %in% tolower(existing_names)) {
      i <- i + 1
      out <- paste0(name, "_", i)
    }
    out
  }

  .clean_reduction_name <- function(x) {
    x <- gsub("^X_", "", x)
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x
  }

  .make_seurat_key <- function(red_name) {
    key <- toupper(gsub("[^A-Za-z0-9]", "", red_name))
    paste0(key, "_")
  }

  .matches_any_pattern <- function(x, patterns) {
    any(vapply(patterns, function(p) grepl(p, x, ignore.case = TRUE), logical(1)))
  }

  .is_reduction_like <- function(name, patterns) {
    .matches_any_pattern(name, patterns)
  }

  .should_skip_obsm <- function(name, patterns) {
    .matches_any_pattern(name, patterns)
  }

  .safe_get_keys <- function(mapping_obj) {
    tryCatch({
      builtins <- reticulate::import_builtins(convert = FALSE)
      keys_py <- builtins$list(mapping_obj$keys())
      keys_r <- reticulate::py_to_r(keys_py)
      unlist(keys_r)
    }, error = function(e) character(0))
  }

  # ============================================================================
  # 2. \u8BFB\u53D6 h5ad
  # ============================================================================
  .msg("Reading H5AD file: ", file)

  ad <- reticulate::import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(file)

  # ============================================================================
  # 3. \u9009\u62E9\u8868\u8FBE\u77E9\u9635\u6765\u6E90
  # ============================================================================
  matrix_data <- NULL
  var_names_py <- NULL
  data_source_used <- NULL

  if (identical(data_type, "raw")) {
    if (.check_raw_exists(adata)) {
      .msg("Using raw counts from adata.raw.X")
      matrix_data <- adata$raw$X
      var_names_py <- adata$raw$var_names
      data_source_used <- "raw"
    } else {
      .warn("adata.raw not found or empty. Falling back to adata.X")
      data_type <- "X"
    }
  }

  if (!identical(data_type, "X") && !identical(data_type, "raw")) {
    layer_keys <- .safe_get_keys(adata$layers)

    if (data_type %in% layer_keys) {
      .msg("Using data from adata.layers['", data_type, "']")
      matrix_data <- adata$layers$get(data_type)
      var_names_py <- adata$var_names
      data_source_used <- paste0("layer:", data_type)
    } else {
      avail_keys <- if (length(layer_keys) == 0) "None" else paste(layer_keys, collapse = ", ")
      .warn("Layer '", data_type, "' not found. Available: [", avail_keys, "]. Falling back to adata.X")
      data_type <- "X"
    }
  }

  if (identical(data_type, "X")) {
    .msg("Using matrix from adata.X")
    matrix_data <- adata$X
    var_names_py <- adata$var_names
    data_source_used <- "X"
  }

  if (is.null(matrix_data)) {
    stop("Failed to retrieve matrix data from h5ad.")
  }

  # ============================================================================
  # 4. \u8868\u8FBE\u77E9\u9635\u8F6C\u6362
  # ============================================================================
  .msg("Converting matrix to sparse format...")
  counts <- .to_sparse_feature_by_cell(matrix_data)

  # obs / var names
  obs_names <- .safe_py_to_r(adata$obs_names$to_list())
  var_names <- .safe_py_to_r(var_names_py$to_list())

  if (is.null(obs_names) || is.null(var_names)) {
    stop("Failed to retrieve obs_names or var_names from AnnData object.")
  }

  obs_names <- as.character(obs_names)
  var_names <- as.character(var_names)

  if (length(obs_names) != ncol(counts)) {
    stop("Cell names length does not match matrix columns: ",
         length(obs_names), " vs ", ncol(counts))
  }

  if (length(var_names) != nrow(counts)) {
    stop("Feature names length does not match matrix rows: ",
         length(var_names), " vs ", nrow(counts))
  }

  if (anyDuplicated(var_names) > 0) {
    if (isTRUE(make_unique_features)) {
      .warn("Duplicate gene names found. Making them unique.")
      var_names <- make.unique(var_names)
    } else {
      .warn("Duplicate gene names found. Seurat may fail unless names are unique.")
    }
  }

  rownames(counts) <- var_names
  colnames(counts) <- obs_names

  # ============================================================================
  # 5. metadata
  # ============================================================================
  .msg("Processing metadata...")
  meta_data <- .safe_py_to_r(adata$obs, default = NULL)

  if (is.null(meta_data)) {
    meta_data <- data.frame(row.names = colnames(counts))
  } else {
    meta_data <- as.data.frame(meta_data)
    if (nrow(meta_data) == ncol(counts)) {
      rownames(meta_data) <- colnames(counts)
    } else {
      .warn("obs metadata row count does not match number of cells. Creating empty metadata instead.")
      meta_data <- data.frame(row.names = colnames(counts))
    }
  }

  # ============================================================================
  # 6. \u521B\u5EFA Seurat \u5BF9\u8C61
  # ============================================================================
  .msg("Creating Seurat object...")
  if (is.null(project)) {
    project <- tools::file_path_sans_ext(basename(file))
  }

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = meta_data,
    assay = assay,
    project = project
  )

  # \u8BB0\u5F55\u6570\u636E\u6765\u6E90
  seurat_obj@misc$h5ad_info <- list(
    file = normalizePath(file, winslash = "/", mustWork = FALSE),
    data_source = data_source_used
  )

  rm(matrix_data, counts)
  gc()

  # ============================================================================
  # 7. \u5BFC\u5165 obsm
  # ============================================================================
  if (isTRUE(import_obsm) && !reticulate::py_is_null_xptr(adata$obsm)) {

    .msg("Importing obsm entries...")
    obsm_keys <- .safe_get_keys(adata$obsm)

    if (length(obsm_keys) > 0) {
      for (key in obsm_keys) {

        clean_key <- .clean_reduction_name(key)
        .msg("Processing obsm: ", key)

        emb <- tryCatch({
          reticulate::py_to_r(adata$obsm$get(key))
        }, error = function(e) {
          .warn("Failed to read obsm '", key, "': ", conditionMessage(e))
          NULL
        })

        if (is.null(emb)) next

        emb_mat <- .to_numeric_matrix(emb)
        if (is.null(emb_mat)) {
          .warn("Skipping obsm '", key, "': cannot convert to a clean numeric matrix.")
          next
        }

        # \u9884\u671F obsm \u662F cells x dims
        if (nrow(emb_mat) != ncol(seurat_obj)) {
          .warn("Skipping obsm '", key, "': row count (", nrow(emb_mat),
                ") does not match number of cells (", ncol(seurat_obj), ").")
          next
        }

        # \u5148\u5224\u65AD\u662F\u5426\u8BE5\u8DF3\u8FC7 reduction \u5BFC\u5165
        if (.should_skip_obsm(clean_key, skip_obsm_patterns)) {
          if (isTRUE(store_non_reduction_obsm_as_meta)) {
            meta_df <- as.data.frame(emb_mat)
            rownames(meta_df) <- colnames(seurat_obj)

            safe_prefix <- gsub("[^A-Za-z0-9_]", "_", clean_key)
            colnames(meta_df) <- paste0(safe_prefix, "_", seq_len(ncol(meta_df)))

            seurat_obj@meta.data <- cbind(
              seurat_obj@meta.data,
              meta_df[rownames(seurat_obj@meta.data), , drop = FALSE]
            )
            .msg("Stored obsm '", key, "' into meta.data instead of reductions.")
          } else {
            .msg("Skipped obsm '", key, "'.")
          }
          next
        }

        # \u5982\u679C\u4E0D\u50CF\u6807\u51C6 reduction\uFF0C\u4E5F\u53EF\u6309\u9009\u9879\u5B58 meta.data
        if (!.is_reduction_like(clean_key, reduction_patterns)) {
          if (isTRUE(store_non_reduction_obsm_as_meta)) {
            meta_df <- as.data.frame(emb_mat)
            rownames(meta_df) <- colnames(seurat_obj)

            safe_prefix <- gsub("[^A-Za-z0-9_]", "_", clean_key)
            colnames(meta_df) <- paste0(safe_prefix, "_", seq_len(ncol(meta_df)))

            seurat_obj@meta.data <- cbind(
              seurat_obj@meta.data,
              meta_df[rownames(seurat_obj@meta.data), , drop = FALSE]
            )
            .msg("Stored non-reduction obsm '", key, "' into meta.data.")
          } else {
            .msg("Skipped non-reduction obsm '", key, "'.")
          }
          next
        }

        # reduction \u547D\u540D\u53BB\u91CD
        existing_reds <- names(seurat_obj@reductions)
        red_name <- .make_unique_case_insensitive(clean_key, existing_reds)

        seurat_key <- .make_seurat_key(red_name)

        rownames(emb_mat) <- colnames(seurat_obj)
        colnames(emb_mat) <- paste0(seurat_key, seq_len(ncol(emb_mat)))

        seurat_obj[[red_name]] <- Seurat::CreateDimReducObject(
          embeddings = emb_mat,
          key = seurat_key,
          assay = assay,
          global = TRUE
        )

        .msg("Added reduction: ", red_name)
      }
    }
  }

  # ============================================================================
  # 8. \u53EF\u9009\u5BFC\u5165 uns
  # ============================================================================
  if (isTRUE(import_uns) && !reticulate::py_is_null_xptr(adata$uns)) {
    .msg("Importing uns into seurat_obj@misc$uns ...")
    uns_keys <- .safe_get_keys(adata$uns)

    uns_list <- list()
    if (length(uns_keys) > 0) {
      for (uk in uns_keys) {
        uns_list[[uk]] <- tryCatch(
          reticulate::py_to_r(adata$uns$get(uk)),
          error = function(e) NULL
        )
      }
    }

    seurat_obj@misc$uns <- uns_list
  }

  .msg("Done.")
  return(seurat_obj)
}

# =============== writeH5AD ================
# =============== 2.writeH5AD  ================
#' @title Export Seurat Object to H5AD (Scanpy)
#' @description Converts a Seurat object into a .h5ad file compatible with Scanpy.
#' Fully compatible with both Seurat V4 and V5 structures.
#'
#' @param object A Seurat object.
#' @param file Output file path (must end in .h5ad).
#' @param assay Name of the assay to export. Default is `DefaultAssay(object)`.
#' @param layer_counts Layer/Slot name for raw counts. Default "counts".
#' @param layer_data Layer/Slot name for normalized data. Default "data".
#' @param export_embeddings Logical. Whether to export dimensionality reductions. Default TRUE.
#' @param verbose Logical. Print progress. Default TRUE.
#'
#' @examples
#' \dontrun{
#'   library(Seurat)
#'   # Use a small Seurat example dataset.
#'   data("pbmc_small")
#'
#'   # Define an output path for demonstration.
#'   out_file <- tempfile(fileext = ".h5ad")
#'
#'   # 1. Basic export using the default assay.
#'   writeH5AD(pbmc_small, file = out_file)
#'
#'   # 2. Custom export with explicit assay and layer names.
#'   # writeH5AD(object = pbmc_small,
#'   #           file = "output_sct.h5ad",
#'   #           assay = "SCT",
#'   #           layer_counts = "counts",
#'   #           layer_data = "data")
#' }
#'
#' @export
#' @importFrom Seurat DefaultAssay GetAssayData GetTissueCoordinates Embeddings Key
#' @importFrom Matrix t
#' @importFrom methods as
writeH5AD <- function(object,
                      file,
                      assay = NULL,
                      layer_counts = "counts",
                      layer_data = "data",
                      export_embeddings = TRUE,
                      verbose = TRUE) {

  # --- 0. \u5185\u90E8\u8F85\u52A9\u51FD\u6570\uFF1AV4/V5 \u517C\u5BB9\u7684\u6570\u636E\u63D0\u53D6 ---
  .get_matrix_data <- function(obj, assay, type) {
    # \u52A8\u6001\u68C0\u67E5\u73AF\u5883\u662F\u5426\u652F\u6301 V5
    # \u6CE8\u610F\uFF1A\u8FD9\u91CC\u4E0D\u4F7F\u7528 importFrom\uFF0C\u800C\u662F\u8FD0\u884C\u65F6\u68C0\u67E5\uFF0C\u907F\u514D V4 \u73AF\u5883\u7F16\u8BD1\u62A5\u9519
    is_v5 <- tryCatch("LayerData" %in% getNamespaceExports("Seurat"), error = function(e) FALSE)
    mat <- NULL

    if (is_v5) {
      # Seurat V5 \u903B\u8F91 (\u4F7F\u7528 :: \u52A8\u6001\u8C03\u7528\uFF0C\u907F\u514D\u9759\u6001\u4F9D\u8D56)
      # \u8FD9\u91CC\u7684 Seurat::Layers \u548C Seurat::LayerData \u5728 V4 \u73AF\u5883\u4E0B\u4E0D\u4F1A\u6267\u884C\u5230\uFF0C
      # \u6240\u4EE5\u4E0D\u4F1A\u62A5\u9519\uFF0C\u4E14\u4E0D\u9700\u8981\u5728 importFrom \u4E2D\u58F0\u660E
      layers_fun <- getExportedValue("Seurat", "Layers")
      layer_data_fun <- getExportedValue("Seurat", "LayerData")
      avail_layers <- layers_fun(obj, assay = assay)

      if (type %in% avail_layers) {
        if(verbose) message(paste("  [V5] Extracting layer:", type))
        mat <- layer_data_fun(obj, assay = assay, layer = type)
      } else if (type == "counts" && "counts" %in% names(obj@assays[[assay]])) {
        # V5 \u5BF9\u8C61\u4F46\u53EF\u80FD\u662F V3/V4 \u8F6C\u6362\u8FC7\u6765\u7684\uFF0C\u6CA1\u6709 layer \u53EA\u6709 slot
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
      } else if (type == "data" && "data" %in% names(obj@assays[[assay]])) {
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
      }
    } else {
      # Seurat V4 \u903B\u8F91
      if(verbose) message(paste("  [V4] Extracting slot:", type))
      try({
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = type)
      }, silent = TRUE)
    }

    # \u5982\u679C\u662F BPCells (V5 on-disk)\uFF0C\u9700\u8981\u8F6C\u4E3A\u5185\u5B58\u77E9\u9635
    if (!is.null(mat) && inherits(mat, "IterableMatrix")) {
      if(verbose) message("  [BPCells] Converting on-disk matrix to in-memory dgCMatrix...")
      mat <- as(mat, "dgCMatrix")
    }

    return(mat)
  }

  # \u83B7\u53D6\u7A7A\u95F4\u5750\u6807 (\u517C\u5BB9 V4/V5)
  .get_spatial_coords <- function(obj, img_key) {
    coords <- Seurat::GetTissueCoordinates(obj, image = img_key)
    cell_ids <- rownames(coords)
    if ("cell" %in% colnames(coords)) cell_ids <- coords$cell

    # \u5C1D\u8BD5\u5339\u914D\u5217\u540D
    target_cols <- c("imagecol", "imagerow")
    if (!all(target_cols %in% colnames(coords))) {
      if (all(c("x", "y") %in% colnames(coords))) {
        target_cols <- c("x", "y")
      } else {
        # \u6700\u540E\u7684\u5C1D\u8BD5\uFF1A\u53D6\u524D\u4E24\u5217\u6570\u503C
        num_cols <- sapply(coords, is.numeric)
        target_cols <- colnames(coords)[num_cols][1:2]
      }
    }
    final_coords <- as.matrix(coords[, target_cols])
    rownames(final_coords) <- cell_ids
    colnames(final_coords) <- c("x", "y")
    return(final_coords)
  }

  # --- 1. \u73AF\u5883\u68C0\u67E5 ---
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Install 'reticulate' first.")

  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (verbose) message(paste("Exporting Assay:", assay))

  # \u52A0\u8F7D Python \u5E93
  ad <- reticulate::import("anndata", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  # --- 2. \u51C6\u5907\u6570\u636E ---

  # 2.1 \u51C6\u5907 X (Layer Data / Normalized)
  mat_X <- .get_matrix_data(object, assay, layer_data)
  if (is.null(mat_X)) {
    if (verbose) message("Layer 'data' not found. Using 'counts' as X.")
    mat_X <- .get_matrix_data(object, assay, layer_counts)
  }
  if (is.null(mat_X)) stop("Could not find expression data (counts or data).")

  mat_X <- Matrix::t(mat_X)
  mat_X <- as(mat_X, "dgCMatrix")

  # 2.2 \u51C6\u5907 Layers (Raw Counts)
  # \u63D0\u524D\u63D0\u53D6 counts \u5E76\u653E\u5165 list\uFF0C\u4F20\u7ED9\u6784\u9020\u51FD\u6570
  layers_list <- list()
  mat_counts <- .get_matrix_data(object, assay, layer_counts)

  if (!is.null(mat_counts)) {
    if (verbose) message(paste0("Preparing layer '", layer_counts, "'..."))
    mat_counts <- Matrix::t(mat_counts)
    mat_counts <- as(mat_counts, "dgCMatrix")
    layers_list[[layer_counts]] <- mat_counts
  }

  # 2.3 \u51C6\u5907 MetaData
  meta_df <- object@meta.data
  i <- sapply(meta_df, is.factor)
  meta_df[i] <- lapply(meta_df[i], as.character) # Factor -> Char

  # --- 3. \u521B\u5EFA AnnData ---
  if (verbose) message("Creating AnnData object...")

  # \u76F4\u63A5\u5728\u6784\u9020\u51FD\u6570\u4E2D\u4F20\u5165 layers
  adata <- ad$AnnData(
    X = mat_X,
    obs = pd$DataFrame(meta_df),
    var = pd$DataFrame(index = colnames(mat_X)),
    layers = if(length(layers_list) > 0) layers_list else NULL
  )

  # --- 4. \u5BFC\u51FA\u964D\u7EF4 (Embeddings) ---
  if (export_embeddings && length(object@reductions) > 0) {
    if (verbose) message("Exporting reductions...")
    for (red in names(object@reductions)) {
      emb <- Seurat::Embeddings(object, reduction = red)
      key_name <- paste0("X_", red)
      if (grepl("^X_", red)) key_name <- red
      adata$obsm[[key_name]] <- np$array(emb)
    }
  }

  # --- 5. \u5BFC\u51FA\u7A7A\u95F4\u4FE1\u606F (Spatial) ---
  if (length(object@images) > 0) {
    if (verbose) message("Exporting spatial data (Images & Coordinates)...")
    img_key <- names(object@images)[1]
    image_obj <- object@images[[img_key]]

    spatial_coords <- .get_spatial_coords(object, img_key)
    cell_order <- unlist(adata$obs_names$to_list())

    spatial_ordered <- matrix(0, nrow = length(cell_order), ncol = 2)
    rownames(spatial_ordered) <- cell_order
    colnames(spatial_ordered) <- c("x", "y")

    common <- intersect(cell_order, rownames(spatial_coords))
    if(length(common) > 0) {
      spatial_ordered[common, ] <- spatial_coords[common, ]
    }

    adata$obsm[["spatial"]] <- np$array(spatial_ordered)

    scale_factors <- image_obj@scale.factors
    spatial_dict <- list()
    spatial_dict[[img_key]] <- list(
      images = list(hires = NULL, lowres = NULL),
      scalefactors = list(
        tissue_hires_scalef = scale_factors$hires,
        tissue_lowres_scalef = scale_factors$lowres,
        fiducial_diameter_fullres = scale_factors$fiducial,
        spot_diameter_fullres = scale_factors$spot
      ),
      metadata = list(chemistry_description = "Visium")
    )
    adata$uns[["spatial"]] <- reticulate::r_to_py(spatial_dict)
  }

  # --- 6. \u5199\u5165\u6587\u4EF6 ---
  if (verbose) message(paste("Writing H5AD to:", file))
  adata$write_h5ad(file)
  if (verbose) message("Done.")
}
