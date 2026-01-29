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
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return A Seurat object.
#' @export
readH5AD <- function(file, data_type = "X", assay = "RNA", verbose = TRUE) {

  # 1. 检查依赖与文件
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. Please install it.")
  }
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # 加载 Python 模块
  ad <- reticulate::import("anndata", convert = FALSE)
  builtins <- reticulate::import_builtins(convert = FALSE) # 加载 Python 内置函数

  if (verbose) message("Reading H5AD file: ", file)
  adata <- ad$read_h5ad(file)

  # ============================================================================
  # 2. 确定读取的数据源 (X / raw / layers)
  # ============================================================================
  matrix_data <- NULL
  var_names_py <- NULL

  # --- 辅助函数：安全检查 raw 是否存在 ---
  check_raw_exists <- function(a) {
    tryCatch({
      !is.null(a$raw) && !is.null(a$raw$X)
    }, error = function(e) FALSE)
  }

  # --- 情况 A: 请求 "raw" ---
  if (data_type == "raw") {
    if (check_raw_exists(adata)) {
      if (verbose) message("Using raw counts from adata.raw.X")
      matrix_data <- adata$raw$X
      var_names_py <- adata$raw$var_names
    } else {
      warning("adata.raw not found or empty. Falling back to adata.X")
      data_type <- "X" # 降级
    }
  }

  # --- 情况 B: 请求 Layers (非 X 非 raw) ---
  if (data_type != "X" && data_type != "raw") {
    # 【修复重点】确保 layers_keys 是一个 R 向量
    layers_keys <- tryCatch({
      # 1. 先用 Python 的 list() 转换 dict_keys
      keys_py <- builtins$list(adata$layers$keys())
      # 2. 转为 R 对象
      keys_r <- reticulate::py_to_r(keys_py)
      # 3. 确保打平为向量 (unlist)
      unlist(keys_r)
    }, error = function(e) character(0))

    if (data_type %in% layers_keys) {
      if (verbose) message(sprintf("Using data from adata.layers['%s']", data_type))
      matrix_data <- adata$layers$get(data_type)
      var_names_py <- adata$var_names
    } else {
      # 如果 layers_keys 是 NULL，使用空字符避免打印报错
      avail_keys <- if (is.null(layers_keys)) "None" else paste(layers_keys, collapse=", ")
      warning(sprintf("Layer '%s' not found. Available: [%s]. Falling back to adata.X",
                      data_type, avail_keys))
      data_type <- "X" # 降级
    }
  }

  # --- 情况 C: 请求 "X" (或回退到 X) ---
  if (data_type == "X") {
    if (verbose) message("Using matrix from adata.X")
    matrix_data <- adata$X
    var_names_py <- adata$var_names
  }

  # ============================================================================
  # 3. 矩阵转换 (内存安全)
  # ============================================================================
  if (verbose) message("Converting matrix to sparse format...")
  counts <- reticulate::py_to_r(matrix_data)

  # 强制转为 dgCMatrix (R 标准稀疏格式)
  if (inherits(counts, "dgRMatrix")) {
    counts <- as(counts, "CsparseMatrix")
  } else if (!inherits(counts, "sparseMatrix")) {
    counts <- as(as.matrix(counts), "CsparseMatrix")
  }

  # 转置: AnnData (Cells x Genes) -> Seurat (Genes x Cells)
  counts <- Matrix::t(counts)

  # 4. 处理行名和列名
  obs_names <- reticulate::py_to_r(adata$obs_names$to_list())
  var_names <- reticulate::py_to_r(var_names_py$to_list())

  # 基因名去重
  if (any(duplicated(var_names))) {
    warning("Duplicate gene names found. Making them unique.")
    var_names <- make.unique(as.character(var_names))
  }

  rownames(counts) <- var_names
  colnames(counts) <- obs_names

  # 5. 元数据
  if (verbose) message("Processing metadata...")
  meta_data <- reticulate::py_to_r(adata$obs)
  if(nrow(meta_data) == ncol(counts)) {
    rownames(meta_data) <- colnames(counts)
  }

  # 6. 创建 Seurat 对象
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta_data, assay = assay)

  # 清理内存
  rm(counts, matrix_data)
  gc()

  # 7. 处理降维信息 (Embeddings)
  if (!reticulate::py_is_null_xptr(adata$obsm)) {

    obsm_keys_py <- builtins$list(adata$obsm$keys())
    obsm_keys <- tryCatch(reticulate::py_to_r(obsm_keys_py), error = function(e) NULL)
    # 同样确保它是向量
    obsm_keys <- unlist(obsm_keys)

    if (length(obsm_keys) > 0) {
      for (key in obsm_keys) {
        if (key == "spatial") next

        emb_matrix <- reticulate::py_to_r(adata$obsm$get(key))

        if (nrow(emb_matrix) != ncol(seurat_obj)) {
          if (verbose) warning(paste("Skipping", key, ": Dimensions do not match cell count."))
          next
        }

        rownames(emb_matrix) <- colnames(seurat_obj)
        clean_key <- gsub("^X_", "", key)
        colnames(emb_matrix) <- paste0(clean_key, "_", 1:ncol(emb_matrix))
        seurat_key <- paste0(toupper(clean_key), "_")

        if (verbose) message("Adding reduction: ", clean_key)

        seurat_obj[[clean_key]] <- Seurat::CreateDimReducObject(
          embeddings = emb_matrix,
          key = seurat_key,
          assay = assay,
          global = TRUE
        )
      }
    }
  }

  if (verbose) message("Done.")
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
#'   # 使用 Seurat 自带的小型示例数据集
#'   data("pbmc_small")
#'
#'   # 定义输出路径 (这里使用临时文件演示)
#'   out_file <- tempfile(fileext = ".h5ad")
#'
#'   # 1. 基础导出 (默认导出 DefaultAssay)
#'   writeH5AD(pbmc_small, file = out_file)
#'
#'   # 2. 自定义参数导出 (指定 Assay 和 Layer 名称)
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

  # --- 0. 内部辅助函数：V4/V5 兼容的数据提取 ---
  .get_matrix_data <- function(obj, assay, type) {
    # 动态检查环境是否支持 V5
    # 注意：这里不使用 importFrom，而是运行时检查，避免 V4 环境编译报错
    is_v5 <- tryCatch("LayerData" %in% getNamespaceExports("Seurat"), error = function(e) FALSE)
    mat <- NULL

    if (is_v5) {
      # Seurat V5 逻辑 (使用 :: 动态调用，避免静态依赖)
      # 这里的 Seurat::Layers 和 Seurat::LayerData 在 V4 环境下不会执行到，
      # 所以不会报错，且不需要在 importFrom 中声明
      avail_layers <- Seurat::Layers(obj, assay = assay)

      if (type %in% avail_layers) {
        if(verbose) message(paste("  [V5] Extracting layer:", type))
        mat <- Seurat::LayerData(obj, assay = assay, layer = type)
      } else if (type == "counts" && "counts" %in% names(obj@assays[[assay]])) {
        # V5 对象但可能是 V3/V4 转换过来的，没有 layer 只有 slot
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
      } else if (type == "data" && "data" %in% names(obj@assays[[assay]])) {
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
      }
    } else {
      # Seurat V4 逻辑
      if(verbose) message(paste("  [V4] Extracting slot:", type))
      try({
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = type)
      }, silent = TRUE)
    }

    # 如果是 BPCells (V5 on-disk)，需要转为内存矩阵
    if (!is.null(mat) && inherits(mat, "IterableMatrix")) {
      if(verbose) message("  [BPCells] Converting on-disk matrix to in-memory dgCMatrix...")
      mat <- as(mat, "dgCMatrix")
    }

    return(mat)
  }

  # 获取空间坐标 (兼容 V4/V5)
  .get_spatial_coords <- function(obj, img_key) {
    coords <- Seurat::GetTissueCoordinates(obj, image = img_key)
    cell_ids <- rownames(coords)
    if ("cell" %in% colnames(coords)) cell_ids <- coords$cell

    # 尝试匹配列名
    target_cols <- c("imagecol", "imagerow")
    if (!all(target_cols %in% colnames(coords))) {
      if (all(c("x", "y") %in% colnames(coords))) {
        target_cols <- c("x", "y")
      } else {
        # 最后的尝试：取前两列数值
        num_cols <- sapply(coords, is.numeric)
        target_cols <- colnames(coords)[num_cols][1:2]
      }
    }
    final_coords <- as.matrix(coords[, target_cols])
    rownames(final_coords) <- cell_ids
    colnames(final_coords) <- c("x", "y")
    return(final_coords)
  }

  # --- 1. 环境检查 ---
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Install 'reticulate' first.")

  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (verbose) message(paste("Exporting Assay:", assay))

  # 加载 Python 库
  ad <- reticulate::import("anndata", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  # --- 2. 准备数据 ---

  # 2.1 准备 X (Layer Data / Normalized)
  mat_X <- .get_matrix_data(object, assay, layer_data)
  if (is.null(mat_X)) {
    if (verbose) message("Layer 'data' not found. Using 'counts' as X.")
    mat_X <- .get_matrix_data(object, assay, layer_counts)
  }
  if (is.null(mat_X)) stop("Could not find expression data (counts or data).")

  mat_X <- Matrix::t(mat_X)
  mat_X <- as(mat_X, "dgCMatrix")

  # 2.2 准备 Layers (Raw Counts)
  # 提前提取 counts 并放入 list，传给构造函数
  layers_list <- list()
  mat_counts <- .get_matrix_data(object, assay, layer_counts)

  if (!is.null(mat_counts)) {
    if (verbose) message(paste0("Preparing layer '", layer_counts, "'..."))
    mat_counts <- Matrix::t(mat_counts)
    mat_counts <- as(mat_counts, "dgCMatrix")
    layers_list[[layer_counts]] <- mat_counts
  }

  # 2.3 准备 MetaData
  meta_df <- object@meta.data
  i <- sapply(meta_df, is.factor)
  meta_df[i] <- lapply(meta_df[i], as.character) # Factor -> Char

  # --- 3. 创建 AnnData ---
  if (verbose) message("Creating AnnData object...")

  # 直接在构造函数中传入 layers
  adata <- ad$AnnData(
    X = mat_X,
    obs = pd$DataFrame(meta_df),
    var = pd$DataFrame(index = colnames(mat_X)),
    layers = if(length(layers_list) > 0) layers_list else NULL
  )

  # --- 4. 导出降维 (Embeddings) ---
  if (export_embeddings && length(object@reductions) > 0) {
    if (verbose) message("Exporting reductions...")
    for (red in names(object@reductions)) {
      emb <- Seurat::Embeddings(object, reduction = red)
      key_name <- paste0("X_", red)
      if (grepl("^X_", red)) key_name <- red
      adata$obsm[[key_name]] <- np$array(emb)
    }
  }

  # --- 5. 导出空间信息 (Spatial) ---
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

  # --- 6. 写入文件 ---
  if (verbose) message(paste("Writing H5AD to:", file))
  adata$write_h5ad(file)
  if (verbose) message("Done.")
}
