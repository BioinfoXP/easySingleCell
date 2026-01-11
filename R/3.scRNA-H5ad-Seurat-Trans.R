# =============== readH5AD ================
# =============== 1.readH5AD  ================
#' @title Read H5AD file to Seurat Object
#' @description Imports an AnnData object (.h5ad) and converts it to a Seurat object,
#' preserving counts, metadata, and dimensionality reductions.
#'
#' @param file Path to the .h5ad file.
#' @param use_raw Logical. If TRUE, tries to use `adata.raw.X` (raw counts) as the expression matrix.
#' If `adata.raw` is not present, falls back to `adata.X`. Default is TRUE.
#' @param assay Name of the assay to create in Seurat object. Default is "RNA".
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return A Seurat object.
#' @export
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject
#' @importFrom Matrix t
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#'   library(Seurat)
#'
#'   # 1. 基础用法：读取 H5AD 文件
#'   # 默认会自动寻找 raw counts
#'   seu <- readH5AD("input_data/dataset.h5ad")
#'   print(seu)
#'
#'   # 2. 强制读取 Normalized Data (adata.X)
#'   # 如果你不需要 raw data，或者 h5ad 里没有 raw
#'   seu_norm <- readH5AD("input_data/dataset.h5ad", use_raw = FALSE)
#'
#'   # 3. 指定 Assay 名称 (例如 Spatial)
#'   seu_spatial <- readH5AD("input_data/spatial.h5ad", assay = "Spatial")
#'
#'   # 4. 验证降维结果
#'   # 读取后可以直接画图
#'   DimPlot(seu, reduction = "umap")
#' }
readH5AD <- function(file, use_raw = TRUE, assay = "RNA", verbose = TRUE) {

  # 1. 检查依赖与文件
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. Please install it.")
  }
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # 加载 Python 模块
  ad <- reticulate::import("anndata", convert = FALSE)
  builtins <- reticulate::import_builtins(convert = FALSE) # 加载 Python 内置函数 list()

  if (verbose) message("Reading H5AD file: ", file)
  adata <- ad$read_h5ad(file)

  # 2. 确定使用 Raw 还是 X
  has_raw <- !reticulate::py_is_null_xptr(adata$raw) && !is.null(adata$raw$X)

  if (use_raw && has_raw) {
    if (verbose) message("Using raw counts from adata.raw.X")
    matrix_data <- adata$raw$X
    var_names <- reticulate::py_to_r(adata$raw$var_names$to_list())
  } else {
    if (use_raw && verbose) message("adata.raw not found or empty. Falling back to adata.X")
    if (!use_raw && verbose) message("Using matrix from adata.X (as requested)")
    matrix_data <- adata$X
    var_names <- reticulate::py_to_r(adata$var_names$to_list())
  }

  # 3. 矩阵转换 (内存安全)
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

    # === 关键修复点 ===
    # adata.obsm.keys() 是 Python 的 dict_keys 对象，不能直接在 R 中循环
    # 必须调用 Python 的 list() 显式转换为列表，再转为 R 向量
    obsm_keys_py <- builtins$list(adata$obsm$keys())
    obsm_keys <- reticulate::py_to_r(obsm_keys_py)

    if (length(obsm_keys) > 0) {
      for (key in obsm_keys) {
        # 跳过 spatial 坐标 (通常不作为 DimReduc 处理，或者需要特殊处理)
        if (key == "spatial") next

        # 获取矩阵
        emb_matrix <- reticulate::py_to_r(adata$obsm$get(key))

        # 简单检查行数匹配
        if (nrow(emb_matrix) != ncol(seurat_obj)) {
          if (verbose) warning(paste("Skipping", key, ": Dimensions do not match cell count."))
          next
        }

        rownames(emb_matrix) <- colnames(seurat_obj)

        # 清洗 Key 名称 (例如 X_umap -> umap)
        clean_key <- gsub("^X_", "", key)

        # 设置列名 (umap_1, umap_2)
        colnames(emb_matrix) <- paste0(clean_key, "_", 1:ncol(emb_matrix))

        # Seurat Key (UMAP_)
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

  # 获取矩阵数据 (兼容 V4 slot 和 V5 layer)
  .get_matrix_data <- function(obj, assay, type) {
    # 尝试方法 A: V5 Layers
    # 检查是否为 V5 对象且有 layers
    is_v5 <- tryCatch("LayerData" %in% getNamespaceExports("Seurat"), error = function(e) FALSE)

    mat <- NULL

    if (is_v5) {
      # Seurat V5 逻辑
      # 检查 layer 是否存在
      avail_layers <- Seurat::Layers(obj, assay = assay)
      if (type %in% avail_layers) {
        if(verbose) message(paste("  [V5] Extracting layer:", type))
        mat <- Seurat::LayerData(obj, assay = assay, layer = type)
      } else if (type == "counts" && "counts" %in% names(obj@assays[[assay]])) {
        # 回退兼容
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
      } else if (type == "data" && "data" %in% names(obj@assays[[assay]])) {
        mat <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
      }
    } else {
      # Seurat V4 逻辑
      if(verbose) message(paste("  [V4] Extracting slot:", type))
      # V4 通常 slot 叫 counts 或 data
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

  # 获取空间坐标 (兼容 V4/V5 返回格式差异)
  .get_spatial_coords <- function(obj, img_key) {
    # V5 可能返回 columns: x, y, cell
    # V4 可能返回 columns: imagecol, imagerow (rownames 是 cell)
    coords <- Seurat::GetTissueCoordinates(obj, image = img_key)

    # 统一获取 Cell IDs
    cell_ids <- rownames(coords)
    if ("cell" %in% colnames(coords)) {
      cell_ids <- coords$cell
    }

    # 统一获取 X, Y (Scanpy 标准: columns=[x, y])
    # 注意: imagecol 是 x, imagerow 是 y
    target_cols <- c("imagecol", "imagerow")

    if (!all(target_cols %in% colnames(coords))) {
      # 尝试 fallback 到 x, y
      if (all(c("x", "y") %in% colnames(coords))) {
        # 注意: Seurat V5 有时 x=row(y), y=col(x)，这很混乱。
        # 通常 Visium 标准是: imagecol=x, imagerow=y.
        # 如果只有 x, y，通常假设对应 col, row
        target_cols <- c("x", "y")
      } else {
        # 最后的尝试：取前两列数值
        num_cols <- sapply(coords, is.numeric)
        target_cols <- colnames(coords)[num_cols][1:2]
      }
    }

    # 提取并重组
    final_coords <- as.matrix(coords[, target_cols])
    rownames(final_coords) <- cell_ids
    colnames(final_coords) <- c("x", "y") # 强制重命名以便 Python 识别

    return(final_coords)
  }

  # --- 1. 环境检查 ---

  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Install 'reticulate' first.")

  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (verbose) message(paste("Exporting Assay:", assay))

  # 加载 Python
  ad <- reticulate::import("anndata", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  # --- 2. 准备矩阵 (X 和 Layers) ---

  # 2.1 获取 X (Data/Normalized)
  mat_X <- .get_matrix_data(object, assay, layer_data)

  # 如果没有 data 层，尝试用 counts 层顶替 X
  if (is.null(mat_X)) {
    if (verbose) message("Layer 'data' not found. Using 'counts' as X.")
    mat_X <- .get_matrix_data(object, assay, layer_counts)
  }

  if (is.null(mat_X)) stop("Could not find expression data (counts or data).")

  # 转置 + 稀疏化
  mat_X <- Matrix::t(mat_X)
  mat_X <- as(mat_X, "dgCMatrix")

  # 2.2 准备 MetaData (类型安全转换)
  meta_df <- object@meta.data
  i <- sapply(meta_df, is.factor)
  meta_df[i] <- lapply(meta_df[i], as.character) # Factor -> Char

  # 2.3 创建 AnnData
  if (verbose) message("Creating AnnData object...")
  adata <- ad$AnnData(
    X = mat_X,
    obs = pd$DataFrame(meta_df),
    var = pd$DataFrame(index = colnames(mat_X))
  )

  # 2.4 添加 Raw Counts (如果存在且不同于 X)
  mat_counts <- .get_matrix_data(object, assay, layer_counts)
  if (!is.null(mat_counts)) {
    if (verbose) message("Adding 'counts' to layers...")
    mat_counts <- Matrix::t(mat_counts)
    mat_counts <- as(mat_counts, "dgCMatrix")
    adata$layers[[layer_counts]] <- mat_counts
  }

  # --- 3. 导出降维 (Embeddings) ---

  if (export_embeddings && length(object@reductions) > 0) {
    if (verbose) message("Exporting reductions...")
    for (red in names(object@reductions)) {
      emb <- Seurat::Embeddings(object, reduction = red)

      # 命名修正: umap -> X_umap
      key_name <- paste0("X_", red)
      if (grepl("^X_", red)) key_name <- red

      adata$obsm[[key_name]] <- np$array(emb)
    }
  }

  # --- 4. 导出空间信息 (Spatial) ---

  if (length(object@images) > 0) {
    if (verbose) message("Exporting spatial data (Images & Coordinates)...")

    # 默认取第一张图
    img_key <- names(object@images)[1]
    image_obj <- object@images[[img_key]]

    # 4.1 获取坐标 (使用兼容函数)
    spatial_coords <- .get_spatial_coords(object, img_key)

    # 确保坐标顺序与 adata.obs (Cells) 一致
    # Seurat 可能会过滤掉一些没有坐标的细胞，或者顺序不同
    common_cells <- intersect(rownames(adata$obs_names$to_list()), rownames(spatial_coords))

    if (length(common_cells) < nrow(spatial_coords)) {
      warning("Some cells in object do not have spatial coordinates.")
    }

    # Scanpy 要求 obsm['spatial'] 必须和 obs 行数一一对应且顺序一致
    # 这里的策略：先创建一个全是 NaN 的矩阵，然后填入有坐标的值
    # 或者让 Scanpy 处理。更安全的方法是只存匹配的。
    # 但由于 h5ad 结构限制，obsm 必须和 shape[0] 一致。

    # 简单处理：我们假设 Seurat 对象里的细胞都有坐标 (如果是空间对象)
    # 按照 adata.obs_names 重排坐标
    cell_order <- unlist(adata$obs_names$to_list())
    spatial_ordered <- spatial_coords[cell_order, , drop=FALSE]

    # 替换 NA 为 0 (防止 Python 报错，虽然 Ideally 不该有 NA)
    spatial_ordered[is.na(spatial_ordered)] <- 0

    adata$obsm[["spatial"]] <- np$array(spatial_ordered)

    # 4.2 Scale Factors
    scale_factors <- image_obj@scale.factors

    spatial_dict <- list()
    spatial_dict[[img_key]] <- list(
      images = list(hires = NULL, lowres = NULL), # 不导出图片内容，防文件过大
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

  # --- 5. 写入 ---
  if (verbose) message(paste("Writing H5AD to:", file))
  adata$write_h5ad(file)
  if (verbose) message("Done.")
}
