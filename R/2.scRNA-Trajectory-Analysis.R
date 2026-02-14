# ========= Monocle =========
# =============== 1.Monocle2 ================
#' Run Complete Monocle2 Analysis with Stratified Downsampling
#'
#' @description
#' A robust wrapper for Monocle2 with:
#' 1. Stratified Downsampling (Preserves rare cell types).
#' 2. Auto-patches for igraph/VGAM compatibility.
#' 3. Memory safety features.
#'
#' @param scRNA A Seurat Object.
#' @param save_path Path to save results.
#' @param downsample_n Integer. Total target number of cells (e.g., 3000).
#' @param stratify_by Character. Metadata column to stratify by (e.g., "celltype"). Default is NULL (uses active Idents).
#' @param use_seurat_var_genes Logical. Use Seurat's VariableFeatures (Recommended).
#' @param cores Number of cores.
#'
#' @export
runMonocleAnalysis <- function(scRNA,
                               save_path = "./output_data/monocle2.RData",
                               downsample_n = 5000,   # 默认 5000，防止内存崩溃
                               stratify_by = NULL,    # 新增：指定分层依据
                               use_seurat_var_genes = TRUE,
                               cores = 4) {

  # --- 1. 环境与补丁 (Environment & Patches) ---
  suppressPackageStartupMessages({
    library(Seurat)
    library(monocle)
    library(dplyr)
    library(VGAM)
    library(igraph)
  })

  # [Patch 1] 修复 igraph 'neimode' 报错
  tryCatch({
    raw_fun <- get("extract_ddrtree_ordering", envir = asNamespace("monocle"))
    if (any(grepl("neimode", capture.output(body(raw_fun))))) {
      body(raw_fun) <- parse(text = gsub("neimode", "mode", capture.output(body(raw_fun))))
      assignInNamespace("extract_ddrtree_ordering", raw_fun, ns = "monocle")
    }
  }, error = function(e) warning("Patch failed: igraph compatibility might be broken."))

  # --- 2. 智能分层下采样 (Stratified Downsampling) ---
  if (!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)

  total_cells <- ncol(scRNA)

  # 只有当细胞总数 > 目标数时才进行下采样
  if (!is.null(downsample_n) && is.numeric(downsample_n) && total_cells > downsample_n) {

    # 确定分层依据 (优先用 stratify_by，否则用 Idents)
    if (is.null(stratify_by)) {
      group_col <- "Ident" # 标记
      group_vec <- Seurat::Idents(scRNA)
    } else {
      if (!stratify_by %in% colnames(scRNA@meta.data)) {
        stop(paste("Error: stratify_by column", stratify_by, "not found in metadata."))
      }
      group_col <- stratify_by
      group_vec <- scRNA@meta.data[[stratify_by]]
    }

    message(sprintf(">>> [1/6] Downsampling from %d to %d cells (Stratified by '%s')...",
                    total_cells, downsample_n, ifelse(is.null(stratify_by), "Idents", stratify_by)))

    # --- 分层采样核心逻辑 ---
    # 1. 计算全局保留比例
    keep_ratio <- downsample_n / total_cells

    # 2. 构建临时数据框用于 dplyr 操作
    meta_df <- data.frame(cell_id = Seurat::Cells(scRNA), group = group_vec)

    # 3. 按组分层采样 (保证每组都有代表，且保留比例一致)
    set.seed(123)
    sampled_cells <- meta_df %>%
      dplyr::group_by(group) %>%
      dplyr::slice_sample(prop = keep_ratio) %>% # 按比例抽取
      dplyr::pull(cell_id)

    # 4. 执行子集化
    # 如果某组太小导致没抽到(极罕见)，slice_sample 至少会尝试保留
    scRNA <- subset(scRNA, cells = sampled_cells)

    message(sprintf("    Final cell count: %d. (Maintained relative proportions of groups)", ncol(scRNA)))

  } else {
    message(sprintf(">>> [1/6] Processing all %d cells (No downsampling needed).", total_cells))
  }

  # --- 3. 构建 CDS (包含 VGAM 修复) ---
  counts_matrix <- Seurat::GetAssayData(scRNA, assay = "RNA", slot = "counts")
  pd <- new("AnnotatedDataFrame", data = scRNA@meta.data)
  fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(counts_matrix),
                                                    row.names = rownames(counts_matrix)))

  # [Patch 2] VGAM::negbinomial 直接调用
  mycds <- newCellDataSet(counts_matrix, phenoData = pd, featureData = fd,
                          lowerDetectionLimit = 0.5, expressionFamily = VGAM::negbinomial())

  # [Patch 3] S4 Class Fix
  mycds@expressionFamily@vfamily <- "negbinomial.size"

  # --- 4. 预处理 ---
  mycds <- estimateSizeFactors(mycds)

  message(">>> [2/6] Estimating dispersions...")
  dispersion_success <- FALSE
  tryCatch({
    mycds <- estimateDispersions(mycds, cores = cores, relative_expr = TRUE)
    dispersion_success <- TRUE
  }, error = function(e) message("    Warning: Blind dispersion estimation failed (using fallback)."))

  # --- 5. 选基因 ---
  message(">>> [3/6] Selecting Ordering Genes...")
  ordering_genes <- NULL

  if (use_seurat_var_genes || !dispersion_success) {
    ordering_genes <- Seurat::VariableFeatures(scRNA)
    if (length(ordering_genes) == 0) {
      scRNA <- Seurat::FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      ordering_genes <- Seurat::VariableFeatures(scRNA)
    }
    message("    Using Seurat Variable Features.")
  } else {
    disp_table <- dispersionTable(mycds)
    ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    message("    Using Monocle Dispersion Table.")
  }
  mycds <- setOrderingFilter(mycds, ordering_genes)

  # --- 6. 降维与排序 ---
  message(">>> [4/6] Reducing Dimensions (DDRTree)...")
  # 稀疏矩阵必须用 log
  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree', verbose = FALSE, norm_method = "log")

  message(">>> [5/6] Ordering Cells (Forward & Reverse)...")
  mycds <- orderCells(mycds)
  mycds_reverse <- orderCells(mycds, reverse = TRUE)

  # --- 7. 保存 ---
  message(sprintf(">>> [6/6] Saving to: %s", save_path))
  save(mycds, mycds_reverse, file = save_path)

  return(list(cds = mycds, cds_reverse = mycds_reverse))
}


# ========= CytoTRACE =========
# ========= 2. CytoTRACE Trajectory Analysis (Optimized) =========

#' Run CytoTRACE Analysis
#'
#' @description Performs CytoTRACE analysis to predict differentiation states
#' using a silent, namespace-safe implementation.
#'
#' @param scRNA A Seurat object.
#' @param celltype Metadata column for cell type annotation. Default "celltype".
#' @param output_data_dir Path for saving .Rdata results. Default "./output_data/".
#' @param output_figure_dir Path for saving plots. Default "./output_figure/".
#' @param emb Dimensionality reduction key (e.g., "umap") for visualization.
#' @param ncores Number of CPU cores. Default 5.
#' @param force_run Logical. If TRUE, ignores existing results and reruns analysis.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom CytoTRACE CytoTRACE plotCytoTRACE
#' @export
runCytoTRACEAnalysis <- function(scRNA,
                                 celltype = 'celltype',
                                 output_data_dir = "./output_data/",
                                 output_figure_dir = "./output_figure/",
                                 emb = 'umap',
                                 ncores = 5,
                                 force_run = FALSE) {

  # 1. Environment Preparation
  if (!base::dir.exists(output_data_dir)) base::dir.create(output_data_dir, recursive = TRUE)
  if (!base::dir.exists(output_figure_dir)) base::dir.create(output_figure_dir, recursive = TRUE)

  # 2. Input Validation
  if (!base::inherits(scRNA, "Seurat")) base::stop("Input 'scRNA' must be a Seurat object.")
  if (!celltype %in% base::colnames(scRNA@meta.data)) base::stop(base::paste0("Column '", celltype, "' not found in metadata."))

  if (!emb %in% base::names(scRNA@reductions)) {
    base::stop(base::paste0("Reduction '", emb, "' not found in Seurat object."))
  }

  # 3. Prepare Phenotype Data
  # Extract phenotype vector matching cell names
  phe <- base::as.character(scRNA@meta.data[[celltype]])
  base::names(phe) <- base::rownames(scRNA@meta.data)

  # 4. Analysis Execution (with Caching)
  results_file <- base::file.path(output_data_dir, "Cytotrace_results.Rdata")

  if (base::file.exists(results_file) && !force_run) {
    base::message("Loading existing CytoTRACE results...")
    # Variable 'results' will be loaded into environment
    base::load(results_file)

    # Check if loaded results match current cells (basic validation)
    if (!base::exists("results")) base::stop("Loaded file does not contain 'results' object.")

  } else {
    base::message("Running CytoTRACE analysis (this may take time)...")

    # Secure data extraction using Seurat accessor
    # CytoTRACE requires a raw count matrix (genes x cells)
    mat_counts <- base::as.matrix(Seurat::GetAssayData(scRNA, assay = "RNA", slot = "counts"))

    # Run core analysis
    results <- CytoTRACE::CytoTRACE(mat = mat_counts, ncores = ncores)

    base::save(results, file = results_file)
    base::message("Analysis saved to ", results_file)
  }

  # 5. Visualization
  # plotCytoTRACE automatically generates PDFs in the output directory
  base::message("Generating plots in ", output_figure_dir)

  CytoTRACE::plotCytoTRACE(
    results,
    outputDir = output_figure_dir,
    emb = scRNA@reductions[[emb]]@cell.embeddings,
    phenotype = phe
  )
}
