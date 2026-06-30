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
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param save_path Path to save results.
#' @param downsample_n Integer. Total target number of cells (e.g., 3000).
#' @param stratify_by Character. Metadata column to stratify by (e.g., "celltype"). Default is NULL (uses active Idents).
#' @param use_seurat_var_genes Logical. Use Seurat's VariableFeatures (Recommended).
#' @param cores Number of cores.
#' @param scRNA Backward-compatible alias for `object`.
#'
#' @export
runMonocleAnalysis <- function(object = NULL,
                               save_path = "./output_data/monocle2.RData",
                               downsample_n = 5000,   # \u9ED8\u8BA4 5000\uFF0C\u9632\u6B62\u5185\u5B58\u5D29\u6E83
                               stratify_by = NULL,    # \u65B0\u589E\uFF1A\u6307\u5B9A\u5206\u5C42\u4F9D\u636E
                               use_seurat_var_genes = TRUE,
                               cores = 4,
                               scRNA = object) {
  if (is.null(scRNA)) stop("A Seurat object must be provided via 'object' or 'scRNA'.")

  # --- 1. \u73AF\u5883\u4E0E\u8865\u4E01 (Environment & Patches) ---
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("monocle", quietly = TRUE)) stop("Package 'monocle' is required.")
  if (!requireNamespace("VGAM", quietly = TRUE)) stop("Package 'VGAM' is required.")
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Package 'Biobase' is required.")
  monocle_fun <- function(name) get(name, envir = asNamespace("monocle"))
  monocle_newCellDataSet <- monocle_fun("newCellDataSet")
  monocle_estimateSizeFactors <- monocle_fun("estimateSizeFactors")
  monocle_estimateDispersions <- monocle_fun("estimateDispersions")
  monocle_dispersionTable <- monocle_fun("dispersionTable")
  monocle_setOrderingFilter <- monocle_fun("setOrderingFilter")
  monocle_reduceDimension <- monocle_fun("reduceDimension")
  monocle_orderCells <- monocle_fun("orderCells")

  # [Patch 1] \u4FEE\u590D igraph 'neimode' \u62A5\u9519
  tryCatch({
    raw_fun <- get("extract_ddrtree_ordering", envir = asNamespace("monocle"))
    if (any(grepl("neimode", utils::capture.output(body(raw_fun))))) {
      body(raw_fun) <- parse(text = gsub("neimode", "mode", utils::capture.output(body(raw_fun))))
      utils::assignInNamespace("extract_ddrtree_ordering", raw_fun, ns = "monocle")
    }
  }, error = function(e) warning("Patch failed: igraph compatibility might be broken."))

  # --- 2. \u667A\u80FD\u5206\u5C42\u4E0B\u91C7\u6837 (Stratified Downsampling) ---
  if (!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)

  total_cells <- ncol(scRNA)

  # \u53EA\u6709\u5F53\u7EC6\u80DE\u603B\u6570 > \u76EE\u6807\u6570\u65F6\u624D\u8FDB\u884C\u4E0B\u91C7\u6837
  if (!is.null(downsample_n) && is.numeric(downsample_n) && total_cells > downsample_n) {

    # \u786E\u5B9A\u5206\u5C42\u4F9D\u636E (\u4F18\u5148\u7528 stratify_by\uFF0C\u5426\u5219\u7528 Idents)
    if (is.null(stratify_by)) {
      group_col <- "Ident" # \u6807\u8BB0
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

    # --- \u5206\u5C42\u91C7\u6837\u6838\u5FC3\u903B\u8F91 ---
    # 1. \u8BA1\u7B97\u5168\u5C40\u4FDD\u7559\u6BD4\u4F8B
    keep_ratio <- downsample_n / total_cells

    # 2. \u6784\u5EFA\u4E34\u65F6\u6570\u636E\u6846\u7528\u4E8E dplyr \u64CD\u4F5C
    meta_df <- data.frame(cell_id = Seurat::Cells(scRNA), group = group_vec)

    # 3. \u6309\u7EC4\u5206\u5C42\u91C7\u6837 (\u4FDD\u8BC1\u6BCF\u7EC4\u90FD\u6709\u4EE3\u8868\uFF0C\u4E14\u4FDD\u7559\u6BD4\u4F8B\u4E00\u81F4)
    set.seed(123)
    sampled_cells <- unlist(lapply(split(meta_df$cell_id, meta_df$group), function(cells) {
      n_keep <- max(1L, floor(length(cells) * keep_ratio))
      sample(cells, min(length(cells), n_keep))
    }), use.names = FALSE)

    # 4. \u6267\u884C\u5B50\u96C6\u5316
    # \u5982\u679C\u67D0\u7EC4\u592A\u5C0F\u5BFC\u81F4\u6CA1\u62BD\u5230(\u6781\u7F55\u89C1)\uFF0Cslice_sample \u81F3\u5C11\u4F1A\u5C1D\u8BD5\u4FDD\u7559
    scRNA <- subset(scRNA, cells = sampled_cells)

    message(sprintf("    Final cell count: %d. (Maintained relative proportions of groups)", ncol(scRNA)))

  } else {
    message(sprintf(">>> [1/6] Processing all %d cells (No downsampling needed).", total_cells))
  }

  # --- 3. \u6784\u5EFA CDS (\u5305\u542B VGAM \u4FEE\u590D) ---
  counts_matrix <- Seurat::GetAssayData(scRNA, assay = "RNA", slot = "counts")
  pd <- Biobase::AnnotatedDataFrame(data = scRNA@meta.data)
  fd <- Biobase::AnnotatedDataFrame(data = data.frame(gene_short_name = rownames(counts_matrix),
                                                      row.names = rownames(counts_matrix)))

  # [Patch 2] VGAM::negbinomial \u76F4\u63A5\u8C03\u7528
  mycds <- monocle_newCellDataSet(counts_matrix, phenoData = pd, featureData = fd,
                                  lowerDetectionLimit = 0.5, expressionFamily = VGAM::negbinomial())

  # [Patch 3] S4 Class Fix
  mycds@expressionFamily@vfamily <- "negbinomial.size"

  # --- 4. \u9884\u5904\u7406 ---
  mycds <- monocle_estimateSizeFactors(mycds)

  message(">>> [2/6] Estimating dispersions...")
  dispersion_success <- FALSE
  tryCatch({
    mycds <- monocle_estimateDispersions(mycds, cores = cores, relative_expr = TRUE)
    dispersion_success <- TRUE
  }, error = function(e) message("    Warning: Blind dispersion estimation failed (using fallback)."))

  # --- 5. \u9009\u57FA\u56E0 ---
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
    disp_table <- monocle_dispersionTable(mycds)
    ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    message("    Using Monocle Dispersion Table.")
  }
  mycds <- monocle_setOrderingFilter(mycds, ordering_genes)

  # --- 6. \u964D\u7EF4\u4E0E\u6392\u5E8F ---
  message(">>> [4/6] Reducing Dimensions (DDRTree)...")
  # \u7A00\u758F\u77E9\u9635\u5FC5\u987B\u7528 log
  mycds <- monocle_reduceDimension(mycds, max_components = 2, method = 'DDRTree', verbose = FALSE, norm_method = "log")

  message(">>> [5/6] Ordering Cells (Forward & Reverse)...")
  mycds <- monocle_orderCells(mycds)
  mycds_reverse <- monocle_orderCells(mycds, reverse = TRUE)

  # --- 7. \u4FDD\u5B58 ---
  message(sprintf(">>> [6/6] Saving to: %s", save_path))
  save(mycds, mycds_reverse, file = save_path)

  return(list(cds = mycds, cds_reverse = mycds_reverse))
}


# ========= CytoTRACE =========
# ========= 2. CytoTRACE Trajectory Analysis (Optimized) =========
#' Run CytoTRACE Analysis (Simplified)
#'
#' @description A streamlined wrapper for CytoTRACE analysis on Seurat objects.
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param group.by Metadata column used as the phenotype annotation. Default `"celltype"`.
#' @param outdir Directory for CytoTRACE result files and plots.
#' @param ncores Number of cores passed to CytoTRACE.
#' @param scRNA Backward-compatible alias for `object`.
#' @export
runCytoTRACEAnalysis <- function(object = NULL,
                                 group.by = 'celltype',
                                 outdir = "./CytoTRACE_results",
                                 ncores = 8,
                                 scRNA = object) {
  if (is.null(scRNA)) stop("A Seurat object must be provided via 'object' or 'scRNA'.")

  # 0. \u4F9D\u8D56\u68C0\u67E5
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("CytoTRACE", quietly = TRUE)) stop("Package 'CytoTRACE' is required.")

  # 1. \u8DEF\u5F84\u51C6\u5907 (\u5408\u5E76\u6570\u636E\u548C\u56FE\u8868\u5230\u4E00\u4E2A\u76EE\u5F55\uFF0C\u51CF\u5C11\u6587\u4EF6\u5939\u5C42\u7EA7)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  result_file <- file.path(outdir, "CytoTRACE_obj.rds")

  # 2. \u6838\u5FC3\u903B\u8F91\uFF1A\u6709\u7F13\u5B58\u8BFB\u7F13\u5B58\uFF0C\u6CA1\u7F13\u5B58\u8DD1\u5206\u6790
  if (file.exists(result_file)) {
    message(">>> Loading cached CytoTRACE results...")
    results <- readRDS(result_file)
  } else {
    message(">>> Running CytoTRACE (this may take a while)...")

    # \u83B7\u53D6\u8868\u8FBE\u77E9\u9635 (\u81EA\u52A8\u5904\u7406 Assay\uFF0C\u65E0\u9700\u786C\u7F16\u7801 "RNA")
    # \u6CE8\u610F\uFF1ACytoTRACE \u5B98\u65B9\u5EFA\u8BAE\u4F7F\u7528 raw counts\uFF0C\u4F46\u5F88\u591A\u573A\u666F\u4E0B matrix \u8F6C\u6362\u6781\u5176\u6D88\u8017\u5185\u5B58
    # \u8FD9\u91CC\u4E3A\u4E86\u517C\u5BB9\u6027\u8F6C\u4E3A matrix\uFF0C\u4F46\u5EFA\u8BAE\u5927\u56FE\u6570\u636E\u6CE8\u610F\u5185\u5B58
    mat <- as.matrix(Seurat::GetAssayData(scRNA, slot = "counts"))

    # \u8FD0\u884C\u5206\u6790
    results <- CytoTRACE::CytoTRACE(mat = mat, ncores = ncores, subsamplesize = 1000) # \u53EF\u9009\uFF1Asubsamplesize \u52A0\u901F
    saveRDS(results, file = result_file)
  }

  # 3. \u53EF\u89C6\u5316
  message(">>> Generating plots...")

  # \u51C6\u5907\u8868\u578B\u6570\u636E
  phenotype <- as.character(scRNA[[group.by]][,1])
  names(phenotype) <- colnames(scRNA)

  # \u51C6\u5907\u964D\u7EF4\u6570\u636E (\u81EA\u52A8\u67E5\u627E umap \u6216 tsne)
  emb <- if ("umap" %in% names(scRNA@reductions)) Seurat::Embeddings(scRNA, "umap") else Seurat::Embeddings(scRNA, "tsne")

  # \u7ED8\u56FE
  CytoTRACE::plotCytoTRACE(results, phenotype = phenotype, emb = emb, outputDir = outdir)

  message(">>> Done! Check: ", outdir)
  return(results)
}
