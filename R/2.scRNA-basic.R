# =============== \u8D28\u63A7 ================
# =============== 1.QC  ================

#' @title scVis: Pure QC (No Plotting)
#' @description Calculates mitochondrial percentage and filters cells based on thresholds.
#' Does NOT generate plots or save files. Fast and robust for Seurat v4/v5.
#'
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param sce Backward-compatible alias for `object`.
#' @param minGene Min detected genes (nFeature_RNA). Default 200.
#' @param maxGene Max detected genes. Default 6000.
#' @param pctMT Max mitochondrial percentage. Default 20.
#' @param maxCounts Max UMI counts (nCount_RNA). Default 20000.
#' @param species "human" (starts with MT-) or "mouse" (starts with mt-).
#' @param mt_pattern Custom regex for mitochondrial genes. If NULL, infers from species.
#'
#' @return A filtered Seurat object.
#' @export
#' @importFrom Seurat PercentageFeatureSet
#'
runScRNAQC <- function(object = NULL,
                       minGene = 200,
                       maxGene = 6000,
                       pctMT = 20,
                       maxCounts = 20000,
                       species = 'human',
                       mt_pattern = NULL,
                       sce = object) {

  # --- 1. \u57FA\u7840\u68C0\u67E5 ---
  if (is.null(sce)) stop("A Seurat object must be provided via 'object' or 'sce'.")

  # --- 2. \u8BA1\u7B97\u7EBF\u7C92\u4F53\u6BD4\u4F8B ---
  if (is.null(mt_pattern)) {
    mt_pattern <- switch(tolower(species),
                         "human" = "^MT-",
                         "mouse" = "^mt-",
                         stop("\u274C Unsupported species. Please provide 'mt_pattern' manually."))
  }

  # \u5982\u679C percent.mt \u5C1A\u672A\u8BA1\u7B97\uFF0C\u5219\u8BA1\u7B97\u5B83
  if (!"percent.mt" %in% colnames(sce@meta.data)) {
    message(paste0("\u2139\uFE0F Calculating mitochondrial percentage using pattern: ", mt_pattern))
    sce[["percent.mt"]] <- Seurat::PercentageFeatureSet(sce, pattern = mt_pattern)
  }

  # --- 3. \u6267\u884C\u8FC7\u6EE4 (\u4F7F\u7528 Metadata \u7D22\u5F15\uFF0C\u6700\u7A33\u5065\u7684\u65B9\u5F0F) ---
  n_before <- ncol(sce)
  meta <- sce@meta.data

  # \u627E\u51FA\u7B26\u5408\u6761\u4EF6\u7684\u7EC6\u80DE ID
  keep_cells <- rownames(meta)[
    meta$nFeature_RNA > minGene &
      meta$nFeature_RNA < maxGene &
      meta$percent.mt < pctMT &
      meta$nCount_RNA < maxCounts
  ]

  n_after <- length(keep_cells)
  n_removed <- n_before - n_after

  # --- 4. \u6253\u5370\u62A5\u544A ---
  cat(sprintf("\n====== scRNA-seq QC Report ======\n"))
  cat(sprintf("Input Cells:      %d\n", n_before))
  cat(sprintf("Filtering Criteria:\n"))
  cat(sprintf("  - Gene Count:   %d < nFeature < %d\n", minGene, maxGene))
  cat(sprintf("  - Mito Percent: < %s%%\n", pctMT))
  cat(sprintf("  - Max UMI:      < %d\n", maxCounts))
  cat(sprintf("---------------------------------\n"))
  cat(sprintf("Remaining Cells:  %d\n", n_after))
  cat(sprintf("Removed Cells:    %d (%.2f%%)\n", n_removed, (n_removed/n_before)*100))
  cat(sprintf("=================================\n\n"))

  if (n_after == 0) stop("\u274C Error: All cells were filtered out! Please adjust thresholds.")

  # --- 5. \u8FD4\u56DE\u5BF9\u8C61 ---
  # \u76F4\u63A5\u5207\u7247\uFF0C\u517C\u5BB9 v4/v5
  return(sce[, keep_cells])
}

# =============== \u9884\u5904\u7406 ================
# =============== 2.run_preprocess ================
#' Run Standard Preprocessing on a Seurat Object
#'
#' @description Normalization, Scaling, PCA, Harmony, and UMAP.
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param dims Dimensions for PCA/neighbor graph/UMAP. Default `1:50`.
#' @param group.by Seurat-style alias for `batch_var`.
#' @param batch_var Metadata column used by Harmony. Default `"orig.ident"`.
#' @param n_features Number of variable features to select. Default `3000`.
#' @param seurat_obj Backward-compatible alias for `object`.
#' @export
run_preprocess <- function(object = NULL,
                           dims = 1:50,
                           group.by = NULL,
                           batch_var = "orig.ident",
                           n_features = 3000,
                           seurat_obj = object) {
  if (is.null(seurat_obj)) stop("A Seurat object must be provided via 'object' or 'seurat_obj'.")
  if (!is.null(group.by)) batch_var <- group.by
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = n_features)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- harmony::RunHarmony(seurat_obj, group.by.vars = batch_var, verbose = FALSE)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = dims, reduction = "harmony", verbose = FALSE)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = dims, reduction = "harmony", verbose = FALSE)
  base::return(seurat_obj)
}

# ========= \u53CC\u7EC6\u80DE\u53BB\u9664 =========
# ========= 3.Doublet Detection (DoubletFinder) =========

#' Run DoubletFinder Analysis (Pure Multi-threaded)
#'
#' @description A streamlined, multi-threaded implementation for Seurat V4/V5.
#' No fallback logic. If cores fail, it fails.
#'
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param split_by Metadata column to split samples (e.g. "orig.ident").
#' @param split.by Seurat-style alias for `split_by`.
#' @param pcs PCs to use (default 1:15).
#' @param ncores Number of cores for paramSweep.
#' @param sct Logical, whether to use SCTransform.
#' @param sce Backward-compatible alias for `object`.
#' @export
runDoubletFinderAnalysis <- function(object = NULL,
                                     split_by = "orig.ident",
                                     split.by = NULL,
                                     pcs = 1:15,
                                     ncores = 4,
                                     sct = FALSE,
                                     sce = object) {
  if (is.null(sce)) stop("A Seurat object must be provided via 'object' or 'sce'.")
  if (!is.null(split.by)) split_by <- split.by

  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) stop("Package 'DoubletFinder' is required.")

  message(paste0("Splitting object by ", split_by, "..."))
  sce_list <- Seurat::SplitObject(sce, split.by = split_by)

  # \u904D\u5386\u6BCF\u4E2A\u6837\u672C
  results_list <- lapply(names(sce_list), function(x) {
    sub_sce <- sce_list[[x]]
    n_cells <- ncol(sub_sce)

    # 1. \u8FC7\u6EE4\u6781\u5C0F\u6837\u672C
    if (n_cells < 50) {
      warning(paste("Skipping", x, ": too few cells."))
      return(NULL)
    }

    message(paste0("Processing: ", x, " (", n_cells, " cells) | Cores: ", ncores))

    # 2. \u5FC5\u987B\u5728\u5B50\u96C6\u4E0A\u91CD\u8DD1 PCA (DoubletFinder \u4F9D\u8D56\u9879)
    # \u4E0D\u7528 tryCatch\uFF0C\u76F4\u63A5\u8DD1\uFF0C\u51FA\u9519\u5373\u505C\u6B62
    if (sct) {
      sub_sce <- Seurat::SCTransform(sub_sce, verbose = FALSE)
      sub_sce <- Seurat::RunPCA(sub_sce, verbose = FALSE)
    } else {
      sub_sce <- Seurat::NormalizeData(sub_sce, verbose = FALSE)
      sub_sce <- Seurat::FindVariableFeatures(sub_sce, verbose = FALSE)
      sub_sce <- Seurat::ScaleData(sub_sce, verbose = FALSE)
      sub_sce <- Seurat::RunPCA(sub_sce, verbose = FALSE)
    }

    # 3. \u8FD0\u884C ParamSweep (\u76F4\u63A5\u5E76\u884C)
    # \u4F7F\u7528 _v3 \u7248\u672C\u4EE5\u652F\u6301 num.cores \u548C sct \u53C2\u6570
    sweep.res <- DoubletFinder::paramSweep_v3(sub_sce, PCs = pcs, sct = sct, num.cores = ncores)

    # 4. \u8BA1\u7B97 pK
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)
    pK_val <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # 5. \u8FD0\u884C\u4E3B\u7A0B\u5E8F
    est_rate <- min(n_cells * 0.000008, 0.15) # 10x \u4F30\u7B97\u516C\u5F0F
    nExp <- round(est_rate * n_cells)

    sub_sce <- DoubletFinder::doubletFinder_v3(sub_sce, PCs = pcs, pN = 0.25, pK = pK_val,
                                               nExp = nExp, reuse.pANN = FALSE, sct = sct)

    # 6. \u63D0\u53D6\u7ED3\u679C
    meta <- sub_sce@meta.data
    pann_col <- grep("^pANN_", colnames(meta), value = TRUE)
    class_col <- grep("^DF.classifications_", colnames(meta), value = TRUE)

    data.frame(
      row.names = rownames(meta),
      Doublet_Score = meta[[utils::tail(pann_col, 1)]],
      Doublet_Class = meta[[utils::tail(class_col, 1)]]
    )
  })

  # \u5408\u5E76\u7ED3\u679C
  message("Merging results...")
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) stop("No samples processed successfully.")

  all_res <- do.call(rbind, results_list)
  sce <- Seurat::AddMetaData(sce, metadata = all_res)

  message("Done.")
  return(sce)
}

# ========= \u53BB\u9664RNA\u6C61\u67D3 =========
# ========= 4.decontX =========
#' Run decontX to estimate ambient RNA contamination
#'
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param assay The assay to use for decontamination. Default is "RNA".
#' @param seed Random seed for reproducibility. Default is 123.
#' @param sce Backward-compatible alias for `object`.
#'
#' @return A Seurat object with 'contamination' added to meta.data.
#' @export
#'
#' @importFrom Seurat GetAssayData AddMetaData
#' @importFrom methods is
run_decontX <- function(object = NULL, assay = "RNA", seed = 123, sce = object) {
  if (is.null(sce)) stop("A Seurat object must be provided via 'object' or 'sce'.")
  # 1. \u68C0\u67E5\u4F9D\u8D56\u5305\u662F\u5426\u5B89\u88C5 (celda \u662F Bioconductor \u5305\uFF0C\u901A\u5E38\u653E\u5728 Suggests \u91CC)
  if (!requireNamespace("celda", quietly = TRUE)) {
    stop("Package 'celda' is required for this function. Please install it using: BiocManager::install('celda')")
  }

  message("Extracting counts matrix and running decontX...")

  # 2. \u8BBE\u7F6E\u968F\u673A\u79CD\u5B50
  set.seed(seed)

  # 3. \u83B7\u53D6 Count \u77E9\u9635 (\u517C\u5BB9 Seurat V4/V5)
  # \u6CE8\u610F\uFF1AdecontX \u9700\u8981 raw counts
  counts_matrix <- Seurat::GetAssayData(sce, slot = "counts", assay = assay)

  # 4. \u8FD0\u884C decontX
  # \u4F7F\u7528 tryCatch \u9632\u6B62\u8BA1\u7B97\u8FC7\u7A0B\u4E2D\u51FA\u9519\u5BFC\u81F4\u6574\u4E2A\u6D41\u7A0B\u5D29\u6E83
  tryCatch({
    decontX_res <- celda::decontX(counts_matrix)

    # 5. \u5C06\u6C61\u67D3\u8BC4\u5206\u6DFB\u52A0\u5230 Seurat \u5143\u6570\u636E
    sce <- Seurat::AddMetaData(sce, metadata = decontX_res$contamination, col.name = "contamination")

    # \u53EF\u9009\uFF1A\u4E5F\u53EF\u4EE5\u628A\u53BB\u6C61\u67D3\u540E\u7684\u77E9\u9635\u52A0\u56DE\u53BB\uFF08\u5982\u679C\u4F60\u9700\u8981\u7528\u53BB\u6C61\u67D3\u7684\u77E9\u9635\u505A\u4E0B\u6E38\u5206\u6790\uFF09
    # decont_counts <- decontX_res$decontXcounts
    # sce[["decontX"]] <- Seurat::CreateAssayObject(counts = decont_counts)

    message("Decontamination score calculated successfully.")
    return(sce)

  }, error = function(e) {
    stop("Error running decontX: ", e$message)
  })
}
