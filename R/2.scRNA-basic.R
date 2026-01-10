# =============== 质控 ================
# =============== 1.QC  ================

#' Filter and Plot scRNA-seq Data
#'
#' @description Filters scRNA-seq data and saves QC violin plots.
#' @param sce A Seurat object.
#' @param minGene Minimum genes. Default 200.
#' @param maxGene Maximum genes. Default 5000.
#' @param pctMT Max mitochondrial %. Default 15.
#' @param maxCounts Max UMI counts. Default 20000.
#' @param species "human" or "mouse".
#' @param pal Numeric ID (1-100) or color vector.
#' @param output_dir Output path.
#' @param width,height PDF dimensions.
#' @return Filtered Seurat object.
#' @importFrom Seurat PercentageFeatureSet VlnPlot NoLegend
#' @importFrom ggplot2 ggsave
#' @export
runScRNAQC <- function(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000,
                       species = 'human', pal = 100, output_dir = "./output_figure/",
                       width = 6, height = 4) {

  if (!base::dir.exists(output_dir)) base::dir.create(output_dir, recursive = TRUE)
  actual_pal <- .get_pal(pal)

  mt_pattern <- base::switch(base::tolower(species), "human" = "^MT-", "mouse" = "^mt-",
                             base::stop("Unsupported species."))

  sce[["percent.mt"]] <- Seurat::PercentageFeatureSet(sce, pattern = mt_pattern)

  .save_qc <- function(obj, stage) {
    for (f in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
      p <- Seurat::VlnPlot(obj, features = f, cols = actual_pal, pt.size = 0) + Seurat::NoLegend()
      ggplot2::ggsave(base::file.path(output_dir, base::paste0(stage, "-", f, ".pdf")),
                      plot = p, width = width, height = height)
    }
  }

  .save_qc(sce, "BeforeQC")
  scRNA_filtered <- base::subset(sce, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                                   percent.mt < pctMT & nCount_RNA < maxCounts)
  .save_qc(scRNA_filtered, "AfterQC")

  base::return(scRNA_filtered)
}


# =============== 预处理 ================
# =============== 2.run_preprocess ================
#' Run Standard Preprocessing on a Seurat Object
#'
#' @description Normalization, Scaling, PCA, Harmony, and UMAP.
#' @param seurat_obj Seurat object.
#' @param dims Dimensions for UMAP/Neighbors. Default 1:50.
#' @param batch_var Column for Harmony.
#' @export
run_preprocess <- function(seurat_obj, dims = 1:50, batch_var = "orig.ident", n_features = 3000) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = n_features)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- harmony::RunHarmony(seurat_obj, group.by.vars = batch_var, verbose = FALSE)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = dims, reduction = "harmony", verbose = FALSE)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = dims, reduction = "harmony", verbose = FALSE)
  base::return(seurat_obj)
}

# ========= 双细胞去除 =========
# ========= 3.Doublet Detection (DoubletFinder) =========

#' Run DoubletFinder Analysis (Pure Multi-threaded)
#'
#' @description A streamlined, multi-threaded implementation for Seurat V4/V5.
#' No fallback logic. If cores fail, it fails.
#'
#' @param sce Seurat object.
#' @param split_by Metadata column to split samples (e.g. "orig.ident").
#' @param pcs PCs to use (default 1:15).
#' @param ncores Number of cores for paramSweep.
#' @param sct Logical, whether to use SCTransform.
#' @export
runDoubletFinderAnalysis <- function(sce,
                                     split_by = "orig.ident",
                                     pcs = 1:15,
                                     ncores = 4,
                                     sct = FALSE) {

  require(Seurat)
  require(DoubletFinder)
  require(dplyr)

  message(paste0("Splitting object by ", split_by, "..."))
  sce_list <- SplitObject(sce, split.by = split_by)

  # 遍历每个样本
  results_list <- lapply(names(sce_list), function(x) {
    sub_sce <- sce_list[[x]]
    n_cells <- ncol(sub_sce)

    # 1. 过滤极小样本
    if (n_cells < 50) {
      warning(paste("Skipping", x, ": too few cells."))
      return(NULL)
    }

    message(paste0("Processing: ", x, " (", n_cells, " cells) | Cores: ", ncores))

    # 2. 必须在子集上重跑 PCA (DoubletFinder 依赖项)
    # 不用 tryCatch，直接跑，出错即停止
    if (sct) {
      sub_sce <- SCTransform(sub_sce, verbose = FALSE)
      sub_sce <- RunPCA(sub_sce, verbose = FALSE)
    } else {
      sub_sce <- NormalizeData(sub_sce, verbose = FALSE)
      sub_sce <- FindVariableFeatures(sub_sce, verbose = FALSE)
      sub_sce <- ScaleData(sub_sce, verbose = FALSE)
      sub_sce <- RunPCA(sub_sce, verbose = FALSE)
    }

    # 3. 运行 ParamSweep (直接并行)
    # 使用 _v3 版本以支持 num.cores 和 sct 参数
    sweep.res <- paramSweep_v3(sub_sce, PCs = pcs, sct = sct, num.cores = ncores)

    # 4. 计算 pK
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_val <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # 5. 运行主程序
    est_rate <- min(n_cells * 0.000008, 0.15) # 10x 估算公式
    nExp <- round(est_rate * n_cells)

    sub_sce <- doubletFinder_v3(sub_sce, PCs = pcs, pN = 0.25, pK = pK_val,
                                nExp = nExp, reuse.pANN = FALSE, sct = sct)

    # 6. 提取结果
    meta <- sub_sce@meta.data
    pann_col <- grep("^pANN_", colnames(meta), value = TRUE)
    class_col <- grep("^DF.classifications_", colnames(meta), value = TRUE)

    data.frame(
      row.names = rownames(meta),
      Doublet_Score = meta[[tail(pann_col, 1)]],
      Doublet_Class = meta[[tail(class_col, 1)]]
    )
  })

  # 合并结果
  message("Merging results...")
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) stop("No samples processed successfully.")

  all_res <- do.call(rbind, results_list)
  sce <- AddMetaData(sce, metadata = all_res)

  message("Done.")
  return(sce)
}

# ========= 去除RNA污染 =========
# ========= 4.decontX =========
#' Run decontX to estimate ambient RNA contamination
#'
#' @param sce A Seurat object.
#' @param assay The assay to use for decontamination. Default is "RNA".
#' @param seed Random seed for reproducibility. Default is 123.
#'
#' @return A Seurat object with 'contamination' added to meta.data.
#' @export
#'
#' @importFrom Seurat GetAssayData AddMetaData
#' @importFrom methods is
run_decontX <- function(sce, assay = "RNA", seed = 123) {
  # 1. 检查依赖包是否安装 (celda 是 Bioconductor 包，通常放在 Suggests 里)
  if (!requireNamespace("celda", quietly = TRUE)) {
    stop("Package 'celda' is required for this function. Please install it using: BiocManager::install('celda')")
  }

  message("Extracting counts matrix and running decontX...")

  # 2. 设置随机种子
  set.seed(seed)

  # 3. 获取 Count 矩阵 (兼容 Seurat V4/V5)
  # 注意：decontX 需要 raw counts
  counts_matrix <- Seurat::GetAssayData(sce, slot = "counts", assay = assay)

  # 4. 运行 decontX
  # 使用 tryCatch 防止计算过程中出错导致整个流程崩溃
  tryCatch({
    decontX_res <- celda::decontX(counts_matrix)

    # 5. 将污染评分添加到 Seurat 元数据
    sce <- Seurat::AddMetaData(sce, metadata = decontX_res$contamination, col.name = "contamination")

    # 可选：也可以把去污染后的矩阵加回去（如果你需要用去污染的矩阵做下游分析）
    # decont_counts <- decontX_res$decontXcounts
    # sce[["decontX"]] <- Seurat::CreateAssayObject(counts = decont_counts)

    message("Decontamination score calculated successfully.")
    return(sce)

  }, error = function(e) {
    stop("Error running decontX: ", e$message)
  })
}
