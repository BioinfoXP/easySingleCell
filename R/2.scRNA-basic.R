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

#' Run DoubletFinder Analysis
#'
#' @description Automatically detects doublets in scRNA-seq data. It splits the object by sample,
#' calculates the optimal pK for each sample, predicts doublets, standardizes the metadata column names,
#' and merges the object back together.
#'
#' @param sce A Seurat object.
#' @param split_by Metadata column to split samples by. Default is "orig.ident".
#' @param pcs Numeric vector of PCs to use. Default is 1:15.
#' @param ncores Number of cores for parallel processing (paramSweep). Default is 4.
#' @param double_rate_mode Strategy for doublet rate:
#' \itemize{
#'   \item "auto": Estimates rate based on cell count (linear fit to 10x Genomics standards: ~0.8% per 1000 cells).
#'   \item Numeric value (e.g., 0.05): Applies a fixed rate to all samples.
#' }
#' @param sct Logical, whether sctransform was used. Default is FALSE.
#' @return A Seurat object with 'Doublet_Score' and 'Doublet_Class' added to metadata.
#' @importFrom Seurat SplitObject Idents
#' @importFrom DoubletFinder paramSweep summarizeSweep find.pK modelHomotypic doubletFinder
#' @importFrom stats approxfun
#' @importFrom dplyr %>%
#' @export
runDoubletFinderAnalysis <- function(sce,
                                     split_by = "orig.ident",
                                     pcs = 1:15,
                                     ncores = 4,
                                     double_rate_mode = "auto",
                                     sct = FALSE) {

  # 1. Input Validation
  if (!base::inherits(sce, "Seurat")) base::stop("sce must be a Seurat object")
  if (!split_by %in% base::colnames(sce@meta.data)) base::stop(base::paste0("Column ", split_by, " not found."))

  # Ensure Idents are set to clusters for homotypic calculation
  # If seurat_clusters doesn't exist, we fallback to a placeholder or stop
  if (!"seurat_clusters" %in% base::colnames(sce@meta.data)) {
    base::warning("'seurat_clusters' not found. Homotypic doublet adjustment will assume random distribution (prop=0).")
    has_clusters <- FALSE
  } else {
    has_clusters <- TRUE
  }

  base::message("Splitting object by ", split_by, "...")
  sce_list <- Seurat::SplitObject(sce, split.by = split_by)

  # 2. Process each sample
  processed_list <- base::lapply(base::names(sce_list), function(sample_name) {
    sub_sce <- sce_list[[sample_name]]
    n_cells <- base::ncol(sub_sce)

    # Skip very small samples
    if (n_cells < 50) {
      base::warning(base::paste0("Sample ", sample_name, " has too few cells (<50). Skipping DoubletFinder."))
      sub_sce$Doublet_Score <- NA
      sub_sce$Doublet_Class <- "Uncertain"
      return(sub_sce)
    }

    base::message(base::paste0("Processing sample: ", sample_name, " (", n_cells, " cells)"))

    # 2.1 Determine Doublet Rate
    if (base::is.numeric(double_rate_mode)) {
      est_rate <- double_rate_mode
    } else if (double_rate_mode == "auto") {
      # 10x Genomics standard: ~0.8% per 1000 cells
      # Rate ~= n_cells * 8e-6
      est_rate <- n_cells * 0.000008
      # Cap reasonably (e.g., max 10% for extremely huge datasets unless specified)
      if (est_rate > 0.15) est_rate <- 0.15
    } else {
      est_rate <- 0.075 # Default fallback
    }

    # 2.2 Parameter Sweep (finding optimal pK)
    # capture output to reduce noise
    sweep_res <- DoubletFinder::paramSweep(sub_sce, PCs = pcs, sct = sct, num.cores = ncores)
    sweep_stats <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep_stats)

    # Extract optimal pK
    pK_bcmvn <- bcmvn$pK[base::which.max(bcmvn$BCmetric)] %>% base::as.character() %>% base::as.numeric()

    # 2.3 Calculate nExp (Expected doublets)
    if (has_clusters) {
      annotations <- sub_sce@meta.data$seurat_clusters
      homotypic_prop <- DoubletFinder::modelHomotypic(annotations)
    } else {
      homotypic_prop <- 0
    }

    nExp_poi <- base::round(est_rate * n_cells)
    nExp_adj <- base::round(nExp_poi * (1 - homotypic_prop))

    # 2.4 Run DoubletFinder
    # Force reuse.pANN to FALSE to ensure clean run
    sub_sce <- DoubletFinder::doubletFinder(
      sub_sce,
      PCs = pcs,
      pN = 0.25,
      pK = pK_bcmvn,
      nExp = nExp_adj,
      reuse.pANN = FALSE,
      sct = sct
    )

    # 2.5 Standardize Metadata Column Names
    # DoubletFinder produces columns like "pANN_0.25_0.09_913" and "DF.classifications_0.25_0.09_913"
    meta_cols <- base::colnames(sub_sce@meta.data)

    # Find the pANN column
    pann_col <- meta_cols[base::grepl("^pANN_", meta_cols)]
    # Find the classification column
    class_col <- meta_cols[base::grepl("^DF.classifications_", meta_cols)]

    # Rename to fixed names for clean merging later
    if (base::length(pann_col) > 0) {
      sub_sce$Doublet_Score <- sub_sce@meta.data[[pann_col[1]]]
      sub_sce[[pann_col[1]]] <- NULL # Remove original messy column
    }

    if (base::length(class_col) > 0) {
      sub_sce$Doublet_Class <- sub_sce@meta.data[[class_col[1]]]
      sub_sce[[class_col[1]]] <- NULL # Remove original messy column
    }

    return(sub_sce)
  })

  # 3. Merge back into one object
  base::message("Merging samples back...")

  if (base::length(processed_list) == 1) {
    sce_final <- processed_list[[1]]
  } else {
    sce_final <- base::merge(processed_list[[1]], y = processed_list[2:base::length(processed_list)])
    # Re-join layers if using Seurat V5 to keep it clean (optional but recommended)
    if (base::grepl("^5", base::packageVersion("Seurat"))) {
      sce_final <- Seurat::JoinLayers(sce_final)
    }
  }

  # Add parameters to misc for reproducibility
  sce_final@misc$doublet_params <- base::list(
    mode = double_rate_mode,
    PCs = pcs,
    split_by = split_by
  )

  base::return(sce_final)
}

# ================= Example Usage =================

# \dontrun{
#   # 假设 filtered_sce 是经过预处理和初步聚类（拥有 seurat_clusters）的对象
#
#   # 模式 1: 自动估算比率（最推荐）
#   # 自动根据细胞数计算：1000个细胞~0.8%双细胞率
#   sce_doublet <- runDoubletFinderAnalysis(
#     sce = filtered_sce,
#     split_by = "orig.ident",
#     ncores = 8,
#     double_rate_mode = "auto"
#   )
#
#   # 模式 2: 手动指定固定比率
#   sce_doublet_manual <- runDoubletFinderAnalysis(
#     sce = filtered_sce,
#     double_rate_mode = 0.05  # 强制 5% 双细胞率
#   )
#
#   # 可视化结果
#   Seurat::DimPlot(sce_doublet, group.by = "Doublet_Class", cols = c("red", "grey"))
# }

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
