#' Add SCENIC AUC to Seurat with Beautified Logging
#'
#' @description
#' Integrates SCENIC AUC scores into a Seurat object.
#' Features robust format detection, dynamic parameter adjustment, and professional logging.
#'
#' @param sce A Seurat object.
#' @param auc_df A data.frame or matrix of AUC scores.
#' @param reduction.name Name of the new reduction (default: "umap_ras").
#' @param assay.name Name of the new assay (default: "RAS").
#' @param dims Integer vector for UMAP dimensions (default: 1:30).
#' @param n.neighbors Integer for UMAP neighbors (default: 30).
#'
#' @return A Seurat object with "RAS" assay and UMAP reduction.
#' @export
add_Umap_RAS <- function(sce, 
                         auc_df, 
                         reduction.name = "umap_ras",
                         assay.name = "RAS",
                         dims = 1:30,
                         n.neighbors = 30L) {
  
  # -- 日志辅助函数 --
  log_header <- function(msg) base::message(base::sprintf("\n── %s %s", msg, base::paste(base::rep("─", 40), collapse = "")))
  log_info   <- function(msg) base::message(base::sprintf(" ℹ %s", msg))
  log_ok     <- function(msg) base::message(base::sprintf(" ✔ %s", msg))
  log_step   <- function(msg) base::message(base::sprintf(" ➤ %s", msg))
  
  log_header("SCENIC Integration Pipeline")
  
  # -- 1. 格式转换与方向检测 --
  auc_mtx <- base::as.matrix(auc_df)
  
  sce_cells <- base::colnames(sce)
  n_total_cells <- base::length(sce_cells)
  
  # 检测交集
  cols_in_sce <- base::sum(base::colnames(auc_mtx) %in% sce_cells)
  rows_in_sce <- base::sum(base::rownames(auc_mtx) %in% sce_cells)
  
  if (cols_in_sce > rows_in_sce) {
    # 格式: [Regulons x Cells] (无需转置)
    log_info(base::sprintf("Format detected: [Regulons x Cells] (%d cols match)", cols_in_sce))
  } else if (rows_in_sce > 0) {
    # 格式: [Cells x Regulons] (需要转置)
    log_info(base::sprintf("Format detected: [Cells x Regulons] (%d rows match)", rows_in_sce))
    log_step("Transposing matrix to match Seurat requirements...")
    auc_mtx <- base::t(auc_mtx)
  } else {
    base::stop("❌ Error: Cannot match cell IDs! Check row/col names.")
  }
  
  # -- 2. 清洗 Regulon 名称 --
  if (base::any(base::grepl("\\.\\.\\.", base::rownames(auc_mtx)))) {
    n_fixed <- base::sum(base::grepl("\\.\\.\\.", base::rownames(auc_mtx)))
    base::rownames(auc_mtx) <- base::gsub("\\.\\.\\.", "(+)", base::rownames(auc_mtx))
    log_ok(base::sprintf("Fixed %d regulon names (replaced '...' with '(+)')", n_fixed))
  }
  
  # -- 3. 匹配细胞 --
  common_cells <- base::intersect(base::colnames(sce), base::colnames(auc_mtx))
  n_common <- base::length(common_cells)
  
  if (n_common < 2) base::stop("❌ Error: Too few common cells (<2).")
  
  # 取子集
  sce_sub <- sce[, common_cells]
  auc_sub <- auc_mtx[, common_cells]
  
  log_ok(base::sprintf("Aligned %d common cells (discarded %d)", n_common, n_total_cells - n_common))
  
  # -- 4. 添加 Assay (详细日志) --
  new_assay <- Seurat::CreateAssayObject(data = auc_sub)
  sce_sub[[assay.name]] <- new_assay
  Seurat::DefaultAssay(sce_sub) <- assay.name
  
  log_step(base::sprintf("Creating Assay '%s':", assay.name))
  base::message(base::sprintf("   ● Features (Regulons): %d", base::nrow(auc_sub)))
  base::message(base::sprintf("   ● Cells: %d", base::ncol(auc_sub)))
  base::message(base::sprintf("   ● Data Range: %.2f - %.2f", base::min(auc_sub), base::max(auc_sub)))
  
  # -- 5. 降维分析 --
  # 动态计算 PC
  n_features <- base::nrow(auc_sub)
  n_pcs_max <- base::min(n_features - 1, 50)
  
  if (n_pcs_max < 2) base::stop("❌ Error: Not enough regulons for PCA.")
  
  log_step("Running Dimensionality Reduction:")
  
  # Scale & PCA
  sce_sub <- Seurat::ScaleData(sce_sub, features = base::rownames(sce_sub), verbose = FALSE)
  sce_sub <- Seurat::RunPCA(sce_sub, features = base::rownames(sce_sub), npcs = n_pcs_max, verbose = FALSE)
  base::message(base::sprintf("   ● PCA: Computed top %d PCs (limited by %d features)", n_pcs_max, n_features))
  
  # UMAP
  valid_dims <- dims[dims <= n_pcs_max]
  if (length(valid_dims) == 0) valid_dims <- 1:n_pcs_max
  
  # 动态调整 neighbors
  actual_neighbors <- base::min(n.neighbors, base::ncol(sce_sub) - 1)
  
  sce_sub <- Seurat::RunUMAP(sce_sub, 
                             dims = valid_dims, 
                             reduction.name = reduction.name, 
                             reduction.key = paste0(assay.name, "UMAP_"), 
                             n.neighbors = actual_neighbors,
                             verbose = FALSE)
  
  base::message(base::sprintf("   ● UMAP: Using dims 1:%d | neighbors=%d", base::max(valid_dims), actual_neighbors))
  
  log_header("Completed")
  base::message(base::sprintf("✅ Visualize with: DimPlot(object, reduction = '%s')", reduction.name))
  base::message("──────────────────────────────────────────\n")
  
  return(sce_sub)
}
