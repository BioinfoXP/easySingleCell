#' Add SCENIC AUC to Seurat with Beautified Logging
#'
#' @description
#' Integrates SCENIC AUC scores into a Seurat object.
#' Features robust format detection, dynamic parameter adjustment, and professional logging.
#'
#' @param object A Seurat object. This Seurat-style argument is preferred.
#' @param auc_df A data.frame or matrix of AUC scores.
#' @param reduction.name Name of the new reduction (default: "umap_ras").
#' @param assay.name Name of the new assay (default: "RAS").
#' @param dims Integer vector for UMAP dimensions (default: 1:30).
#' @param n.neighbors Integer for UMAP neighbors (default: 30).
#' @param sce Backward-compatible alias for `object`.
#' @param reduction Seurat-style alias for `reduction.name`.
#' @param assay Seurat-style alias for `assay.name`.
#'
#' @return A Seurat object with "RAS" assay and UMAP reduction.
#' @export
add_Umap_RAS <- function(object = NULL,
                         auc_df = NULL, 
                         reduction.name = "umap_ras",
                         assay.name = "RAS",
                         dims = 1:30,
                         n.neighbors = 30L,
                         sce = object,
                         reduction = NULL,
                         assay = NULL) {
  if (is.null(sce)) stop("A Seurat object must be provided via 'object' or 'sce'.")
  if (is.null(auc_df)) stop("'auc_df' is required.")
  if (!is.null(reduction)) reduction.name <- reduction
  if (!is.null(assay)) assay.name <- assay
  
  # -- \u65E5\u5FD7\u8F85\u52A9\u51FD\u6570 --
  log_header <- function(msg) base::message(base::sprintf("\n\u2500\u2500 %s %s", msg, base::paste(base::rep("\u2500", 40), collapse = "")))
  log_info   <- function(msg) base::message(base::sprintf(" \u2139 %s", msg))
  log_ok     <- function(msg) base::message(base::sprintf(" \u2714 %s", msg))
  log_step   <- function(msg) base::message(base::sprintf(" \u27A4 %s", msg))
  
  log_header("SCENIC Integration Pipeline")
  
  # -- 1. \u683C\u5F0F\u8F6C\u6362\u4E0E\u65B9\u5411\u68C0\u6D4B --
  auc_mtx <- base::as.matrix(auc_df)
  
  sce_cells <- base::colnames(sce)
  n_total_cells <- base::length(sce_cells)
  
  # \u68C0\u6D4B\u4EA4\u96C6
  cols_in_sce <- base::sum(base::colnames(auc_mtx) %in% sce_cells)
  rows_in_sce <- base::sum(base::rownames(auc_mtx) %in% sce_cells)
  
  if (cols_in_sce > rows_in_sce) {
    # \u683C\u5F0F: [Regulons x Cells] (\u65E0\u9700\u8F6C\u7F6E)
    log_info(base::sprintf("Format detected: [Regulons x Cells] (%d cols match)", cols_in_sce))
  } else if (rows_in_sce > 0) {
    # \u683C\u5F0F: [Cells x Regulons] (\u9700\u8981\u8F6C\u7F6E)
    log_info(base::sprintf("Format detected: [Cells x Regulons] (%d rows match)", rows_in_sce))
    log_step("Transposing matrix to match Seurat requirements...")
    auc_mtx <- base::t(auc_mtx)
  } else {
    base::stop("\u274C Error: Cannot match cell IDs! Check row/col names.")
  }
  
  # -- 2. \u6E05\u6D17 Regulon \u540D\u79F0 --
  if (base::any(base::grepl("\\.\\.\\.", base::rownames(auc_mtx)))) {
    n_fixed <- base::sum(base::grepl("\\.\\.\\.", base::rownames(auc_mtx)))
    base::rownames(auc_mtx) <- base::gsub("\\.\\.\\.", "(+)", base::rownames(auc_mtx))
    log_ok(base::sprintf("Fixed %d regulon names (replaced '...' with '(+)')", n_fixed))
  }
  
  # -- 3. \u5339\u914D\u7EC6\u80DE --
  common_cells <- base::intersect(base::colnames(sce), base::colnames(auc_mtx))
  n_common <- base::length(common_cells)
  
  if (n_common < 2) base::stop("\u274C Error: Too few common cells (<2).")
  
  # \u53D6\u5B50\u96C6
  sce_sub <- sce[, common_cells]
  auc_sub <- auc_mtx[, common_cells]
  
  log_ok(base::sprintf("Aligned %d common cells (discarded %d)", n_common, n_total_cells - n_common))
  
  # -- 4. \u6DFB\u52A0 Assay (\u8BE6\u7EC6\u65E5\u5FD7) --
  new_assay <- Seurat::CreateAssayObject(data = auc_sub)
  sce_sub[[assay.name]] <- new_assay
  Seurat::DefaultAssay(sce_sub) <- assay.name
  
  log_step(base::sprintf("Creating Assay '%s':", assay.name))
  base::message(base::sprintf("   \u25CF Features (Regulons): %d", base::nrow(auc_sub)))
  base::message(base::sprintf("   \u25CF Cells: %d", base::ncol(auc_sub)))
  base::message(base::sprintf("   \u25CF Data Range: %.2f - %.2f", base::min(auc_sub), base::max(auc_sub)))
  
  # -- 5. \u964D\u7EF4\u5206\u6790 --
  # \u52A8\u6001\u8BA1\u7B97 PC
  n_features <- base::nrow(auc_sub)
  n_pcs_max <- base::min(n_features - 1, 50)
  
  if (n_pcs_max < 2) base::stop("\u274C Error: Not enough regulons for PCA.")
  
  log_step("Running Dimensionality Reduction:")
  
  # Scale & PCA
  sce_sub <- Seurat::ScaleData(sce_sub, features = base::rownames(sce_sub), verbose = FALSE)
  sce_sub <- Seurat::RunPCA(sce_sub, features = base::rownames(sce_sub), npcs = n_pcs_max, verbose = FALSE)
  base::message(base::sprintf("   \u25CF PCA: Computed top %d PCs (limited by %d features)", n_pcs_max, n_features))
  
  # UMAP
  valid_dims <- dims[dims <= n_pcs_max]
  if (length(valid_dims) == 0) valid_dims <- 1:n_pcs_max
  
  # \u52A8\u6001\u8C03\u6574 neighbors
  actual_neighbors <- base::min(n.neighbors, base::ncol(sce_sub) - 1)
  
  sce_sub <- Seurat::RunUMAP(sce_sub, 
                             dims = valid_dims, 
                             reduction.name = reduction.name, 
                             reduction.key = paste0(assay.name, "UMAP_"), 
                             n.neighbors = actual_neighbors,
                             verbose = FALSE)
  
  base::message(base::sprintf("   \u25CF UMAP: Using dims 1:%d | neighbors=%d", base::max(valid_dims), actual_neighbors))
  
  log_header("Completed")
  base::message(base::sprintf("\u2705 Visualize with: DimPlot(object, reduction = '%s')", reduction.name))
  base::message("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  
  return(sce_sub)
}
