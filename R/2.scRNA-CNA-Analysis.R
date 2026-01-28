#' @title Run CopyKAT by Sample
#' @description
#' Runs CopyKAT on a Seurat object by splitting it into samples.
#' @param sce A Seurat object containing both tumor and normal/immune cells.
#' @param target_class String. The cell type to investigate (e.g., "Epithelial").
#' @param ref_classes Vector of strings. The cell types to serve as normal references (e.g., c("T_cell", "Myeloid")).
#' @param sample_col String. Column name for sample IDs. Default: "orig.ident".
#' @param celltype_col String. Column name for cell types. Default: "celltype".
#' @param output_dir String. Directory to save results. Default: "./CopyKAT".
#' @param overwrite Logical. If TRUE, re-runs CopyKAT even if output files exist. Default: FALSE.
#' @param min_cells Integer. Minimum number of cells required per sample. Default: 25.
#' @param n.cores Integer. Number of cores. Default: 20.
#' @param KS.cut Numeric. Segmentation sensitivity (0-1). Default: 0.1.
#' @param ngene.chr Integer. Min genes per chromosome. Default: 5.
#' @param win.size Integer. Window size. Default: 25.
#' @param ... Additional arguments passed directly to copykat().
#'
#' @return A Seurat object with added 'copykat.pred' metadata.
#' @export
#'
#' @examples
#' \dontrun{
#'   # 1. Basic Run (Resume mode by default)
#'   # It will skip samples that already have output files in ./CopyKAT
#'   sce_res <- runCopyKAT(
#'     sce = sce_full,
#'     target_class = "Epithelial",
#'     ref_classes = c("T_cell", "Myeloid")
#'   )
#'
#'   # 2. Force Re-run (Overwrite existing results)
#'   # Useful if you changed parameters (e.g., KS.cut) and need to update results
#'   sce_res <- runCopyKAT(
#'     sce = sce_full,
#'     target_class = "Epithelial",
#'     ref_classes = c("T_cell"),
#'     overwrite = TRUE,
#'     n.cores = 30
#'   )
#'
#'   # 3. Advanced Tuning
#'   sce_res <- runCopyKAT(
#'     sce = sce_full,
#'     target_class = "Epithelial",
#'     ref_classes = c("T_cell", "Macrophage"),
#'     KS.cut = 0.15,      # Less sensitive segmentation
#'     win.size = 20,      # Smaller window size
#'     output_dir = "./Analysis/CopyKAT_V2"
#'   )
#'
#'   # 4. Visualization
#'   Seurat::DimPlot(sce_res, group.by = "copykat.pred")
#' }
runCopyKAT <- function(sce,
                       target_class,
                       ref_classes,
                       sample_col = "orig.ident",
                       celltype_col = "celltype",
                       output_dir = "./CopyKAT",
                       overwrite = FALSE,
                       min_cells = 25,
                       # --- Explicit CopyKAT Parameters ---
                       n.cores = 20,
                       KS.cut = 0.1,
                       ngene.chr = 5,
                       win.size = 25,
                       ...) {

  # 1. Dependency Check
  # CopyKAT MUST be attached to access internal data (full.anno).
  if (!"package:copykat" %in% search()) {
    message(">> Loading copykat package...")
    library(copykat)
  }
  if (getRversion() < "4.1.0") stop("Requires R >= 4.1.0 for |> pipe.")

  required_pkgs <- c("Seurat", "dplyr", "tibble")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) stop("Missing: ", paste(missing_pkgs, collapse = ", "))

  # 2. Setup
  message(">> Checking inputs...")
  options(future.globals.maxSize = 100 * 1024^3)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (!all(c(sample_col, celltype_col) %in% colnames(sce@meta.data))) {
    stop("Column not found in metadata.")
  }

  # 3. Subset
  message(sprintf(">> Subsetting: Target ['%s'] vs Ref [%s]", target_class, paste(head(ref_classes, 2), collapse=",")))
  sce_target <- sce[, sce[[celltype_col, drop=TRUE]] %in% target_class]
  sce_ref <- sce[, sce[[celltype_col, drop=TRUE]] %in% ref_classes]

  sample_ids <- unique(sce_target[[sample_col, drop=TRUE]])
  prediction_results <- list()

  original_wd <- getwd()
  on.exit(setwd(original_wd))

  # 4. Loop
  for (sid in sample_ids) {
    message(sprintf("\n>>> Processing Sample: %s", sid))

    # --- Resume Logic ---
    sample_subdir <- file.path(output_dir, sid)
    # Check for possible file patterns (some versions add S_, some don't)
    expected_file_1 <- file.path(sample_subdir, paste0("S_", sid, "_copykat_prediction.txt"))
    expected_file_2 <- file.path(sample_subdir, paste0(sid, "_copykat_prediction.txt")) # Fallback

    target_file <- if(file.exists(expected_file_1)) expected_file_1 else if(file.exists(expected_file_2)) expected_file_2 else NULL

    if (!overwrite && !is.null(target_file)) {
      message(sprintf("   [Found] %s. Skipping computation.", basename(target_file)))
      tryCatch({
        pred <- read.table(target_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
        pred$sample_id <- sid
        prediction_results[[sid]] <- pred
        next
      }, error = function(e) message("   [Warning] Read failed, re-running."))
    }

    # --- Run Logic ---
    tryCatch({
      cells_tgt <- colnames(sce_target)[sce_target[[sample_col]] == sid]
      cells_ref <- colnames(sce_ref)[sce_ref[[sample_col]] == sid]

      if (length(cells_tgt) < min_cells || length(cells_ref) < min_cells) {
        message("   [Skip] Not enough cells.")
        next
      }

      mat_tgt <- Seurat::GetAssayData(sce_target, slot="counts", assay="RNA")[, cells_tgt, drop=FALSE]
      mat_ref <- Seurat::GetAssayData(sce_ref, slot="counts", assay="RNA")[, cells_ref, drop=FALSE]

      genes <- intersect(rownames(mat_tgt), rownames(mat_ref))
      mat_comb <- as.matrix(cbind(mat_tgt[genes,], mat_ref[genes,]))
      mat_comb <- mat_comb[rowSums(mat_comb > 0) >= 5, , drop=FALSE]

      if (!dir.exists(sample_subdir)) dir.create(sample_subdir)
      setwd(sample_subdir)

      res <- copykat(
        rawmat = mat_comb, norm.cell.names = cells_ref,
        sam.name = paste0("S_", sid), id.type = "S", output.seg = FALSE,
        n.cores = n.cores, KS.cut = KS.cut, ngene.chr = ngene.chr, win.size = win.size, ...
      )

      if (is.data.frame(res$prediction)) {
        pred <- res$prediction
        pred$sample_id <- sid
        prediction_results[[sid]] <- pred
      }
      rm(mat_comb, mat_tgt, mat_ref, res); gc(verbose=F)

    }, error = function(e) message(sprintf("   [Error] %s", e$message)))
    setwd(original_wd)
  }

  # 5. Merge & Auto-Fix IDs
  if (length(prediction_results) == 0) {
    warning(">> No results generated.")
    return(sce_target)
  }

  message("\n>> Merging results...")
  final_pred_df <- do.call(rbind, prediction_results)

  # =======================================================
  # CRITICAL FIX: CLEAN PREFIXES
  # =======================================================
  # Problem: CopyKAT adds "SampleID." prefix (e.g., "GSM123.GSE_GSM123_Cell1")
  # Target: Seurat has "GSE_GSM123_Cell1"

  message(">> Checking and cleaning Cell IDs...")
  raw_ids <- rownames(final_pred_df)
  target_ids <- colnames(sce_target)

  # 1. First Pass: Try Exact Match
  match_rate_raw <- sum(raw_ids %in% target_ids) / length(raw_ids)

  if (match_rate_raw < 0.1) {
    message("   Low match rate detected. Attempting to strip 'SampleID.' prefixes...")

    clean_ids <- raw_ids

    # Iterate through each sample to safely remove its specific prefix
    if (!is.null(final_pred_df$sample_id)) {
      for (sid in unique(final_pred_df$sample_id)) {
        # Define Pattern: Start of string (^), followed by SampleID, followed by dot (\.)
        # Example pattern: ^GSM3516662\.
        prefix_pattern <- paste0("^", sid, "\\.")

        # Apply replacement
        clean_ids <- sub(prefix_pattern, "", clean_ids)
      }
    }

    # Check match rate after cleaning
    match_rate_clean <- sum(clean_ids %in% target_ids) / length(clean_ids)
    message(sprintf("   Match rate improved: %.1f%% -> %.1f%%", match_rate_raw*100, match_rate_clean*100))

    if (match_rate_clean > match_rate_raw) {
      rownames(final_pred_df) <- clean_ids
    }
  }
  # =======================================================

  write.csv(final_pred_df, file = file.path(output_dir, "AllSamples_CopyKAT_prediction.csv"))

  # Annotate
  meta_df <- final_pred_df |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::filter(cell_id %in% colnames(sce_target)) |>
    dplyr::select(cell_id, dplyr::any_of(c("copykat.pred", "copykat.test"))) |>
    dplyr::distinct(cell_id, .keep_all = TRUE) |>
    tibble::column_to_rownames("cell_id")

  sce_target <- Seurat::AddMetaData(sce_target, metadata = meta_df)

  if ("copykat.pred" %in% colnames(sce_target@meta.data)) {
    sce_target@meta.data$copykat.pred[is.na(sce_target@meta.data$copykat.pred)] <- "Not_Run"
    print(table(sce_target@meta.data$copykat.pred))
  }

  message(">> Done.")
  return(sce_target)
}
