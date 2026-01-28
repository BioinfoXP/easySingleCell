#' @title Run CopyKAT by Sample
#' @description
#' Runs CopyKAT on a Seurat object sample by sample.
#'
#' @param sce A Seurat object.
#' @param target_class String. The cell type to investigate (e.g., "Epithelial").
#' @param ref_classes Vector of strings. The normal reference cell types (e.g., c("T_cell", "Myeloid")).
#' @param sample_col String. Metadata column for sample IDs. Default: "orig.ident".
#' @param celltype_col String. Metadata column for cell types. Default: "celltype".
#' @param output_dir String. Directory to save results. Default: "./CopyKAT".
#' @param overwrite Logical. If TRUE, re-runs existing samples. Default: FALSE.
#' @param min_cells Integer. Minimum cells required per sample (target + ref). Default: 25.
#' @param n.cores Integer. Number of cores. Default: 20.
#' @param ... Additional arguments passed to copykat() (e.g., KS.cut, ngene.chr).
#'
#' @return Invisible data.frame of the combined results. Main output is the CSV file.
#' @export
#'
#' @examples
#' \dontrun{
#'   # 1. Basic Run
#'   # Results will be saved to "./CopyKAT/AllSamples_CopyKAT_prediction.csv"
#'   runCopyKAT(
#'     sce = sce_full,
#'     target_class = "Epithelial",
#'     ref_classes = c("T_cell", "Macrophage"),
#'     n.cores = 10
#'   )
#'
#'   # 2. Advanced: Manual Annotation (How to use the output)
#'   # Since the function doesn't return a Seurat object, do this afterwards:
#'
#'   # a) Read the generated CSV
#'   res <- read.csv("./CopyKAT/AllSamples_CopyKAT_prediction.csv")
#'
#'   # b) Check match rate (using the auto-cleaned ID)
#'   print(mean(res$cleaned_id %in% colnames(sce_full)))
#'
#'   # c) Add to Seurat (Safe Join)
#'   library(tidyverse)
#'   meta_to_add <- res %>%
#'      distinct(cleaned_id, .keep_all = TRUE) %>%
#'      column_to_rownames("cleaned_id") %>%
#'      select(copykat.pred)
#'
#'   sce_full <- AddMetaData(sce_full, meta_to_add)
#'   DimPlot(sce_full, group.by="copykat.pred")
#' }
runCopyKAT <- function(sce,
                       target_class,
                       ref_classes,
                       sample_col = "orig.ident",
                       celltype_col = "celltype",
                       output_dir = "./CopyKAT",
                       overwrite = FALSE,
                       min_cells = 25,
                       n.cores = 20,
                       KS.cut = 0.1,
                       ngene.chr = 5,
                       win.size = 25,
                       ...) {

  # --- 1. Environment Setup ---
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  if (!"package:copykat" %in% search()) {
    message(">> Loading copykat package...")
    library(copykat)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Safe working directory handling
  original_wd <- getwd()
  on.exit(setwd(original_wd))

  # --- 2. Input Validation ---
  message(">> Validating inputs...")
  if (!all(c(sample_col, celltype_col) %in% colnames(sce@meta.data))) {
    stop(sprintf("Columns '%s' or '%s' not found in metadata.", sample_col, celltype_col))
  }

  # Identify Target and Reference Cells
  message(sprintf(">> Subsetting: Target ['%s'] vs Ref [%s]", target_class, paste(head(ref_classes, 2), collapse=",")))

  target_cells <- rownames(sce@meta.data)[sce@meta.data[[celltype_col]] %in% target_class]
  ref_cells <- rownames(sce@meta.data)[sce@meta.data[[celltype_col]] %in% ref_classes]

  if(length(target_cells) == 0) stop("No target cells found!")
  if(length(ref_cells) == 0) stop("No reference cells found!")

  # Get Sample IDs that actually contain target cells
  sample_ids <- unique(sce@meta.data[target_cells, sample_col])
  prediction_results <- list()

  # --- 3. Main Loop ---
  for (sid in sample_ids) {
    message(sprintf("\n>>> Processing Sample: %s", sid))

    sample_subdir <- file.path(output_dir, sid)
    # Define possible output filenames (CopyKAT acts differently depending on version/inputs)
    expected_file_1 <- file.path(sample_subdir, paste0("S_", sid, "_copykat_prediction.txt"))
    expected_file_2 <- file.path(sample_subdir, paste0(sid, "_copykat_prediction.txt"))

    # 3.1 Check for Existing Results (Resume)
    found_file <- if(file.exists(expected_file_1)) expected_file_1 else if(file.exists(expected_file_2)) expected_file_2 else NULL

    if (!overwrite && !is.null(found_file)) {
      message(sprintf("   [Found] %s. Loading...", basename(found_file)))
      tryCatch({
        pred <- read.table(found_file, header = TRUE, check.names = FALSE, sep = "\t")
        if(!"cell_id" %in% colnames(pred)) pred <- tibble::rownames_to_column(pred, "cell_id")
        pred$sample_id <- sid
        prediction_results[[sid]] <- pred
        next
      }, error = function(e) message("   [Warning] File corrupt. Re-running."))
    }

    # 3.2 Run CopyKAT
    tryCatch({
      # Subset Data for current sample
      curr_cells <- rownames(sce@meta.data)[sce@meta.data[[sample_col]] == sid]
      curr_tgt <- intersect(curr_cells, target_cells)
      curr_ref <- intersect(curr_cells, ref_cells)

      if (length(curr_tgt) < 5 || length(curr_ref) < min_cells) {
        message(sprintf("   [Skip] Not enough cells (Target: %d, Ref: %d).", length(curr_tgt), length(curr_ref)))
        next
      }

      # Create Matrix (Memory Safe)
      mat_tgt <- Seurat::GetAssayData(sce, slot="counts", assay="RNA")[, curr_tgt, drop=FALSE]
      mat_ref <- Seurat::GetAssayData(sce, slot="counts", assay="RNA")[, curr_ref, drop=FALSE]

      genes <- intersect(rownames(mat_tgt), rownames(mat_ref))
      mat_comb <- as.matrix(cbind(mat_tgt[genes,], mat_ref[genes,]))
      mat_comb <- mat_comb[rowSums(mat_comb > 0) >= 5, , drop=FALSE]

      if (!dir.exists(sample_subdir)) dir.create(sample_subdir, recursive = TRUE)
      setwd(sample_subdir)

      res <- copykat(
        rawmat = mat_comb, norm.cell.names = curr_ref,
        sam.name = paste0("S_", sid), id.type = "S",
        output.seg = FALSE, n.cores = n.cores,
        KS.cut = KS.cut, ngene.chr = ngene.chr, win.size = win.size, ...
      )

      if (is.data.frame(res$prediction)) {
        pred <- res$prediction
        pred <- tibble::rownames_to_column(pred, "cell_id")
        pred$sample_id <- sid
        prediction_results[[sid]] <- pred
        message("   [Success] Done.")
      }

      rm(mat_comb, mat_tgt, mat_ref, res); gc(verbose=FALSE)

    }, error = function(e) message(sprintf("   [Error] %s", e$message)))

    setwd(original_wd)
  }

  # --- 4. Export Results ---
  if (length(prediction_results) == 0) {
    warning(">> No results generated.")
    return(NULL)
  }

  message("\n>> Aggregating results...")
  final_df <- dplyr::bind_rows(prediction_results)

  # Auto-Cleaning IDs (Remove prefixes like "S_Sample1." to match Seurat)
  final_df <- final_df %>%
    dplyr::mutate(
      raw_id = cell_id,
      cleaned_id = mapply(function(id, samp) {
        # Pattern: Starts with S_SampleID. OR Just SampleID.
        gsub(paste0("^(S_)?", samp, "\\."), "", id)
      }, raw_id, sample_id)
    ) %>%
    dplyr::select(sample_id, raw_id, cleaned_id, copykat.pred, everything())

  out_file <- file.path(output_dir, "AllSamples_CopyKAT_prediction.csv")
  write.csv(final_df, file = out_file, row.names = FALSE)

  message(sprintf(">> SUCCESS. Results saved to: %s", out_file))
  message(">> Use read.csv() to load this file and AddMetaData() to annotate your object.")

  return(invisible(final_df))
}
