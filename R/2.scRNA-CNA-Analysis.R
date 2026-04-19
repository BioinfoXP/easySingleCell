#' @title Run CopyKAT by Sample or with Shared Reference
#' @description
#' Runs CopyKAT on a Seurat object sample by sample.
#' Supports either sample-specific references or a shared global reference pool.
#'
#' @param sce A Seurat object.
#' @param target_class String. The cell type to investigate.
#' @param ref_classes Vector of strings. The normal reference cell types.
#' @param sample_col String. Metadata column for sample IDs. Default: "orig.ident".
#' @param celltype_col String. Metadata column for cell types. Default: "celltype".
#' @param assay String. Assay to use. Default: "RNA".
#' @param output_dir String. Directory to save results. Default: "./CopyKAT".
#' @param overwrite Logical. If TRUE, re-runs existing samples. Default: FALSE.
#' @param ref_mode String. "sample" for per-sample references; "shared" for global shared reference. Default: "sample".
#' @param min_cells Integer. Minimum reference cells required in sample mode. Default: 25.
#' @param min_target_cells Integer. Minimum target cells required per sample. Default: 5.
#' @param min_ref_cells Integer. Minimum reference cells required. In shared mode, refers to global ref size. Default: 25.
#' @param max_ref_cells Integer or NULL. In shared mode, optionally downsample global refs to this size. Default: 3000.
#' @param n.cores Integer. Number of cores. Default: 20.
#' @param seed Integer. Random seed for shared reference downsampling. Default: 123.
#' @param KS.cut Numeric. Passed to copykat().
#' @param ngene.chr Integer. Passed to copykat().
#' @param win.size Integer. Passed to copykat().
#' @param ... Additional arguments passed to copykat().
#'
#' @return Invisible data.frame of the combined results.
#' @export
runCopyKAT <- function(sce,
                            target_class,
                            ref_classes,
                            sample_col = "orig.ident",
                            celltype_col = "celltype",
                            assay = "RNA",
                            output_dir = "./CopyKAT",
                            overwrite = FALSE,
                            ref_mode = c("sample", "shared"),
                            min_cells = 25,
                            min_target_cells = 5,
                            min_ref_cells = 25,
                            max_ref_cells = 3000,
                            n.cores = 20,
                            seed = 123,
                            KS.cut = 0.1,
                            ngene.chr = 5,
                            win.size = 25,
                            ...) {

  ref_mode <- match.arg(ref_mode)

  # --- 1. Environment Setup ---
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  if (!"package:copykat" %in% search()) {
    message(">> Loading copykat package...")
    library(copykat)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)

  # --- 2. Input Validation ---
  message(">> Validating inputs...")
  if (!inherits(sce, "Seurat")) stop("'sce' must be a Seurat object.")
  if (!all(c(sample_col, celltype_col) %in% colnames(sce@meta.data))) {
    stop(sprintf("Columns '%s' or '%s' not found in metadata.", sample_col, celltype_col))
  }
  if (!assay %in% names(sce@assays)) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }

  meta <- sce@meta.data

  # --- 3. Identify Target and Reference Cells ---
  message(sprintf(">> Target class: '%s'", target_class))
  message(sprintf(">> Reference mode: '%s'", ref_mode))
  message(sprintf(">> Candidate reference classes: [%s]", paste(ref_classes, collapse = ", ")))

  target_cells <- rownames(meta)[meta[[celltype_col]] %in% target_class]
  ref_cells_all <- rownames(meta)[meta[[celltype_col]] %in% ref_classes]

  if (length(target_cells) == 0) stop("No target cells found!")
  if (length(ref_cells_all) == 0) stop("No reference cells found!")

  # Remove overlap just in case
  ref_cells_all <- setdiff(ref_cells_all, target_cells)

  # Shared reference preparation
  ref_cells_use <- NULL
  if (ref_mode == "shared") {
    if (length(ref_cells_all) < min_ref_cells) {
      stop(sprintf("Global reference cells are too few: %d (< %d).",
                   length(ref_cells_all), min_ref_cells))
    }

    if (!is.null(max_ref_cells) && length(ref_cells_all) > max_ref_cells) {
      set.seed(seed)
      ref_cells_use <- sample(ref_cells_all, max_ref_cells)
      message(sprintf(">> Downsampled shared references: %d -> %d",
                      length(ref_cells_all), length(ref_cells_use)))
    } else {
      ref_cells_use <- ref_cells_all
      message(sprintf(">> Using all shared references: %d", length(ref_cells_use)))
    }
  }

  sample_ids <- unique(meta[target_cells, sample_col])
  prediction_results <- list()

  message(">> Loading count matrix...")
  counts <- Seurat::GetAssayData(sce, slot = "counts", assay = assay)

  # --- 4. Main Loop ---
  for (sid in sample_ids) {
    message(sprintf("\n>>> Processing Sample: %s", sid))

    sample_subdir <- file.path(output_dir, sid)
    expected_file_1 <- file.path(sample_subdir, paste0("S_", sid, "_copykat_prediction.txt"))
    expected_file_2 <- file.path(sample_subdir, paste0(sid, "_copykat_prediction.txt"))

    found_file <- if (file.exists(expected_file_1)) {
      expected_file_1
    } else if (file.exists(expected_file_2)) {
      expected_file_2
    } else {
      NULL
    }

    # 4.1 Resume
    if (!overwrite && !is.null(found_file)) {
      message(sprintf("   [Found] %s. Loading...", basename(found_file)))
      tryCatch({
        pred <- read.table(found_file, header = TRUE, check.names = FALSE, sep = "\t")
        if (!"cell_id" %in% colnames(pred)) pred <- tibble::rownames_to_column(pred, "cell_id")
        pred$sample_id <- sid
        prediction_results[[sid]] <- pred
        next
      }, error = function(e) {
        message("   [Warning] Existing file corrupt. Re-running.")
      })
    }

    # 4.2 Run CopyKAT
    tryCatch({
      curr_cells <- rownames(meta)[meta[[sample_col]] == sid]
      curr_tgt <- intersect(curr_cells, target_cells)

      if (length(curr_tgt) < min_target_cells) {
        message(sprintf("   [Skip] Not enough target cells: %d", length(curr_tgt)))
        next
      }

      if (ref_mode == "sample") {
        curr_ref <- intersect(curr_cells, ref_cells_all)
        if (length(curr_ref) < min_cells) {
          message(sprintf("   [Skip] Not enough cells (Target: %d, Ref: %d).",
                          length(curr_tgt), length(curr_ref)))
          next
        }
      } else {
        curr_ref <- ref_cells_use
        if (length(curr_ref) < min_ref_cells) {
          message(sprintf("   [Skip] Shared references too few: %d", length(curr_ref)))
          next
        }
      }

      message(sprintf("   [Info] Target cells: %d | Ref cells: %d",
                      length(curr_tgt), length(curr_ref)))

      mat_tgt <- counts[, curr_tgt, drop = FALSE]
      mat_ref <- counts[, curr_ref, drop = FALSE]

      genes <- intersect(rownames(mat_tgt), rownames(mat_ref))
      mat_comb <- as.matrix(cbind(
        mat_tgt[genes, , drop = FALSE],
        mat_ref[genes, , drop = FALSE]
      ))
      mat_comb <- mat_comb[rowSums(mat_comb > 0) >= 5, , drop = FALSE]

      if (!dir.exists(sample_subdir)) dir.create(sample_subdir, recursive = TRUE)
      setwd(sample_subdir)

      res <- copykat(
        rawmat = mat_comb,
        norm.cell.names = curr_ref,
        sam.name = paste0("S_", sid),
        id.type = "S",
        output.seg = FALSE,
        n.cores = n.cores,
        KS.cut = KS.cut,
        ngene.chr = ngene.chr,
        win.size = win.size,
        ...
      )

      if (is.data.frame(res$prediction)) {
        pred <- res$prediction
        pred <- tibble::rownames_to_column(pred, "cell_id")
        pred$sample_id <- sid
        prediction_results[[sid]] <- pred
        message("   [Success] Done.")
      } else {
        message("   [Warning] No valid prediction returned.")
      }

      rm(mat_tgt, mat_ref, mat_comb, res)
      gc(verbose = FALSE)

    }, error = function(e) {
      message(sprintf("   [Error] %s", e$message))
    })

    setwd(original_wd)
  }

  # --- 5. Export Results ---
  if (length(prediction_results) == 0) {
    warning(">> No results generated.")
    return(NULL)
  }

  message("\n>> Aggregating results...")
  final_df <- dplyr::bind_rows(prediction_results)

  final_df <- final_df %>%
    dplyr::mutate(
      raw_id = cell_id,
      cleaned_id = mapply(function(id, samp) {
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
