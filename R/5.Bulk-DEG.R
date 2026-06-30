# =============== Limma \u5DEE\u5F02\u5206\u6790 ================
# =============== 1.Limma  ================
#' @title Run Standard Limma Differential Expression Analysis (Auto-Aligned)
#' @description A robust wrapper for limma. Automatically aligns expression matrix and metadata
#' based on sample IDs before analysis to prevent dimension mismatch errors.
#'
#' @param exp_mat Normalized expression matrix (genes x samples).
#' @param metadata Data frame containing sample annotations. Rownames must match exp_mat colnames.
#' @param group_col Column name in metadata defining the groups.
#' @param case_group Name of the Case group.
#' @param control_group Name of the Control group.
#' @param adj_method P-value adjustment method. Default "BH".
#' @param p_cutoff Adjusted P-value cutoff. Default 0.05.
#' @param logfc_cutoff Log2 Fold Change cutoff. Default 1.
#'
#' @return A data frame with DEG results and 'change' column.
#' @export
#' @importFrom limma modelMatrix lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' \dontrun{
#'   # 1. Create a mock expression matrix (100 genes x 6 samples).
#'   set.seed(123)
#'   exp_data <- matrix(rnorm(600), nrow = 100, ncol = 6)
#'   colnames(exp_data) <- c("S1", "S2", "S3", "S4", "S5", "S6")
#'   rownames(exp_data) <- paste0("Gene", 1:100)
#'
#'   # Simulate differential expression for the first five genes in S4-S6.
#'   exp_data[1:5, 4:6] <- exp_data[1:5, 4:6] + 3
#'
#'   # 2. Create metadata with one extra unused sample to test alignment.
#'   meta_data <- data.frame(
#'     row.names = c("S1", "S2", "S3", "S4", "S5", "S6", "S7_Unused"),
#'     condition = c("Ctrl", "Ctrl", "Ctrl", "Treat", "Treat", "Treat", "Treat"),
#'     batch = c(1, 1, 1, 2, 2, 2, 3)
#'   )
#'
#'   # 3. Run the analysis. The function keeps only matched samples.
#'   deg_res <- BulkLimma(exp_mat = exp_data,
#'                        metadata = meta_data,
#'                        group_col = "condition",
#'                        case_group = "Treat",
#'                        control_group = "Ctrl")
#'
#'   # 4. Inspect results.
#'   head(deg_res)
#'   table(deg_res$change)
#' }
BulkLimma <- function(exp_mat,
                      metadata,
                      group_col,
                      case_group,
                      control_group,
                      adj_method = "BH",
                      p_cutoff = 0.05,
                      logfc_cutoff = 1) {

  # --- 1. \u81EA\u52A8\u5BF9\u9F50\u4E0E\u68C0\u67E5 (\u5173\u952E\u4FEE\u590D) ---

  # \u627E\u51FA\u5171\u6709\u7684\u6837\u672C ID
  common_samples <- intersect(colnames(exp_mat), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("Error: No common samples found between expression matrix columns and metadata rownames. Please check your IDs.")
  }

  message(paste0("Matching samples: ", length(common_samples), " shared samples found."))

  # \u5982\u679C\u6570\u91CF\u4E0D\u4E00\u81F4\uFF0C\u7ED9\u51FA\u63D0\u793A
  if (length(common_samples) != ncol(exp_mat) || length(common_samples) != nrow(metadata)) {
    message("Subsetting data to matched samples only...")
  }

  # \u5BF9\u9F50\u6570\u636E\uFF1A\u53EA\u4FDD\u7559\u5171\u6709\u6837\u672C\uFF0C\u4E14\u987A\u5E8F\u4E00\u81F4
  exp_mat <- exp_mat[, common_samples]
  metadata <- metadata[common_samples, , drop = FALSE]

  # \u63D0\u53D6\u5206\u7EC4\u5411\u91CF
  group_vec <- metadata[[group_col]]

  # --- 2. \u68C0\u67E5\u5206\u7EC4\u6709\u6548\u6027 ---
  if (!case_group %in% group_vec || !control_group %in% group_vec) {
    stop(paste0("Group names '", case_group, "' or '", control_group, "' not found in column '", group_col, "'."))
  }

  # --- 3. \u6784\u5EFA\u8BBE\u8BA1\u77E9\u9635 ---
  # \u5904\u7406\u7EC4\u540D\u4E2D\u7684\u7279\u6B8A\u5B57\u7B26 (\u5982\u7A7A\u683C, -)
  clean_groups <- make.names(group_vec)
  clean_case <- make.names(case_group)
  clean_control <- make.names(control_group)

  # \u786E\u4FDD Control \u5728\u524D (\u4F5C\u4E3A\u5206\u6BCD)
  group_factor <- factor(clean_groups, levels = c(clean_control, clean_case))

  # \u8FD9\u91CC\u7684 design \u884C\u6570\u73B0\u5728\u4E00\u5B9A\u7B49\u4E8E exp_mat \u5217\u6570
  design <- stats::model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)

  # --- 4. \u7EBF\u6027\u62DF\u5408 ---
  fit <- limma::lmFit(exp_mat, design)

  # --- 5. \u6784\u5EFA\u5BF9\u6BD4 ---
  contrast_str <- paste0(clean_case, " - ", clean_control)
  message(paste0("Running Contrast: ", contrast_str))

  cont.matrix <- limma::makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2)

  # --- 6. \u63D0\u53D6\u7ED3\u679C ---
  DEG <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = adj_method)
  DEG <- tibble::rownames_to_column(DEG, var = "symbol")

  # \u6DFB\u52A0 UP/DOWN \u6807\u8BB0
  DEG$change <- "NOT"
  DEG$change[DEG$adj.P.Val < p_cutoff & DEG$logFC > logfc_cutoff] <- "UP"
  DEG$change[DEG$adj.P.Val < p_cutoff & DEG$logFC < -logfc_cutoff] <- "DOWN"

  return(DEG)
}



# =============== GSEA\u5BCC\u96C6\u5206\u6790  ================
# =============== 2. GSEA  ================
#' @title Run GSEA Analysis (GO) with Pre-ranked List
#' @description A clean wrapper for clusterProfiler::gseGO.
#' Note: The input 'gene_rank' must be a named numeric vector sorted in decreasing order.
#'
#' @param gene_rank A named numeric vector (names are genes, values are logFC or other metrics).
#' @param org_db Organism database. Default org.Hs.eg.db.
#' @param ont Ontology: "BP", "MF", or "CC". Default "BP".
#' @param key_type Key type of input gene names. Default "SYMBOL".
#' @param p_cutoff Adjusted P-value cutoff. Default 0.05.
#' @param min_size Minimum gene set size. Default 10.
#' @param max_size Maximum gene set size. Default 500.
#' @param seed Random seed. Default 123.
#'
#' @return A gseaResult object.
#' @export
#' @importFrom clusterProfiler gseGO
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#'   # 1. Prepare the ranked list yourself
#'   gene_rank <- sort(setNames(deg$logFC, deg$symbol), decreasing = TRUE)
#'
#'   # 2. Run GSEA
#'   gse_res <- BulkGseGO(gene_rank, ont = "BP")
#' }
BulkGseGO <- function(gene_rank,
                      org_db = org.Hs.eg.db,
                      ont = "BP",
                      key_type = "SYMBOL",
                      p_cutoff = 0.05,
                      min_size = 10,
                      max_size = 500,
                      seed = 123) {

  # 1. Basic Input Validation
  if (!is.numeric(gene_rank) || is.null(names(gene_rank))) {
    stop("Error: 'gene_rank' must be a named numeric vector (e.g., setNames(logFC, symbol)).")
  }

  # 2. Ensure Sort Order (Safety measure)
  # Although the user prepares the list, GSEA strictly requires decreasing order.
  # Re-sorting ensures the function doesn't crash if the user forgot.
  gene_rank <- sort(gene_rank, decreasing = TRUE)

  message(paste0("Running GSEA (", ont, ")... Input gene count: ", length(gene_rank)))

  # 3. Run gseGO
  set.seed(seed)

  gse_res <- clusterProfiler::gseGO(
    geneList     = gene_rank,
    OrgDb        = org_db,
    ont          = ont,
    keyType      = key_type,
    minGSSize    = min_size,
    maxGSSize    = max_size,
    pvalueCutoff = p_cutoff,
    verbose      = TRUE,
    seed         = TRUE
  )

  # 4. Result Check
  if (is.null(gse_res) || nrow(gse_res) == 0) {
    warning("No significant GSEA terms found.")
    return(NULL)
  }

  message(paste("GSEA Done. Significant terms:", nrow(gse_res)))
  return(gse_res)
}
