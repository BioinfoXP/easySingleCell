#' Run GEO Data Miner (CLI Version)
#'
#' @param keyword String. The main topic (e.g., "HCC immunotherapy").
#' @param context String. Detailed criteria for scoring (e.g., "Must be clinical tissue").
#' @param limit Integer. Max number of datasets to analyze. Default 10.
#' @param output_file String. Path to save Excel file. If NULL, returns data frame only.
#' @param cores Integer. Number of CPU cores. Default is auto-detect.
#' @param api_key String. LLM API Key. Defaults to Sys.getenv("OPENAI_API_KEY").
#' @param base_url String. API Endpoint. Default "https://api.openai.com/v1".
#' @param model String. Model name. Default "gpt-4o-mini".
#'
#' @return A data frame containing the mining results.
#' @export
#' @importFrom openxlsx createWorkbook addWorksheet writeData createStyle addStyle setColWidths saveWorkbook
#' @importFrom progressr with_progress
#'
#' @examples
#' \dontrun{
#' # 1. Simple search for HCC immunotherapy
#' results <- geo_mine(
#'   keyword = "HCC immunotherapy",
#'   limit = 5,
#'   api_key = "sk-..."
#' )
#'
#' # 2. Advanced search with context and DeepSeek API
#' # This searches for lung cancer data specifically looking for EGFR mutations
#' # and excluding cell lines.
#' results_deepseek <- geo_mine(
#'   keyword = "Lung Cancer",
#'   context = "Must be EGFR mutation positive. Tissue samples only. No cell lines.",
#'   limit = 20,
#'   output_file = "lung_egfr_results.xlsx",
#'   base_url = "https://api.deepseek.com/v1",
#'   model = "deepseek-chat",
#'   api_key = "sk-..."
#' )
#' }
geo_mine <- function(keyword,
                     context = "",
                     limit = 10,
                     output_file = "geo_results.xlsx",
                     cores = parallelly::availableCores() - 1,
                     api_key = Sys.getenv("OPENAI_API_KEY"),
                     base_url = "https://api.gpt.ge/v1",
                     model = "gpt-4o-mini") {

  if (api_key == "") stop("API Key is missing. Provide it as argument or set OPENAI_API_KEY env var.")
  if (cores < 1) cores <- 1

  # 包装 progressr 以在终端显示进度条
  results <- progressr::with_progress({
    geo_mine_engine(keyword, context, limit, cores, api_key, base_url, model)
  })

  if (is.null(results) || nrow(results) == 0) {
    message("No results found.")
    return(invisible(NULL))
  }

  message(sprintf("Done! Found %d datasets.", nrow(results)))

  # --- Excel 导出逻辑 (移植自 server.R) ---
  if (!is.null(output_file)) {
    message(sprintf("Saving to %s...", output_file))

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "AI Analysis Results")
    openxlsx::writeData(wb, "AI Analysis Results", results, startRow = 1, startCol = 1)

    # 定义样式
    header_style <- openxlsx::createStyle(fontSize = 11, fontColour = "#FFFFFF", fgFill = "#2c3e50",
                                          halign = "center", textDecoration = "bold")
    style_high <- openxlsx::createStyle(fgFill = "#27ae60", fontColour = "white")
    style_med  <- openxlsx::createStyle(fgFill = "#f39c12", fontColour = "white")
    style_low  <- openxlsx::createStyle(fgFill = "#95a5a6", fontColour = "white")

    # 应用样式
    openxlsx::addStyle(wb, "AI Analysis Results", header_style, rows = 1, cols = 1:ncol(results), gridExpand = TRUE)

    # 分数条件格式化
    if ("Match_Score" %in% names(results)) {
      scores <- results$Match_Score
      rows_high <- which(scores >= 80) + 1
      rows_med  <- which(scores >= 50 & scores < 80) + 1
      rows_low  <- which(scores < 50) + 1
      col_idx <- which(names(results) == "Match_Score")

      if(length(rows_high) > 0) openxlsx::addStyle(wb, "AI Analysis Results", style_high, rows=rows_high, cols=col_idx)
      if(length(rows_med) > 0) openxlsx::addStyle(wb, "AI Analysis Results", style_med, rows=rows_med, cols=col_idx)
      if(length(rows_low) > 0) openxlsx::addStyle(wb, "AI Analysis Results", style_low, rows=rows_low, cols=col_idx)
    }

    # 设置列宽
    openxlsx::setColWidths(wb, "AI Analysis Results", cols = 1:ncol(results), widths = "auto")
    openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  }

  return(results)
}
