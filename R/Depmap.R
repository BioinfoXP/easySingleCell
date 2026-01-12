# =============== DepMap Core ================
# =============== 1.DepmapPrepare ================

#' @title Prepare DepMap Data (Manual Download Guide)
#' @description Checks for DepMap CSV files. If missing, provides specific download links
#' (including fast mirror and official sources) and instructions.
#' If present, processes data into an .Rdata file.
#'
#' @param target_genes Character vector. Genes of interest to extract.
#' @param data_dir String. Directory to store data. Default "./Depmap".
#' @param force_process Logical. If TRUE, re-processes the CSVs even if RData exists. Default FALSE.
#'
#' @export
#' @importFrom data.table fread
#' @importFrom dplyr distinct
DepmapPrepare <- function(target_genes,
                          data_dir = "./Depmap",
                          force_process = FALSE) {

  # 1. 建立目录
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }

  # 定义必须存在的 3 个文件名
  required_files <- c(
    "CRISPRGeneEffect.csv",
    "CRISPRGeneDependency.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
  )

  # 检查文件是否缺失
  missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]

  # === 场景 A: 文件缺失 -> 打印详细指引并停止 ===
  if (length(missing_files) > 0) {
    abs_path <- normalizePath(data_dir, mustWork = FALSE)

    message("\n=========================================================================")
    message("❌ MISSING DATA FILES / 数据文件缺失")
    message("-------------------------------------------------------------------------")
    message("The following files are missing in your data directory:")
    for (f in missing_files) {
      message(paste("   -", f))
    }

    message("\nPLEASE DOWNLOAD MANUALLY / 请手动下载:")
    message("-------------------------------------------------------------------------")

    message("🚀 Option 1: Fast Mirror (Recommended / 推荐高速下载):")
    message("   🔗 https://www.123865.com/s/nnlSTd-JUj0h?pwd=DPMP")
    message("   🔑 Password: DPMP")

    message("\n🌐 Option 2: Official Direct Links (DepMap 25Q3):")
    message("   (Requires DepMap Login / 需要登录官网)")
    message("   1️⃣ CRISPRGeneEffect.csv:")
    message("      👉 https://depmap.org/portal/api/download/file?releasename=DepMap%20Public%2025Q3&filename=CRISPRGeneEffect.csv")
    message("   2️⃣ CRISPRGeneDependency.csv:")
    message("      👉 https://depmap.org/portal/api/download/file?releasename=DepMap%20Public%2025Q3&filename=CRISPRGeneDependency.csv")
    message("   3️⃣ OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv:")
    message("      👉 https://depmap.org/portal/api/download/file?releasename=DepMap%20Public%2025Q3&filename=OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")

    message("\n📂 ACTION REQUIRED:")
    message("   Download the files and move them into this folder:")
    message(paste("   📂", abs_path))
    message("=========================================================================\n")

    stop("Process halted. Please download the files and try again.")
  }

  # === 场景 B: 文件存在 -> 开始处理 ===
  rdata_path <- file.path(data_dir, "DepMap_Processed.Rdata")

  # 如果 RData 已存在且不强制更新，直接结束
  if (file.exists(rdata_path) && !force_process) {
    message("✅ Processed data already exists: ", rdata_path)
    message("   Skipping processing. Use force_process=TRUE to overwrite.")
    return(invisible(NULL))
  }

  message("✅ Files found! Starting data processing...")
  message(">>> Step 1/3: Reading CSVs (This may take 1-2 minutes)...")

  # 内部读取函数
  read_clean <- function(fname) {
    fpath <- file.path(data_dir, fname)
    df <- data.table::fread(fpath, data.table = FALSE)
    colnames(df)[1] <- "ModelID"
    # 清洗列名: "A1BG (1)" -> "A1BG"
    colnames(df) <- c("ModelID", sub("\\s\\(.*\\)", "", colnames(df)[-1]))
    df <- df[!duplicated(df$ModelID), ]
    rownames(df) <- df$ModelID
    df$ModelID <- NULL
    return(df)
  }

  df_eff <- read_clean(required_files[1])
  df_dep <- read_clean(required_files[2])
  df_exp <- read_clean(required_files[3])

  message(">>> Step 2/3: Intersecting cell lines and genes...")
  common_cells <- Reduce(intersect, list(rownames(df_eff), rownames(df_dep), rownames(df_exp)))

  valid_genes <- intersect(target_genes, colnames(df_eff))
  if (length(valid_genes) == 0) stop("No valid genes found in the datasets! Check gene symbols.")

  missing_genes <- setdiff(target_genes, valid_genes)
  if (length(missing_genes) > 0) {
    warning("Some genes were not found: ", paste(missing_genes, collapse = ", "))
  }

  depmap_list <- list(
    effect = df_eff[common_cells, valid_genes, drop = FALSE],
    dependency = df_dep[common_cells, valid_genes, drop = FALSE],
    expression = df_exp[common_cells, valid_genes, drop = FALSE],
    genes = valid_genes,
    timestamp = Sys.time()
  )

  message(paste(">>> Step 3/3: Saving RData to", rdata_path))
  save(depmap_list, file = rdata_path)
  message("🎉 Done! You can now run DepmapBox(), DepmapHeatmap(), etc.")
}

# =============== DepMap Viz ================
# =============== 2.DepmapBox ================

#' @title Visualize DepMap Gene Effect (With Flexible Metadata Filtering)
#' @description Plots Gene Effect/Dependency scores. Supports flexible sample filtering via
#' dplyr-style expressions on metadata or explicit IDs.
#'
#' @param genes A character vector of gene symbols.
#' @param rdata_path Path to processed RData.
#' @param data_type "effect" or "dependency".
#' @param cell_filter (Optional) A dplyr-style filtering expression.
#' Example: \code{cell_filter = OncotreeLineage == "Lung"}.
#' @param subset_ids (Optional) Character vector of ModelIDs to include.
#' @param threshold Threshold line value.
#' @param add_dots Show jitter points?
#' @param jitter.alpha Transparency (0-1).
#' @param jitter.color Color of points.
#' @param pal Custom palette.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_hline labs theme element_text unit scale_fill_manual
#' @importFrom cowplot theme_cowplot
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter
#' @importFrom rlang enquo quo_is_null
#'
#' @examples
#' \dontrun{
#'   # 1. Basic Usage
#'   DepmapBox(genes = "TP53")
#'
#'   # 2. Metadata Filtering (Standard Column Names)
#'   DepmapBox(genes = "TP53", cell_filter = OncotreeLineage == "Lung" & Sex == "Female")
#'
#'   # 3. Combine with AI Selection (Unified Naming)
#'   selected_ids <- DepmapMetaSelect("Lung Adenocarcinoma", "OncotreePrimaryDisease", api_key="...")
#'   DepmapBox(genes = "KRAS", subset_ids = selected_ids)
#' }
DepmapBox <- function(genes,
                      rdata_path = "./Depmap/DepMap_Processed.Rdata",
                      data_type = c("effect", "dependency"),
                      cell_filter = NULL,
                      subset_ids = NULL,
                      threshold = NULL,
                      add_dots = TRUE,
                      jitter.alpha = 0.1,
                      jitter.color = "grey40",
                      pal = NULL) {

  # 1. Check Arguments
  data_type <- match.arg(data_type)
  if (missing(genes)) stop("Error: Please specify 'genes'.")
  if (!file.exists(rdata_path)) stop(paste("Error: RData not found at", rdata_path))

  # 2. Config
  if (data_type == "effect") {
    target_obj <- "depmap_effect"
    title_str <- "CRISPR Gene Effect"
    ylab_str <- "Gene Effect (Chronos)"
    if (is.null(threshold)) threshold <- -0.5
  } else {
    target_obj <- "depmap_dependency"
    title_str <- "CRISPR Gene Dependency"
    ylab_str <- "Gene Dependency"
    if (is.null(threshold)) threshold <- 0.5
  }

  # 3. Load Data
  env <- new.env()
  load(rdata_path, envir = env)
  if (!exists(target_obj, envir = env)) stop("Invalid RData.")
  dat_matrix <- get(target_obj, envir = env)

  has_meta <- exists("depmap_metadata", envir = env)
  meta_df <- if (has_meta) env$depmap_metadata else NULL

  # 4. === Filtering Logic ===
  filter_quo <- rlang::enquo(cell_filter)

  # 4.1 Apply cell_filter
  if (!rlang::quo_is_null(filter_quo)) {
    if (!has_meta) stop("Metadata required for 'cell_filter'.")
    message("ℹ️ Applying metadata filter...")
    filtered_meta <- dplyr::filter(meta_df, !!filter_quo)
    target_ids <- rownames(filtered_meta)
  } else {
    target_ids <- rownames(dat_matrix)
  }

  # 4.2 Apply subset_ids
  if (!is.null(subset_ids)) {
    target_ids <- intersect(target_ids, subset_ids)
    if (length(target_ids) == 0) stop("Error: Intersection of filters is empty.")
  }

  # 4.3 Final Intersection with Matrix
  final_ids <- intersect(rownames(dat_matrix), target_ids)
  if (length(final_ids) == 0) stop("Error: No cell lines matched.")

  if (length(final_ids) < nrow(dat_matrix)) {
    message(sprintf("ℹ️ Subset: Plotting %d / %d cell lines.", length(final_ids), nrow(dat_matrix)))
  }

  dat_matrix <- dat_matrix[final_ids, , drop = FALSE]

  # 5. Process & Plot
  valid_genes <- intersect(genes, colnames(dat_matrix))
  if (length(valid_genes) == 0) stop("Genes not found.")

  df_sub <- data.frame(dat_matrix[, valid_genes, drop=FALSE])
  df_sub$ModelID <- rownames(dat_matrix)
  df_long <- tidyr::pivot_longer(df_sub, -ModelID, names_to = "Gene", values_to = "Score")

  if (is.null(pal)) pal <- c("#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Gene, y = Score, fill = Gene)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6, size = 0.4) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "grey30", size = 0.5) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")
    ) +
    ggplot2::labs(title = title_str, x = "", y = ylab_str)

  if (add_dots) {
    p <- p + ggplot2::geom_jitter(width = 0.2, color = jitter.color, size = 1.0, alpha = jitter.alpha)
  }

  return(p)
}



# =============== DepMap Viz ================
# =============== 3.DepmapScatter ================

#' @title Visualize Expression vs. Effect Correlation (With Filtering)
#' @description Plots the correlation between Gene Expression and Gene Effect.
#' Supports metadata filtering and manual stats calculation.
#'
#' @param genes A character vector of gene symbols.
#' @param rdata_path Path to processed RData.
#' @param data_type "effect" or "dependency".
#' @param cell_filter (Optional) dplyr-style metadata filter.
#' @param subset_ids (Optional) Character vector of ModelIDs.
#' @param cor_method "pearson" or "spearman".
#' @param point.color Color of points.
#' @param point.alpha Transparency (0-1).
#' @param line.color Color of regression line.
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme element_text annotate
#' @importFrom ggpubr ggarrange
#' @importFrom cowplot theme_cowplot
#' @importFrom stats cor.test complete.cases
#' @importFrom dplyr filter
#' @importFrom rlang enquo quo_is_null
#'
#' @examples
#' \dontrun{
#'   # 1. Basic Usage
#'   DepmapScatter(genes = "TP53")
#'
#'   # 2. Metadata Filtering
#'   DepmapScatter(genes = "TP53", cell_filter = OncotreeLineage == "Lung")
#'
#'   # 3. Combine with AI Selection (Unified Naming)
#'   selected_ids <- DepmapMetaSelect("Lung Adenocarcinoma", "OncotreePrimaryDisease", api_key="...")
#'   DepmapScatter(genes = "KRAS", subset_ids = selected_ids)
#' }
DepmapScatter <- function(genes,
                          rdata_path = "./Depmap/DepMap_Processed.Rdata",
                          data_type = c("effect", "dependency"),
                          cell_filter = NULL,
                          subset_ids = NULL,
                          cor_method = "pearson",
                          point.color = "navy",
                          point.alpha = 0.5,
                          line.color = "orange") {

  # 1. Check Arguments
  data_type <- match.arg(data_type)
  if (missing(genes)) stop("Error: Please specify 'genes'.")
  if (!file.exists(rdata_path)) stop(paste("Error: RData not found at", rdata_path))

  # 2. Config
  if (data_type == "effect") {
    target_obj_name <- "depmap_effect"
    y_lab_base <- "Gene Effect (Chronos)"
  } else {
    target_obj_name <- "depmap_dependency"
    y_lab_base <- "Gene Dependency"
  }

  # 3. Load Data
  env <- new.env()
  load(rdata_path, envir = env)
  if (!exists("depmap_expression", envir = env)) stop("Expression data not found.")
  if (!exists(target_obj_name, envir = env)) stop("Target data not found.")

  dat_exp <- env$depmap_expression
  dat_tar <- get(target_obj_name, envir = env)

  has_meta <- exists("depmap_metadata", envir = env)
  meta_df <- if (has_meta) env$depmap_metadata else NULL

  # 4. === Filtering Logic ===
  filter_quo <- rlang::enquo(cell_filter)

  if (!rlang::quo_is_null(filter_quo)) {
    if (!has_meta) stop("Metadata required for 'cell_filter'.")
    message("ℹ️ Applying metadata filter...")
    filtered_meta <- dplyr::filter(meta_df, !!filter_quo)
    target_ids <- rownames(filtered_meta)
  } else {
    target_ids <- rownames(dat_exp)
  }

  if (!is.null(subset_ids)) {
    target_ids <- intersect(target_ids, subset_ids)
    if (length(target_ids) == 0) stop("Error: No cells left after 'subset_ids'.")
  }

  # Intersection of Exp + Target + Filter
  common_ids <- Reduce(intersect, list(rownames(dat_exp), rownames(dat_tar), target_ids))

  if (length(common_ids) == 0) stop("Error: 0 common cell lines found.")
  if (length(common_ids) < nrow(dat_exp)) message(sprintf("ℹ️ Subset: Analyzing %d cell lines.", length(common_ids)))

  # 5. Process & Plot
  valid_genes <- intersect(genes, intersect(colnames(dat_exp), colnames(dat_tar)))
  if (length(valid_genes) == 0) stop("Genes not found.")

  plot_list <- lapply(valid_genes, function(g) {
    x_val <- dat_exp[common_ids, g]
    y_val <- dat_tar[common_ids, g]

    clean_idx <- stats::complete.cases(x_val, y_val)
    x_val <- x_val[clean_idx]
    y_val <- y_val[clean_idx]

    ct <- stats::cor.test(x_val, y_val, method = cor_method)
    r_val <- ct$estimate
    p_val <- ct$p.value

    p_str <- if(p_val < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.2e", p_val)
    stats_label <- paste0("R = ", sprintf("%.2f", r_val), ", ", p_str)

    df_plot <- data.frame(Expr = x_val, YVal = y_val)

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Expr, y = YVal)) +
      ggplot2::geom_point(color = point.color, alpha = point.alpha, size = 1.5) +
      ggplot2::geom_smooth(method = "lm", color = line.color, fill = "grey85", size = 1.0) +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = stats_label,
                        hjust = -0.1, vjust = 1.5, size = 5, fontface = "italic") +
      ggplot2::labs(title = g, x = "Gene Expression (log2 TPM+1)", y = y_lab_base) +
      cowplot::theme_cowplot() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16))

    return(p)
  })

  if (length(plot_list) == 1) return(plot_list[[1]])
  else return(ggpubr::ggarrange(plotlist = plot_list, ncol = min(3, length(plot_list)), nrow = ceiling(length(plot_list)/3)))
}



# =============== DepMap AI ================
# =============== 4.DepmapMetaSelect ================

#' @title AI-Powered Fuzzy Metadata Search
#' @description Performs semantic/fuzzy search within a specific metadata column.
#' Returns ModelIDs matching the query.
#'
#' @param query String. Fuzzy search term (e.g., "Lung Cancer").
#' @param col_name String. Metadata column (e.g., "OncotreePrimaryDisease").
#' @param rdata_path Path to processed RData.
#' @param model LLM model name.
#' @param api_key OpenAI API Key.
#' @param base_url API Base URL.
#' @param verbose Print debug info.
#'
#' @return A character vector of ModelIDs.
#' @export
#' @importFrom openai OpenAI
#' @importFrom jsonlite fromJSON
#' @importFrom glue glue
#' @importFrom stringr str_remove_all
#'
#' @examples
#' \dontrun{
#'   # ========================================================
#'   # Example: Search and Visualize (Unified Naming)
#'   # ========================================================
#'
#'   # 1. Select IDs using AI
#'   selected_ids <- DepmapMetaSelect(
#'     query = "Lung Cancer",
#'     col_name = "OncotreePrimaryDisease",
#'     api_key = "sk-xxxxxx"
#'   )
#'
#'   # 2. Pass selected_ids to Boxplot
#'   DepmapBox(genes = "TP53", subset_ids = selected_ids)
#'
#'   # 3. Pass selected_ids to Scatterplot
#'   DepmapScatter(genes = "KRAS", subset_ids = selected_ids)
#' }
DepmapMetaSelect <- function(query,
                             col_name,
                             rdata_path = "./Depmap/DepMap_Processed.Rdata",
                             model = "gpt-4o",
                             api_key = NULL,
                             base_url = "https://api.gpt.ge/v1",
                             verbose = TRUE) {

  # 1. Setup
  pkgs <- c("openai", "jsonlite", "glue", "stringr", "dplyr")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse=", "))

  if (is.null(api_key)) api_key <- Sys.getenv("OPENAI_API_KEY")
  if (api_key == "") stop("Please provide 'api_key'.")

  # 2. Load Metadata
  if (verbose) message(">>> Loading Metadata...")
  env <- new.env()
  load(rdata_path, envir = env)
  if (!exists("depmap_metadata", envir = env)) stop("Metadata not found.")
  meta <- env$depmap_metadata

  if (!col_name %in% colnames(meta)) stop("Invalid column name: ", col_name)

  # 3. Prepare Context
  valid_values <- unique(na.omit(meta[[col_name]]))
  values_str <- paste(valid_values, collapse = ", ")
  if (verbose) message(glue::glue("ℹ️  Search Space: {length(valid_values)} values in '{col_name}'."))

  # 4. Prompt
  sys_prompt <- glue::glue("
    You are a Semantic Search Engine for biomedical metadata.
    --- TASK ---
    Identify categories from the 'Valid Categories' list that match the User's Query.
    --- RULES ---
    1. **Fuzzy Matching**: Match semantically related terms.
    2. **Exact Output**: Return strings EXACTLY as listed.
    3. **JSON Only**: {{ \"matches\": [\"String1\", \"String2\"] }}
    --- VALID CATEGORIES ---
    [{values_str}]
  ")

  user_msg <- glue::glue("User Query: '{query}'")

  # 5. API Call
  if (verbose) message(glue::glue(">>> Searching for '{query}' in '{col_name}'..."))
  client <- openai::OpenAI(api_key = api_key, base_url = base_url)

  resp <- tryCatch({
    client$chat$completions$create(
      model = model,
      messages = list(list(role = "system", content = sys_prompt), list(role = "user", content = user_msg)),
      temperature = 0,
      response_format = list(type = "json_object")
    )
  }, error = function(e) stop("API Error: ", e$message))

  # 6. Parse
  raw_text <- resp$choices[[1]]$message$content
  clean_text <- raw_text %>%
    stringr::str_remove_all("^```json") %>%
    stringr::str_remove_all("^```") %>%
    stringr::str_remove_all("```$") %>%
    trimws()

  parsed <- tryCatch({ jsonlite::fromJSON(clean_text) }, error = function(e) stop("JSON Parsing Failed."))
  matched_vals <- unlist(parsed$matches)
  verified_vals <- intersect(matched_vals, valid_values)

  if (length(verified_vals) == 0) {
    warning("⚠️ AI found no matching values.")
    return(character(0))
  }

  if (verbose) {
    message("✅ AI Matches Found:")
    print(verified_vals)
  }

  # 7. Filter & Return
  result_ids <- rownames(meta[meta[[col_name]] %in% verified_vals, , drop=FALSE])
  if (verbose) message(glue::glue(">>> Found {length(result_ids)} cell lines."))

  return(result_ids)
}
