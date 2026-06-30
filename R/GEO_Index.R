# =============== GEO AI Tool ================
# =============== GEO_Classify_AI ================

#' @title AI-Powered GEO Data Classifier (Batch Processing)
#' @description Automatically classifies GEO datasets into 5 standardized types:
#' **scRNA-seq**, **stRNA-seq**, **Bulk RNA-seq**, **Microarray**, or **Other**.
#' Uses batch processing to handle large datasets efficiently.
#'
#' @section Input Data Requirements:
#' The input \code{input_data} must contain at least:
#' \itemize{
#'   \item \strong{ID Column}: Unique identifier (e.g., "GSE", "Accession").
#'   \item \strong{Description Column}: Text describing the study (e.g., "Title", "Summary", "Description").
#'   \item \strong{Platform Column (Optional)}: Sequencing platform info (e.g., "Platform", "GPL").
#' }
#'
#' @param input_data Data.frame. The raw metadata table.
#' @param batch_size Numeric. Number of rows to process per API call (Default 50).
#' @param api_key OpenAI-compatible API key. Defaults to `OPENAI_API_KEY`.
#' @param model Character or NULL. LLM model. If \code{NULL}, uses the
#'   internal general AI model from \code{R/AI_config.R}.
#' @param base_url Character or NULL. OpenAI-compatible API base URL. If
#'   \code{NULL}, uses the internal base URL from \code{R/AI_config.R}.
#' @param endpoint Character or NULL. One of \code{"auto"}, \code{"chat"}, or
#'   \code{"responses"}. \code{NULL} uses the package default \code{"auto"}.
#' @param verbose Logical. Show progress bar.
#'
#' @return A data.frame with a new column \code{AI_SeqType} appended.
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom glue glue
#' @importFrom dplyr bind_rows
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#'   # 1. Read raw GEO metadata.
#'   # raw_df <- read.csv("geo_metadata.csv")
#'
#'   # 2. Run classification.
#'   classified_df <- GEO_Classify_AI(
#'     input_data = raw_df,
#'     api_key = "sk-xxxxxx"
#'   )
#'
#'   # 3. Inspect classification summary.
#'   table(classified_df$AI_SeqType)
#'
#'   # 4. Optionally save results for downstream screening.
#'   # save(classified_df, file = "geo_classified.Rdata")
#' }
GEO_Classify_AI <- function(input_data,
                            batch_size = 50,
                            api_key = NULL,
                            model = NULL,
                            base_url = NULL,
                            endpoint = NULL,
                            verbose = TRUE) {

  # 1. Setup
  provider <- .easyAI_build_provider(
    api_key = api_key,
    model = model,
    base_url = base_url,
    endpoint = endpoint,
    task = "general"
  )
  model <- provider$model
  api_key <- provider$api_key

  if (!is.data.frame(input_data) || nrow(input_data) == 0) {
    warning("Input data is empty.")
    return(input_data)
  }

  # 2. Smart Column Detection
  cols <- names(input_data)

  # ID (GSE)
  id_col <- grep("\u6570\u636E\u96C6\u7F16\u53F7|GSE|Accession|dataset_id", cols, ignore.case=T, value=T)[1]
  if(is.na(id_col)) id_col <- cols[1]

  # Platform
  plat_col <- grep("\u6D4B\u5E8F\u5E73\u53F0|Platform|GPL|Instrument|sequencing_platform", cols, ignore.case=T, value=T)[1]

  # Description/Design
  desc_col <- grep("\u5B9E\u9A8C\u8BBE\u8BA1|Title|Summary|Description|experimental_design", cols, ignore.case=T, value=T)[1]

  if(is.na(desc_col)) stop("Could not find a 'Description/Title' column. AI needs text to classify.")

  if (verbose) message(glue::glue("\u2139\uFE0F  Columns Identified: ID='{id_col}', Platform='{plat_col}', Desc='{desc_col}'"))

  # 3. Batch Processing Setup
  n_rows <- nrow(input_data)
  n_batches <- ceiling(n_rows / batch_size)

  results_list <- list()

  if (verbose) {
    message(glue::glue("\U0001F680 Starting Classification ({model}): {n_rows} datasets in {n_batches} batches..."))
    pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
  }

  # 4. Loop Batches
  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_rows)
    batch_df <- input_data[start_idx:end_idx, ]

    # Prepare Text for AI
    batch_text <- apply(batch_df, 1, function(row) {
      p_txt <- if(!is.na(plat_col)) paste0("[Plat: ", row[plat_col], "]") else ""
      d_txt <- as.character(row[desc_col])
      # Truncate to save tokens (keep enough context)
      if (nchar(d_txt) > 300) d_txt <- paste0(substr(d_txt, 1, 297), "...")
      paste0(row[id_col], " | ", p_txt, " ", d_txt)
    })

    # Construct Prompt
    sys_prompt <- glue::glue("
      You are a Bioinformatics Classifier.

      --- TASK ---
      Classify each dataset into ONE category based on Platform and Description.

      --- CATEGORIES ---
      1. **scRNA-seq**: Single cell RNA sequencing (Keywords: 10x, Chromium, Single cell, scRNA, Drop-seq, Smart-seq).
      2. **stRNA-seq**: Spatial Transcriptomics (Keywords: Visium, Spatial, GeoMx, Slide-seq).
      3. **Microarray**: Chip/Array based (Keywords: GPL, Affymetrix, Agilent, Array, Expression profiling by array).
      4. **Bulk RNA-seq**: High throughput sequencing of tissue/bulk (Keywords: Illumina HiSeq/NextSeq, RNA-Seq) AND NOT Single Cell.
      5. **Other**: Methylation, ChIP-seq, miRNA, or Unknown.

      --- INPUT FORMAT ---
      ID | [Plat: ...] Description

      --- OUTPUT FORMAT ---
      JSON Object: {{ \"classifications\": [ {{ \"id\": \"GSE...\", \"type\": \"Category\" }}, ... ] }}

      --- DATA BATCH ---
      {paste(batch_text, collapse = '\n')}
    ")

    # API Call with Retry Logic
    success <- FALSE
    retry_count <- 0

    while (!success && retry_count < 3) {
      tryCatch({
        content <- .easyAI_call_provider_messages(
          messages = list(list(role = "system", content = sys_prompt)),
          provider = provider,
          api_key = api_key,
          temperature = 0,
          response_format = list(type = "json_object")
        )

        # Parse JSON
        parsed <- jsonlite::fromJSON(content)
        res_df <- as.data.frame(parsed$classifications)

        # Mapping back to batch_df robustly
        batch_results <- data.frame(ID = batch_df[[id_col]], AI_SeqType = "Unknown", stringsAsFactors = FALSE)

        # Match IDs (Robust merge to handle potential AI shuffling)
        match_idx <- match(batch_results$ID, res_df$id)
        valid_idx <- !is.na(match_idx)
        batch_results$AI_SeqType[valid_idx] <- res_df$type[match_idx[valid_idx]]

        results_list[[i]] <- batch_results
        success <- TRUE

      }, error = function(e) {
        retry_count <<- retry_count + 1
        if(verbose) message(glue::glue("\n\u26A0\uFE0F Batch {i} failed (Attempt {retry_count}): {e$message}"))
        Sys.sleep(1) # Wait before retry
      })
    }

    if (!success) {
      results_list[[i]] <- data.frame(
        ID = batch_df[[id_col]],
        AI_SeqType = "Unknown",
        stringsAsFactors = FALSE
      )
    }

    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # 5. Combine & Merge Results
  all_classifications <- dplyr::bind_rows(results_list)

  # Merge back to original dataframe (Preserve original order)
  final_df <- input_data
  # Add/Update the column
  final_df$AI_SeqType <- all_classifications$AI_SeqType[match(final_df[[id_col]], all_classifications$ID)]

  # 6. Summary Report
  if (verbose) {
    message("\n\u2705 Classification Complete! Summary:")
    print(table(final_df$AI_SeqType))
  }

  return(final_df)
}


# =============== GEO AI Screener ================
# =============== GEO_Screen_AI (Pure AI) ================

#' @title AI-Powered GEO Dataset Screener (Pure AI Edition)
#' @description Rely exclusively on LLM to screen datasets based on semantic understanding.
#' No regex/keyword matching is used. The AI reads every row's context (Disease, Type, Species, Description).
#' Supports unlimited rows via automatic batch processing.
#'
#' @section Input Data Requirements:
#' The \code{input_data} must be a \strong{data.frame} or \strong{tibble}.
#' While the function attempts to auto-detect columns, for best results, your table should contain:
#' \itemize{
#'   \item \strong{ID Column}: Unique identifier (e.g., named "GSE", "Accession", "dataset_id").
#'   \item \strong{Description Column}: Text describing the study (e.g., named "Title", "Summary", "Description", "experimental_design").
#'   \item \strong{Context Columns (Recommended)}:
#'     \itemize{
#'       \item "Species" or "Organism" (e.g., "Homo sapiens").
#'       \item "Platform" or "GPL" (e.g., "GPL570").
#'       \item "AI_SeqType" or "DataType" (e.g., "scRNA-seq", generated by \code{GEO_Classify_AI}).
#'       \item "disease_type" (e.g., "Liver Cancer").
#'     }
#'   \item \strong{Sample Size}: Numeric or string (e.g., "Sample_Count", "N", "sample_size") for sorting.
#' }
#'
#' @param input_data Data.frame. The metadata table. See "Input Data Requirements".
#' @param disease String. Disease criterion (e.g. "pancreatic cancer"). AI checks `disease_type` or Description.
#' @param data_type String. Sequencing type (e.g. "spatial transcriptomics"). AI checks `AI_SeqType` or Platform.
#' @param target_species String. Species criterion. AI checks `species` column.
#' @param other_req String. Any additional requirements (e.g., "Exclude cell lines").
#' @param batch_size Numeric. Rows per AI call (Default 500).
#' @param api_key OpenAI-compatible API key. Defaults to `OPENAI_API_KEY`.
#' @param model Character or NULL. LLM model. If \code{NULL}, uses the
#'   internal general AI model from \code{R/AI_config.R}.
#' @param base_url Character or NULL. OpenAI-compatible API base URL. If
#'   \code{NULL}, uses the internal base URL from \code{R/AI_config.R}.
#' @param endpoint Character or NULL. One of \code{"auto"}, \code{"chat"}, or
#'   \code{"responses"}. \code{NULL} uses the package default \code{"auto"}.
#' @param verbose Print debug info.
#'
#' @return A filtered, sorted, and deduplicated data.frame.
#' @export
#' @importFrom dplyr filter select arrange desc distinct mutate
#' @importFrom glue glue
#' @importFrom jsonlite fromJSON
#' @importFrom stringr str_extract
#' @importFrom rlang sym
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#'   # ========================================================
#'   # 1. Build mock data.
#'   # ========================================================
#'   # Assume this data was processed by GEO_Classify_AI.
#'   mock_df <- data.frame(
#'     dataset_id = c("GSE001", "GSE002", "GSE003", "GSE004"),
#'     species = c("Homo sapiens", "Mus musculus", "Homo sapiens", "Homo sapiens"),
#'     AI_SeqType = c("scRNA-seq", "Bulk RNA-seq", "stRNA-seq", "Microarray"),
#'     disease_type = c("Pancreatic Cancer", "Liver Cancer", "Pancreatic Ductal Adenocarcinoma", "Lung Cancer"),
#'     sample_size = c(10, 50, 5, 100),
#'     experimental_design = c(
#'       "Single cell analysis of 10 PDAC patients.",
#'       "Bulk RNA-seq of mouse liver tumor models.",
#'       "Visium spatial transcriptomics of pancreatic tumor sections.",
#'       "Gene expression profiling of lung cancer cell lines."
#'     )
#'   )
#'
#'   # ========================================================
#'   # 2. Basic screening: find pancreatic cancer scRNA-seq data.
#'   # ========================================================
#'   # AI can match related concepts such as PDAC and Pancreatic Cancer.
#'   res1 <- GEO_Screen_AI(
#'     input_data = mock_df,
#'     disease = "pancreatic cancer",
#'     data_type = "scRNA-seq",
#'     api_key = "sk-xxxxxx"
#'   )
#'   # Expected result: GSE001.
#'
#'   # ========================================================
#'   # 3. Semantic screening: find spatial transcriptomics data.
#'   # ========================================================
#'   # AI can match "spatial transcriptomics" to stRNA-seq or Visium.
#'   # It can also recognize PDAC as a pancreatic cancer context.
#'   res2 <- GEO_Screen_AI(
#'     input_data = mock_df,
#'     disease = "Pancreatic Cancer",
#'     data_type = "spatial transcriptomics",
#'     api_key = "sk-xxxxxx"
#'   )
#'   # Expected result: GSE003.
#'
#'   # ========================================================
#'   # 4. Additional logic: exclude cell lines.
#'   # ========================================================
#'   res3 <- GEO_Screen_AI(
#'     input_data = mock_df,
#'     other_req = "Exclude cell lines",
#'     api_key = "sk-xxxxxx"
#'   )
#'   # Expected result: exclude GSE004.
#' }
GEO_Screen_AI <- function(input_data,
                          disease = NULL,
                          data_type = NULL,
                          target_species = NULL,
                          other_req = NULL,
                          batch_size = 500,
                          api_key = NULL,
                          model = NULL,
                          base_url = NULL,
                          endpoint = NULL,
                          verbose = TRUE) {

  # 1. Environment & Input Check
  provider <- .easyAI_build_provider(
    api_key = api_key,
    model = model,
    base_url = base_url,
    endpoint = endpoint,
    task = "general"
  )
  model <- provider$model
  api_key <- provider$api_key
  pkgs <- c("dplyr", "glue", "jsonlite", "stringr")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse=", "))

  if (!is.data.frame(input_data) || nrow(input_data) == 0) {
    warning("Input data is empty.")
    return(input_data)
  }

  # 2. Smart Column Detection (Locate Context)
  cols <- names(input_data)

  # ID
  id_col <- grep("\u6570\u636E\u96C6\u7F16\u53F7|GSE|Accession|dataset_id", cols, ignore.case=T, value=T)[1]
  if(is.na(id_col)) id_col <- cols[1]

  # Context Columns (Match user provided str structure)
  sp_col <- grep("\u7269\u79CD|Species|Organism|species", cols, ignore.case=T, value=T)[1]
  type_col <- grep("AI_SeqType|SequencingType|DataType", cols, ignore.case=T, value=T)[1]
  dis_col <- grep("disease_type|Disease|Condition", cols, ignore.case=T, value=T)[1]
  desc_col <- grep("\u5B9E\u9A8C\u8BBE\u8BA1|Title|Summary|Description|experimental_design", cols, ignore.case=T, value=T)[1]
  plat_col <- grep("\u6D4B\u5E8F\u5E73\u53F0|Platform|GPL|sequencing_platform", cols, ignore.case=T, value=T)[1]

  # Sample Size (for sorting only)
  num_cols <- names(input_data)[sapply(input_data, is.numeric)]
  samp_col_num <- grep("\u6837\u672C\u91CF|Sample|Count|N_|sample_size", num_cols, ignore.case=T, value=T)[1]
  samp_col_any <- grep("\u6837\u672C\u91CF|Sample|Count|N_|sample_size", cols, ignore.case=T, value=T)[1]
  samp_col <- if(!is.na(samp_col_num)) samp_col_num else samp_col_any

  if (verbose) {
    message(glue::glue("\u2139\uFE0F  Context Columns Found:\n   - ID: {id_col}\n   - Species: {sp_col}\n   - Type: {type_col}\n   - Disease: {dis_col}\n   - Desc: {desc_col}"))
  }

  # 3. Construct AI Criteria Prompt
  criteria_list <- c()
  if (!is.null(disease)) criteria_list <- c(criteria_list, glue::glue("- Target Disease: {disease} (Check [Disease] tag or Description)"))
  if (!is.null(data_type)) criteria_list <- c(criteria_list, glue::glue("- Data Type: {data_type} (Check [Type] tag)"))
  if (!is.null(target_species)) criteria_list <- c(criteria_list, glue::glue("- Species: {target_species} (Check [Sp] tag)"))
  if (!is.null(other_req)) criteria_list <- c(criteria_list, glue::glue("- Other Req: {other_req}"))

  if (length(criteria_list) == 0) {
    warning("No screening criteria provided. Returning original data.")
    return(input_data)
  }
  criteria_str <- paste(criteria_list, collapse = "\n")

  # 4. Batch AI Processing
  n_total <- nrow(input_data)
  n_batches <- ceiling(n_total / batch_size)
  all_selected_ids <- c()

  if (verbose) {
    message(glue::glue("\U0001F9E0 Starting Pure AI Scan ({model}): {n_total} datasets in {n_batches} batches..."))
    pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
  }

  for (i in 1:n_batches) {
    # Slice Batch
    start_i <- (i - 1) * batch_size + 1
    end_i <- min(i * batch_size, n_total)
    batch_df <- input_data[start_i:end_i, ]

    # Construct Context String per Row
    data_summary <- apply(batch_df, 1, function(row) {
      # \u663E\u5F0F\u5730\u5C06\u5173\u952E\u5217\u7684\u5185\u5BB9\u6807\u8BB0\u51FA\u6765\uFF0C\u5582\u7ED9 AI
      info_parts <- c()
      if(!is.na(sp_col)) info_parts <- c(info_parts, paste0("[Sp: ", row[sp_col], "]"))
      if(!is.na(type_col)) info_parts <- c(info_parts, paste0("[Type: ", row[type_col], "]"))
      if(!is.na(dis_col)) info_parts <- c(info_parts, paste0("[Disease: ", row[dis_col], "]"))

      # \u8865\u5145\u63CF\u8FF0\u4FE1\u606F
      extra_desc <- ""
      if(!is.na(desc_col)) extra_desc <- as.character(row[desc_col])
      # \u5982\u679C\u63CF\u8FF0\u592A\u957F\uFF0C\u622A\u65AD\u4EE5\u8282\u7701 Token
      if (nchar(extra_desc) > 300) extra_desc <- paste0(substr(extra_desc, 1, 297), "...")

      info_str <- paste(c(info_parts, extra_desc), collapse = " ")
      paste0(row[id_col], ": ", info_str)
    })

    # System Prompt
    sys_prompt <- glue::glue("
      You are a Bioinformatics Expert.

      --- TASK ---
      Select dataset IDs that match the User Criteria.

      --- CRITERIA ---
      {criteria_str}

      --- RULES ---
      1. **Semantic Matching**:
         - User '\u80F0\u817A\u764C' matches data 'Pancreatic Cancer', 'PDAC' or context implying pancreas.
         - User '\u7A7A\u95F4\u8F6C\u5F55\u7EC4' matches data 'stRNA-seq', 'Visium', 'Spatial'.
         - User '\u5355\u7EC6\u80DE' matches 'scRNA-seq', 'Single Cell', '10x'.
      2. **Output**: Valid JSON object {{ \"selected_ids\": [\"ID1\", \"ID2\"] }}

      --- DATA BATCH ---
      {paste(data_summary, collapse = '\n')}
    ")

    # API Call
    tryCatch({
      content <- .easyAI_call_provider_messages(
        messages = list(list(role = "system", content = sys_prompt)),
        provider = provider,
        api_key = api_key,
        temperature = 0,
        response_format = list(type = "json_object")
      )

      # Parse IDs
      ids <- jsonlite::fromJSON(content)$selected_ids
      all_selected_ids <- c(all_selected_ids, ids)

    }, error = function(e) {
      if(verbose) message(glue::glue("\n\u26A0\uFE0F Batch {i} Error: {e$message}"))
    })

    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)

  # 5. Result Filtering
  if (length(all_selected_ids) > 0) {
    if (verbose) message(glue::glue("\u2705 AI Selected Total: {length(all_selected_ids)} IDs."))
    final_df <- input_data[input_data[[id_col]] %in% all_selected_ids, ]
  } else {
    message("\u26A0\uFE0F  No datasets matched criteria.")
    return(input_data[0, ])
  }

  # 6. Post-Processing (Sort & Deduplicate)
  final_df <- dplyr::distinct(final_df, !!rlang::sym(id_col), .keep_all = TRUE)

  if (!is.na(samp_col)) {
    # \u63D0\u53D6\u6570\u5B57\u6392\u5E8F
    sample_nums <- stringr::str_extract(as.character(final_df[[samp_col]]), "^\\d+")
    if (all(is.na(sample_nums))) sample_nums <- final_df[[samp_col]]

    # \u521B\u5EFA\u4E34\u65F6\u5217\u6392\u5E8F\uFF0C\u9632\u6B62\u6C61\u67D3\u539F\u6570\u636E
    final_df$tmp_sort_N <- suppressWarnings(as.numeric(sample_nums))
    final_df <- final_df[order(final_df$tmp_sort_N, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    final_df$tmp_sort_N <- NULL
  }

  # Clean Columns (Put context cols first)
  prio_cols <- c(id_col, sp_col, type_col, dis_col, samp_col, desc_col)
  prio_cols <- prio_cols[!is.na(prio_cols) & prio_cols %in% names(final_df)]
  other_cols <- setdiff(names(final_df), prio_cols)
  final_df <- dplyr::select(final_df, dplyr::all_of(c(prio_cols, other_cols)))

  return(final_df)
}
