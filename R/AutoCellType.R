#' Automated Cell Type Annotation with GPT Models (Stable Version)
#'
#' @title AutoCellType
#' @description
#' A production-ready function for cell type annotation with GPT models.
#'
#' Key Features:
#' 1. **Robust Parsing**: Uses Regex to extract JSON, handling cases where the model adds conversational text.
#' 2. **Context Aware**: Explicitly handles 'major' vs 'subtype' annotation levels in the prompt.
#' 3. **Quality Control**: Automatically flags mitochondrial/stress signatures.
#'
#' @param input Data frame (must contain: cluster, gene, avg_log2FC) or a named list of gene vectors.
#' @param tissuename Character. Tissue origin (e.g., "Liver").
#' @param cellname Character (Optional). Parent cell type for sub-clustering (e.g., "T cells").
#' @param background Character (Optional). Experimental context (e.g., "Sorafenib treated HCC").
#' @param annotation_level Character. "major" (lineage) or "subtype" (functional state).
#' @param model Character. Model name (default: "gpt-5-nano").
#' @param topgenenumber Integer. Top N markers per cluster (default: 20).
#' @param base_url Character. API endpoint (default: "https://api.gpt.ge/v1").
#' @param api_key Character. API Key (if NULL, checks OPENAI_API_KEY env var).
#' @param retries Integer. Max API retry attempts (default: 3).
#' @param verbose Logical. Print progress logs (default: TRUE).
#'
#' @return A tibble with columns: Cluster, Prediction, Confidence, Reasoning.
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- 1. Prepare Mock Data (Simulated Liver T cell subtypes) ---
#'   library(tibble)
#'   test_markers <- tibble::tribble(
#'     ~cluster, ~gene,     ~avg_log2FC,
#'     "0",      "CD3D",    2.5,
#'     "0",      "CD8A",    2.3,
#'     "0",      "GZMB",    1.8,
#'     "1",      "CD3E",    2.4,
#'     "1",      "FOXP3",   2.5, # Treg
#'     "1",      "IL2RA",   1.8,
#'     "2",      "CD3D",    2.1,
#'     "2",      "PDCD1",   1.9, # Exhausted
#'     "2",      "HAVCR2",  1.6
#'   )
#'
#'   # --- 2. Run Annotation ---
#'   # Replace with your actual API Key
#'   my_key <- "sk-your_actual_api_key_here"
#'
#'   result <- AutoCellType(
#'     input = test_markers,
#'     tissuename = "Liver Cancer",
#'     cellname = "T cells",
#'     annotation_level = "subtype",
#'     background = "Hepatocellular Carcinoma immunotherapy treated",
#'     api_key = my_key
#'   )
#'
#'   # --- 3. View Results ---
#'   print(result)
#' }
AutoCellType <- function(
    input,
    tissuename,
    cellname = NULL,
    background = "",
    annotation_level = "major",
    model = "gpt-5-nano",
    topgenenumber = 20,
    base_url = "https://api.gpt.ge/v1",
    api_key = NULL,
    retries = 3,
    verbose = TRUE
) {

  # ----------------------------------------------------------------------------
  # 1. Dependency Check
  # ----------------------------------------------------------------------------
  check_deps <- function() {
    # Added 'stringr' for robust JSON extraction
    pkgs <- c("openai", "dplyr", "glue", "tibble", "jsonlite", "purrr", "stringr")
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "))
  }

  # ----------------------------------------------------------------------------
  # 2. Data Formatting
  # ----------------------------------------------------------------------------
  format_markers <- function(dat, n) {
    if (inherits(dat, "data.frame")) {
      req_cols <- c("cluster", "gene", "avg_log2FC")
      if (!all(req_cols %in% colnames(dat))) stop("Input must contain columns: ", paste(req_cols, collapse=", "))
      dat |>
        dplyr::filter(avg_log2FC > 0) |>
        dplyr::group_by(cluster) |>
        dplyr::slice_max(order_by = avg_log2FC, n = n) |>
        dplyr::summarise(genes = paste0(gene, collapse = ", "), .groups = "drop") |>
        tibble::deframe()
    } else if (inherits(dat, "list")) {
      sapply(dat, function(x) paste0(head(x, n), collapse = ", "))
    } else {
      stop("Input must be a data.frame or a named list.")
    }
  }

  # ----------------------------------------------------------------------------
  # 3. Prompt Engineering (Level Aware)
  # ----------------------------------------------------------------------------
  make_prompt <- function() {
    role <- "You are a distinguished expert in single-cell transcriptomics and immunology."

    # Context
    ctx <- glue::glue("Tissue Origin: {tissuename}.")
    if (!is.null(cellname)) ctx <- paste(ctx, "\nParent Population:", cellname)
    if (nzchar(background)) ctx <- paste(ctx, "\nExperimental Background:", background)

    # Dynamic Task Instruction based on annotation_level
    level_instruction <- if (annotation_level == "major") {
      "LEVEL: MAJOR LINEAGE. Identify broad categories only (e.g., T cells, B cells, Myeloid, Fibroblasts). Do NOT sub-cluster."
    } else {
      "LEVEL: DETAILED SUBTYPE. Identify specific functional states (e.g., CD8+ Exhausted, Treg, Kupffer cells, Inflammatory Macs)."
    }

    # Core Prompt
    glue::glue("
      {role}

      --- CONTEXT ---
      {ctx}

      --- TASK ---
      Target Level: {annotation_level}
      {level_instruction}

      I will provide clusters with their top marker genes. You must annotate them based on biological evidence.

      --- RULES ---
      1. **Strict Granularity**: Adhere strictly to the requested '{annotation_level}' level.
      2. **Tissue Logic**: Only predict cell types valid for {tissuename}.
      3. **QC**: If markers are predominantly Mitochondrial (MT-) or Ribosomal (RPS/RPL), label as 'Low Quality'.
      4. **Format**: OUTPUT RAW JSON ONLY. Do not output markdown text or code blocks.

      --- JSON STRUCTURE ---
      [
        {{
          \"cluster_id\": \"Identifier\",
          \"cell_type\": \"Annotation\",
          \"confidence\": \"High/Medium/Low\",
          \"reasoning\": \"Brief explanation referencing key genes.\"
        }}
      ]
    ")
  }

  # ----------------------------------------------------------------------------
  # 4. API Execution & Robust Parsing
  # ----------------------------------------------------------------------------
  run_batch <- function(batch_ids, batch_genes, client, sys_prompt) {
    input_str <- paste(batch_ids, ": ", batch_genes, collapse = "\n")
    user_msg <- glue::glue("Annotate these clusters:\n{input_str}")

    res_df <- NULL

    for (i in 1:retries) {
      if (verbose) message(sprintf("    -> Calling API (Attempt %d/%d)...", i, retries))

      out <- tryCatch({
        client$chat$completions$create(
          model = model,
          messages = list(
            list(role = "system", content = sys_prompt),
            list(role = "user", content = user_msg)
          ),
          temperature = 0.1, # Low temp for consistency
          response_format = list(type = "json_object")
        )
      }, error = function(e) e)

      if (!inherits(out, "error")) {
        # --- ROBUST EXTRACTION ---
        parse_res <- tryCatch({
          raw_text <- out$choices[[1]]$message$content

          # 1. Regex Extraction: Find the array [...] inside the text
          # (?s) allows dot to match newlines
          json_match <- stringr::str_extract(raw_text, "(?s)\\[.*\\]")

          # Use extracted match if found, otherwise raw text
          target_json <- if (!is.na(json_match)) json_match else raw_text

          # 2. Clean Markdown artifacts
          clean_json <- gsub("```json|```", "", target_json)

          # 3. Parse JSON
          json_dat <- jsonlite::fromJSON(clean_json)

          # 4. Validate Structure (handle single object vs list)
          if (!is.data.frame(json_dat)) {
            if (is.list(json_dat) && length(json_dat) == 1) json_dat <- json_dat[[1]]
          }

          if (is.data.frame(json_dat)) {
            colnames(json_dat) <- tolower(colnames(json_dat))
            json_dat
          } else {
            stop("JSON parsed but structure is invalid.")
          }

        }, error = function(e) e)

        if (!inherits(parse_res, "error")) {
          res_df <- parse_res
          break # Success
        } else {
          if (verbose) message("    !! Parse Error (Invalid JSON). Retrying...")
        }
      } else {
        if (verbose) message("    !! Network Error: ", out$message)
      }
      Sys.sleep(1) # Rate limit buffer
    }
    return(res_df)
  }

  # ----------------------------------------------------------------------------
  # 5. Main Execution
  # ----------------------------------------------------------------------------
  tryCatch({
    check_deps()

    key <- if (is.null(api_key)) Sys.getenv("OPENAI_API_KEY") else api_key
    if (!nzchar(key)) stop("API Key missing. Please provide 'api_key' or set OPENAI_API_KEY.")

    client <- openai::OpenAI(api_key = key, base_url = base_url)

    if (verbose) message(glue::glue("=== AutoCellType (Stable): {tissuename} [{annotation_level}] ==="))

    markers_vec <- format_markers(input, topgenenumber)
    sys_prompt <- make_prompt()

    # Batching (Size 5 is safe for context window)
    batch_size <- 5
    batches <- split(names(markers_vec), ceiling(seq_along(markers_vec)/batch_size))

    results <- purrr::map_dfr(names(batches), function(b) {
      ids <- batches[[b]]
      if (verbose) message(glue::glue("Processing Batch {b}..."))

      res <- run_batch(ids, markers_vec[ids], client, sys_prompt)

      if (is.null(res)) {
        return(tibble::tibble(
          cluster_id = ids,
          cell_type = "Failed",
          confidence = "None",
          reasoning = "API/Parse Error"
        ))
      }
      return(res)
    })

    final <- results |>
      dplyr::rename(
        Cluster = cluster_id,
        Prediction = cell_type,
        Confidence = confidence,
        Reasoning = reasoning
      ) |>
      dplyr::mutate(Annotation_Level = annotation_level) |>
      dplyr::select(Cluster, Prediction, Confidence, Reasoning, Annotation_Level)

    if (verbose) message("=== Done ===")
    return(final)

  }, error = function(e) stop("Fatal Error in AutoCellType: ", e$message))
}
