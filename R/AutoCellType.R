#' Automated Cell Type Annotation with GPT Models
#'
#' @title AutoCellType
#' @description
#' Leverages GPT models to annotate cell types based on marker genes.
#' This function enforces STRICT JSON output to ensure stability and provides
#' "Reasoning" and "Confidence" scores to make the results interpretable.
#'
#' @param input Data frame (must contain columns: cluster, gene, avg_log2FC) or a named list of gene vectors.
#' @param tissuename Character string. Tissue origin (e.g., "Liver", "PBMC").
#' @param cellname Optional character string. Parent cell type (e.g., "T cells", used for sub-clustering).
#' @param background Optional character string. **Crucial for smart annotation**. Experimental context (e.g., "Tumor microenvironment, Sorafenib resistance").
#' @param annotation_level Character string. "major" (lineage level) or "subtype" (specific state).
#' @param model Character string. Model name (default: "gpt-4o").
#' @param topgenenumber Integer. Number of top marker genes to use per cluster (default: 20).
#' @param base_url Character string. API endpoint (default: "https://api.gpt.ge/v1").
#' @param api_key Character string. API Key (if NULL, looks for OPENAI_API_KEY env variable).
#' @param retries Integer. Maximum retry attempts for API calls (default: 3).
#' @param verbose Logical. Whether to print progress logs (default: TRUE).
#'
#' @return A tibble with columns: Cluster, Cell_Type, Confidence, Reasoning, Annotation_Level.
#' @export
#'
#' @examples
#' \dontrun{
#'   # --- 1. Prepare Test Data (Simulated Liver T cell subtypes) ---
#'   library(tibble)
#'   test_markers <- tibble::tribble(
#'     ~cluster, ~gene,     ~avg_log2FC,
#'     "0",      "CD3D",    2.5,
#'     "0",      "CD8A",    2.3,
#'     "0",      "GZMB",    1.8,
#'     "1",      "CD3E",    2.4,
#'     "1",      "FOXP3",   2.5, # Treg marker
#'     "1",      "IL2RA",   1.8,
#'     "2",      "CD3D",    2.1,
#'     "2",      "PDCD1",   1.9, # Exhaustion marker
#'     "2",      "HAVCR2",  1.6
#'   )
#'
#'   # --- 2. Run Annotation ---
#'   # Replace with your actual API Key
#'   my_key <- "sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
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
#'   # Expected output should include a 'Reasoning' column explaining
#'   # why Cluster 1 is Treg (based on FOXP3) etc.
#' }
AutoCellType <- function(
    input,
    tissuename,
    cellname = NULL,
    background = "",
    annotation_level = "major",
    model = "gpt-4o",
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
    pkgs <- c("openai", "dplyr", "glue", "tibble", "jsonlite", "purrr")
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "))
  }

  # ----------------------------------------------------------------------------
  # 2. Data Formatting
  # ----------------------------------------------------------------------------
  format_markers <- function(dat, n) {
    if (inherits(dat, "data.frame")) {
      req_cols <- c("cluster", "gene", "avg_log2FC")
      if (!all(req_cols %in% colnames(dat))) {
        stop("Input data.frame must contain columns: ", paste(req_cols, collapse=", "))
      }
      # Extract top N genes and convert to named vector
      dat %>%
        dplyr::filter(avg_log2FC > 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice_max(order_by = avg_log2FC, n = n) %>%
        dplyr::summarise(genes = paste0(gene, collapse = ", "), .groups = "drop") %>%
        tibble::deframe()
    } else if (inherits(dat, "list")) {
      sapply(dat, function(x) paste0(head(x, n), collapse = ", "))
    } else {
      stop("Input must be a data.frame or a named list.")
    }
  }

  # ----------------------------------------------------------------------------
  # 3. Prompt Construction
  # ----------------------------------------------------------------------------
  make_prompt <- function() {
    role <- "You are an expert biologist and single-cell data analyst."

    # Build Context
    ctx <- glue::glue("Tissue: {tissuename}.")
    if (!is.null(cellname)) ctx <- paste(ctx, "Parent Cell Type:", cellname)
    if (nzchar(background)) ctx <- paste(ctx, "Experimental Context:", background)

    # Build Task Description
    task <- if (annotation_level == "major") {
      "Identify the major cell lineage (e.g., T cells, B cells, Epithelial)."
    } else {
      "Identify the specific subtype or functional state (e.g., CD8+ Exhausted T cells, Tregs)."
    }

    # Core Instructions (Enforcing JSON)
    glue::glue("
      {role}

      CONTEXT:
      {ctx}

      TASK:
      {task}

      INSTRUCTIONS:
      1. I will provide clusters with their top marker genes.
      2. Analyze the markers to determine the cell type.
      3. Provide 'Confidence' (High/Medium/Low).
      4. Provide 'Reasoning' (Which genes justify the decision?).
      5. OUTPUT STRICT JSON ONLY. No markdown blocks.

      JSON EXAMPLE:
      [
        {{
          \"cluster_id\": \"0\",
          \"cell_type\": \"Macrophages\",
          \"confidence\": \"High\",
          \"reasoning\": \"High expression of CD68 and CD163.\"
        }}
      ]
    ")
  }

  # ----------------------------------------------------------------------------
  # 4. API Execution Logic
  # ----------------------------------------------------------------------------
  run_batch <- function(batch_ids, batch_genes, client, sys_prompt) {
    # Construct User Input: "ClusterID: Genes..."
    input_str <- paste(batch_ids, ": ", batch_genes, collapse = "\n")
    user_msg <- glue::glue("Annotate these:\n{input_str}")

    res_df <- NULL

    # Retry Loop
    for (i in 1:retries) {
      if (verbose) message(sprintf("    -> Calling API (Attempt %d/%d)...", i, retries))

      out <- tryCatch({
        client$chat$completions$create(
          model = model,
          messages = list(
            list(role = "system", content = sys_prompt),
            list(role = "user", content = user_msg)
          ),
          temperature = 0.2, # Low temperature for consistent formatting
          response_format = list(type = "json_object")
        )
      }, error = function(e) e)

      # Handle Response
      if (!inherits(out, "error")) {
        # Parse JSON
        parse_res <- tryCatch({
          raw_json <- out$choices[[1]]$message$content
          clean_json <- gsub("```json|```", "", raw_json) # Clean Markdown if present
          json_data <- jsonlite::fromJSON(clean_json)

          # Handle nested lists (e.g., {"clusters": [...]})
          if (!is.data.frame(json_data)) {
            if (is.list(json_data) && length(json_data) == 1) json_data <- json_data[[1]]
          }
          # Validate structure
          if (is.data.frame(json_data)) {
            colnames(json_data) <- tolower(colnames(json_data))
            json_data
          } else {
            stop("JSON structure mismatch")
          }
        }, error = function(e) e)

        if (!inherits(parse_res, "error")) {
          res_df <- parse_res
          break # Success, exit loop
        } else {
          if (verbose) message("    !! JSON Parse Error, retrying...")
        }
      } else {
        if (verbose) message("    !! Network Error: ", out$message)
      }
      Sys.sleep(2) # Avoid rate limits
    }

    return(res_df)
  }

  # ----------------------------------------------------------------------------
  # 5. Main Execution
  # ----------------------------------------------------------------------------
  tryCatch({
    # A. Initialization
    check_deps()

    # API Key Handling
    key <- if (is.null(api_key)) Sys.getenv("OPENAI_API_KEY") else api_key
    if (!nzchar(key)) stop("API Key not found. Please provide 'api_key' argument or set OPENAI_API_KEY env variable.")

    client <- openai::OpenAI(api_key = key, base_url = base_url)

    if (verbose) message(glue::glue("=== Annotation Started: {tissuename} ({annotation_level}) ==="))

    # B. Process Inputs
    markers_vec <- format_markers(input, topgenenumber)
    sys_prompt <- make_prompt()
    cluster_ids <- names(markers_vec)

    # C. Batch Processing
    # Batch size limited to 5 to prevent token overflow or JSON truncation
    batch_size <- 5
    batches <- split(cluster_ids, ceiling(seq_along(cluster_ids)/batch_size))

    results_list <- list()

    for (b_idx in names(batches)) {
      curr_ids <- batches[[b_idx]]
      if (verbose) message(glue::glue("Processing Batch {b_idx} (Clusters: {paste(curr_ids, collapse=', ')})..."))

      batch_res <- run_batch(curr_ids, markers_vec[curr_ids], client, sys_prompt)

      # Fallback if batch fails completely
      if (is.null(batch_res)) {
        batch_res <- tibble::tibble(
          cluster_id = curr_ids,
          cell_type = "Failed",
          confidence = "None",
          reasoning = "API or Parse Error"
        )
      }
      results_list[[b_idx]] <- batch_res
    }

    # D. Final Assembly
    final_df <- dplyr::bind_rows(results_list) %>%
      dplyr::rename(
        Cluster = cluster_id,
        Prediction = cell_type,
        Confidence = confidence,
        Reasoning = reasoning
      ) %>%
      dplyr::mutate(Annotation_Level = annotation_level) %>%
      dplyr::select(Cluster, Prediction, Confidence, Reasoning, Annotation_Level)

    if (verbose) message("=== Annotation Complete ===")
    return(final_df)

  }, error = function(e) {
    stop("Fatal Error in AutoCellType: ", e$message)
  })
}
