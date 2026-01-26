#' @title AutoCellType
#' @description
#' A high-precision, parallelized single-cell annotation engine (v7.0).
#'
#' **Architecture Shift**:
#' 1. **Single-Shot Strategy**: Processes clusters strictly **one by one**. This eliminates "Model Laziness" (skipping clusters in a batch) and ensures maximum reasoning quality per cluster.
#' 2. **User-Controlled Parallelism**: Uses multi-threading (`n_cores`) to execute single-shot requests concurrently, maintaining high throughput.
#' 3. **Seurat Adapter**: Auto-detects V3/V4/V5 formats and fixes sorting bugs.
#' 4. **Zero-NA Guarantee**: Automatically fills missing subtype fields with lineage information.
#'
#' @param input **Required**. The input marker data. Supports two formats:
#'   \itemize{
#'     \item \strong{Data Frame}: Output from Seurat's \code{FindAllMarkers()}. Must contain \code{cluster} and \code{gene} (or rownames) columns, plus \code{avg_log2FC} (or \code{avg_logFC}).
#'     \item \strong{List}: A named list of marker vectors, e.g., \code{list("Cluster0" = c("GeneA", "GeneB"), ...)}.
#'   }
#' @param tissuename **Character**. The tissue of origin (e.g., "Human Liver", "Mouse Brain"). Serves as the primary biological context.
#' @param n_cores **Integer**. Number of parallel threads.
#'   \itemize{
#'     \item \code{1} (Default): Sequential processing (safest for debugging).
#'     \item \code{>1}: Parallel processing (Recommended: 4-8). Requires high API rate limits.
#'   }
#' @param cell_ontology **Constraint Mode**. Controls the standardization of the 'primary_lineage' field:
#'   \itemize{
#'     \item \code{NULL} (Default): Uses the built-in **Universal Ontology** (~40 types).
#'     \item \code{Character Vector}: **Recommended for Sub-clustering**. A user-defined list of allowed subtypes.
#'   }
#' @param prior_info **Character** (Optional). External biological context. E.g., "T cell sub-clustering".
#' @param model **Character**. LLM model name. Default \code{"gpt-4o-mini"}. \code{"gpt-4o"} is recommended for sub-clustering.
#' @param topgenenumber **Integer**. The number of top markers (sorted by LogFC) to send per cluster. Default: 20.
#' @param p_val_thresh **Numeric**. Significance threshold. Default: 0.05.
#' @param base_url **Character**. API endpoint. Default: "https://api.gpt.ge/v1".
#' @param api_key **Character**. OpenAI-compatible API Key. Defaults to \code{Sys.getenv("OPENAI_API_KEY")}.
#' @param retries **Integer**. Max retry attempts per cluster. Default: 3.
#' @param verbose **Logical**. Print progress bar. Default: \code{TRUE}.
#'
#' @return A \code{tibble} containing standardized columns: \code{cluster_id}, \code{primary_lineage}, \code{detailed_subtype}, \code{functional_state}, \code{confidence}, \code{reasoning}.
#'
#' @examples
#' \dontrun{
#'   # =========================================================================
#'   # Example 1: Standard Annotation (Sequential)
#'   # =========================================================================
#'   library(dplyr)
#'
#'   # 1. Create Mock Seurat Output
#'   # Cluster 0: T cells (CD3D, CD8A)
#'   # Cluster 1: B cells (MS4A1, CD79A)
#'   # Cluster 2: Fibroblasts (COL1A1, DCN)
#'   mock_markers <- data.frame(
#'     cluster = c(rep("0", 3), rep("1", 3), rep("2", 3)),
#'     gene = c("CD3D", "CD8A", "GZMK",      # C0
#'              "MS4A1", "CD79A", "CD19",    # C1
#'              "COL1A1", "DCN", "PDGFRA"),  # C2
#'     avg_log2FC = c(2.5, 2.3, 1.8, 2.4, 2.2, 2.0, 3.1, 2.8, 2.5),
#'     p_val_adj = rep(0.001, 9)
#'   )
#'
#'   # 2. Run
#'   res <- AutoCellType(
#'     input = mock_markers,
#'     tissuename = "Human PBMC",
#'     n_cores = 1, # Sequential
#'     api_key = "sk-..."
#'   )
#'   print(res)
#'
#'   # =========================================================================
#'   # Example 2: High-Performance Parallel Mode
#'   # =========================================================================
#'   # If you have 20+ clusters, use n_cores = 4 or 8 to speed up by 4x-8x.
#'
#'   res_parallel <- AutoCellType(
#'     input = mock_markers,
#'     tissuename = "Human Tissue",
#'     n_cores = 4, # Use 4 threads
#'     model = "gpt-4o-mini",
#'     api_key = "sk-..."
#'   )
#'
#'   # =========================================================================
#'   # Example 3: Sub-clustering (Custom Ontology + Context)
#'   # =========================================================================
#'   # Scenario: You have re-clustered T cells and need to distinguish subtypes.
#'
#'   t_ontology <- c("CD8 Naive", "CD8 Effector", "CD8 Exhausted", "Treg")
#'
#'   res_sub <- AutoCellType(
#'     input = mock_markers, # Assume these are T cell sub-clusters
#'     tissuename = "Human Liver",
#'     prior_info = "Sub-clustering of CD3+ T cells.", # Inject context
#'     cell_ontology = t_ontology, # Enforce specific names
#'     model = "gpt-4o", # Better model for subtle differences
#'     n_cores = 4
#'   )
#' }
#' @export
AutoCellType <- function(
    input,
    tissuename,
    cell_ontology = NULL,
    prior_info = NULL,
    n_cores = 1,
    model = "gpt-4o-mini",
    topgenenumber = 20,
    p_val_thresh = 0.05,
    base_url = "https://api.gpt.ge/v1",
    api_key = NULL,
    retries = 3,
    verbose = TRUE
) {

  # ============================================================================
  # 1. Configuration & Ontology
  # ============================================================================
  STANDARD_COLS <- c("cluster_id", "primary_lineage", "detailed_subtype", "functional_state", "confidence", "reasoning")

  default_ontology <- c(
    # Immune
    "T cell", "B cell", "Plasma cell", "NK cell", "Myeloid cell", "Macrophage", "Monocyte",
    "Dendritic cell", "Neutrophil", "Mast cell", "Eosinophil", "Basophil",
    # Stromal
    "Fibroblast", "Mesenchymal cell", "Smooth muscle cell", "Pericyte", "Adipocyte",
    "Mesothelial cell", "Chondrocyte", "Osteoblast", "Stromal cell",
    # Epithelial/Endothelial
    "Epithelial cell", "Endothelial cell", "Lymphatic endothelial cell", "Keratinocyte",
    "Hepatocyte", "Cholangiocyte", "Pneumocyte", "Podocyte", "Mesangial cell",
    # Neural
    "Neuron", "Glial cell", "Schwann cell", "Cardiomyocyte",
    # Embryonic
    "Stem cell", "Progenitor cell", "Embryonic stem cell", "Trophoblast",
    "Germ cell", "Erythrocyte", "Platelet", "Melanocyte"
  )

  current_ontology <- if (is.character(cell_ontology)) cell_ontology else default_ontology
  should_enforce <- !isFALSE(cell_ontology)

  # ============================================================================
  # 2. Dependency Check
  # ============================================================================
  check_deps <- function() {
    pkgs <- c("openai", "dplyr", "glue", "tibble", "jsonlite", "stringr", "future", "future.apply", "progressr")
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) stop("AutoCellType requires missing packages: ", paste(missing, collapse = ", "))
  }

  # ============================================================================
  # 3. Data Formatting
  # ============================================================================
  format_markers <- function(dat, n, p_thresh) {
    if (inherits(dat, "data.frame")) {
      cols <- colnames(dat)
      if (!"gene" %in% cols) {
        if ("feature" %in% cols) dat <- dplyr::rename(dat, gene = feature)
        else dat$gene <- rownames(dat)
      }
      if (!"avg_log2FC" %in% cols && "avg_logFC" %in% cols) dat <- dplyr::rename(dat, avg_log2FC = avg_logFC)

      # Strict numeric conversion
      dat$avg_log2FC <- as.numeric(dat$avg_log2FC)
      if("p_val_adj" %in% cols) dat$p_val_adj <- as.numeric(dat$p_val_adj)

      dat_clean <- dat |>
        dplyr::filter(if("p_val_adj" %in% colnames(dat)) p_val_adj < p_thresh else TRUE) |>
        dplyr::filter(avg_log2FC > 0) |>
        dplyr::group_by(cluster) |>
        dplyr::slice_max(order_by = avg_log2FC, n = n) |>
        dplyr::summarise(genes = paste0(gene, collapse = ", "), .groups = "drop") |>
        dplyr::mutate(cluster = as.character(cluster))

      return(tibble::deframe(dat_clean))
    } else if (inherits(dat, "list")) {
      sapply(dat, function(x) paste0(head(x, n), collapse = ", "))
    } else stop("Input must be a Data Frame or a List.")
  }

  # ============================================================================
  # 4. Prompt Engineering (Single Shot)
  # ============================================================================
  make_system_prompt <- function() {
    ctx <- glue::glue("Tissue Origin: {tissuename}")
    if (!is.null(prior_info)) ctx <- paste(ctx, glue::glue("Context: {prior_info}"), sep = "\n")

    ontology_str <- if(should_enforce) paste(paste0("- ", current_ontology), collapse = "\n") else "No strict list."

    instruction <- if(should_enforce) {
      glue::glue("For 'primary_lineage', SELECT EXACTLY from:\n{ontology_str}")
    } else "Use standard scientific names."

    glue::glue("
      You are an expert Cell Biologist. Annotate the SINGLE cluster provided.

      --- CONTEXT ---
      {ctx}
      {instruction}

      --- RULES ---
      1. Analyze the markers carefully.
      2. 'detailed_subtype' is MANDATORY. Do not leave it empty.
      3. Return a JSON Object (not array) for this single cluster.

      --- OUTPUT FORMAT ---
      {{
        \"cluster_id\": \"id\",
        \"primary_lineage\": \"Standard Name\",
        \"detailed_subtype\": \"Specific Subtype\",
        \"functional_state\": \"State\",
        \"confidence\": \"High/Medium/Low\",
        \"reasoning\": \"Logic\"
      }}
    ")
  }

  # ============================================================================
  # 5. Utilities
  # ============================================================================
  standardize_one_row <- function(df, expected_cols) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    colnames(df) <- tolower(colnames(df))

    missing <- setdiff(expected_cols, colnames(df))
    if (length(missing) > 0) for (col in missing) df[[col]] <- NA

    # N/A Killer
    if ("detailed_subtype" %in% colnames(df) && "primary_lineage" %in% colnames(df)) {
      df$detailed_subtype <- ifelse(
        is.na(df$detailed_subtype) | df$detailed_subtype == "" | df$detailed_subtype == "N/A",
        df$primary_lineage, df$detailed_subtype
      )
    }
    return(df[, expected_cols, drop = FALSE])
  }

  align_lineage <- function(df, ontology_list) {
    if (!should_enforce || is.null(df) || nrow(df) == 0) return(df)

    val <- df$primary_lineage[1]
    if (is.na(val)) { df$primary_lineage <- "Unknown"; return(df) }

    if (val %in% ontology_list) return(df)

    idx <- match(gsub("s$", "", tolower(val)), gsub("s$", "", tolower(ontology_list)))
    if (!is.na(idx)) df$primary_lineage <- ontology_list[idx]

    return(df)
  }

  # ============================================================================
  # 6. Single Cluster Runner
  # ============================================================================
  run_single_cluster <- function(cid, cgenes, client, sys_prompt, retry_limit) {
    user_msg <- glue::glue("Cluster ID: {cid}\nTop Markers: {cgenes}")
    last_error <- "Unknown Error"

    for (i in 1:retry_limit) {
      if (i > 1) Sys.sleep(2^(i-1))

      api_res <- tryCatch({
        client$chat$completions$create(
          model = model,
          messages = list(list(role = "system", content = sys_prompt), list(role = "user", content = user_msg)),
          temperature = 0.1, max_tokens = 1000, response_format = list(type = "json_object")
        )
      }, error = function(e) list(error = e$message))

      if (!is.null(api_res$error)) { last_error <- paste("API Error:", api_res$error); next }

      df <- tryCatch({
        raw <- api_res$choices[[1]]$message$content
        json_match <- stringr::str_extract(raw, "(?s)\\{.*\\}")
        target <- if (!is.na(json_match)) json_match else raw
        parsed <- jsonlite::fromJSON(target)

        if (is.list(parsed) && !is.data.frame(parsed)) dplyr::bind_rows(parsed)
        else if (is.data.frame(parsed)) parsed
        else NULL
      }, error = function(e) {
        if(i == retry_limit) last_error <<- paste("JSON Error:", e$message)
        NULL
      })

      if (!is.null(df) && nrow(df) > 0) {
        df <- standardize_one_row(df, STANDARD_COLS)
        df$cluster_id <- as.character(cid)
        df <- align_lineage(df, current_ontology)
        return(df)
      }
    }

    return(tibble::tibble(
      cluster_id = as.character(cid),
      primary_lineage = "Failed",
      detailed_subtype = last_error,
      functional_state = NA, confidence = "None", reasoning = NA
    ))
  }

  # ============================================================================
  # 7. Main Execution (Parallel/Sequential)
  # ============================================================================
  tryCatch({
    check_deps()
    key <- if (is.null(api_key)) Sys.getenv("OPENAI_API_KEY") else api_key
    if (!nzchar(key)) stop("API Key not found.")
    client <- openai::OpenAI(api_key = key, base_url = base_url)

    markers_vec <- format_markers(input, topgenenumber, p_val_thresh)
    if (length(markers_vec) == 0) stop("No valid markers.")

    sys_prompt <- make_system_prompt()
    cluster_ids <- names(markers_vec)
    n_total <- length(cluster_ids)

    mode_msg <- if (n_cores > 1) glue::glue("Parallel ({n_cores} cores)") else "Sequential"
    if (verbose) message(glue::glue("=== AutoCellType v7.0 | {n_total} Clusters | {mode_msg} ==="))

    worker_fun <- function(cid) {
      run_single_cluster(cid, markers_vec[cid], client, sys_prompt, retries)
    }

    if (n_cores > 1) {
      future::plan(future::multisession, workers = n_cores)
      res_list <- progressr::with_progress({
        p <- progressr::progressor(steps = n_total)
        future.apply::future_lapply(cluster_ids, function(cid) {
          p(sprintf("Cluster %s", cid))
          worker_fun(cid)
        }, future.seed = TRUE)
      })
    } else {
      res_list <- purrr::map(cluster_ids, function(cid) {
        if(verbose) message(glue::glue("Processing Cluster {cid}..."))
        worker_fun(cid)
      })
    }

    final_df <- dplyr::bind_rows(res_list) |> dplyr::as_tibble()
    if (verbose) message("=== Done ===")
    return(final_df)

  }, error = function(e) stop("Error: ", e$message))
}
