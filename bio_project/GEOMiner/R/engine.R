.log <- function(msg) {
  cat(file = stderr(), sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}


#' @importFrom openai OpenAI
#' @importFrom httr RETRY status_code content
#' @importFrom rvest read_html html_elements html_text2
#' @importFrom xml2 xml_remove
#' @importFrom jsonlite fromJSON
#' @importFrom glue glue
#' @importFrom stringr str_match
#' @noRd

.process_single_gse <- function(gse_id, user_query, user_context, api_key, base_url, model, known_n = "?") {

  # 1. 爬取网页
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", trimws(gse_id))

  page_html <- tryCatch({
    resp <- httr::RETRY("GET", url, times = 3, quiet = TRUE, pause_base = 1)
    if (httr::status_code(resp) != 200) stop("HTTP Error")
    rvest::read_html(httr::content(resp, as = "text", encoding = "UTF-8"))
  }, error = function(e) return(NULL))

  if (is.null(page_html)) return(NULL)

  # 2. 预处理文本
  tryCatch({ xml2::xml_remove(rvest::html_elements(page_html, "head|script|style|nav|footer")) }, error = function(e) NULL)

  # 获取完整文本用于正则备用
  full_text <- gsub("\\s+", " ", rvest::html_text2(page_html))
  # 截取用于 LLM
  raw_text_llm <- substr(full_text, 1, 12000)

  # 3. 确定 N_Samples (优先使用传入的 known_n，其次正则，最后放弃 LLM 提取数字)
  final_n_samples <- as.character(known_n)

  # 如果 API 没传回来数字 (即 "?")，或者为 0，尝试用正则从网页里扣
  if (final_n_samples == "?" || final_n_samples == "0") {
    # GEO 网页通常包含 "Samples (24)" 这样的格式
    regex_match <- stringr::str_match(full_text, "Samples \\(\\s*(\\d+)\\s*\\)")
    if (!is.na(regex_match[1,2])) {
      final_n_samples <- regex_match[1,2]
    }
  }

  # 4. LLM 解析
  client <- openai::OpenAI(api_key = api_key, base_url = base_url)
  ctx_prompt <- if (nchar(user_context) > 2) glue::glue("USER CONTEXT constraints: '{user_context}'. \nSCORING: If contradicts Context, score < 50. If supports Context, score > 80.") else "Score based on keyword relevance."

  # Prompt
  sys_prompt <- glue::glue("
    TASK: Extract metadata for {gse_id}.
    KEYWORDS: {user_query}
    {ctx_prompt}
    RAW METADATA: {raw_text_llm}

    RULES:
    -Always use professional English not other language.
    - Return ONLY valid JSON. No markdown, no extra text.
    - All fields must be present. If unknown, use \"\".
    - match_score must be an integer 0-100.
    - score_reason must be <= 100 words and justify the score (relevance, sample size, clinical detail, organism).
    - scoring fully considering {user_query} and {ctx_prompt}
    MATCH SCORE GUIDELINES:
    - 90-100 (High): Perfect/near-perfect match; rich clinical/phenotype detail and/or large sample size; correct organism.
    - 70-89 (Medium): Topic matches but missing key clinical details, limited sample size, or partial mismatch.
    - 0-69 (Low): Irrelevant dataset, wrong organism, wrong modality, or insufficient info.

    OUTPUT JSON FORMAT:
    {{
      \"match_score\": 0,
      \"score_reason\": \"\",
      \"organism\": \"\",
      \"groups\": \"\",
      \"seq_strategy\": \"\",
      \"description\": \"\",
      \"abstract\": \"\",
      \"year\": \"\"
    }}
")

  meta <- NULL
  for (attempt in 1:2) {
    tryCatch({
      resp <- client$chat$completions$create(
        model = model,
        messages = list(list(role = "user", content = sys_prompt)),
        response_format = list(type = "json_object"),
        temperature = 0.2
      )
      meta <- jsonlite::fromJSON(resp$choices[[1]]$message$content)
      if (!is.null(meta$match_score)) break
    }, error = function(e) Sys.sleep(1))
  }

  if (is.null(meta)) {
    return(data.frame(GSE = gse_id, Match_Score = 0, Score_Reason = "LLM Failed", Organism = "-", Groups = "-", N_Samples = final_n_samples, Seq_Strategy = "-", Description = "-", Abstract = "-", Year = "-", URL = url, stringsAsFactors = FALSE))
  }

  # 5. 返回结果 (使用我们确定的 final_n_samples)
  data.frame(
    GSE = gse_id,
    Match_Score = as.integer(meta$match_score),
    Score_Reason = as.character(meta$score_reason),
    Organism = as.character(meta$organism),
    Groups = as.character(meta$groups),
    N_Samples = final_n_samples,
    Seq_Strategy = as.character(meta$seq_strategy),
    Description = as.character(meta$description),
    Abstract = as.character(meta$abstract),
    Year = as.character(meta$year),
    URL = url,
    stringsAsFactors = FALSE
  )
}


#' @noRd
.ai_generate_query <- function(keyword, api_key, base_url, model) {
  .log(paste("Generating query for:", keyword))
  client <- openai::OpenAI(api_key = api_key, base_url = base_url)

  prompt <- glue::glue("
    ROLE: Senior Data Curator for NCBI GEO.
    TASK: Translate User Input into a precise NCBI Entrez Search Query.
    USER INPUT: '{keyword}'
    RULES:
        1. Identify core entities.
        2. Expand with boolean OR.
        3. Combine with boolean AND.
        4. MANDATORY: Append 'AND gse[Entry Type]'.
        5. OUTPUT: Raw query string only.
        6. Always use professional English not other language.
        7. Parse all the keywords of the user provide, especially when user
        provide the underline/space/comma to split.
    OUTPUT: Raw query string only.
  ")

  tryCatch({
    resp <- client$chat$completions$create(
      model = model,
      messages = list(list(role = "user", content = prompt)),
      temperature = 0.3
    )
    q <- gsub("^\"|\"$|`", "", resp$choices[[1]]$message$content)
    clean_q <- trimws(q)
    .log(paste("Query generated:", clean_q))
    return(clean_q)
  }, error = function(e) {
    .log(paste("Query gen failed:", e$message))
    return(NULL)
  })
}



#' 核心挖掘引擎
#' @importFrom rentrez entrez_search entrez_summary
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor
#' @importFrom parallelly availableCores
#' @import dplyr
#' @noRd
geo_mine_engine <- function(keyword, context, limit, n_cores, api_key, base_url, model, update_fn = NULL) {

  # 1. Gen Query
  if (!is.null(update_fn)) update_fn("Generating AI Query...")
  final_query <- .ai_generate_query(keyword, api_key, base_url, model)
  if (is.null(final_query)) stop("Query generation failed.")

  # 2. Search NCBI
  if (!is.null(update_fn)) update_fn("Searching NCBI Database...")
  search_res <- tryCatch(
    rentrez::entrez_search(db="gds", term=final_query, retmax=max(100, limit*5), use_history=TRUE),
    error = function(e) stop(paste("NCBI Error:", e$message))
  )
  if (length(search_res$ids) == 0) return(NULL)

  # 3. Filter IDs AND Extract N_Samples
  if (!is.null(update_fn)) update_fn("Filtering Valid GSE IDs...")

  target_gses <- c()
  gse_n_map <- list() # 用来存储 ID -> 样本数 的映射

  idx <- 1; chunk_size <- 50 # 稍微调小 chunk 以防超时
  while(length(target_gses) < limit && idx <= length(search_res$ids)) {
    end_idx <- min(idx + chunk_size - 1, length(search_res$ids))
    ids_batch <- search_res$ids[idx:end_idx]

    summ <- tryCatch(rentrez::entrez_summary(db="gds", id=ids_batch), error=function(e) NULL)

    if (!is.null(summ)) {
      # rentrez 返回单个结果时不是 list，需要转换
      s_list <- if(!is.null(summ$uid)) list(summ) else summ

      for (item in s_list) {
        acc <- item$accession
        # 只取 GSE
        if (grepl("^GSE", acc)) {
          if (!(acc %in% target_gses)) {
            target_gses <- c(target_gses, acc)

            # 提取样本数 (n_samples)
            ns <- "?"
            if (!is.null(item$n_samples)) ns <- as.character(item$n_samples)
            gse_n_map[[acc]] <- ns
          }
        }
        if (length(target_gses) >= limit) break
      }
    }
    idx <- idx + chunk_size
    # 安全中断
    if(idx/chunk_size > 10) break
  }

  if(length(target_gses) == 0) return(NULL)

  # 4. 执行分析
  allowed_cores <- parallelly::availableCores()
  safe_cores <- min(as.integer(n_cores), as.integer(allowed_cores))
  results_list <- list()

  if (safe_cores > 1) {
    # --- 多核模式 ---
    .log(paste("Mode: Parallel (", safe_cores, "cores)"))
    p <- progressr::progressor(along = target_gses)
    plan(multisession, workers = safe_cores)

    results_list <- future_lapply(target_gses, function(id) {
      p(message = sprintf("Analyzing %s", id))

      # 获取该 ID 对应的样本数
      known_n <- if(!is.null(gse_n_map[[id]])) gse_n_map[[id]] else "?"

      # 调用函数 (必须传 known_n)
      .process_single_gse(id, keyword, context, api_key, base_url, model, known_n)
    },
    future.seed = TRUE,
    future.packages = c("httr", "rvest", "openai", "jsonlite", "glue", "xml2", "stringr"),
    future.chunk.size = 1
    )
    plan(sequential)

  } else {
    # --- 单核模式 (带 UI 回调) ---
    .log("Mode: Single Core")
    for (i in seq_along(target_gses)) {
      id <- target_gses[i]

      # 更新 UI
      msg <- sprintf("Analyzing %s (%d/%d)...", id, i, length(target_gses))
      if (!is.null(update_fn)) {
        update_fn(msg)
        Sys.sleep(0.1) # 让 Shiny 喘口气
      }

      # 获取该 ID 对应的样本数
      known_n <- if(!is.null(gse_n_map[[id]])) gse_n_map[[id]] else "?"

      # 调用函数 (必须传 known_n)
      res <- .process_single_gse(id, keyword, context, api_key, base_url, model, known_n)
      results_list[[i]] <- res
    }
  }

  # 5. 合并
  final_df <- do.call(rbind, Filter(Negate(is.null), results_list))
  if (!is.null(final_df)) {
    final_df <- final_df %>% arrange(desc(Match_Score))
    attr(final_df, "generated_query") <- final_query
  }
  return(final_df)
}
