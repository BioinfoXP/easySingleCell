# =============== easySingleCell AI defaults ================

# Internal helpers for maintaining OpenAI-compatible model settings in one
# place. Update this file when the package default model, API environment
# variable, or compatible base URL changes.
.esc_ai_defaults <- function() {
  list(
    api_env = "OPENAI_API_KEY",
    base_url = "https://api.gpt.ge/v1",
    endpoint = "auto",
    models = list(
      general = "gpt-5.5",
      annotation = "gpt-4o-mini"
    )
  )
}

.esc_ai_api_key <- function(api_key = NULL) {
  if (!is.null(api_key) && length(api_key) == 1 && nzchar(api_key)) {
    return(api_key)
  }
  Sys.getenv(.esc_ai_defaults()$api_env)
}

.esc_ai_base_url <- function(base_url = NULL) {
  .easyAI_normalize_base_url(base_url)
}

.esc_ai_model <- function(task = "general", model = NULL) {
  if (!is.null(model) && length(model) == 1 && nzchar(model)) {
    return(model)
  }
  defaults <- .esc_ai_defaults()$models
  task <- as.character(task)[1]
  if (task %in% names(defaults)) {
    return(defaults[[task]])
  }
  defaults$general
}

.easyAI_has_value <- function(x) {
  !is.null(x) && length(x) == 1 && !is.na(x) && nzchar(trimws(as.character(x)))
}

.easyAI_normalize_base_url <- function(base_url = NULL) {
  if (!.easyAI_has_value(base_url)) {
    return(.esc_ai_defaults()$base_url)
  }

  url <- trimws(as.character(base_url)[1])
  url <- sub("/+$", "", url)
  url_no_query <- sub("[?#].*$", "", url)
  if (grepl("^https?://[^/]+$", url_no_query, ignore.case = TRUE)) {
    return(paste0(url_no_query, "/v1"))
  }
  url_no_query
}

.easyAI_normalize_endpoint <- function(endpoint = NULL) {
  if (!.easyAI_has_value(endpoint)) {
    return(.esc_ai_defaults()$endpoint)
  }
  endpoint <- tolower(trimws(as.character(endpoint)[1]))
  aliases <- list(
    auto = c("auto", "default"),
    chat = c("chat", "completion", "completions", "chat.completions", "chat_completion", "chat-completions"),
    responses = c("responses", "response", "responses-api", "response-api")
  )
  for (name in names(aliases)) {
    if (endpoint %in% aliases[[name]]) {
      return(name)
    }
  }
  stop("'endpoint' must be one of 'auto', 'chat', or 'responses'.", call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.easyAI_resolve_settings <- function(api_key = NULL,
                                     model = NULL,
                                     base_url = NULL,
                                     endpoint = NULL,
                                     task = "general",
                                     config = list()) {
  defaults <- .esc_ai_defaults()
  env_key <- Sys.getenv(defaults$api_env)
  config <- list()

  key_from_config <- FALSE
  key_source <- "missing"
  if (.easyAI_has_value(api_key)) {
    resolved_key <- trimws(as.character(api_key)[1])
    key_source <- "argument"
  } else if (.easyAI_has_value(env_key)) {
    resolved_key <- trimws(as.character(env_key)[1])
    key_source <- defaults$api_env
  } else {
    resolved_key <- ""
  }

  base_from_config <- FALSE
  if (.easyAI_has_value(base_url)) {
    resolved_base <- .easyAI_normalize_base_url(base_url)
    base_source <- "argument"
  } else {
    resolved_base <- .easyAI_normalize_base_url(defaults$base_url)
    base_source <- "package default"
  }

  model_from_config <- FALSE
  if (.easyAI_has_value(model)) {
    resolved_model <- trimws(as.character(model)[1])
    model_source <- "argument"
  } else {
    resolved_model <- .esc_ai_model(task, NULL)
    model_source <- "package default"
  }

  endpoint_from_config <- FALSE
  if (.easyAI_has_value(endpoint)) {
    resolved_endpoint <- .easyAI_normalize_endpoint(endpoint)
    endpoint_source <- "argument"
  } else {
    resolved_endpoint <- defaults$endpoint
    endpoint_source <- "package default"
  }

  list(
    api_key = resolved_key,
    base_url = resolved_base,
    model = resolved_model,
    endpoint = resolved_endpoint,
    source = list(
      api_key = key_source,
      base_url = base_source,
      model = model_source,
      endpoint = endpoint_source
    ),
    from_config = list(
      api_key = key_from_config,
      base_url = base_from_config,
      model = model_from_config,
      endpoint = endpoint_from_config
    )
  )
}

.easyAI_build_provider <- function(api_key = NULL,
                                   model = NULL,
                                   base_url = NULL,
                                   endpoint = NULL,
                                   task = "general",
                                   require_key = TRUE,
                                   config = list()) {
  settings <- .easyAI_resolve_settings(
    api_key = api_key,
    model = model,
    base_url = base_url,
    endpoint = endpoint,
    task = task,
    config = config
  )
  provider <- .easyAI_provider(settings)
  if (isTRUE(require_key) && !.easyAI_has_value(provider$api_key)) {
    stop(
      "Please provide 'api_key' or set OPENAI_API_KEY.",
      call. = FALSE
    )
  }
  provider
}

.easyAI_provider <- function(settings, name = "openai-compatible") {
  if (is.null(settings) || !is.list(settings)) {
    stop("'settings' must be a list returned by .easyAI_resolve_settings().", call. = FALSE)
  }

  base_url <- if (.easyAI_has_value(settings$base_url)) {
    .easyAI_normalize_base_url(settings$base_url)
  } else {
    .easyAI_normalize_base_url(.esc_ai_defaults()$base_url)
  }
  uses_responses <- .easyAI_uses_responses_endpoint(base_url)
  endpoint_preference <- .easyAI_normalize_endpoint(settings$endpoint)
  endpoint <- if (endpoint_preference == "auto" && uses_responses) {
    "responses"
  } else if (endpoint_preference == "auto") {
    "auto"
  } else {
    endpoint_preference
  }

  structure(
    list(
      name = name,
      model = if (.easyAI_has_value(settings$model)) trimws(as.character(settings$model)[1]) else .esc_ai_defaults()$models$general,
      base_url = base_url,
      api_key = if (.easyAI_has_value(settings$api_key)) trimws(as.character(settings$api_key)[1]) else "",
      endpoint = endpoint,
      endpoint_preference = endpoint_preference,
      uses_responses = uses_responses,
      temperature = 0.1,
      source = settings$source %||% list(),
      from_config = settings$from_config %||% list()
    ),
    class = "easyAI_provider"
  )
}

.easyAI_provider_label <- function(provider) {
  if (is.null(provider) || !inherits(provider, "easyAI_provider")) {
    return("provider: unknown")
  }
  key_status <- if (.easyAI_has_value(provider$api_key)) "set" else "missing"
  paste0(
    provider$name,
    " | model = ", provider$model,
    " | endpoint = ", provider$endpoint,
    " (", provider$endpoint_preference %||% "auto", ")",
    " | base_url = ", provider$base_url,
    " | key = ", key_status
  )
}

.easyAI_endpoint_order <- function(provider) {
  if (is.null(provider) || !inherits(provider, "easyAI_provider")) {
    stop("'provider' must be an easyAI_provider object.", call. = FALSE)
  }

  preference <- .easyAI_normalize_endpoint(provider$endpoint_preference %||% provider$endpoint)
  if (identical(preference, "chat")) {
    return("chat")
  }
  if (identical(preference, "responses")) {
    return("responses")
  }
  if (isTRUE(provider$uses_responses) || identical(provider$endpoint, "responses")) {
    return(c("responses", "chat"))
  }
  c("chat", "responses")
}
