# =============== OpenAI-compatible AI provider internals ================

.easyAI_uses_responses_endpoint <- function(base_url) {
  if (!.easyAI_has_value(base_url)) {
    return(FALSE)
  }
  url <- sub("[?#].*$", "", trimws(as.character(base_url)[1]))
  grepl("/responses/?$", url, ignore.case = TRUE)
}

.easyAI_strip_endpoint_url <- function(base_url) {
  url <- trimws(as.character(base_url)[1])
  url <- sub("[?#].*$", "", url)
  url <- sub("/+$", "", url)
  url <- sub("/chat/completions$", "", url, ignore.case = TRUE)
  url <- sub("/responses$", "", url, ignore.case = TRUE)
  url
}

.easyAI_chat_url <- function(base_url) {
  url <- trimws(as.character(base_url)[1])
  url <- sub("[?#].*$", "", url)
  url <- sub("/+$", "", url)
  if (grepl("/chat/completions$", url, ignore.case = TRUE)) {
    return(url)
  }
  paste0(.easyAI_strip_endpoint_url(url), "/chat/completions")
}

.easyAI_responses_url <- function(base_url) {
  url <- trimws(as.character(base_url)[1])
  url <- sub("[?#].*$", "", url)
  url <- sub("/+$", "", url)
  if (grepl("/responses$", url, ignore.case = TRUE)) {
    return(url)
  }
  paste0(.easyAI_strip_endpoint_url(url), "/responses")
}

.easyAI_call_provider_messages <- function(messages,
                                           provider,
                                           api_key,
                                           temperature = 0.1,
                                           response_format = NULL,
                                           max_tokens = NULL,
                                           call_chat = .easyAI_call_chat_messages,
                                           call_responses = .easyAI_call_responses_messages) {
  if (is.null(provider) || !inherits(provider, "easyAI_provider")) {
    stop("'provider' must be an easyAI_provider object.", call. = FALSE)
  }
  if (is.null(messages) || !is.list(messages) || length(messages) == 0) {
    stop("'messages' must be a non-empty list.", call. = FALSE)
  }

  order <- .easyAI_endpoint_order(provider)
  errors <- list()
  for (endpoint in order) {
    routed_provider <- provider
    routed_provider$endpoint <- endpoint
    value <- tryCatch(
      {
        if (identical(endpoint, "chat")) {
          call_chat(
            messages = messages,
            provider = routed_provider,
            api_key = api_key,
            temperature = temperature,
            response_format = response_format,
            max_tokens = max_tokens
          )
        } else {
          call_responses(
            messages = messages,
            provider = routed_provider,
            api_key = api_key,
            temperature = temperature,
            response_format = response_format,
            max_tokens = max_tokens
          )
        }
      },
      error = function(e) e
    )
    if (!inherits(value, "condition")) {
      return(value)
    }
    errors[[endpoint]] <- value
  }

  stop(.easyAI_endpoint_error_message(errors, order), call. = FALSE)
}

.easyAI_call_chat_messages <- function(messages,
                                       provider,
                                       api_key,
                                       temperature = 0.1,
                                       response_format = NULL,
                                       max_tokens = NULL) {
  body <- list(
    model = provider$model,
    messages = messages
  )
  body <- .easyAI_add_generation_options(
    body = body,
    provider = provider,
    endpoint = "chat",
    temperature = temperature,
    max_tokens = max_tokens
  )
  if (!is.null(response_format)) {
    body$response_format <- response_format
  }

  req <- httr2::request(.easyAI_chat_url(provider$base_url))
  req <- httr2::req_method(req, "POST")
  req <- httr2::req_headers(
    req,
    Authorization = paste("Bearer", api_key),
    Accept = "application/json"
  )
  req <- httr2::req_timeout(req, 60)
  req <- httr2::req_body_json(req, body, auto_unbox = TRUE)
  resp <- httr2::req_perform(req)
  parsed <- httr2::resp_body_json(resp, simplifyVector = FALSE)
  .easyAI_extract_response_text(parsed)
}

.easyAI_call_responses_messages <- function(messages,
                                            provider,
                                            api_key,
                                            temperature = 0.1,
                                            response_format = NULL,
                                            max_tokens = NULL) {
  formatted <- .easyAI_responses_messages(messages)
  body <- list(
    model = provider$model,
    input = formatted$input
  )
  if (.easyAI_has_value(formatted$instructions)) {
    body$instructions <- formatted$instructions
  }
  body <- .easyAI_add_generation_options(
    body = body,
    provider = provider,
    endpoint = "responses",
    temperature = temperature,
    max_tokens = max_tokens
  )
  if (!is.null(response_format)) {
    body$text <- .easyAI_responses_text_format(response_format)
  }

  req <- httr2::request(.easyAI_responses_url(provider$base_url))
  req <- httr2::req_method(req, "POST")
  req <- httr2::req_headers(
    req,
    Authorization = paste("Bearer", api_key),
    Accept = "application/json"
  )
  req <- httr2::req_timeout(req, 60)
  req <- httr2::req_body_json(req, body, auto_unbox = TRUE)
  resp <- httr2::req_perform(req)
  parsed <- httr2::resp_body_json(resp, simplifyVector = FALSE)
  .easyAI_extract_response_text(parsed)
}

.easyAI_is_reasoning_model <- function(model) {
  if (!.easyAI_has_value(model)) {
    return(FALSE)
  }
  model <- tolower(trimws(as.character(model)[1]))
  grepl("^(o[0-9]|gpt-5)", model, perl = TRUE)
}

.easyAI_add_generation_options <- function(body,
                                           provider,
                                           endpoint = c("chat", "responses"),
                                           temperature = NULL,
                                           max_tokens = NULL) {
  endpoint <- .easyAI_normalize_endpoint(endpoint[1])
  if (identical(endpoint, "auto")) {
    endpoint <- provider$endpoint %||% "chat"
  }
  is_reasoning <- .easyAI_is_reasoning_model(provider$model)

  if (identical(endpoint, "chat")) {
    if (!is.null(temperature) && !is_reasoning) {
      body$temperature <- temperature
    }
    if (!is.null(max_tokens)) {
      if (is_reasoning) {
        body$max_completion_tokens <- max_tokens
      } else {
        body$max_tokens <- max_tokens
      }
    }
  } else {
    if (!is.null(temperature)) {
      body$temperature <- temperature
    }
    if (!is.null(max_tokens)) {
      body$max_output_tokens <- max_tokens
    }
  }

  body
}

.easyAI_responses_messages <- function(messages) {
  roles <- vapply(messages, function(message) {
    role <- message$role %||% "user"
    tolower(as.character(role)[1])
  }, character(1))
  contents <- vapply(messages, function(message) {
    paste(as.character(message$content %||% ""), collapse = "\n")
  }, character(1))

  instructions <- paste(contents[roles %in% c("system", "developer") & nzchar(contents)], collapse = "\n\n")
  keep <- !(roles %in% c("system", "developer")) & nzchar(contents)
  input <- lapply(which(keep), function(i) {
    role <- roles[[i]]
    if (!role %in% c("user", "assistant")) {
      role <- "user"
    }
    list(role = role, content = contents[[i]])
  })
  if (length(input) == 0) {
    input <- list(list(
      role = "user",
      content = "Return the requested result using the instructions."
    ))
  }
  list(instructions = instructions, input = input)
}

.easyAI_responses_text_format <- function(response_format) {
  if (is.null(response_format)) {
    return(NULL)
  }
  fmt_type <- response_format$type %||% NULL
  if (.easyAI_has_value(fmt_type) && identical(tolower(as.character(fmt_type)[1]), "json_object")) {
    return(list(format = list(type = "json_object")))
  }
  if (.easyAI_has_value(fmt_type) && identical(tolower(as.character(fmt_type)[1]), "json_schema")) {
    out <- response_format
    names(out)[names(out) == "type"] <- "type"
    return(list(format = out))
  }
  list(format = response_format)
}

.easyAI_endpoint_error_message <- function(errors, order = names(errors)) {
  if (length(errors) == 0) {
    return("Model call failed before an endpoint was attempted.")
  }
  pieces <- vapply(order, function(endpoint) {
    err <- errors[[endpoint]]
    if (is.null(err)) {
      return(NULL)
    }
    status <- .easyAI_condition_status(err)
    prefix <- if (is.na(status)) {
      paste0(endpoint, " failed")
    } else {
      paste0(endpoint, " failed with HTTP ", status)
    }
    paste0(prefix, ": ", conditionMessage(err))
  }, character(1))
  paste(pieces[nzchar(pieces)], collapse = "\n")
}

.easyAI_condition_status <- function(error) {
  if (!inherits(error, "condition")) {
    return(NA_integer_)
  }
  resp <- error$response %||% error$resp
  if (!is.null(resp)) {
    status <- tryCatch(httr2::resp_status(resp), error = function(e) NA_integer_)
    if (!is.na(status)) {
      return(as.integer(status))
    }
  }
  msg <- conditionMessage(error)
  hit <- regmatches(msg, regexpr("\\bHTTP\\s+([0-9]{3})\\b", msg, perl = TRUE))
  if (length(hit) && nzchar(hit)) {
    return(as.integer(sub(".*HTTP\\s+([0-9]{3}).*", "\\1", hit, perl = TRUE)))
  }
  hit <- regmatches(msg, regexpr("\\b([0-9]{3})\\b", msg, perl = TRUE))
  if (length(hit) && nzchar(hit)) {
    return(as.integer(hit))
  }
  NA_integer_
}

.easyAI_extract_response_text <- function(resp) {
  if (is.null(resp)) {
    return("")
  }
  if (is.character(resp)) {
    return(paste(resp, collapse = "\n"))
  }
  if (!is.null(resp$output_text)) {
    return(paste(unlist(resp$output_text), collapse = "\n"))
  }
  if (!is.null(resp$choices) && length(resp$choices) > 0) {
    first <- resp$choices[[1]]
    if (!is.null(first$message$content)) {
      return(paste(unlist(first$message$content), collapse = "\n"))
    }
    if (!is.null(first$text)) {
      return(paste(unlist(first$text), collapse = "\n"))
    }
  }
  if (!is.null(resp$output) && length(resp$output) > 0) {
    pieces <- unlist(lapply(resp$output, function(item) {
      if (is.null(item$content)) {
        return(NULL)
      }
      unlist(lapply(item$content, function(content_item) {
        content_item$text %||% content_item$output_text %||% content_item$summary %||% NULL
      }))
    }))
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces) > 0) {
      return(paste(pieces, collapse = "\n"))
    }
  }
  jsonlite::toJSON(resp, auto_unbox = TRUE)
}
