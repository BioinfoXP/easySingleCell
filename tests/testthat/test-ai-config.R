test_that("AI defaults are maintained through one internal configuration layer", {
  defaults <- .esc_ai_defaults()

  expect_equal(defaults$api_env, "OPENAI_API_KEY")
  expect_equal(defaults$models$general, "gpt-5.5")
  expect_equal(defaults$endpoint, "auto")
  expect_equal(.esc_ai_base_url(), defaults$base_url)
  expect_equal(.esc_ai_model("general"), defaults$models$general)
  expect_equal(.esc_ai_model("annotation"), defaults$models$annotation)
})

test_that("AI API key resolver prefers explicit keys and then environment keys", {
  old_key <- Sys.getenv("OPENAI_API_KEY", unset = NA)
  on.exit({
    if (is.na(old_key)) {
      Sys.unsetenv("OPENAI_API_KEY")
    } else {
      Sys.setenv(OPENAI_API_KEY = old_key)
    }
  }, add = TRUE)

  Sys.setenv(OPENAI_API_KEY = "sk-env-test")

  expect_equal(.esc_ai_api_key("sk-explicit-test"), "sk-explicit-test")
  expect_equal(.esc_ai_api_key(NULL), "sk-env-test")
})

test_that("AI settings resolver uses explicit values, environment key, and package defaults", {
  with_no_openai_key({
    resolved_default <- .easyAI_resolve_settings(
      api_key = NULL,
      model = NULL,
      base_url = NULL,
      endpoint = NULL
    )
    expect_equal(resolved_default$api_key, "")
    expect_equal(resolved_default$base_url, .esc_ai_defaults()$base_url)
    expect_equal(resolved_default$model, .esc_ai_defaults()$models$general)
    expect_equal(resolved_default$endpoint, "auto")
    expect_false(any(unlist(resolved_default$from_config)))

    Sys.setenv(OPENAI_API_KEY = "sk-env")
    resolved_env <- .easyAI_resolve_settings(
      api_key = NULL,
      model = NULL,
      base_url = NULL,
      endpoint = NULL
    )
    expect_equal(resolved_env$api_key, "sk-env")
    expect_equal(resolved_env$source$api_key, "OPENAI_API_KEY")

    resolved_explicit <- .easyAI_resolve_settings(
      api_key = "sk-explicit",
      model = "gpt-5.5",
      base_url = "https://explicit.example",
      endpoint = "chat"
    )
    expect_equal(resolved_explicit$api_key, "sk-explicit")
    expect_equal(resolved_explicit$model, "gpt-5.5")
    expect_equal(resolved_explicit$base_url, "https://explicit.example/v1")
    expect_equal(resolved_explicit$endpoint, "chat")
    expect_equal(resolved_explicit$source$model, "argument")
  })
})

test_that("easyAI recognizes chat and responses style base URLs", {
  expect_equal(.easyAI_normalize_endpoint(NULL), "auto")
  expect_equal(.easyAI_normalize_endpoint("completion"), "chat")
  expect_equal(.easyAI_normalize_endpoint("chat.completions"), "chat")
  expect_equal(.easyAI_normalize_endpoint("response"), "responses")
  expect_error(.easyAI_normalize_endpoint("bad-endpoint"), "endpoint")

  expect_equal(.easyAI_normalize_base_url("https://api.gpt.ge"), "https://api.gpt.ge/v1")
  expect_equal(.easyAI_normalize_base_url("https://api.gpt.ge/"), "https://api.gpt.ge/v1")
  expect_equal(.easyAI_normalize_base_url("https://api.gpt.ge/v1/"), "https://api.gpt.ge/v1")
  expect_equal(
    .easyAI_normalize_base_url("https://lucen.cc/v1/responses/"),
    "https://lucen.cc/v1/responses"
  )
  expect_false(.easyAI_uses_responses_endpoint("https://api.example/v1"))
  expect_true(.easyAI_uses_responses_endpoint("https://api.example/v1/responses"))
  expect_true(.easyAI_uses_responses_endpoint("https://api.example/v1/responses/"))

  expect_equal(
    .easyAI_responses_url("https://api.example/v1/responses/"),
    "https://api.example/v1/responses"
  )
  expect_equal(
    .easyAI_responses_url("https://api.example/v1"),
    "https://api.example/v1/responses"
  )
})

test_that("AI provider resolves OpenAI-compatible metadata", {
  with_no_openai_key({
    settings <- .easyAI_resolve_settings(
      api_key = "sk-provider",
      model = "gpt-5.5",
      base_url = "https://provider.example/v1/responses",
      endpoint = "auto"
    )
    provider <- .easyAI_provider(settings)

    expect_s3_class(provider, "easyAI_provider")
    expect_equal(provider$name, "openai-compatible")
    expect_equal(provider$model, "gpt-5.5")
    expect_equal(provider$api_key, "sk-provider")
    expect_equal(provider$base_url, "https://provider.example/v1/responses")
    expect_equal(provider$endpoint, "responses")
    expect_equal(provider$endpoint_preference, "auto")
    expect_true(provider$uses_responses)
    expect_equal(provider$source$model, "argument")
    expect_match(.easyAI_provider_label(provider), "openai-compatible")
    expect_match(.easyAI_provider_label(provider), "responses")
  })
})

test_that("easyAI endpoint ordering follows aisdk-style provider preferences", {
  expect_equal(
    .easyAI_endpoint_order(.easyAI_provider(list(
      api_key = "sk-test",
      base_url = "https://api.example/v1",
      model = "gpt-test",
      endpoint = "auto"
    ))),
    c("chat", "responses")
  )

  expect_equal(
    .easyAI_endpoint_order(.easyAI_provider(list(
      api_key = "sk-test",
      base_url = "https://api.example/v1/responses",
      model = "gpt-test",
      endpoint = "auto"
    ))),
    c("responses", "chat")
  )

  expect_equal(
    .easyAI_endpoint_order(.easyAI_provider(list(
      api_key = "sk-test",
      base_url = "https://api.example/v1",
      model = "gpt-test",
      endpoint = "responses"
    ))),
    "responses"
  )
})

test_that("easyAI generation options adapt to reasoning-style model APIs", {
  chat_provider <- .easyAI_provider(list(
    api_key = "sk-test",
    base_url = "https://api.example/v1",
    model = "gpt-4o-mini",
    endpoint = "chat"
  ))
  reasoning_provider <- .easyAI_provider(list(
    api_key = "sk-test",
    base_url = "https://api.example/v1",
    model = "gpt-5.5",
    endpoint = "chat"
  ))

  chat_body <- .easyAI_add_generation_options(
    list(model = chat_provider$model),
    provider = chat_provider,
    endpoint = "chat",
    temperature = 0.2,
    max_tokens = 500
  )
  reasoning_body <- .easyAI_add_generation_options(
    list(model = reasoning_provider$model),
    provider = reasoning_provider,
    endpoint = "chat",
    temperature = 0.2,
    max_tokens = 500
  )
  responses_body <- .easyAI_add_generation_options(
    list(model = reasoning_provider$model),
    provider = reasoning_provider,
    endpoint = "responses",
    temperature = 0.2,
    max_tokens = 500
  )

  expect_equal(chat_body$temperature, 0.2)
  expect_equal(chat_body$max_tokens, 500)
  expect_null(chat_body$max_completion_tokens)
  expect_null(reasoning_body$temperature)
  expect_null(reasoning_body$max_tokens)
  expect_equal(reasoning_body$max_completion_tokens, 500)
  expect_equal(responses_body$temperature, 0.2)
  expect_null(responses_body$max_tokens)
  expect_equal(responses_body$max_output_tokens, 500)
})

test_that("OpenAI-compatible message dispatcher falls back across endpoints", {
  provider <- .easyAI_provider(list(
    api_key = "sk-test",
    base_url = "https://api.example/v1",
    model = "gpt-test",
    endpoint = "auto"
  ))
  called <- character()

  ans <- .easyAI_call_provider_messages(
    messages = list(list(role = "user", content = "Return JSON.")),
    provider = provider,
    api_key = "sk-test",
    temperature = 0,
    response_format = list(type = "json_object"),
    call_chat = function(...) {
      called <<- c(called, "chat")
      stop("HTTP 404 Not Found.", call. = FALSE)
    },
    call_responses = function(...) {
      called <<- c(called, "responses")
      "{\"ok\":true}"
    }
  )

  expect_equal(ans, "{\"ok\":true}")
  expect_equal(called, c("chat", "responses"))
})
