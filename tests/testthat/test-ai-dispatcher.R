test_that("model-backed AI helpers use the shared easyAI provider dispatcher", {
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")

  calls <- list()
  fake_dispatch <- function(messages,
                            provider,
                            api_key,
                            temperature = 0.1,
                            response_format = NULL,
                            max_tokens = NULL,
                            ...) {
    calls[[length(calls) + 1]] <<- list(
      messages = messages,
      provider = provider,
      api_key = api_key,
      temperature = temperature,
      response_format = response_format,
      max_tokens = max_tokens
    )
    text <- paste(vapply(messages, `[[`, character(1), "content"), collapse = "\n")
    if (grepl("cluster_id", text, fixed = TRUE)) {
      return('{"cluster_id":"0","primary_lineage":"T cell","detailed_subtype":"CD8 T cell","functional_state":"cytotoxic","confidence":"High","reasoning":"CD3D and CD8A markers"}')
    }
    if (grepl("classifications", text, fixed = TRUE)) {
      return('{"classifications":[{"id":"GSE001","type":"scRNA-seq"}]}')
    }
    if (grepl("selected_ids", text, fixed = TRUE)) {
      return('{"selected_ids":["GSE001"]}')
    }
    '{"matches":["Lung Adenocarcinoma"]}'
  }

  testthat::local_mocked_bindings(
    .easyAI_call_provider_messages = fake_dispatch,
    .package = "easySingleCell"
  )

  markers <- data.frame(
    cluster = rep("0", 3),
    gene = c("CD3D", "CD8A", "GZMK"),
    avg_log2FC = c(2.5, 2.1, 1.8),
    p_val_adj = c(0.001, 0.001, 0.002),
    stringsAsFactors = FALSE
  )
  celltypes <- AutoCellType(
    input = markers,
    tissuename = "Human PBMC",
    n_cores = 1,
    api_key = "sk-test",
    model = "gpt-test",
    base_url = "https://example.test/v1/responses",
    endpoint = "responses",
    verbose = FALSE
  )
  expect_equal(celltypes$primary_lineage, "T cell")

  geo <- data.frame(
    dataset_id = "GSE001",
    species = "Homo sapiens",
    AI_SeqType = "scRNA-seq",
    disease_type = "Pancreatic Cancer",
    sample_size = 12,
    experimental_design = "10x single-cell RNA-seq from pancreatic cancer tissue.",
    stringsAsFactors = FALSE
  )

  classified <- GEO_Classify_AI(
    input_data = geo,
    batch_size = 1,
    api_key = "sk-test",
    model = "gpt-test",
    base_url = "https://example.test/v1/responses",
    endpoint = "responses",
    verbose = FALSE
  )
  expect_equal(classified$AI_SeqType, "scRNA-seq")

  screened <- GEO_Screen_AI(
    input_data = geo,
    disease = "pancreatic cancer",
    data_type = "single cell",
    batch_size = 1,
    api_key = "sk-test",
    model = "gpt-test",
    base_url = "https://example.test/v1/responses",
    endpoint = "responses",
    verbose = FALSE
  )
  expect_equal(screened$dataset_id, "GSE001")

  rdata_path <- tempfile(fileext = ".Rdata")
  depmap_metadata <- data.frame(
    OncotreePrimaryDisease = c("Lung Adenocarcinoma", "Melanoma"),
    row.names = c("ACH-001", "ACH-002"),
    stringsAsFactors = FALSE
  )
  save(depmap_metadata, file = rdata_path)

  ids <- DepmapMetaSelect(
    query = "lung cancer",
    col_name = "OncotreePrimaryDisease",
    rdata_path = rdata_path,
    api_key = "sk-test",
    model = "gpt-test",
    base_url = "https://example.test/v1/responses",
    endpoint = "responses",
    verbose = FALSE
  )
  expect_equal(ids, "ACH-001")

  expect_length(calls, 4)
  expect_true(all(vapply(calls, function(x) inherits(x$provider, "easyAI_provider"), logical(1))))
  expect_true(all(vapply(calls, function(x) identical(x$provider$endpoint_preference, "responses"), logical(1))))
  expect_true(all(vapply(calls, function(x) identical(x$provider$endpoint, "responses"), logical(1))))
  expect_true(all(vapply(calls, function(x) identical(x$response_format$type, "json_object"), logical(1))))
})

test_that("model-backed AI helpers reuse environment key and package provider defaults", {
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")

  with_no_openai_key({
    Sys.setenv(OPENAI_API_KEY = "sk-env")

    calls <- list()
    fake_dispatch <- function(messages,
                              provider,
                              api_key,
                              temperature = 0.1,
                              response_format = NULL,
                              max_tokens = NULL,
                              ...) {
      calls[[length(calls) + 1]] <<- list(provider = provider, api_key = api_key)
      text <- paste(vapply(messages, `[[`, character(1), "content"), collapse = "\n")
      if (grepl("cluster_id", text, fixed = TRUE)) {
        return('{"cluster_id":"0","primary_lineage":"T cell","detailed_subtype":"CD8 T cell","functional_state":"cytotoxic","confidence":"High","reasoning":"CD3D and CD8A markers"}')
      }
      if (grepl("classifications", text, fixed = TRUE)) {
        return('{"classifications":[{"id":"GSE001","type":"scRNA-seq"}]}')
      }
      if (grepl("selected_ids", text, fixed = TRUE)) {
        return('{"selected_ids":["GSE001"]}')
      }
      '{"matches":["Lung Adenocarcinoma"]}'
    }

    testthat::local_mocked_bindings(
      .easyAI_call_provider_messages = fake_dispatch,
      .package = "easySingleCell"
    )

    markers <- data.frame(
      cluster = rep("0", 2),
      gene = c("CD3D", "CD8A"),
      avg_log2FC = c(2.5, 2.1),
      p_val_adj = c(0.001, 0.001),
      stringsAsFactors = FALSE
    )
    invisible(AutoCellType(
      input = markers,
      tissuename = "Human PBMC",
      n_cores = 1,
      verbose = FALSE
    ))

    geo <- data.frame(
      dataset_id = "GSE001",
      experimental_design = "10x single-cell RNA-seq from pancreatic cancer tissue.",
      stringsAsFactors = FALSE
    )
    invisible(GEO_Classify_AI(input_data = geo, batch_size = 1, verbose = FALSE))

    geo$AI_SeqType <- "scRNA-seq"
    invisible(GEO_Screen_AI(
      input_data = geo,
      disease = "pancreatic cancer",
      batch_size = 1,
      verbose = FALSE
    ))

    rdata_path <- tempfile(fileext = ".Rdata")
    depmap_metadata <- data.frame(
      OncotreePrimaryDisease = "Lung Adenocarcinoma",
      row.names = "ACH-001",
      stringsAsFactors = FALSE
    )
    save(depmap_metadata, file = rdata_path)
    invisible(DepmapMetaSelect(
      query = "lung cancer",
      col_name = "OncotreePrimaryDisease",
      rdata_path = rdata_path,
      verbose = FALSE
    ))

    expect_length(calls, 4)
    expect_true(all(vapply(calls, function(x) identical(x$api_key, "sk-env"), logical(1))))
    expect_equal(
      vapply(calls, function(x) x$provider$model, character(1)),
      c("gpt-4o-mini", "gpt-5.5", "gpt-5.5", "gpt-5.5")
    )
    expect_true(all(vapply(calls, function(x) identical(x$provider$base_url, "https://api.gpt.ge/v1"), logical(1))))
    expect_true(all(vapply(calls, function(x) identical(x$provider$endpoint_preference, "auto"), logical(1))))
  })
})
