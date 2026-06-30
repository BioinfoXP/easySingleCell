test_that("Ro/e symbol labels are centered around expected value 1", {
  roe_mat <- matrix(
    c(2.2, 2.0, 1.7, 1.5, 1.2, 1.1, 1.0, 0.9, 0.8, 0.67, 0.55, 0.5, 0.3),
    nrow = 1,
    dimnames = list("celltype_a", paste0("group_", seq_len(13)))
  )

  labels <- easySingleCell:::.scVisRoeLabels(roe_mat, "symbol")

  expect_equal(
    unname(as.character(labels[1, ])),
    c("+++", "+++", "++", "++", "+", "", "", "", "-", "-", "--", "--", "---")
  )
})

test_that("Ro/e numeric labels keep two decimal places", {
  roe_mat <- matrix(
    c(1, 1.234, 0.996),
    nrow = 1,
    dimnames = list("celltype_a", c("group_a", "group_b", "group_c"))
  )

  labels <- easySingleCell:::.scVisRoeLabels(roe_mat, "numeric")

  expect_equal(
    unname(as.character(labels[1, ])),
    c("1.00", "1.23", "1.00")
  )
})

test_that("Ro/e labels can be disabled", {
  roe_mat <- matrix(1, nrow = 1, dimnames = list("celltype_a", "group_a"))

  expect_false(easySingleCell:::.scVisRoeLabels(roe_mat, "none"))
})

test_that("Ro/e heatmap breaks stay non-negative and keep 1 as visual center", {
  roe_mat <- matrix(
    c(3.4, 0.4, 1),
    nrow = 1,
    dimnames = list("celltype_a", c("high", "low", "expected"))
  )

  breaks <- easySingleCell:::.scVisRoeBreaks(roe_mat)

  expect_equal(length(breaks), 101)
  expect_true(min(breaks) >= 0)
  expect_equal(breaks[51], 1)
  expect_lte(max(roe_mat), max(breaks))
  expect_gte(min(roe_mat), min(breaks))
})

test_that("scVisRoePlot exposes separate row and column font controls", {
  formals <- formals(scVisRoePlot)

  expect_true("font.size.row" %in% names(formals))
  expect_true("font.size.col" %in% names(formals))
  expect_true("symbol.cutoffs" %in% names(formals))
  expect_identical(as.character(formals$font.size.row), "font.size")
  expect_identical(as.character(formals$font.size.col), "font.size")
})

test_that("Ro/e heatmap returns pheatmap object with Ro/e matrix attached", {
  roe_mat <- matrix(
    c(1.2, 0.8, 2.1, 0.4),
    nrow = 2,
    dimnames = list(c("type_a", "type_b"), c("group_a", "group_b"))
  )
  pheatmap_object <- structure(list(tree_row = NULL), class = "pheatmap")

  local_mocked_bindings(
    .package = "pheatmap",
    pheatmap = function(mat, ...) {
      expect_equal(mat, roe_mat)
      pheatmap_object
    }
  )

  plotted <- easySingleCell:::.scVisRoeHeatmap(
    roe_mat,
    title = "Ro/e",
    display.mode = "symbol",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    font.size = 10,
    font.size.row = 9,
    font.size.col = 8
  )

  expect_s3_class(plotted, "pheatmap")
  expect_equal(attr(plotted, "roe_mat"), roe_mat)
})
