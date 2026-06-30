test_that("package source does not ship bundled data directory", {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(root, "DESCRIPTION")) && dirname(root) != root) {
    root <- dirname(root)
  }
  if (!file.exists(file.path(root, "DESCRIPTION"))) {
    skip("Package root could not be located.")
  }

  expect_false(file.exists(file.path(root, "data")))
  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
  expect_false("LazyData" %in% names(desc))
})

test_that("interactive easyAI assistant is not exported or documented", {
  expect_false("easyAI" %in% getNamespaceExports("easySingleCell"))
  expect_false(exists("easyAI", envir = asNamespace("easySingleCell"), inherits = FALSE))
  expect_false(file.exists(test_path("../../man/easyAI.Rd")))
  expect_false(dir.exists(test_path("../../inst/skills/easySingleCell-ai")))
})

test_that("workflow vignettes are formal package vignettes", {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(root, "DESCRIPTION")) && dirname(root) != root) {
    root <- dirname(root)
  }
  if (!file.exists(file.path(root, "DESCRIPTION"))) {
    skip("Package root could not be located.")
  }

  rbuildignore <- readLines(file.path(root, ".Rbuildignore"), warn = FALSE)
  expect_false(any(grepl("^\\^vignettes/?", rbuildignore)))
  expect_true(file.exists(file.path(root, "vignettes", "SingleCell_Workflow.Rmd")))
  expect_true(file.exists(file.path(root, "vignettes", "Bulk_GEO_DepMap_Workflow.Rmd")))
})
