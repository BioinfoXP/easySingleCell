test_that("runScRNAQC accepts Seurat-style object argument", {
  counts <- Matrix::Matrix(
    c(
      20, 0, 0,
      20, 0, 0,
      20, 0, 0,
      20, 0, 0,
      20, 0, 0,
      0, 100, 100
    ),
    nrow = 6,
    sparse = TRUE,
    dimnames = list(
      c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "MT-ND1"),
      c("cell1", "cell2", "cell3")
    )
  )

  object <- Seurat::CreateSeuratObject(counts = counts, project = "compat")
  filtered <- runScRNAQC(
    object = object,
    minGene = 1,
    maxGene = 10,
    pctMT = 50,
    maxCounts = 200,
    species = "human"
  )

  expect_s4_class(filtered, "Seurat")
  expect_equal(colnames(filtered), c("cell1", "cell2"))
})

test_that("convert_id supports AnnotationDbi keytypes beyond symbol and ensembl", {
  converted <- convert_id(
    c("7157", "1956"),
    from_type = "entrezid",
    to_type = "symbol",
    verbose = FALSE
  )

  expect_true(all(c("TP53", "EGFR") %in% converted$Converted_ID))
  expect_equal(unique(converted$FromType), "ENTREZID")
  expect_equal(unique(converted$ToType), "SYMBOL")
})
