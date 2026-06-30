test_that("scGSEA marker ranking uses Seurat logFC columns safely", {
  markers <- data.frame(
    gene = c("GeneA", "GeneB", "GeneC", "GeneB"),
    avg_log2FC = c(0.5, 2, NA, 1),
    p_val_adj = c(0.2, 0.01, 0.8, 0.04)
  )

  ranked <- .scGSEA_rank_markers(markers, gene_col = "gene")

  expect_named(ranked, c("GeneB", "GeneA"))
  expect_equal(unname(ranked), c(2, 0.5))
})

test_that("scGSEA marker ranking can use an explicit ranking column", {
  markers <- data.frame(
    symbol = c("GeneA", "GeneB", "GeneC"),
    stat = c(-1, 3, 0.5)
  )

  ranked <- .scGSEA_rank_markers(markers, rank.by = "stat", gene_col = "symbol")

  expect_named(ranked, c("GeneB", "GeneC", "GeneA"))
  expect_equal(unname(ranked), c(3, 0.5, -1))
})

test_that("scGSEA validates species and optional OrgDb packages", {
  expect_identical(.scGSEA_orgdb("human", NULL), org.Hs.eg.db::org.Hs.eg.db)
  expect_error(.scGSEA_orgdb("zebrafish", NULL), "Unsupported species")
})

test_that("scGSEA allows custom OrgDb without forcing built-in species", {
  markers <- data.frame(
    gene = c("GeneA", "GeneB"),
    avg_log2FC = c(1, -1)
  )

  expect_error(
    scGSEA(
      markers = markers,
      gene_col = "gene",
      species = "zebrafish",
      OrgDb = list(custom = TRUE),
      min_size = 99,
      verbose = FALSE
    ),
    "Not enough ranked genes"
  )
})

test_that("scGSEA marker mode returns ranked genes and GSEA result", {
  markers <- data.frame(
    gene = c("GeneA", "GeneB", "GeneC", "GeneD"),
    avg_log2FC = c(1.2, -0.5, 2.5, 0.1)
  )
  captured <- list()
  fake_gsea <- function(geneList, OrgDb, ont, keyType, minGSSize, maxGSSize,
                        pvalueCutoff, verbose, seed) {
    captured$geneList <<- geneList
    captured$ont <<- ont
    data.frame(ID = "GO:0000001", Description = "mock term")
  }

  res <- .scGSEA_from_markers(
    markers = markers,
    gene_col = "gene",
    min_size = 2,
    verbose = FALSE,
    gsea_fun = fake_gsea
  )

  expect_s3_class(res, "easy_scGSEA")
  expect_named(res$gene_rank, c("GeneC", "GeneA", "GeneD", "GeneB"))
  expect_identical(captured$geneList, res$gene_rank)
  expect_equal(captured$ont, "BP")
  expect_equal(nrow(res$gsea), 1)
})
