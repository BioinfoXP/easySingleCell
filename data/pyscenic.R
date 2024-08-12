makeMetaCells <- function(seu, min.cells = 10, reduction = "umap", dims = 1:2, k.param = 10, cores = 20) {
  seu <- seu %>%
    FindNeighbors(reduction = reduction, dims = dims, k.param = k.param) %>%
    FindClusters(res = 50)
  metadata <- seu@meta.data
  metadata$METACELL_ID <- factor(metadata$seurat_clusters)
  dge_mat <- seu[["RNA"]]@counts

  dge_mat_mc <- parallel::mclapply(levels(metadata$METACELL_ID), function(xx) {
    cells <- rownames(subset(metadata, METACELL_ID == xx))
    Matrix::rowSums(dge_mat[, cells])
  }, mc.cores = cores)
  dge_mat_mc <- do.call(cbind, dge_mat_mc)

  metacell_metadata <- metadata[["METACELL_ID"]] %>% table() %>% as.data.frame()
  colnames(metacell_metadata) <- c("METACELL_ID", "CELL_COUNT")
  rownames(metacell_metadata) <- metacell_metadata[["METACELL_ID"]]

  kept.cells <- subset(metacell_metadata, CELL_COUNT >= min.cells)[["METACELL_ID"]]
  metacells <- list(
    mat = dge_mat_mc[, kept.cells],
    metadata = metacell_metadata[kept.cells, ]
  )
  colnames(metacells$mat) <- paste0(seu@project.name, ".METACELL_", kept.cells)
  rownames(metacells$metadata) <- colnames(metacells$mat)
  return(metacells)
}








