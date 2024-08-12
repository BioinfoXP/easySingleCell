#' @title Run GSEA Analysis For ScRNA
#' @description This function performs Gene Set Enrichment Analysis (GSEA) on Seurat object data.
#' @param sce A Seurat object containing the single-cell RNA-seq data.
#' @param ident1 The identity class to be used as the reference group.
#' @param ident2 The identity class to be compared against the reference group.
#' @param logfc_threshold The log fold change threshold for identifying differentially expressed genes. Default is 0.1.
#' @param gmt_paths A list of paths to GMT files for GSEA analysis. Default includes Hallmark, GO BP, and KEGG pathways for human.
#' @param pvalue_cutoff The p-value cutoff for GSEA results. Default is 0.05.
#' @param p_adjust_method The method for p-value adjustment. Default is 'none'.
#' @return A list of GSEA results for each GMT file.
#' @export
#' @import Seurat
#' @import dplyr
#' @import clusterProfiler
#' @import fgsea
#' @examples
#' \dontrun{
#' sce <- readRDS('./path_to_your_seurat_object.rds')
#' gsea_results <- run_gsea_analysis(
#'   sce = sce,
#'   ident1 = "cell_type_1",
#'   ident2 = "cell_type_2"
#' )
#' }

runScGSEA <- function(sce, ident1, ident2, logfc_threshold = 0.1,
                              gmt_paths = list("./data/h.all.v7.4.symbols.gmt",
                                               "./data/c5.go.bp.v7.4.symbols.gmt",
                                               "./data/c2.cp.kegg.v7.4.symbols.gmt"),
                              pvalue_cutoff = 0.05, p_adjust_method = 'none') {

  # Load necessary libraries
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(fgsea)

  # Ensure the identity classes are set correctly
  Idents(sce) <- ident1

  # Find differentially expressed genes
  sce_edg <- FindMarkers(sce, ident.1 = ident1, ident.2 = ident2, logfc.threshold = logfc_threshold)

  # Add gene symbols
  sce_edg$SYMBOL <- rownames(sce_edg)
  names(sce_edg)[names(sce_edg) == "avg_log2FC"] <- "logFC"
  sce_edg <- sce_edg %>% arrange(desc(logFC))

  # Create gene list
  geneList <- sce_edg$logFC
  names(geneList) <- sce_edg$SYMBOL

  # Read GMT files
  gmt_list <- lapply(gmt_paths, function(path) {
    if (file.exists(path)) {
      return(read.gmt(path))
    } else {
      stop(paste("File not found:", path))
    }
  })

  # Run GSEA analysis
  gsea_results <- lapply(gmt_list, function(gmt) {
    GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_adjust_method)
  })

  # Get GSEA results
  gsea_results_list <- lapply(gsea_results, function(gsea) gsea@result)

  return(gsea_results_list)
}

# Example call
# sce <- readRDS('./path_to_your_seurat_object.rds')
# gsea_results <- run_gsea_analysis(
#   sce = sce,
#   ident1 = "cell_type_1",
#   ident2 = "cell_type_2"
# )
