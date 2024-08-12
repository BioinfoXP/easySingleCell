#' @title Perform Metabolism Analysis and Generate DotPlot
#' @description This function performs metabolism analysis on Seurat object data and generates a DotPlot of metabolic pathways.
#' @param seurat_obj A Seurat object containing the spatial transcriptomics data.
#' @param output_file Path to the output Rdata file for saving the Seurat object with metabolism scores. Default is './output_data/Figure4/scMetabolism.Rdata'.
#' @param pdf_file Path to the output PDF file for saving the DotPlot. Default is './output_figure/Figure2/scMetabolism.pdf'.
#' @param pathways Number of top pathways to plot. Default is 20.
#' @param cell_type_column Column name in the Seurat object metadata representing cell types. Default is "cell_type".
#' @param width Width of the output PDF. Default is 6.
#' @param height Height of the output PDF. Default is 10.
#' @return None
#' @export
#' @import scMetabolism
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' seurat_obj <- readRDS('./path_to_your_seurat_object.rds')
#' performMetabolismAnalysis(
#'   seurat_obj = seurat_obj,
#'   pathways = 20,
#'   cell_type_column = "cell_type"
#' )
#' }





runScMetabolismAnalysis <- function(seurat_obj,
                                      output_file = './output_data/Figure4/scMetabolism.Rdata',
                                      pdf_file = './output_figure/Figure2/scMetabolism.pdf',
                                      pathways = 20,
                                      cell_type_column = "cell_type",
                                      width = 6,
                                      height = 10) {

  library(scMetabolism)
  library(Seurat)
  library(ggplot2)

  # 计算代谢通路得分
  human_countexp_Seurat <- sc.metabolism.Seurat(
    obj = seurat_obj,
    method = "VISON",
    imputation = FALSE,
    ncores = 10,
    metabolism.type = "KEGG"
  )

  # 保存Seurat对象
  save(human_countexp_Seurat, file = output_file)

  # 绘制代谢通路的DotPlot
  dotplot <- DotPlot.metabolism(
    obj = human_countexp_Seurat,
    pathway = rownames(human_countexp_Seurat@assays[["METABOLISM"]][["score"]])[1:pathways],
    phenotype = cell_type_column,
    norm = "y"
  ) + xlab('')

  # 保存图像为PDF
  ggsave(pdf_file, plot = dotplot, width = width, height = height)
}

# Example call
# seurat_obj <- readRDS('./path_to_your_seurat_object.rds')
# performMetabolismAnalysis(
#   seurat_obj = seurat_obj,
#   pathways = 20,
#   cell_type_column = "cell_type"
# )
