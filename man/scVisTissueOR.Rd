% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scVisual_panel.R
\name{scVisTissueOR}
\alias{scVisTissueOR}
\title{Visualize Tissue Odds Ratio (OR) Analysis}
\usage{
scVisTissueOR(
  scRNA,
  group = "orig.ident",
  celltype = "celltype",
  output_dir = "./output_figure",
  output_file = "tissue_OR.pdf",
  width = 5,
  height = 4
)
}
\arguments{
\item{scRNA}{A Seurat object containing scRNA-seq data.}

\item{group}{A character string specifying the grouping variable in the metadata. Default is 'orig.ident'.}

\item{celltype}{A character string specifying the cell type variable in the metadata. Default is 'celltype'.}

\item{output_dir}{A character string specifying the prefix for output files. Default is "./output_figure".}

\item{output_file}{A character string specifying the output file name. Default is "tissue_OR.pdf".}

\item{width}{Numeric value specifying the width of the output plot. Default is 5.}

\item{height}{Numeric value specifying the height of the output plot. Default is 4.}
}
\value{
A list with tissue OR analysis results and a heatmap plot.
}
\description{
This function performs tissue odds ratio analysis and plots the results.
}
\examples{
\dontrun{
# Assuming `sce` is a pre-existing Seurat object
scVisTissueOR(sce, group = "group", celltype = "celltype", output_dir = "./output_figure", output_file = "tissue_OR.pdf", width = 5, height = 4)
}
}
