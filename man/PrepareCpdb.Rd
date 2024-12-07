% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_panel.R
\name{PrepareCpdb}
\alias{PrepareCpdb}
\title{Preprocess and Save Single-cell Data for CellPhoneDB Analysis}
\usage{
PrepareCpdb(
  scRNA,
  output_dir = "./output_data/",
  celltype_column = "celltype",
  normalization_threshold = 0.001,
  prefix = ""
)
}
\arguments{
\item{scRNA}{A Seurat object containing the single-cell RNA data.}

\item{output_dir}{A character string specifying the output directory where results will be saved. Default is './output_data/'.}

\item{celltype_column}{A character string specifying the column name in the metadata that contains cell type information. Default is 'celltype'.}

\item{normalization_threshold}{A numeric value specifying the threshold for filtering genes based on mean expression. Default is 0.001.}

\item{prefix}{A character string specifying the prefix to be added to the output file names. Default is "".}
}
\value{
None. The function saves the normalized counts and metadata files in the specified output directory.
}
\description{
This function preprocesses single-cell data, normalizes it, and saves the results for CellPhoneDB analysis.
}
\examples{
\dontrun{
PrepareCpdb(
  scRNA = sce,
  output_dir = './output_data/',
  celltype_column = 'celltype',
  normalization_threshold = 0.001,
  prefix = ''
)
}
}
