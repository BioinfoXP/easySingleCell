% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_panel.R
\name{RunCopyKAT}
\alias{RunCopyKAT}
\title{Run CopyKAT Analysis}
\usage{
RunCopyKAT(
  seurat_obj,
  output_dir = "./output_data/copykat",
  sam_name = "sample",
  id_type = "S",
  cell_line = "no",
  ngene_chr = 1,
  win_size = 25,
  KS_cut = 0.15,
  distance = "euclidean",
  n_cores = 1,
  ...
)
}
\arguments{
\item{seurat_obj}{A Seurat object containing the RNA assay with count data.}

\item{output_dir}{A character string specifying the directory to save CopyKAT output. Default is './output_data/copykat'.}

\item{sam_name}{A character string specifying the sample name for CopyKAT. Default is 'sample'.}

\item{id_type}{A character string specifying the ID type for CopyKAT. Default is 'S'.}

\item{cell_line}{A character string specifying whether the data is from a cell line. Default is 'no'.}

\item{ngene_chr}{An integer specifying the minimum number of genes per chromosome. Default is 1.}

\item{win_size}{An integer specifying the window size for CopyKAT. Default is 25.}

\item{KS_cut}{A numeric value specifying the KS cutoff for CopyKAT. Default is 0.15.}

\item{distance}{A character string specifying the distance metric for CopyKAT. Default is 'euclidean'.}

\item{n_cores}{An integer specifying the number of cores to use for CopyKAT. Default is 1.}

\item{...}{Additional parameters passed to the \code{copykat} function.}
}
\value{
A Seurat object with CopyKAT prediction results added to its metadata.
}
\description{
This function performs CopyKAT analysis on a given Seurat object and saves the results.
}
\examples{
\dontrun{
sce.tumor = RunCopyKAT(
  seurat_obj = sce.tumor,
  output_dir = './output_data/copykat',
  sam_name = 'hepatocyte',
  n_cores = 30
)
}
}
