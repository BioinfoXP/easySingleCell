% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_panel.R
\name{MarkDoubletsWrapper}
\alias{MarkDoubletsWrapper}
\title{Mark Doublets in Seurat Object}
\usage{
MarkDoubletsWrapper(seu, PCs = 1:10, split.by = NULL)
}
\arguments{
\item{seu}{A Seurat object.}

\item{PCs}{A numeric vector specifying the principal components to use for doublet detection. Default is 1:10.}

\item{split.by}{A character string specifying the column in meta.data to split the analysis by. Default is "orig.ident".}
}
\value{
A Seurat object with doublet information added.
}
\description{
This function wraps the MarkDoublets function from the DoubletFinder package to identify doublets in a Seurat object.
}
\examples{
\dontrun{
seu <- MarkDoubletsWrapper(seu = seu, PCs = 1:10, split.by = "orig.ident")
}
}
