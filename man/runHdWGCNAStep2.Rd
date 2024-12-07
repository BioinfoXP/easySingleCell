% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_panel.R
\name{runHdWGCNAStep2}
\alias{runHdWGCNAStep2}
\title{Perform HdWGCNA Analysis Step 2}
\usage{
runHdWGCNAStep2(
  sce,
  soft_power = NULL,
  name = "HdWGCNA",
  group_by = NULL,
  n_hubs = 50,
  output_dir = "./output_data/"
)
}
\arguments{
\item{sce}{A Seurat object after hdWGCNA step 1.}

\item{soft_power}{An integer specifying the soft power value for constructing the co-expression network.}

\item{name}{A character string specifying the name of the experiment, used for naming the topological overlap matrix written to disk. Default is "HdWGCNA".}

\item{group_by}{A character string specifying the metadata column to group by in the DotPlot.}

\item{n_hubs}{An integer specifying the number of hub genes to extract. Default is 50.}

\item{output_dir}{A character string specifying the directory to save the output. Default is "./output_data/".}
}
\value{
A list containing the Seurat object with updated metadata and the hub genes data frame.
}
\description{
This function performs the second step of HdWGCNA analysis on a Seurat object, including constructing the co-expression network, plotting KME, extracting hub genes, and generating a DotPlot for hMEs.
}
\examples{
\dontrun{
sce <- qread('./output_data/sce_hdWGCNA_step1.qs')
results <- runHdWGCNAStep2(sce, soft_power = 9, name = "MyExperiment", group_by = "cell_type", n_hubs = 50)
seurat_obj <- results$seurat_obj
hub_df <- results$hub_df
}
}
