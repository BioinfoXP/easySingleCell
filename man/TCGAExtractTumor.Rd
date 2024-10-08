% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGA_panel.R
\name{TCGAExtractTumor}
\alias{TCGAExtractTumor}
\title{Data Preprocessing for Tumor Samples}
\usage{
TCGAExtractTumor(
  expr,
  surv,
  cli,
  output_file = "./output_data/TCGAExtractTumor.RData"
)
}
\arguments{
\item{expr}{A matrix or data frame containing TPM expression data with genes as rows and samples as columns.}

\item{surv}{A data frame containing survival data with columns 'sample', 'OS', and 'OS.time'.}

\item{cli}{A data frame containing clinical data with a column 'submitter_id.samples'.}

\item{output_file}{A character string specifying the file path to save the preprocessed data as an RData file. Default is "./output_data/preprocessed_data.RData".}
}
\value{
A list containing the preprocessed expression, survival, and clinical data for tumor samples.
}
\description{
This function preprocesses the expression, survival, and clinical data for tumor samples.
}
\examples{
\dontrun{
expr <- readRDS('./path_to_your_expression_data.rds')
surv <- read.csv('./path_to_your_survival_data.csv')
cli <- read.csv('./path_to_your_clinical_data.csv')
TCGA_Tumor <- TCGAExtractTumor(expr, surv, cli, output_file = "./output_data/TCGAExtractTumor.RData")
}
}
