#' @title Analyze Single Cell With Spatial Transcriptomics Data
#' @description This function performs MistyR analysis on a single spatial transcriptomics sample.
#' @param spatial_data Seurat object of spatial data.
#' @param spot_mixture Data frame of spot mixture results (e.g., from cell2location or RCTD).
#' @param output_dir Directory to save the output results.
#' @param out_prefix Prefix for the output files.
#' @return None
#' @export
#' @import tidyverse
#' @import Seurat
#' @import mistyR
#' @import future
#' @examples
#' \dontrun{
#' cell2loc <- read.csv('st_cell2location_res_SRR20330030.csv', row.names = 1)
#' row.names(cell2loc) <- sapply(strsplit(row.names(cell2loc), '_'), '[', 2)
#' colnames(cell2loc) <- sapply(strsplit(colnames(cell2loc), 'w_sf_'), '[', 2)
#' runMistyRAnalysis(spatial_data = spatial_data, spot_mixture = cell2loc, output_dir = "./output_data/", out_prefix = "sample1")
#' }
#'

runMistyRAnalysis <- function(spatial_data, spot_mixture, output_dir, out_prefix) {
  library(tidyverse)
  library(Seurat)
  library(mistyR)
  library(future)
  source(system.file("data", "misty_utilities.R", package = "easySingleCell"))
  source(system.file("data", "run_colocalization.R", package = "easySingleCell"))

  # Check and handle unmatched barcodes
  unmatched_barcodes <- setdiff(colnames(spatial_data), rownames(spot_mixture))
  if (length(unmatched_barcodes) > 0) {
    spatial_data$barcode <- colnames(spatial_data)
    spatial_data <- subset(spatial_data, barcode %in% setdiff(colnames(spatial_data), unmatched_barcodes))
  }

  # Create a new Assay object and set it as the default Assay
  spatial_data[["predictions"]] <- CreateAssayObject(t(spot_mixture))
  DefaultAssay(spatial_data) <- 'predictions'

  # Define useful features
  useful_features <- rownames(spatial_data)
  useful_features <- useful_features[!useful_features %in% "prolif"]

  # Run colocalization analysis
  misty_output <- run_colocalization(
    spatial_data = spatial_data,
    assay = 'predictions',
    useful_features = useful_features,
    out_label = out_prefix,
    misty_out_alias = output_dir
  )

  # Collect and visualize Misty results
  misty_results <- collect_results(misty_output)
  plot_folder <- file.path(misty_output, "plots")

  # Create directory if it doesn't exist
  if (!dir.exists(plot_folder)) {
    dir.create(plot_folder, recursive = TRUE)
  }

  # Create a PDF file and generate various plots
  pdf(file = file.path(plot_folder, paste0(out_prefix, "_summary_plots.pdf")))

  mistyR::plot_improvement_stats(misty_results)
  mistyR::plot_view_contributions(misty_results)

  mistyR::plot_interaction_heatmap(misty_results, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_results, "intra", cutoff = 0)

  mistyR::plot_interaction_heatmap(misty_results, "juxta_2", cutoff = 0)
  mistyR::plot_interaction_communities(misty_results, "juxta_2", cutoff = 0)

  mistyR::plot_interaction_heatmap(misty_results, "para_5", cutoff = 0)
  mistyR::plot_interaction_communities(misty_results, "para_5", cutoff = 0)

  dev.off()
}

# Example call
# runMistyRAnalysis(
#   spatial_data = spatial_seurat_object,
#   spot_mixture = spot_mixture_results,
#   output_dir = './output_data/stRNA/',
#   out_prefix = 'sample_1'
# )
