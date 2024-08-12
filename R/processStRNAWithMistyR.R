run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results"
                             # folder name for output
) {

  mistyR::clear_cache()

  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )

  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )

  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias
  )
}


# Extracts data from an specific assay from a Seurat object
# and aligns the IDs to the geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  print(assay)
  data <- GetAssayData(visium.slide, assay = assay) %>%
    as.matrix() %>%
    t() %>%
    as_tibble(rownames = NA)

  return(data %>% slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)

  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}

# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {

  mistyR::clear_cache()

  view.data.init <- create_initial_view(data)

  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }

  if (view.type == "intra") {
    data.red <- view.data.init[["intraview"]]$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)

    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )

    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }

  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }

  return(misty.view)
}

# Builds automatic MISTy workflow and runs it
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default") {

  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }

  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)

  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))

  # Create other views
  view.names <- names(view.data.filt)

  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )

  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )


  # Run MISTy
  run_misty(pline.views, out.alias, cached = FALSE)
}

#
# Bug in collecting results
#

collect_results_v2 <- function(folders){
  samples <- R.utils::getAbsolutePath(folders)
  message("\nCollecting improvements")
  improvements <- samples %>% furrr::future_map_dfr(function(sample) {
    performance <- readr::read_table2(paste0(sample, .Platform$file.sep,
                                             "performance.txt"), na = c("", "NA", "NaN"), col_types = readr::cols()) %>%
      dplyr::distinct()
    performance %>% dplyr::mutate(sample = sample, gain.RMSE = 100 *
                                    (.data$intra.RMSE - .data$multi.RMSE)/.data$intra.RMSE,
                                  gain.R2 = 100 * (.data$multi.R2 - .data$intra.R2),
    )
  }, .progress = TRUE) %>% tidyr::pivot_longer(-c(.data$sample,
                                                  .data$target), names_to = "measure")
  message("\nCollecting contributions")
  contributions <- samples %>% furrr::future_map_dfr(function(sample) {
    coefficients <- readr::read_table2(paste0(sample, .Platform$file.sep,
                                              "coefficients.txt"), na = c("", "NA", "NaN"), col_types = readr::cols()) %>%
      dplyr::distinct()
    coefficients %>% dplyr::mutate(sample = sample, .after = "target") %>%
      tidyr::pivot_longer(cols = -c(.data$sample, .data$target),
                          names_to = "view")
  }, .progress = TRUE)
  improvements.stats <- improvements %>% dplyr::filter(!stringr::str_starts(.data$measure,
                                                                            "p\\.")) %>% dplyr::group_by(.data$target, .data$measure) %>%
    dplyr::summarise(mean = mean(.data$value), sd = stats::sd(.data$value),
                     cv = .data$sd/.data$mean, .groups = "drop")
  contributions.stats <- dplyr::inner_join((contributions %>%
                                              dplyr::filter(!stringr::str_starts(.data$view, "p\\.") &
                                                              .data$view != "intercept") %>% dplyr::group_by(.data$target,
                                                                                                             .data$view) %>% dplyr::summarise(mean = mean(.data$value),
                                                                                                                                              .groups = "drop_last") %>% dplyr::mutate(fraction = abs(.data$mean)/sum(abs(.data$mean))) %>%
                                              dplyr::ungroup()), (contributions %>% dplyr::filter(stringr::str_starts(.data$view,
                                                                                                                      "p\\.") & !stringr::str_detect(.data$view, "intercept")) %>%
                                                                    dplyr::group_by(.data$target, .data$view) %>% dplyr::mutate(view = stringr::str_remove(.data$view,
                                                                                                                                                           "^p\\.")) %>% dplyr::summarise(p.mean = mean(.data$value),
                                                                                                                                                                                          p.sd = stats::sd(.data$value), .groups = "drop")), by = c("target",
                                                                                                                                                                                                                                                    "view"))
  message("\nCollecting importances")
  importances <- samples %>% furrr::future_map(function(sample) {
    targets <- contributions.stats %>% dplyr::pull(.data$target) %>%
      unique() %>% sort()
    views <- contributions.stats %>% dplyr::pull(.data$view) %>%
      unique()
    maps <- views %>% furrr::future_map(function(view) {
      all.importances <- targets %>% purrr::map(~readr::read_csv(paste0(sample,
                                                                        .Platform$file.sep, "importances_", .x, "_",
                                                                        view, ".txt"), col_types = readr::cols()) %>%
                                                  dplyr::distinct() %>% dplyr::rename(feature = target))
      features <- all.importances %>% purrr::map(~.x$feature) %>%
        unlist() %>% unique() %>% sort()
      pvalues <- contributions %>% dplyr::filter(sample ==
                                                   !!sample, view == paste0("p.", !!view)) %>% dplyr::mutate(value = 1 -
                                                                                                               .data$value)
      all.importances %>% purrr::imap_dfc(~tibble::tibble(feature = features,
                                                          zero.imp = 0) %>% dplyr::left_join(.x, by = "feature") %>%
                                            dplyr::arrange(.data$feature) %>% dplyr::mutate(imp = scale(.data$imp)[,
                                                                                                                   1], `:=`(!!targets[.y], .data$zero.imp + (.data$imp *
                                                                                                                                                               (pvalues %>% dplyr::filter(target == targets[.y]) %>%
                                                                                                                                                                  dplyr::pull(.data$value))))) %>% dplyr::select(targets[.y])) %>%
        dplyr::mutate(Predictor = features)
    }) %>% `names<-`(views)
  }, .progress = TRUE) %>% `names<-`(samples)
  message("\nAggregating")
  importances.aggregated <- importances %>% purrr::reduce(function(acc,
                                                                   l) {
    acc %>% purrr::map2(l, ~(((.x %>% dplyr::select(-.data$Predictor)) +
                                (.y %>% dplyr::select(-.data$Predictor))) %>% dplyr::mutate(Predictor = .x %>%
                                                                                              dplyr::pull(.data$Predictor))))
  }) %>% purrr::map(~.x %>% dplyr::mutate_if(is.numeric, ~./length(samples)))
  return(list(improvements = improvements, improvements.stats = improvements.stats,
              contributions = contributions, contributions.stats = contributions.stats,
              importances = importances, importances.aggregated = importances.aggregated))
}













# run_colocalization.R

#' Run Colocalization
#'
#' This function runs colocalization analysis using the Misty Seurat analysis.
#'
#' @param slide A Visium slide object.
#' @param assay The assay to use for the analysis.
#' @param useful_features A vector of useful features to include in the analysis.
#' @param out_label A label for the output.
#' @param misty_out_alias An alias for the Misty output.
#'
#' @return A character string representing the output path.
#' @import Seurat

run_colocalization <- function(slide, assay, useful_features, out_label, misty_out_alias) {
  # Define assays for each view
  view_assays <- list(
    "main" = assay,
    "juxta" = assay,
    "para" = assay
  )

  # Define features for each view
  view_features <- list(
    "main" = useful_features,
    "juxta" = useful_features,
    "para" = useful_features
  )

  # Define spatial context for each view
  view_types <- list(
    "main" = "intra",
    "juxta" = "juxta",
    "para" = "para"
  )

  # Define additional parameters for each view
  view_params <- list(
    "main" = NULL,
    "juxta" = 2,
    "para" = 5
  )

  # Define output path
  misty_out <- paste0(misty_out_alias, out_label, "_", assay)

  # Run Misty Seurat analysis
  run_misty_seurat(
    visium.slide = slide,
    view.assays = view_assays,
    view.features = view_features,
    view.types = view_types,
    view.params = view_params,
    spot.ids = NULL,
    out.alias = misty_out
  )

  return(misty_out)
}








# ============== runMistyRAnalysis ===========

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
