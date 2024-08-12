# ================== 1. runHdWGCNAStep1 ===========
#' @title Perform Initial HdWGCNA Analysis to Find Soft Threshold
#' @description This function performs the initial steps of HdWGCNA analysis on a Seurat object to find the soft threshold.
#' @param sce A Seurat object containing scRNA-seq data.
#' @param name A character string specifying the common name for HdWGCNA experiment and group analysis.
#' @param group_by A character vector specifying the columns in sce@meta.data to group by.
#' @param output_dir A character string specifying the directory to save the output. Default is "./output_data/".
#' @param gene_select A character string specifying the gene selection approach. Default is "fraction".
#' @param fraction A numeric value specifying the fraction of cells that a gene needs to be expressed in order to be included. Default is 0.05.
#' @param reduction A character string specifying the dimensionality reduction to perform KNN on. Default is "harmony".
#' @param k An integer specifying the nearest-neighbors parameter. Default is 25.
#' @param max_shared An integer specifying the maximum number of shared cells between two metacells. Default is 10.
#' @param network_type A character string specifying the network type for TestSoftPowers. Default is "signed".
#' @return A Seurat object with updated HdWGCNA analysis.
#' @export
#' @import Seurat
#' @import hdWGCNA
#' @import qs
#' @import patchwork
#' @examples
#' \dontrun{
#' sce <- qread('./output_data/sce.qc.qs')
#' sce <- runHdWGCNAStep1(sce, name = "HdWGCNA", group_by = c("cell_type", "Sample"), output_dir = "./output_data/")
#' }

runHdWGCNAStep1 <- function(sce,
                            name,
                            group_by,
                            output_dir = "./output_data/",
                            gene_select = "fraction",
                            fraction = 0.05,
                            reduction = "harmony",
                            k = 25,
                            max_shared = 10,
                            network_type = "signed") {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(hdWGCNA)
  library(qs)
  library(patchwork)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Step 1: Setup for WGCNA
  sce <- SetupForWGCNA(
    seurat_obj = sce,
    gene_select = gene_select,
    fraction = fraction,
    wgcna_name = name
  )

  # Step 2: Construct metacells in each group
  sce <- MetacellsByGroups(
    seurat_obj = sce,
    group.by = group_by,
    reduction = reduction,
    k = k,
    max_shared = max_shared,
    ident.group = group_by[1]
  )

  # Step 3: Normalize metacell expression matrix
  sce <- NormalizeMetacells(sce)

  # Step 4: Set DatExpr for WGCNA analysis
  sce <- SetDatExpr(
    seurat_obj = sce,
    group_name = unique(sce@meta.data[, group_by[1]]),
    group.by = group_by[1]
  )

  # Step 5: Test different soft powers
  sce <- TestSoftPowers(
    seurat_obj = sce,
    networkType = network_type
  )

  # Step 6: Plot the results
  plot_list <- PlotSoftPowers(sce)

  # Save the plot
  plot_file <- file.path(output_dir, "soft_threshold_plot_step1.pdf")
  ggsave(filename = plot_file, plot = wrap_plots(plot_list, ncol = 2), width = 10, height = 10)

  # Save the Seurat object
  sce_file <- file.path(output_dir, "sce_hdWGCNA_step1.qs")
  qsave(sce, file = sce_file)

  # Return the Seurat object
  return(sce)
}

# Example usage
# sce <- qread('./output_data/sce.qc.qs')
# sce <- runHdWGCNAStep1(sce, name = "HdWGCNA", group_by = c("cell_type", "Sample"), output_dir = "./output_data/")


# ================== 2. runHdWGCNAStep2 ===========
#' @title Perform HdWGCNA Analysis Step 2
#' @description This function performs the second step of HdWGCNA analysis on a Seurat object, including constructing the co-expression network, plotting KME, extracting hub genes, and generating a DotPlot for hMEs.
#' @param sce A Seurat object containing scRNA-seq data.
#' @param soft_power An integer specifying the soft power value for constructing the co-expression network.
#' @param name A character string specifying the name of the experiment, used for naming the topological overlap matrix written to disk. Default is "HdWGCNA".
#' @param group_by A character string specifying the metadata column to group by in the DotPlot.
#' @param n_hubs An integer specifying the number of hub genes to extract. Default is 50.
#' @param output_dir A character string specifying the directory to save the output. Default is "./output_data/".
#' @return A list containing the Seurat object with updated metadata and the hub genes data frame.
#' @export
#' @import Seurat
#' @import hdWGCNA
#' @import qs
#' @import ggplot2
#' @examples
#' \dontrun{
#' sce <- qread('./output_data/sce_hdWGCNA_step1.qs')
#' results <- runHdWGCNAStep2(sce, soft_power = 9, name = "MyExperiment", group_by = "cell_type", n_hubs = 50)
#' seurat_obj <- results$seurat_obj
#' hub_df <- results$hub_df
#' }

runHdWGCNAStep2 <- function(sce,
                            soft_power = NULL,
                            name = "HdWGCNA",
                            group_by = NULL,
                            n_hubs = 50,
                            output_dir = "./output_data/") {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(hdWGCNA)
  library(qs)
  library(ggplot2)

  # Check if necessary parameters are provided
  if (is.null(sce)) {
    stop("Error: 'sce' must be provided.")
  }

  # Check if necessary parameters are provided
  if (is.null(group_by)) {
    stop("Error: 'group_by' must be provided.")
  }


  if (is.null(soft_power)) {
    stop("Error: 'soft_power' must be provided.")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Construct co-expression network
  sce <- ConstructNetwork(
    sce,
    soft_power = soft_power,
    setDatExpr = FALSE,
    overwrite_tom = TRUE,
    tom_name = name
  )

  sce <- ModuleEigengenes(sce)
  sce <- ModuleConnectivity(sce)

  # Save the Seurat object
  sce_file <- file.path(output_dir, "sce_hdWGCNA_step2.qs")
  qsave(sce, file = sce_file)

  # Plot dendrogram
  dendrogram_file <- file.path(output_dir, "hdWGCNA_dendrogram.pdf")
  pdf(dendrogram_file, width = 9, height = 6)
  dendrogram_plot <- PlotDendrogram(sce, main = 'hdWGCNA Dendrogram')
  print(dendrogram_plot)
  dev.off()

  # Plot genes ranked by kME for each module
  kme_file <- file.path(output_dir, "hdWGCNA_kme_plot.pdf")
  pdf(kme_file, width = 12, height = 4)
  kme_plot <- PlotKMEs(sce, ncol = 5)
  print(kme_plot)
  dev.off()

  # Get hub genes
  hub_df <- GetHubGenes(sce, n_hubs = n_hubs)
  hub_file <- file.path(output_dir, "hdWGCNA_hub_genes.csv")
  write.csv(hub_df, file = hub_file, row.names = FALSE)

  # Get hMEs from Seurat object
  MEs <- GetMEs(sce, harmonized = TRUE)
  modules <- GetModules(sce)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']

  # Add hMEs to Seurat meta-data
  sce@meta.data <- cbind(sce@meta.data, MEs)

  # Plot with Seurat's DotPlot function
  dot_plot <- DotPlot(sce, features = mods, group.by = group_by) +
    RotatedAxis() +
    scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue')

  dot_plot_file <- file.path(output_dir, "hdWGCNA_dot_plot.pdf")
  pdf(dot_plot_file, width = 6, height = 5)
  print(dot_plot)
  dev.off()

  # Return the Seurat object and hub genes data frame
  return(list(seurat_obj = sce, hub_df = hub_df))
}

# Example usage
# sce <- qread('./output_data/sce_hdWGCNA_step1.qs')
# results <- runHdWGCNAStep2(sce, soft_power = 9, name = "MyExperiment", group_by = "cell_type", n_hubs = 50)
# seurat_obj <- results$seurat_obj
# hub_df <- results$hub_df
