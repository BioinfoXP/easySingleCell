# ============ 1. scVisFeaturePlot =============

#' @title Feature Plot for Single-cell Data
#' @description This function creates a feature plot for single-cell data using a specified reduction method.
#' @param scRNA A Seurat object containing single-cell data.
#' @param feature A character string specifying the feature to be plotted.
#' @param reduction A character string specifying the reduction method to be used (default is "umap").
#' @param pt.size A numeric value specifying the size of the points (default is 0.0001).
#' @param max.cutoff A numeric value specifying the maximum cutoff for the feature values (default is 1.5).
#' @param cols A vector of colors to be used for the plot.
#' @param title A character string specifying the title of the plot. If NULL, the feature name will be used as the title.
#' @return A ggplot object representing the feature plot.
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(ggplot2)
#'
#' # Example Seurat object
#' seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' seurat_obj <- NormalizeData(seurat_obj)
#' seurat_obj <- FindVariableFeatures(seurat_obj)
#' seurat_obj <- ScaleData(seurat_obj)
#' seurat_obj <- RunPCA(seurat_obj)
#' seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#'
#' # Example feature plot
#' plot <- scVisFeaturePlot(
#'   scRNA = seurat_obj,
#'   feature = "CD3E",
#'   reduction = "umap",
#'   pt.size = 0.5,
#'   max.cutoff = 1.5,
#'   cols =  c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),
#'   title = NULL
#' )
#' print(plot)
#' }

scVisFeaturePlot <- function(scRNA, feature, reduction = "umap", pt.size = 0.0001, max.cutoff = 1.5, cols =  c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"), title = NULL) {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Set title to feature name if title is NULL
  if (is.null(title)) {
    title <- feature
  }

  # Create the feature plot
  plot <- FeaturePlot(
    object = scRNA,
    features = feature,
    reduction = reduction,
    pt.size = pt.size,
    max.cutoff = max.cutoff,
    cols = cols
  ) +
    scale_x_continuous("") +
    scale_y_continuous("") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    ggtitle(title)

  return(plot)
}

# Example usage
# library(Seurat)
# library(ggplot2)
#
# # Example Seurat object
# seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#
# # Example feature plot
# plot <- scVisFeaturePlot(
#   scRNA = seurat_obj,
#   feature = "PC_1",
#   reduction = "umap",
#   pt.size = 0.5,
#   max.cutoff = 1.5,
#   cols = c("blue", "red"),
#   title = NULL
# )
# print(plot)



# ============ 2. scVisDimPlot ==========
#' @title Plot Dimensionality Reduction Visualization
#' @description This function plots a dimensionality reduction visualization (e.g., UMAP, TSNE) of a Seurat object.
#' @param scRNA A Seurat object containing scRNA-seq data.
#' @param reduction A character string specifying the reduction method to use. Default is "umap".
#' @param group.by A character string specifying the metadata column to group by. Default is 'celltype'.
#' @param split.by A character string specifying the metadata column to split the plot by. Default is NULL.
#' @param label A logical value indicating whether to label the clusters. Default is TRUE.
#' @param repel A logical value indicating whether to use label repelling. Default is TRUE.
#' @param pt.size A numeric value specifying the size of points. Default is 0.8.
#' @param cols A vector of colors to use for the plot. Default is a predefined set of colors.
#' @return A ggplot2 object representing the dimensionality reduction plot (UMAP/TSNE).
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Assuming `seurat_obj` is a pre-existing Seurat object
#' scVisDimPlot(seurat_obj, group.by = "celltype")
#' }

scVisDimPlot <- function(scRNA,
                         reduction = "umap",
                         group.by = 'celltype',
                         split.by = NULL,
                         label = TRUE,
                         repel = TRUE,
                         pt.size = 0.8,
                         cols = c("#F1788D", "#54990F", "#E6550D", "#843C39", "#3182BD", "#8C6D31",
                                           "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", "#6BAED6", "#9ECAE1",
                                           "#AD494A", "#A1D99B", "#C7E9C0", "#99600F", "#C3BC3F", "#D6616B",
                                           "#FF7F00", "#1B9E77", "#FDAE6B", "#66A61E", "#E6550D", "#E7969C",
                                           '#53A85F', '#EFD4B4', '#9BCFB9', '#F0C0BD', '#FB6264', '#4E9C99',
                                           '#70C055', '#E98A27', '#FEBC28', "#F1788D", "#54990F", "#E6550D",
                                           "#843C39", "#66A61E", "#E6550D", "#E7969C")) {
                                           # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Plot UMAP/TSNE
  plot <- Seurat::DimPlot(scRNA,
                          reduction = reduction,
                          group.by = group.by,
                          split.by = split.by,
                          label = label,
                          repel = repel,
                          pt.size = pt.size,
                          cols = cols)

  return(plot)
}

# Example usage
# Assuming `seurat_obj` is a pre-existing Seurat object
# scVisDimPlot(seurat_obj, group.by = "celltype")


# ============ 3. scVisCellRatioPlot ==========

#' @title Plot Single Cell Ratio Statistics
#' @description This function plots the proportion of different cell types in a Seurat object.
#' @param scRNA A Seurat object containing scRNA-seq data.
#' @param plot_by A character string specifying the variable to plot by. Default is "celltype".
#' @param meta.include A vector of metadata columns to include in the plot.
#' @param group_by A character string specifying the grouping variable.
#' @param shape_by A character string specifying the variable to use for point shapes. Default is NULL.
#' @param custom_fill_colors A vector of custom colors to use for the plot. Default is NULL.
#' @param group_by.point A character string specifying the variable to use for grouping points. Default is NULL.
#' @param color_by A character string specifying the variable to use for coloring points. Default is NULL.
#' @param pb A logical value indicating whether to use a progress bar. Default is FALSE.
#' @param comparisons A list of comparisons to make in the plot. Default is my_comparisons.
#' @param ncol An integer specifying the number of columns in the facet wrap. Default is NULL.
#' @param label A character string specifying the label format for significance. Options are 'p.format' or 'p.signif'. Default is c("p.format","p.signif").
#' @param label.x A numeric value specifying the x position of the label. Default is NA.
#' @param pt.size A numeric value specifying the size of points. Default is NA.
#' @return A ggplot2 object representing the single cell ratio statistics.
#' @export
#' @import Seurat
#' @import ggplot2
#' @import reshape2
#' @examples
#' \dontrun{
#' # Assuming `scedata` is a pre-existing Seurat object
#' my_comparisons <- list(c("BM", "GM"))
#'
#' # Plotting a bar plot
#' scVisCellRatioPlot(scedata, group_by = "group",
#'                    meta.include = c("group","orig.ident"),
#'                    color_by = 'cell.type')
#'
#' # Plotting a grouped box plot 1
#' scVisCellRatioPlot(scedata, group_by = "group",
#'                    meta.include = c("group","orig.ident"),
#'                    comparisons = my_comparisons, color_by = 'group',
#'                    group_by.point = "orig.ident", label.x = 1, pt.size = 3,
#'                    label = 'p.format', ncol = 3)
#'
#' # Plotting a grouped box plot 2
#' scVisCellRatioPlot(scedata, group_by = "group",
#'                    meta.include = c("group","orig.ident"),
#'                    comparisons = my_comparisons, color_by = 'orig.ident',
#'                    group_by.point = "orig.ident", label.x = 1, pt.size = 3,
#'                    label = 'p.format', ncol = 3)
#'
#' # Plotting with different shapes for points
#' scVisCellRatioPlot(scedata, group_by = "group",
#'                    meta.include = c("group","orig.ident"),
#'                    comparisons = my_comparisons, color_by = 'orig.ident',
#'                    group_by.point = "orig.ident", label.x = 1, pt.size = 3,
#'                    label = 'p.format', ncol = 3, shape_by = 'group')
#' }

scVisCellRatioPlot <- function(scRNA, plot_by = "celltype", meta.include = NULL,
                               group_by = NULL, shape_by = NULL,
                               custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL,
                               pb = FALSE, comparisons = my_comparisons,
                               ncol = NULL, label = c("p.format","p.signif"),
                               label.x = NA, pt.size = NA) {

  plot_by <- match.arg(plot_by)
  if (is.null(group_by)) {
    group_by <- "null.group"
  }

  shapes <- NULL
  if (!is.null(shape_by)) {
    shapes <- c(16, 15, 3, 7, 8, 18, 5, 6, 2, 4, 1, 17)
  }

  fq <- prop.table(table(scRNA@meta.data$celltype, scRNA@meta.data[,"orig.ident"]), margin = 2) * 100
  df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", "orig.ident"))

  uniques <- apply(scRNA@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(scRNA@meta.data[, names(uniques[uniques <= 100])])
  ei <- unique(ei[, colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)

  if (is.null(x = ncol)) {
    ncol <- 3
    if (length(unique(df$celltype)) > 9) {
      ncol <- 4
    }
    if (length(unique(df$celltype)) > 20) {
      ncol <- 5
    }
  }

  custom_fill_colors <- c(RColorBrewer::brewer.pal(9, "Oranges")[2],
                          RColorBrewer::brewer.pal(9, "Reds")[6],
                          RColorBrewer::brewer.pal(9, "Oranges")[5],
                          RColorBrewer::brewer.pal(9, "Blues")[4:9])

  p <- ggplot(df, aes_string(y = "freq", x = group_by)) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(face = "bold", size = 14),
          axis.ticks.x = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black', size = 12),
          axis.text.y = element_text(color = 'black', hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(color = 'black', size = 14))

  if (plot_by == "cell.type" && color_by == "cell.type") {
    p <- p + facet_wrap(group_by, scales = "free_x") +
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), position = "fill", stat = "identity") +
      scale_fill_manual("cell.type", values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
                                                         "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C",
                                                         "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02")) +
                                                           scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
      theme(panel.border = element_blank())
  } else {
    p <- p + facet_wrap("cell.type", scales = "free_y", ncol = ncol) +
      guides(fill = FALSE) +
      geom_boxplot(aes_string(x = group_by), alpha = 0.25, outlier.color = NA) +
      geom_point(size = 4, position = position_jitter(width = 0.25),
                 aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)) +
      scale_shape_manual(values = shapes) +
      theme(panel.grid.major = element_line(color = "grey", size = 0.25)) +
      scale_color_manual(values = custom_fill_colors) +
      scale_fill_manual(values = custom_fill_colors) +
      scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.2))) +
      ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label, method = "t.test")
  }

  return(p)
}



# ============ 4. scVisTissueOR ==========

#' @title Visualize Tissue Odds Ratio (OR) Analysis
#' @description This function performs tissue odds ratio analysis and plots the results.
#' @param scRNA A Seurat object containing scRNA-seq data.
#' @param group A character string specifying the grouping variable in the metadata. Default is 'orig.ident'.
#' @param celltype A character string specifying the cell type variable in the metadata. Default is 'celltype'.
#' @param output_prefix A character string specifying the prefix for output files. Default is "./output_figure".
#' @param output_file A character string specifying the output file name. Default is "/tissue_OR.pdf".
#' @param width Numeric value specifying the width of the output plot. Default is 5.
#' @param height Numeric value specifying the height of the output plot. Default is 4.
#' @return A list with tissue OR analysis results and a heatmap plot.
#' @export
#' @import Seurat
#' @import ggplot2
#' @import export
#' @examples
#' \dontrun{
#' # Assuming `sce` is a pre-existing Seurat object
#' scVisTissueOR(sce, group = "group", celltype = "celltype", output_prefix = "./output_figure")
#' }

scVisTissueOR <- function(scRNA, group = 'orig.ident', celltype = 'celltype', output_prefix = "./output_figure", output_file = "/tissue_OR.pdf", width = 5, height = 4) {
  # Load required source file
  source('./config/tissue_OR.R')
  # Extract metadata
  meta <- scRNA@meta.data

  # Assign group and celltype to new columns
  meta$loc <- meta[, group]
  meta$meta.cluster <- meta[, celltype]

  # Perform tissue OR analysis
  OR_immune_list <- analyze_tissue_dist(meta_data = meta, output_prefix = output_prefix)

  # Plot heatmap
  p <- plot_heatmap(OR_immune_list)

  # Ensure the output directory exists
  dir.create(output_prefix, showWarnings = FALSE, recursive = TRUE)

  # Export the plot to a PDF file
  export::graph2pdf(p, file = paste0(output_prefix, output_file), width = width, height = height)

  return(list(OR_immune_list = OR_immune_list, heatmap_plot = p))
}

# Example usage
# Assuming `sce` is a pre-existing Seurat object
# scVisTissueOR(sce, group = "group", celltype = "celltype", output_prefix = "./output_figure", output_file = "/tissue_OR.pdf", width = 5, height = 4)




# ============ 5. CPDB Visualization =============


#' @title CPDB Visualization
#' @description This function visualizes the results from CellPhoneDB analysis.
#' @param sce A Seurat object containing the single-cell RNA data.
#' @param pvals_path A string specifying the path to the p-values file.
#' @param means_path A string specifying the path to the means file.
#' @param decon_path A string specifying the path to the deconvolution file.
#' @param celltype_key A string specifying the key for cell type in the metadata. Default is "celltype".
#' @param sender A character vector specifying the sender cell types.
#' @param receiver A character vector specifying the receiver cell types.
#' @param output_dir A string specifying the output directory. Default is "./output/".
#' @param output_filename A string specifying the output filename. Default is "cpdb.pdf".
#' @param palette1 A character vector specifying the color palette for dot plot. Default is c("darkblue", "yellow", "red").
#' @param palette2 A character vector specifying the color palette for heatmap. Default is c("navy", "white", "firebrick3").
#' @return None. The function saves the visualizations in a PDF file.
#' @export
#' @import Seurat
#' @import pheatmap
#' @import ktplots
#' @import ggplot2
#' @examples
#' \dontrun{
#' cpdb_visualization(
#'   sce = seurat_obj,
#'   pvals_path = "path/to/pvals.txt",
#'   means_path = "path/to/means.txt",
#'   decon_path = "path/to/decon.txt",
#'   sender = c("CellTypeA"),
#'   receiver = c("CellTypeB"),
#'   output_dir = "./output_data/",
#'   output_filename = "cpdb.pdf"
#' )
#' }

cpdb_visualization <- function(sce, pvals_path, means_path, decon_path,
                               celltype_key = "celltype",
                               sender, receiver,
                               output_dir = "./output_data/",
                               output_filename = "cpdb.pdf",
                               palette1 = c("darkblue", "yellow", "red"),
                               palette2 = c("navy", "white", "firebrick3")) {

  # Ensure necessary libraries are loaded
  library(pheatmap)
  library(ktplots)
  library(Seurat)
  library(ggplot2)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set up the PDF output
  pdf(file.path(output_dir, output_filename))

  # Define color palettes
  my_palette1 <- colorRampPalette(palette1)(n = 100)
  my_palette2 <- colorRampPalette(palette2)(100)

  # Load data
  pvals <- read.delim(pvals_path, check.names = FALSE)
  means <- read.delim(means_path, check.names = FALSE)
  decon <- read.delim(decon_path, check.names = FALSE)

  # Prepare cell types for plotting
  cell_type1 <- paste(sender, collapse = '|')
  cell_type2 <- paste(receiver, collapse = '|')

  # Generate interaction table
  table <- plot_cpdb(scdata = sce, keep_significant_only = TRUE,
                     return_table = TRUE,
                     cell_type1 = cell_type1,
                     cell_type2 = cell_type2,
                     celltype_key = celltype_key,
                     means = means, pvals = pvals)

  # Prepare data for dot plot
  plot_data1 <- table[c(2:5)]
  colnames(plot_data1)[1:4] <- c('interacting_pair', 'cell_pairs', 'means', 'pvals')
  plot_data1$interacting_pair <- str_split(plot_data1$interacting_pair, ">@<>@<>@<", simplify = TRUE)[, 2]

  # Data preprocessing
  plot_data1$means[is.na(plot_data1$means)] <- 0
  plot_data1$pvals[is.na(plot_data1$pvals)] <- 1

  # Dot plot
  plot1 <- ggplot(plot_data1, aes(x = cell_pairs, y = interacting_pair)) +
    geom_point(aes(size = -log10(pvals), color = means)) +
    scale_size_continuous(range = c(0, 3), breaks = c(0, 1.0, 2.0)) +
    scale_color_gradientn('Mean expression', colors = my_palette1, limits = c(0, 1)) +
    theme_bw() +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 7),
          axis.text.y = element_text(size = 7, colour = "black"),
          axis.title = element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  print(plot1)

  # Heatmap
  plot_data2 <- plot_cpdb_heatmap(pvals = pvals, return_tables = TRUE,
                                  cellheight = 20, cellwidth = 20,
                                  cluster_cols = FALSE, cluster_rows = FALSE,
                                  cell_types = unique(sce@meta.data[,celltype_key]))$count_network
  pheatmap(plot_data2, color = my_palette2,
           cluster_rows = FALSE, cluster_cols = FALSE,
           border_color = "white",
           cellwidth = 30, cellheight = 30)

  # Chord diagram
  plot_cpdb3(scdata = as.SingleCellExperiment(sce),
             cell_type1 = cell_type1, cell_type2 = cell_type2,
             celltype_key = celltype_key, means = means, pvals = pvals,
             keep_significant_only = TRUE, deconvoluted = decon)

  dev.off()
}
