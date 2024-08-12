#' Process scRNA Data with Monocle
#'
#' This function processes a Seurat object using Monocle to perform dimensionality reduction and cell ordering.
#'
#' @param scRNA A Seurat Object.
#' @param save_path A path to save your results. Default is "./output_data/monocle2.Rdata".
#'
#' @return A CellDataSet object from Monocle.
#' @export
#'
#' @examples
#' \dontrun{
#' runMonocle2Analysis(scRNA = your_seurat_object, save_path = "./output_data/monocle2.Rdata")
#' }
runMonocle2Analysis <- function(scRNA, save_path = "./output_data/monocle2.Rdata") {
  # Check if necessary packages are installed
  required_packages <- c("monocle", "tidyverse", "Seurat")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Validate inputs
  if (!inherits(scRNA, "Seurat")) {
    stop("The 'scRNA' parameter must be a Seurat object.")
  }
  if (!is.character(save_path) || length(save_path) != 1) {
    stop("The 'save_path' parameter must be a single character string.")
  }

  # Create output directory if it does not exist
  output_dir <- dirname(save_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Load required libraries
  library(monocle)
  library(tidyverse)

  # Convert Seurat object to matrix
  data <- as.matrix(scRNA@assays$RNA@counts)

  # Create phenoData AnnotatedDataFrame
  pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)

  # Create featureData AnnotatedDataFrame
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

  # Create CellDataSet object
  mycds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())

  # Estimate size factors and dispersions
  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds, cores = 14, relative_expr = TRUE)

  # Get highly variable genes
  disp_table <- dispersionTable(mycds)
  disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

  # Set ordering filter
  mycds <- setOrderingFilter(mycds, disp.genes)

  # Plot ordering genes
  plot_ordering_genes(mycds)

  # Reduce dimensions and order cells
  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
  mycds <- orderCells(mycds)

  # Optionally, reverse the order of cells
  mycds_reverse <- orderCells(mycds, reverse = TRUE)

  # Save mycds object
  save(mycds, mycds_reverse, file = save_path)

  # Return the CellDataSet object
  return(mycds)
}

#' Plot scRNA Trajectories
#'
#' This function plots cell trajectories from a CellDataSet object and saves the plots to specified paths.
#'
#' @param mycds A CellDataSet object from Monocle.
#' @param output_dir Directory to save the output plots. Default is "./output_figure/".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plotScRNATrajectories(mycds = your_cds_object, output_dir = "./output_figure/")
#' }
plotScRNATrajectories <- function(mycds, output_dir = "./output_figure/") {
  # Check if necessary packages are installed
  required_packages <- c("monocle", "ggplot2")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Load required libraries
  library(monocle)
  library(ggplot2)

  # Plot cell trajectory by cell type
  plot1 <- plot_cell_trajectory(mycds, color_by = "cell_type", cell_size = 0.5, show_branch_points = FALSE) +
    theme(text = element_text(size = 14), legend.position = "right") +
    scale_color_manual(values = c('#6761A1', '#9BC4DA', '#E3C39F', '#4189B9'))
  ggsave(filename = file.path(output_dir, "monocle_celltype.pdf"), plot = plot1, width = 5, height = 3)

  # Plot cell trajectory by group
  plot2 <- plot_cell_trajectory(mycds, color_by = "group", cell_size = 0.5, show_branch_points = FALSE) +
    theme(text = element_text(size = 14), legend.position = "right") +
    scale_color_manual(values = c('#6761A1', '#9BC4DA', '#E3C39F', '#4189B9'))
  ggsave(filename = file.path(output_dir, "monocle_group.pdf"), plot = plot2, width = 5, height = 3)

  # Plot cell trajectory by pseudotime
  plot3 <- plot_cell_trajectory(mycds, cell_size = 1, show_branch_points = TRUE, color_by = "Pseudotime") +
    theme(legend.position = "top")
  ggsave(filename = file.path(output_dir, "monocle_pseudotime.pdf"), plot = plot3, width = 4, height = 3)
}


#' Generate Branch Enrichment Plot for monocle2
#'
#' This function generates a branch enrichment plot for single-cell RNA sequencing data using monocle2.
#'
#' @param sce A Seurat object.
#' @param mycds A CellDataSet object from monocle2.
#' @param cell_type_col A character string specifying the column name for cell types.
#' @param branch_point An integer specifying the branch point.
#' @param num_clusters An integer specifying the number of clusters. Default is 3.
#' @param top_n_markers An integer specifying the number of top markers to select. Default is 50.
#' @param cores An integer specifying the number of cores to use. Default is 1.
#' @param pvalue_cutoff A numeric value specifying the p-value cutoff for enrichment analysis. Default is 0.05.
#' @param topn_enrich An integer specifying the number of top enriched terms to display. Default is 8.
#' @param seed An integer specifying the random seed. Default is 5201314.
#' @param species A character string specifying the species ('human' or 'mouse'). Default is 'human'.
#' @param num_mark_genes An integer specifying the number of marker genes to randomly select. Default is 25.
#' @param pdf_path A character string specifying the path to save the PDF file. Default is './output_figure/branch-enrich.pdf'.
#' @param pdf_height A numeric value specifying the height of the PDF. Default is 9.
#' @param pdf_width A numeric value specifying the width of the PDF. Default is 16.
#' @param plot_type A character string specifying the type of plot. Default is "both".
#' @param column_names_rot A numeric value specifying the rotation angle for column names. Default is 45.
#' @param show_row_dend A logical value specifying whether to show row dendrogram. Default is FALSE.
#' @param markGenes_side A character string specifying the side to display marker genes. Default is "left".
#' @param go_colors A vector of colors for GO terms. Default is jjAnno::useMyCol("calm", n = 3).
#'
#' @return NULL. The function is called for its side effects.
#' @export
#'
#' @examples
#' \dontrun{
#' generateBranchEnrichmentPlot(sce = sce, mycds = mycds, cell_type_col = 'cell_type', branch_point = 1)
#' }
generateBranchEnrichmentPlot <- function(sce, mycds, cell_type_col, branch_point, num_clusters = 3, top_n_markers = 50,
                                         cores = 1, pvalue_cutoff = 0.05, topn_enrich = 8, seed = 5201314, species = 'human',
                                         num_mark_genes = 25, pdf_path = './output_figure/branch-enrich.pdf',
                                         pdf_height = 9, pdf_width = 16, plot_type = "both",
                                         column_names_rot = 45, show_row_dend = FALSE,
                                         markGenes_side = "left", go_colors = jjAnno::useMyCol("calm", n = 3)) {

  # Load necessary libraries
  required_packages <- c("org.Hs.eg.db", "org.Mm.eg.db", "ClusterGVis", "dplyr", "scutilsR", "monocle", "ggplot2", "jjAnno")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Validate inputs
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The 'sce' parameter must be a SingleCellExperiment object.")
  }
  if (!inherits(mycds, "CellDataSet")) {
    stop("The 'mycds' parameter must be a CellDataSet object from monocle2.")
  }
  if (!is.character(cell_type_col) || length(cell_type_col) != 1) {
    stop("The 'cell_type_col' parameter must be a single character string.")
  }
  if (!is.numeric(branch_point) || length(branch_point) != 1) {
    stop("The 'branch_point' parameter must be a single numeric value.")
  }
  if (!is.numeric(num_clusters) || length(num_clusters) != 1) {
    stop("The 'num_clusters' parameter must be a single numeric value.")
  }
  if (!is.numeric(top_n_markers) || length(top_n_markers) != 1) {
    stop("The 'top_n_markers' parameter must be a single numeric value.")
  }
  if (!is.numeric(cores) || length(cores) != 1) {
    stop("The 'cores' parameter must be a single numeric value.")
  }
  if (!is.numeric(pvalue_cutoff) || length(pvalue_cutoff) != 1) {
    stop("The 'pvalue_cutoff' parameter must be a single numeric value.")
  }
  if (!is.numeric(topn_enrich) || length(topn_enrich) != 1) {
    stop("The 'topn_enrich' parameter must be a single numeric value.")
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("The 'seed' parameter must be a single numeric value.")
  }
  if (!is.character(species) || !species %in% c('human', 'mouse')) {
    stop("The 'species' parameter must be either 'human' or 'mouse'.")
  }
  if (!is.numeric(num_mark_genes) || length(num_mark_genes) != 1) {
    stop("The 'num_mark_genes' parameter must be a single numeric value.")
  }
  if (!is.character(pdf_path) || length(pdf_path) != 1) {
    stop("The 'pdf_path' parameter must be a single character string.")
  }
  if (!is.numeric(pdf_height) || length(pdf_height) != 1) {
    stop("The 'pdf_height' parameter must be a single numeric value.")
  }
  if (!is.numeric(pdf_width) || length(pdf_width) != 1) {
    stop("The 'pdf_width' parameter must be a single numeric value.")
  }
  if (!is.character(plot_type) || length(plot_type) != 1) {
    stop("The 'plot_type' parameter must be a single character string.")
  }
  if (!is.numeric(column_names_rot) || length(column_names_rot) != 1) {
    stop("The 'column_names_rot' parameter must be a single numeric value.")
  }
  if (!is.logical(show_row_dend) || length(show_row_dend) != 1) {
    stop("The 'show_row_dend' parameter must be a single logical value.")
  }
  if (!is.character(markGenes_side) || length(markGenes_side) != 1) {
    stop("The 'markGenes_side' parameter must be a single character string.")
  }
  if (!is.vector(go_colors) || length(go_colors) == 0) {
    stop("The 'go_colors' parameter must be a non-empty vector.")
  }

  # Set cell type
  Idents(sce) <- cell_type_col

  # Find all marker genes
  cell_marker <- scutilsR::mcFindAllMarkers(sce)

  # Select top n marker genes for each cluster
  top <- cell_marker %>%
    group_by(cluster) %>%
    top_n(n = top_n_markers, wt = avg_log2FC)

  # Generate branch heatmap data
  df <- plot_genes_branched_heatmap2(mycds[unique(top$Gene.name.uniq),],
                                     branch_point = branch_point,
                                     num_clusters = num_clusters,
                                     cores = cores,
                                     use_gene_short_name = TRUE,
                                     show_rownames = TRUE)

  # Enrichment analysis
  OrgDb <- if (species == 'human') org.Hs.eg.db else org.Mm.eg.db
  organism <- if (species == 'human') 'hsa' else 'mmu'

  enrich <- enrichCluster(object = df, OrgDb = OrgDb,
                          type = "BP", organism = organism,
                          pvalueCutoff = pvalue_cutoff, topn = topn_enrich,
                          seed = seed)

  # Randomly select marker genes
  markGenes <- sample(unique(df$wide.res$gene), num_mark_genes, replace = FALSE)

  # Plot and save PDF
  pdf(pdf_path, height = pdf_height, width = pdf_width, onefile = FALSE)
  visCluster(object = df, plot.type = plot_type, column_names_rot = column_names_rot,
             show_row_dend = show_row_dend, markGenes = markGenes,
             markGenes.side = markGenes_side, annoTerm.data = enrich,
             go.col = c(rep(go_colors, each = topn_enrich)),
             add.bar = TRUE, line.side = markGenes_side)
  dev.off()

  message("Branch enrichment plot saved to ", pdf_path)
}
