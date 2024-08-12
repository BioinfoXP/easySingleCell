# =============== 1.QC ================

#' @title Filter and Plot scRNA-seq Data
#' @description This function filters scRNA-seq data based on specified criteria, generates violin plots for QC metrics, and saves the plots as PDF files.
#' @param sce A Seurat object containing scRNA-seq data.
#' @param minGene Minimum number of genes detected.
#' @param maxGene Maximum number of genes detected.
#' @param pctMT Maximum percentage of mitochondrial genes.
#' @param maxCounts Maximum number of RNA counts.
#' @param pal A character vector specifying the colors for the violin plots.
#' @param species A character string specifying the species ("human" or "mouse").
#' @param output_dir A string specifying the output directory for the plots. Default is "./output_figure/".
#' @param width Width of the PDF file. Default is 6.
#' @param height Height of the PDF file. Default is 4.
#' @return A filtered Seurat object. The function also saves the violin plots in PDF files.
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' sce <- CreateSeuratObject(counts = your_data)
#' pal <- c("#F1788D", "#54990F")
#' runScRNAQC(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000, species = 'human', pal = pal)
#' }

runScRNAQC <- function(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000,
                       species = 'human', pal = c("#F1788D", "#54990F","#E6550D","#843C39", "#3182BD","#8C6D31",
                                                           "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", "#6BAED6",
                                                           "#9ECAE1", "#AD494A", "#A1D99B", "#C7E9C0", "#99600F",
                                                           "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B",
                                                           "#66A61E","#E6550D", "#E7969C",'#53A85F'),
                                                           output_dir = "./output_figure/", width = 6, height = 4) {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Determine the mitochondrial gene pattern based on species
  if (species == 'human') {
    mt_pattern <- "^MT-"
  } else if (species == 'mouse') {
    mt_pattern <- "^mt-"
  } else {
    stop("Unsupported species. Please specify 'human' or 'mouse'.")
  }

  # Calculate percentage of mitochondrial genes
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = mt_pattern)

  # Generate initial violin plots
  p1 <- VlnPlot(sce, assay = 'RNA', features = c("nCount_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
  p2 <- VlnPlot(sce, assay = 'RNA', features = c("nFeature_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
  p3 <- VlnPlot(sce, assay = 'RNA', features = c("percent.mt"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()

  # Save initial violin plots to PDF files
  ggsave(filename = paste0(output_dir, 'BeforeQC-nCount_RNA.pdf'), plot = p1, width = width, height = height)
  ggsave(filename = paste0(output_dir, 'BeforeQC-nFeature.pdf'), plot = p2, width = width, height = height)
  ggsave(filename = paste0(output_dir, 'BeforeQC-percent-mt.pdf'), plot = p3, width = width, height = height)

  # Filter the data
  scRNA_filtered <- subset(sce, subset = nFeature_RNA > minGene &
                             nFeature_RNA < maxGene &
                             percent.mt < pctMT & nCount_RNA < maxCounts)

  # Generate violin plots after filtering
  p1_filtered <- VlnPlot(scRNA_filtered, assay = 'RNA', features = c("nCount_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
  p2_filtered <- VlnPlot(scRNA_filtered, assay = 'RNA', features = c("nFeature_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
  p3_filtered <- VlnPlot(scRNA_filtered, assay = 'RNA', features = c("percent.mt"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()

  # Save the plots to PDF files
  ggsave(filename = paste0(output_dir, 'AfterQC-nCount_RNA.pdf'), plot = p1_filtered, width = width, height = height)
  ggsave(filename = paste0(output_dir, 'AfterQC-nFeature.pdf'), plot = p2_filtered, width = width, height = height)
  ggsave(filename = paste0(output_dir, 'AfterQC-percent-mt.pdf'), plot = p3_filtered, width = width, height = height)

  return(scRNA_filtered)
}

# Example usage
# sce <- CreateSeuratObject(counts = your_data)
# pal <- c("#F1788D", "#54990F")
# runScRNAQC(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000, species = 'human', pal = pal)


# Example usage
# sce <- CreateSeuratObject(counts = your_data)
# pal <- c("#F1788D", "#54990F")
# runScRNAQC(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000, species = 'human', pal = pal)





# ============== 2. Run standard normalize ==============
#' @title Run Normalization on a Seurat Object
#' @description This function normalizes a Seurat object and performs dimensionality reduction.
#' @param seurat_obj A Seurat object containing the single-cell RNA data.
#' @param dims A numeric vector specifying the dimensions to use in the normalization process. Default is 1:30.
#' @param batch_var A character string specifying the batch variable to use for Harmony integration. Default is "orig.ident".
#' @return A Seurat object with normalized data and dimensionality reduction.
#' @export
#' @import Seurat
#' @importFrom harmony RunHarmony
#' @examples
#' \dontrun{
#' normalized_sce <- run_normalize(seurat_obj = sce, dims = 1:30, batch_var = "batch")
#' }

run_normalize <- function(seurat_obj, dims = 1:30, batch_var = "orig.ident") {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  scale_genes <- VariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = scale_genes)
  seurat_obj <- RunPCA(seurat_obj, features = scale_genes)
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = batch_var)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, reduction = "harmony")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction = "harmony")
  return(seurat_obj)
}

# ============== 3. Resolution by clusterTree ==============
#' @title Generate and Save Clustree Plot
#' @description This function generates and saves a clustree plot for a Seurat object across specified clustering resolutions.
#' @param sce A Seurat object containing scRNA-seq data.
#' @param resolutions A numeric vector specifying the clustering resolutions to analyze. Default is seq(0.1, 1.4, 0.2).
#' @param prefix A character string specifying the prefix for clustering resolution columns. Default is "RNA_snn_res.".
#' @param output_dir A character string specifying the directory to save the output. Default is "./output_figure/".
#' @param output_filename A character string specifying the output filename for the clustree plot. Default is "clustree_plot.pdf".
#' @return A ggplot object.
#' @export
#' @import Seurat
#' @import clustree
#' @import ggplot2
#' @examples
#' \dontrun{
#' sce <- qread('./output_data/sce.qc.qs')
#' clustree_plot(sce)
#' }

clustree_plot <- function(sce, resolutions = seq(0.1, 1.4, 0.2), prefix = "RNA_snn_res.",
                          output_dir = "./output_figure/", output_filename = "clustree_plot.pdf") {
  # Ensure necessary libraries are loaded
  library(clustree)
  library(ggplot2)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Loop through each resolution and perform clustering
  for (i in resolutions) {
    sce <- FindClusters(object = sce, resolution = i)
  }

  # Generate clustree plot
  clustree_plot <- clustree(sce, prefix = prefix)

  # Save the plot
  ggsave(filename = file.path(output_dir, output_filename), plot = clustree_plot, width = 10, height = 10)

  # Return the plot
  return(clustree_plot)
}

# Example usage
# sce <- qread('./output_data/sce.qc.qs')
# clustree_plot(sce)



# ============== 4. Mark Doublets ==============
#' @title Mark Doublets in Seurat Object
#' @description This function wraps the MarkDoublets function from the DoubletFinder package to identify doublets in a Seurat object.
#' @param seu A Seurat object.
#' @param PCs A numeric vector specifying the principal components to use for doublet detection. Default is 1:10.
#' @param split.by A character string specifying the column in meta.data to split the analysis by. Default is "orig.ident".
#' @return A Seurat object with doublet information added.
#' @export
#' @import Seurat
#' @import DoubletFinder
#' @examples
#' \dontrun{
#' seu <- MarkDoubletsWrapper(seu = seu, PCs = 1:10, split.by = "orig.ident")
#' }

MarkDoubletsWrapper <- function(seu, PCs = 1:10, split.by = "orig.ident") {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(DoubletFinder)
  library(scutilsR)

  # Call the MarkDoublets function
  seu <- MarkDoublets(seu = seu, PCs = PCs, split.by = split.by)

  return(seu)
}

# Example usage
# seu <- MarkDoubletsWrapper(seu = seu, PCs = 1:10, split.by = "orig.ident")


# ============ 5.Plot umap/tsne Reduction =====
#' @title Plot UMAP or t-SNE for scRNA-seq Data
#' @description This function plots UMAP or t-SNE for scRNA-seq data using Seurat's DimPlot function.
#' @param sce A Seurat object containing scRNA-seq data.
#' @param reduction A character string specifying the reduction method ("umap" or "tsne"). Default is "umap".
#' @param group.by A character string specifying the grouping variable. Default is 'celltype'.
#' @param label A logical value indicating whether to label the clusters. Default is TRUE.
#' @param repel A logical value indicating whether to use label repelling. Default is TRUE.
#' @param pt.size A numeric value specifying the point size. Default is 0.01.
#' @param cols A character vector specifying the colors for the plot. Default is a predefined color palette.
#' @return A ggplot object.
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' sce <- CreateSeuratObject(counts = your_data)
#' plotScRNAReduction(sce, reduction = "umap", group.by = 'celltype')
#' plotScRNAReduction(sce, reduction = "tsne", group.by = 'celltype')
#' }

plotScRNAReduction <- function(sce, reduction = "umap", group.by = 'celltype', label = TRUE, repel = TRUE, pt.size = 0.01,
                               cols = c("#F1788D", "#54990F", "#E6550D", "#843C39", "#3182BD", "#8C6D31",
                                                 "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", "#6BAED6",
                                                 "#9ECAE1", "#AD494A", "#A1D99B", "#C7E9C0", "#99600F",
                                                 "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B",
                                                 "#66A61E", "#E6550D", "#E7969C", '#53A85F')) {
                                                 # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Check if the reduction method is supported
  if (!reduction %in% c("umap", "tsne")) {
    stop("Unsupported reduction method. Please specify 'umap' or 'tsne'.")
  }

  # Plot UMAP or t-SNE
  plot <- Seurat::DimPlot(sce, reduction = reduction, group.by = group.by, label = label, repel = repel, pt.size = pt.size, cols = cols)

  return(plot)
}

# Example usage
# sce <- CreateSeuratObject(counts = your_data)
# plotScRNAReduction(sce, reduction = "umap", group.by = 'celltype')
# plotScRNAReduction(sce, reduction = "tsne", group.by = 'celltype')

# ============ 6.Plot Feature in umap/tsne =====
#' @title Plot Feature Expression on UMAP or t-SNE for scRNA-seq Data
#' @description This function plots feature expression on UMAP or t-SNE for scRNA-seq data using Seurat's FeaturePlot function.
#' @param sce A Seurat object containing scRNA-seq data.
#' @param feature A character(genes) string specifying the feature to plot.
#' @param reduction A character string specifying the reduction method ("umap", "tsne", etc.). Default is "umap".
#' @param pt.size A numeric value specifying the point size. Default is 0.0001.
#' @param max.cutoff A numeric value specifying the maximum cutoff for the feature. Default is 1.5.
#' @param cols A character vector specifying the colors for the plot. Default is a predefined color palette.
#' @param title A character string specifying the title of the plot.
#' @return A ggplot object.
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' sce <- CreateSeuratObject(counts = your_data)
#' plotScRNAFeaturePlot(sce, feature = "GeneX", reduction = "umap", title = "Feature Plot")
#' plotScRNAFeaturePlot(sce, feature = "GeneY", reduction = "tsne", title = "Feature Plot")
#' }

plotScRNAFeaturePlot <- function(sce, feature, reduction = "umap", pt.size = 0.0001, max.cutoff = 1.5,
                                 cols = c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"), title) {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Plot feature expression
  plot <- FeaturePlot(
    object = sce,
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
# sce <- CreateSeuratObject(counts = your_data)
# plotScRNAFeaturePlot(sce, feature = "GeneX", reduction = "umap", title = "Feature Plot")
# plotScRNAFeaturePlot(sce, feature = "GeneY", reduction = "tsne", title = "Feature Plot")












