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
#' run_clustree(sce)
#' }

run_clustree <- function(sce, resolutions = seq(0.1, 1.4, 0.2), prefix = "RNA_snn_res.",
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
# run_clustree(sce)



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


















#  ================= 5. runScMetabolismAnalysis ===========
#' @title Perform Metabolism Analysis and Generate DotPlot
#' @description This function performs metabolism analysis on Seurat object data and generates a DotPlot of metabolic pathways.
#' @param scRNA A Seurat object containing the single-cell RNA-seq data.
#' @param output_file Path to the output Rdata file for saving the Seurat object with metabolism scores. Default is './output_data/Figure4/scMetabolism.Rdata'.
#' @param pdf_file Path to the output PDF file for saving the DotPlot. Default is './output_figure/Figure2/scMetabolism.pdf'.
#' @param npathways Number of top pathways to plot. Default is 20.
#' @param cell_type_column Column name in the Seurat object metadata representing cell types. Default is "cell_type".
#' @param width Width of the output PDF. Default is 6.
#' @param height Height of the output PDF. Default is 10.
#' @param method Method for metabolism score calculation. Default is "VISION".
#' @return A Seurat object with metabolism scores added.
#' @export
#' @import scMetabolism
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' scRNA <- readRDS('./path_to_your_seurat_object.rds')
#' runScMetabolismAnalysis(
#'   scRNA = scRNA,
#'   npathways = 20,
#'   cell_type_column = "cell_type"
#' )
#' }

runScMetabolismAnalysis <- function(scRNA,
                                    output_file = './output_data/Figure4/scMetabolism.Rdata',
                                    pdf_file = './output_figure/Figure2/scMetabolism.pdf',
                                    npathways = 20,
                                    cell_type_column = "cell_type",
                                    width = 6,
                                    height = 10,
                                    method = "VISION") {

  library(scMetabolism)
  library(Seurat)
  library(ggplot2)

  # Compute metabolism pathway scores
  scRNA <- sc.metabolism.Seurat(
    obj = scRNA,
    method = method,
    imputation = FALSE,
    ncores = 10,
    metabolism.type = "KEGG"
  )

  # Save Seurat object
  save(scRNA, file = output_file)

  # Generate DotPlot for metabolic pathways
  dotplot <- DotPlot.metabolism(
    obj = scRNA,
    pathway = rownames(scRNA@assays[["METABOLISM"]][["score"]])[1:npathways],
    phenotype = cell_type_column,
    norm = "y"
  ) + xlab('')

  # Save plot as PDF
  ggsave(pdf_file, plot = dotplot, width = width, height = height)

  return(scRNA)
}

# ==========  6.runScGSEA  ==========

#' @title Run GSEA Analysis For ScRNA
#' @description This function performs Gene Set Enrichment Analysis (GSEA) on Seurat object data.
#' @param scRNA A Seurat object containing the single-cell RNA-seq data.
#' @param ident1 The identity class to be used as the reference group.
#' @param ident2 The identity class to be compared against the reference group.
#' @param logfc_threshold The log fold change threshold for identifying differentially expressed genes. Default is 0.1.
#' @param gmt_paths A list of paths to GMT files for GSEA analysis. Default includes Hallmark, GO BP, and KEGG pathways for human.
#' @param pvalue_cutoff The p-value cutoff for GSEA results. Default is 0.05.
#' @param p_adjust_method The method for p-value adjustment. Default is 'none'.
#' @param minGSSize Minimum size of each gene set for analyzing. Default is 10.
#' @param maxGSSize Maximum size of each gene set for analyzing. Default is 500.
#' @param nPerm Number of permutations. Default is 1000.
#' @return A list containing GSEA results for each GMT file and the differentially expressed genes.
#' @export
#' @import Seurat
#' @import dplyr
#' @import clusterProfiler
#' @import fgsea
#' @examples
#' \dontrun{
#' scRNA <- readRDS('./path_to_your_seurat_object.rds')
#' gsea_results <- runScGSEA(
#'   scRNA = scRNA,
#'   ident1 = "cell_type_1",
#'   ident2 = "cell_type_2"
#' )
#' }

runScGSEA <- function(scRNA, ident1, ident2, logfc_threshold = 0.1,
                      gmt_paths = list(system.file("data", "h.all.v7.4.symbols.gmt", package = "easySingleCell"),
                                       system.file("data", "c5.go.bp.v7.4.symbols.gmt", package = "easySingleCell"),
                                       system.file("data", "c2.cp.kegg.v7.4.symbols.gmt", package = "easySingleCell")
                      ),
                      pvalue_cutoff = 0.05, p_adjust_method = 'none',
                      minGSSize = 10, maxGSSize = 500, nPerm = 1000) {

  # Load necessary libraries
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(fgsea)

  # Find differentially expressed genes
  sce_edg <- FindMarkers(scRNA, ident.1 = ident1, ident.2 = ident2, logfc.threshold = logfc_threshold)

  # Add gene symbols
  sce_edg$SYMBOL <- rownames(sce_edg)
  names(sce_edg)[names(sce_edg) == "avg_log2FC"] <- "logFC"
  sce_edg <- sce_edg %>% arrange(desc(logFC))

  # Create gene list
  geneList <- sce_edg$logFC
  names(geneList) <- sce_edg$SYMBOL

  # Read GMT files
  gmt_list <- lapply(gmt_paths, function(path) {
    if (file.exists(path)) {
      return(read.gmt(path))
    } else {
      stop(paste("File not found:", path))
    }
  })

  # Run GSEA analysis
  gsea_results <- lapply(gmt_list, function(gmt) {
    GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_adjust_method,
         minGSSize = minGSSize, maxGSSize = maxGSSize, nPerm = nPerm)
  })

  # Get GSEA results
  gsea_results_list <- lapply(gsea_results, function(gsea) gsea@result)

  return(list(gsea_results_list, sce_edg))
}

# Example call
# scRNA <- readRDS('./path_to_your_seurat_object.rds')
# gsea_results <- runScGSEA(
#   scRNA = scRNA,
#   ident1 = "cell_type_1",
#   ident2 = "cell_type_2"
# )



# ==========  7.runCytoTRACEAnalysis  ==========
#' Run CytoTRACE Analysis
#'
#' This function runs CytoTRACE analysis on a Seurat object and saves the results and plots.
#'
#' @param scRNA A Seurat object.
#' @param celltype A character string specifying the column name in `scRNA@meta.data` that contains cell type information. Default is 'celltype'.
#' @param output_data_dir Directory to save the output data. Default is "./output_data/".
#' @param output_figure_dir Directory to save the output figures. Default is "./output_data/".
#' @param ncores Number of cores to use for CytoTRACE analysis. Default is 5.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' runCytoTRACEAnalysis(scRNA = your_seurat_object, celltype = 'celltype', output_data_dir = "./output_data/", output_figure_dir = "./output_data/", ncores = 5)
#' }
runCytoTRACEAnalysis <- function(scRNA, celltype = 'celltype',
                                 output_data_dir = "./output_data/", output_figure_dir = "./output_data/", ncores = 5) {
  # Check if necessary packages are installed
  required_packages <- c("Seurat", "CytoTRACE")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Create output directories if they do not exist
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }
  if (!dir.exists(output_figure_dir)) {
    dir.create(output_figure_dir, recursive = TRUE)
  }

  # Load required libraries
  library(Seurat)
  library(CytoTRACE)

  # Validate inputs
  if (!inherits(scRNA, "Seurat")) {
    stop("The 'scRNA' parameter must be a Seurat object.")
  }
  if (!is.character(celltype) || length(celltype) != 1) {
    stop("The 'celltype' parameter must be a single character string.")
  }
  if (!celltype %in% colnames(scRNA@meta.data)) {
    stop(paste("The specified celltype column '", celltype, "' does not exist in scRNA@meta.data.", sep = ""))
  }
  if (!is.character(output_data_dir) || length(output_data_dir) != 1) {
    stop("The 'output_data_dir' parameter must be a single character string.")
  }
  if (!is.character(output_figure_dir) || length(output_figure_dir) != 1) {
    stop("The 'output_figure_dir' parameter must be a single character string.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0) {
    stop("The 'ncores' parameter must be a single positive number.")
  }

  # Run CytoTRACE analysis
  phe <- scRNA@meta.data[, celltype]
  phe <- as.character(phe)
  names(phe) <- rownames(scRNA@meta.data)

  mat_3k <- as.matrix(scRNA@assays$RNA@counts)

  results <- CytoTRACE(mat = mat_3k, ncores = ncores)
  save(results, file = file.path(output_data_dir, "Cytotrace_results.Rdata"))

  plotCytoTRACE(results, outputDir = output_figure_dir, emb = scRNA@reductions$umap@cell.embeddings)
}


# ========= 8. runCogapsAnalysis =========
#' @title Run CoGAPS Analysis on Single-cell Data
#' @description This function runs CoGAPS analysis on single-cell data for a specified pattern number and saves the results.
#' @param scRNA A Seurat object containing single-cell data.
#' @param nPatterns A numeric value specifying the number of patterns to be used in CoGAPS analysis.
#' @param out_dir A character string specifying the output directory where results will be saved. Default is './output_data/cogaps'.
#' @param nIterations A numeric value specifying the number of iterations for CoGAPS (default is 1000).
#' @param seed A numeric value specifying the random seed for reproducibility (default is 42).
#' @param nThreads A numeric value specifying the number of threads to be used (default is 4).
#' @return The CoGAPS result object.
#' @export
#' @import Seurat
#' @import CoGAPS
#' @import parallel
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(CoGAPS)
#' library(parallel)
#'
#' # Example Seurat object
#' seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' seurat_obj <- NormalizeData(seurat_obj)
#' seurat_obj <- FindVariableFeatures(seurat_obj)
#' seurat_obj <- ScaleData(seurat_obj)
#' seurat_obj <- RunPCA(seurat_obj)
#' seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#'
#' # Run CoGAPS analysis
#' runCogapsAnalysis(
#'   scRNA = seurat_obj,
#'   nPatterns = 4,
#'   out_dir = './output_data/cogaps',
#'   nIterations = 1000,
#'   seed = 42,
#'   nThreads = 5
#' )
#' }

runCogapsAnalysis <- function(scRNA, nPatterns = 3, out_dir = './output_data/cogaps',
                              nIterations = 1000, seed = 42, nThreads = 5) {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(CoGAPS)
  library(parallel)
  library(magrittr)

  # Extract count matrix and normalize
  counts <- GetAssayData(object = scRNA, slot = "counts") %>% as.matrix()
  norm_counts <- log1p(counts)

  # Create output directory if it does not exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Set CoGAPS parameters
  params <- CogapsParams(
    nIterations = nIterations,
    seed = seed,
    nPatterns = nPatterns,
    sparseOptimization = TRUE,
    distributed = "genome-wide"
  )

  # Run CoGAPS analysis
  cogaps_result <- CoGAPS(norm_counts, params, nThreads = nThreads)

  # Save the result
  saveRDS(cogaps_result, file = file.path(out_dir, paste0("cogaps_result_k_", nPatterns, ".rds")))

  # Return the result object
  return(cogaps_result)
}



# ========= 9. runCogapsAnalysis =========

#' Process scRNA Data with Monocle
#'
#' This function processes a Seurat object using Monocle2 to perform dimensionality reduction and cell ordering.
#'
#' @param scRNA A Seurat Object.
#' @param save_path A path to save your results. Default is "./output_data/monocle2.Rdata".
#'
#' @return A CellDataSet object from Monocle.
#' @export
#'
#' @examples
#' \dontrun{
#' runMonocleAnalysis(scRNA = your_seurat_object, save_path = "./output_data/monocle2.Rdata")
#' }
runMonocleAnalysis <- function(scRNA, save_path = "./output_data/monocle2.Rdata") {
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
#' @param celltype_col Column name in `mycds` pData for cell type. Default is "celltype".
#' @param group_col Column name in `mycds` pData for group. Default is "group".
#' @param cell_size Size of the cells in the plot. Default is 0.5.
#' @param show_branch_points Logical, whether to show branch points in the trajectory. Default is FALSE.
#' @param color_palette A vector of colors to use for the plots. Default is a predefined set of colors.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plotMonocleTrajectories(mycds = your_cds_object, output_dir = "./output_figure/")
#' }
plotMonocleTrajectories <- function(mycds, output_dir = "./output_figure/",
                                  celltype_col = "celltype", group_col = "group",
                                  cell_size = 0.5, show_branch_points = FALSE,
                                  color_palette = c('#6761A1', '#9BC4DA', '#E3C39F', '#4189B9')) {
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
  plot1 <- plot_cell_trajectory(mycds, color_by = celltype_col, cell_size = cell_size, show_branch_points = show_branch_points) +
    theme(text = element_text(size = 14), legend.position = "right") +
    scale_color_manual(values = color_palette)
  ggsave(filename = file.path(output_dir, "monocle_celltype.pdf"), plot = plot1, width = 5, height = 3)

  # Plot cell trajectory by group
  plot2 <- plot_cell_trajectory(mycds, color_by = group_col, cell_size = cell_size, show_branch_points = show_branch_points) +
    theme(text = element_text(size = 14), legend.position = "right") +
    scale_color_manual(values = color_palette)
  ggsave(filename = file.path(output_dir, "monocle_group.pdf"), plot = plot2, width = 5, height = 3)

  # Plot cell trajectory by pseudotime
  plot3 <- plot_cell_trajectory(mycds, cell_size = cell_size, show_branch_points = TRUE, color_by = "Pseudotime") +
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
#' @import Seurat
#' @import ClusterGVis
#' @import clusterProfiler
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' plotMonocleBranchEnrich(sce = sce, mycds = mycds, cell_type_col = 'cell_type', branch_point = 1)
#' }
plotMonocleBranchEnrich <- function(sce, mycds, cell_type_col, branch_point, num_clusters = 3, top_n_markers = 50,
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

  library(future);library(future.apply)
  plan(multisession,workers =1)


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

# =================== 10.runCellChatAnalysis ===================
#' Run CellChat Analysis
#'
#' This function runs CellChat analysis on a Seurat object. The analysis can be performed on the entire dataset or on multiple subsets of the data, depending on the specified grouping column. The function normalizes and clusters the data, and saves the results.
#'
#' @param sce A Seurat object.
#' @param celltype A character string specifying the column name in `sce@meta.data` that contains cell type information. Default is 'cell_type'.
#' @param groupby A character string specifying the column name in `sce@meta.data` to split the object by. If NULL, the object will not be split. Default is NULL.
#' @param species A character string specifying the species. Must be either 'human' or 'mouse'.
#' @param output_data_dir Directory to save the output data. Default is "./output_data/".
#' @param ncores Number of cores to use for parallel processing. Default is 1.
#' @param workers Number of workers to use for future processing. Default is 5.
#' @param min_cells Minimum number of cells to filter communication. Default is 10.
#'
#' @export
#' @import CellChat
#' @examples
#' \dontrun{
#' runCellChatAnalysis(sce = your_seurat_object, celltype = 'cell_type', groupby = 'group', species = 'human', output_data_dir = "./output_data/", ncores = 4, workers = 5, min_cells = 10)
#' }
runCellChatAnalysis <- function(sce, celltype = 'cell_type', groupby = NULL, species = 'human',
                                output_data_dir = "./output_data/", ncores = 1, workers = 5, min_cells = 10) {
  # Check if necessary packages are installed
  required_packages <- c("Seurat", "CellChat", "future", "dplyr", "qs")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  })

  # Create output directory if it does not exist
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }

  # Load required libraries
  library(Seurat)
  library(CellChat)
  library(future)
  library(dplyr)
  library(qs)

  # Validate inputs
  if (!inherits(sce, "Seurat")) {
    stop("The 'sce' parameter must be a Seurat object.")
  }
  if (!is.character(celltype) || length(celltype) != 1) {
    stop("The 'celltype' parameter must be a single character string.")
  }
  if (!celltype %in% colnames(sce@meta.data)) {
    stop(paste("The specified celltype column '", celltype, "' does not exist in sce@meta.data.", sep = ""))
  }
  if (!is.null(groupby) && (!is.character(groupby) || length(groupby) != 1)) {
    stop("The 'groupby' parameter must be a single character string or NULL.")
  }
  if (!is.null(groupby) && !groupby %in% colnames(sce@meta.data)) {
    stop(paste("The specified groupby column '", groupby, "' does not exist in sce@meta.data.", sep = ""))
  }
  if (!is.character(species) || length(species) != 1 || !(species %in% c("human", "mouse"))) {
    stop("The 'species' parameter must be a single character string and must be either 'human' or 'mouse'.")
  }
  if (!is.character(output_data_dir) || length(output_data_dir) != 1) {
    stop("The 'output_data_dir' parameter must be a single character string.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0) {
    stop("The 'ncores' parameter must be a single positive number.")
  }
  if (!is.numeric(workers) || length(workers) != 1 || workers <= 0) {
    stop("The 'workers' parameter must be a single positive number.")
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells <= 0) {
    stop("The 'min_cells' parameter must be a single positive number.")
  }

  # Select appropriate CellChat database and PPI data
  if (species == "human") {
    CellChatDB <- CellChatDB.human
    PPI <- PPI.human
  } else if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
    PPI <- PPI.mouse
  }

  # Split object by group if specified
  if (!is.null(groupby)) {
    sce_list <- SplitObject(sce, split.by = groupby)
    names_list <- names(sce_list)
  } else {
    sce_list <- list(sce)
    names_list <- "all"
  }

  # Normalize and cluster
  sce_list <- mclapply(1:length(sce_list), function(i) {
    sce_list[[i]] %>%
      NormalizeData() %>%
      FindClusters()
  }, mc.cores = ncores)

  # Run CellChat analysis
  cellchat_res <- mclapply(1:length(sce_list), function(i) {
    data.input <- sce_list[[i]]@assays$RNA@data
    identity <- sce_list[[i]]@meta.data
    Idents(sce_list[[i]]) <- celltype

    cellchat <- createCellChat(data.input, meta = identity, group.by = celltype)
    cellchat <- addMeta(cellchat, meta = identity)
    cellchat <- setIdent(cellchat, ident.use = celltype)
    groupSize <- as.numeric(table(cellchat@idents))

    # Use the selected CellChat database
    cellchat@DB <- CellChatDB

    # Subset data and identify overexpressed genes and interactions
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = workers)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = min_cells)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    return(cellchat)
  }, mc.cores = ncores)

  # Save results
  names(cellchat_res) <- names_list
  save(sce_list, cellchat_res, file = file.path(output_data_dir, "cellchat_res.Rdata"))
}



# ================== 11. runHdWGCNAStep1 ===========
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


# ================== 12. runHdWGCNAStep2 ===========
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


# ============== 13.runMistyRAnalysis ===========
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


# =============== 14. runInferCNV ==========
#' @title Prepare Data for InferCNV
#'
#' @description  This function prepares the single-cell expression data for InferCNV analysis.
#'
#' @param sce_epi Single-cell expression data for the tumor epithelial cells.
#' @param sce_refer Single-cell expression data for the reference cells (normal epithelial or immune cells).
#' @param celltype Column name in meta.data that contains cell type information.
#' @param infercnv_path Path to save InferCNV results.
#' @param name Name for the InferCNV run.
#' @return A list containing the merged single-cell expression data and the path to saved files.
#' @export
#' @import qs
#' @import dplyr
#' @import Seurat
#' @example
#' \dontrun{
#' # Ensure necessary libraries are loaded
#' library(infercnv)
#' library(qs)
#' # Example usage:
#' # Prepare data using a custom function
#' prepared_data <- PrepareDataForInferCNV(
#'   sce_epi = sce_epi,
#'   sce_refer = sce_refer,
#'   celltype = 'celltype',
#'   infercnv_path = './output_data/inferCNV/',
#'   name = 'infer_run'
#' )
#' }
PrepareDataForInferCNV <- function(sce_epi, sce_refer, celltype = 'celltype',
                                   infercnv_path = './output_data', name = 'infer_run') {
  dir.create(infercnv_path, recursive = TRUE, showWarnings = FALSE)

  # Normalize and cluster reference cells
  sce_refer <- sce_refer %>%
    NormalizeData() %>%
    FindClusters()

  # Merge and normalize data
  sce_infer <- merge(sce_epi, sce_refer) %>%
    NormalizeData() %>%
    FindClusters()

  # Save count matrix and cell type labels
  qsave(as.matrix(sce_infer[["RNA"]]@counts), file = paste0(infercnv_path, name, '.qs'))
  write.table(sce_infer@meta.data[, celltype], file = paste0(infercnv_path, name, '.celltype.label.txt'),
              sep = "\t", quote = FALSE, col.names = FALSE)

  return(list(sce_infer = sce_infer, infercnv_path = infercnv_path, name = name))
}






#' Run InferCNV Analysis
#'
#' This function performs InferCNV analysis on the prepared single-cell expression data.
#'
#' @param infercnv_path Path to the previous prepared saved InferCNV results. eg: './output_data', please don't use './output_data/'.
#' @param gene_order_file Path to the gene order file. It will provided by this package default.
#' @param ref_group_names A vector for reference group names for normal cells.
#' @param name Name for the InferCNV run.
#' @param num_threads Number of threads to use for the analysis.
#' @return The InferCNV object containing the analysis results.
#' @export
#' @import infercnv
#' @import qs
#' @examples
#' \dontrun{
#' # Ensure necessary libraries are loaded
#' library(infercnv)
#' library(qs)
#'
#' # Example usage:
#' # Prepare data using a custom function
#' prepared_data <- PrepareDataForInferCNV(
#'   sce_epi = sce_epi,
#'   sce_refer = sce_refer,
#'   celltype = 'cell_type',
#'   infercnv_path = './output_data/inferCNV/',
#'   name = 'infer_run'
#' )
#'
#' # Run InferCNV analysis
#' infercnv_obj <- RunInferCNVAnalysis(
#'   infercnv_path = prepared_data$infercnv_path,
#'   gene_order_file = system.file("extdata", "hg38_gencode_v27.txt", package = "easySingleCell"),
#'   ref_group_names = c("T/NK", "B Cell"),
#'   name = prepared_data$name,
#'   num_threads = 16
#' )
#'
#' # Please check your input! This function will do the followings:
#' # Load data for InferCNV
#' # matrix_counts: path to the .qs file containing the raw counts matrix
#' # annotations_file: path to the .txt file containing cell type annotations
#' # out_path: directory to save the output of the InferCNV analysis
#' # matrix_counts <- qread(paste0(infercnv_path, '/', name, '.qs'))
#' # annotations_file <- paste0(infercnv_path, '/', name, '.celltype.label.txt')
#' # out_path <- paste0(infercnv_path, '/', name)
#' }
RunInferCNVAnalysis <- function(infercnv_path = './output_data',
                                gene_order_file = system.file("extdata", "hg38_gencode_v27.txt", package = "easySingleCell"),
                                ref_group_names,
                                name = 'infer_run',
                                num_threads = 12) {
  # Load data for InferCNV
  # matrix_counts: path to the .qs file containing the raw counts matrix
  # annotations_file: path to the .txt file containing cell type annotations
  # out_path: directory to save the output of the InferCNV analysis
  matrix_counts <- qread(paste0(infercnv_path, '/', name, '.qs'))
  annotations_file <- paste0(infercnv_path, '/', name, '.celltype.label.txt')
  out_path <- paste0(infercnv_path, '/', name)

  # Create InferCNV object
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = matrix_counts,
                                       annotations_file = annotations_file,
                                       delim = "\t",
                                       gene_order_file = gene_order_file,
                                       ref_group_names = ref_group_names,
                                       chr_exclude = c("chrY", "chrM"))

  # Run InferCNV analysis
  infercnv_obj <- infercnv::run(infercnv_obj, cutoff = 0.1, out_dir = out_path,
                                no_prelim_plot = TRUE, cluster_by_groups = TRUE,
                                denoise = TRUE, HMM = FALSE, min_cells_per_gene = 10,
                                num_threads = num_threads, write_expr_matrix = TRUE)

  # Save the InferCNV object
  saveRDS(infercnv_obj, file = paste0(out_path, '/infercnv_obj.rds'))

  return(infercnv_obj)
}






#' Cluster and Visualize InferCNV Results
#'
#' This function performs clustering on the InferCNV results and generates a heatmap.
#'
#' @param infercnv_obj The InferCNV object containing the analysis results.
#' @param gene_order_file Path to the gene order file.
#' @param ref_cell_name Reference cell names for clustering and visualization.
#' @param obs_cell_name Observed cell names for clustering and visualization.
#' @param ref_group Reference group label for clustering and visualization.
#' @param obs_group Observed group label for clustering and visualization.
#' @param k_clusters Number of clusters for k-means clustering (default is 5).
#' @param heatmap_colors Colors for the heatmap (default is c("#2166ac", "white", "#b2182b")).
#' @param cluster_colors Colors for the clusters (default is c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#843C39")).
#' @param out_path Path to save the heatmap and clustering results.
#' @return A data frame containing the clustering results.
#' @export
#' @import ComplexHeatmap
#' @import dplyr
#' @import grid
#' @import scales
#' @import ggplot2
#' @import ggpubr
#' @import ggsignif
#' @import tidyverse
#' @import ggprism
#' @import RColorBrewer
#' @examples
#' \dontrun{
#' library(infercnv)
#' library(ComplexHeatmap)
#' library(dplyr)
#' library(grid)
#' library(scales)
#' library(ggplot2)
#' library(ggpubr)
#' library(ggsignif)
#' library(tidyverse)
#' library(ggprism)
#' library(vioplot)
#' library(RColorBrewer)
#'
#' Example usage:
#' Assume `infercnv_obj` is the result from `RunInferCNVAnalysis`
#' Kmeas.res <- RunInferCNVCluster(
#'   infercnv_obj = infercnv_obj,
#'   gene_order_file = system.file("extdata", "hg38_gencode_v27.txt", package = "easySingleCell"),
#'   ref_cell_name = c("T/NK","B Cell"),#Or normal epithlial
#'   obs_cell_name = "epithlial",#epithlial in tumor sample
#'   ref_group = "immune",# or normal epi
#'   obs_group = "epithelial",# epithlial in tumor sample
#'   k_clusters = 5,
#'   heatmap_colors = c("#2166ac", "white", "#b2182b"),
#'   cluster_colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#843C39"),
#'   out_path = "./output_data/inferCNV/"
#' )
#' }
RunInferCNVCluster <- function(infercnv_obj, gene_order_file= system.file("extdata", "hg38_gencode_v27.txt", package = "easySingleCell"),
                               ref_cell_name, obs_cell_name, ref_group, obs_group,
                                        k_clusters = 6, heatmap_colors = c("#2166ac", "white", "#b2182b"),
                                        cluster_colors =c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#843C39"),
                                        out_path) {

  # Check if k_clusters equals length of cluster_colors
  if (k_clusters != length(cluster_colors)) {
    stop("The number of clusters (k_clusters) must be equal to the length of cluster_colors.")
  }



   expr <- infercnv_obj@expr.data

  # Load gene position information
  gene_pos <- read.delim(gene_order_file, header = FALSE)
  gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr), ]
  new_cluster <- unique(gene_pos$V2)

  top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                      gp = gpar(col = "white"),
                                                      labels_gp = gpar(cex = 1, col = "black"),
                                                      height = unit(5, "mm")))

  # Extract reference and observed cell positions
  ref_cell <- c(infercnv_obj@reference_grouped_cell_indices[ref_cell_name][[1]])
  obs_cell <- c(infercnv_obj@observation_grouped_cell_indices[obs_cell_name][[1]])
  cell_anno <- data.frame(cell_id = c(colnames(expr)[ref_cell], colnames(expr)[obs_cell]),
                          group = c(rep(ref_group, length(ref_cell)), rep(obs_group, length(obs_cell))))

  # Perform k-means clustering
  set.seed(123)
  kmeans_result <- kmeans(t(expr), k_clusters)
  kmeans_df <- data.frame(kmeans_result$cluster)
  colnames(kmeans_df) <- "k_cluster"
  kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
  kmeans_df <- kmeans_df %>% inner_join(cell_anno, by = "cell_id") %>% arrange(k_cluster)
  kmeans_df$k_cluster <- as.factor(kmeans_df$k_cluster)

  # Create row annotations
  annotation_row <- data.frame(k_cluster = kmeans_df$k_cluster, group = kmeans_df$group)
  row.names(annotation_row) <- kmeans_df$cell_id
  saveRDS(annotation_row, file = paste0(out_path, '/Kmeans.rds'))

  names(cluster_colors) <- as.character(1:k_clusters)

  left_anno <- rowAnnotation(df = annotation_row,
                             col = list(group = eval(parse(text = paste0("c(", ref_group, "='#00A0877F',", obs_group, "='#E64B357F')"))),
                                        k_cluster = cluster_colors),
                             show_annotation_name = FALSE)

  # Plot heatmap
  pdf(paste0(out_path, '/CNV_heatmap.pdf'), width = 9.5, height = 4)
  ht <- Heatmap(t(log2(expr))[rownames(annotation_row), ],
                col = colorRamp2::colorRamp2(c(-0.5, 0, 0.5), heatmap_colors),
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = FALSE, show_row_names = FALSE,
                column_split = factor(gene_pos$V2, new_cluster),
                heatmap_legend_param = list(title = "inferCNV",
                                            direction = "vertical",
                                            title_position = "leftcenter-rot",
                                            legend_height = unit(3, "cm")),
                left_annotation = left_anno,
                row_title = NULL,
                column_title = NULL,
                top_annotation = top_color,
                border = TRUE)
  draw(ht, heatmap_legend_side = "right")
  dev.off()

  return(annotation_row)
}



#' Run Full InferCNV Pipeline
#'
#' This function runs the full pipeline for InferCNV analysis, including data preparation,
#' InferCNV analysis, and clustering/visualization of the results.
#'
#' @param sce_epi Single-cell expression data for the tumor epithelial cells.
#' @param sce_refer Single-cell expression data for the reference cells (normal epithelial or immune cells).
#' @param celltype Column name in meta.data that contains cell type information.
#' @param infercnv_path Path to save InferCNV results.
#' @param name Name for the InferCNV run.
#' @param gene_order_file Path to the gene order file (default is "hg38_gencode_v27.txt" from "easySingleCell" package).
#' @param ref_group_names Reference group names for normal cells.
#' @param ref_cell_name Reference cell names for clustering and visualization.
#' @param obs_cell_name Observed cell names for clustering and visualization.
#' @param ref_group Reference group label for clustering and visualization.
#' @param obs_group Observed group label for clustering and visualization.
#' @param k_clusters Number of clusters for k-means clustering (default is 5).
#' @param heatmap_colors Colors for the heatmap (default is c("#2166ac", "white", "#b2182b")).
#' @param cluster_colors Colors for the clusters (default is c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#843C39")).
#' @param num_threads Number of threads to use for the analysis.
#' @return A list containing the InferCNV object and the clustering results.
#' @export
#' @import infercnv
#' @import qs
#' @import dplyr
#' @import Seurat
#' @import ComplexHeatmap
#' @import grid
#' @import scales
#' @import ggplot2
#' @import ggpubr
#' @import ggsignif
#' @import tidyverse
#' @import ggprism
#' @import vioplot
#' @import RColorBrewer
#' @examples
#' \dontrun{
#' # Ensure necessary libraries are loaded
#' library(infercnv)
#' library(qs)
#' library(dplyr)
#' library(Seurat)
#' library(ComplexHeatmap)
#' library(grid)
#' library(scales)
#' library(ggplot2)
#' library(ggpubr)
#' library(ggsignif)
#' library(tidyverse)
#' library(ggprism)
#' library(vioplot)
#' library(RColorBrewer)
#'
#' # Example usage:
#' # Run the full InferCNV pipeline
#' results <- runInferCNVPipeline(
#'   sce_epi = sce_epi,
#'   sce_refer = sce_refer,
#'   celltype = 'celltype',
#'   infercnv_path = './output_data/inferCNV/',
#'   name = 'infer_run',
#'   ref_group_names = c("T/NK", "B Cell"),
#'   ref_cell_name = c("T/NK","B Cell"),
#'   obs_cell_name = "epithlial",
#'   ref_group = "immune",
#'   obs_group = "epithelial",
#'   k_clusters = 5,
#'   heatmap_colors = c("#2166ac", "white", "#b2182b"),
#'   cluster_colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#843C39"),
#'   num_threads = 16
#' )
#' }
runInferCNVPipeline <- function(sce_epi, sce_refer, celltype = 'celltype',
                                    infercnv_path = './output_data', name = 'infer_run',
                                    gene_order_file = system.file("extdata", "hg38_gencode_v27.txt", package = "easySingleCell"),
                                    ref_group_names,
                                    ref_cell_name, obs_cell_name, ref_group, obs_group,
                                    k_clusters = 5, heatmap_colors = c("#2166ac", "white", "#b2182b"),
                                    cluster_colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#843C39"),
                                    num_threads = 12) {

  # Check if k_clusters equals length of cluster_colors
  if (k_clusters != length(cluster_colors)) {
    stop("The number of clusters (k_clusters) must be equal to the length of cluster_colors.")
  }

  # Step 1: Prepare Data for InferCNV
  prepared_data <- PrepareDataForInferCNV(
    sce_epi = sce_epi,
    sce_refer = sce_refer,
    celltype = celltype,
    infercnv_path = infercnv_path,
    name = name
  )

  # Step 2: Run InferCNV Analysis
  infercnv_obj <- RunInferCNVAnalysis(
    infercnv_path = prepared_data$infercnv_path,
    gene_order_file = gene_order_file,
    ref_group_names = ref_group_names,
    name = prepared_data$name,
    num_threads = num_threads
  )

  # Step 3: Cluster and Visualize InferCNV Results
  clustering_results <- RunInferCNVCluster(
    infercnv_obj = infercnv_obj,
    gene_order_file = gene_order_file,
    ref_cell_name = ref_cell_name,
    obs_cell_name = obs_cell_name,
    ref_group = ref_group,
    obs_group = obs_group,
    k_clusters = k_clusters,
    heatmap_colors = heatmap_colors,
    cluster_colors = cluster_colors,
    out_path = paste0(infercnv_path, '/', name)
  )

  return(list(infercnv_obj = infercnv_obj, clustering_results = clustering_results))
}
