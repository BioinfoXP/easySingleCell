# ========= Monocle =========
# =============== 1.Monocle2 ================

#' Process scRNA Data with Monocle2
#'
#' @description Performs dimensionality reduction and cell ordering using Monocle2.
#' @param scRNA Seurat Object.
#' @param save_path Rdata save path.
#' @param cores CPU cores.
#' @return A list containing Monocle CDS and reversed CDS.
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Seurat GetAssayData
#' @export
runMonocleAnalysis <- function(scRNA, save_path = "./output_data/monocle2.Rdata",
                               cores = 4, mean_expression = 0.1, disp_ratio = 1) {

  if (!base::dir.exists(base::dirname(save_path))) base::dir.create(base::dirname(save_path), recursive = TRUE)

  counts <- Seurat::GetAssayData(scRNA, assay = "RNA", slot = "counts")
  pd <- Biobase::AnnotatedDataFrame(data = scRNA@meta.data)
  fd <- Biobase::AnnotatedDataFrame(data = base::data.frame(gene_short_name = base::rownames(counts),
                                                            row.names = base::rownames(counts)))

  mycds <- monocle::newCellDataSet(base::as.matrix(counts), phenoData = pd, featureData = fd,
                                   expressionFamily = monocle::negbinomial.size())

  mycds <- monocle::estimateSizeFactors(mycds)
  mycds <- monocle::estimateDispersions(mycds, cores = cores, relative_expr = TRUE)

  disp_table <- monocle::dispersionTable(mycds)
  disp_genes <- base::subset(disp_table, mean_expression >= mean_expression &
                               dispersion_empirical >= disp_ratio * dispersion_fit)$gene_id

  mycds <- monocle::setOrderingFilter(mycds, disp_genes)
  mycds <- monocle::reduceDimension(mycds, max_components = 2, method = 'DDRTree', verbose = FALSE)
  mycds <- monocle::orderCells(mycds)
  mycds_reverse <- monocle::orderCells(mycds, reverse = TRUE)

  base::save(mycds, mycds_reverse, file = save_path)
  base::return(base::list(mycds = mycds, mycds_reverse = mycds_reverse))
}

# ========= CytoTRACE =========
# ========= 2. CytoTRACE Trajectory Analysis (Optimized) =========

#' Run CytoTRACE Analysis
#'
#' @description Performs CytoTRACE analysis to predict differentiation states
#' using a silent, namespace-safe implementation.
#'
#' @param scRNA A Seurat object.
#' @param celltype Metadata column for cell type annotation. Default "celltype".
#' @param output_data_dir Path for saving .Rdata results. Default "./output_data/".
#' @param output_figure_dir Path for saving plots. Default "./output_figure/".
#' @param emb Dimensionality reduction key (e.g., "umap") for visualization.
#' @param ncores Number of CPU cores. Default 5.
#' @param force_run Logical. If TRUE, ignores existing results and reruns analysis.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom CytoTRACE CytoTRACE plotCytoTRACE
#' @export
runCytoTRACEAnalysis <- function(scRNA,
                                 celltype = 'celltype',
                                 output_data_dir = "./output_data/",
                                 output_figure_dir = "./output_figure/",
                                 emb = 'umap',
                                 ncores = 5,
                                 force_run = FALSE) {

  # 1. Environment Preparation
  if (!base::dir.exists(output_data_dir)) base::dir.create(output_data_dir, recursive = TRUE)
  if (!base::dir.exists(output_figure_dir)) base::dir.create(output_figure_dir, recursive = TRUE)

  # 2. Input Validation
  if (!base::inherits(scRNA, "Seurat")) base::stop("Input 'scRNA' must be a Seurat object.")
  if (!celltype %in% base::colnames(scRNA@meta.data)) base::stop(base::paste0("Column '", celltype, "' not found in metadata."))

  if (!emb %in% base::names(scRNA@reductions)) {
    base::stop(base::paste0("Reduction '", emb, "' not found in Seurat object."))
  }

  # 3. Prepare Phenotype Data
  # Extract phenotype vector matching cell names
  phe <- base::as.character(scRNA@meta.data[[celltype]])
  base::names(phe) <- base::rownames(scRNA@meta.data)

  # 4. Analysis Execution (with Caching)
  results_file <- base::file.path(output_data_dir, "Cytotrace_results.Rdata")

  if (base::file.exists(results_file) && !force_run) {
    base::message("Loading existing CytoTRACE results...")
    # Variable 'results' will be loaded into environment
    base::load(results_file)

    # Check if loaded results match current cells (basic validation)
    if (!base::exists("results")) base::stop("Loaded file does not contain 'results' object.")

  } else {
    base::message("Running CytoTRACE analysis (this may take time)...")

    # Secure data extraction using Seurat accessor
    # CytoTRACE requires a raw count matrix (genes x cells)
    mat_counts <- base::as.matrix(Seurat::GetAssayData(scRNA, assay = "RNA", slot = "counts"))

    # Run core analysis
    results <- CytoTRACE::CytoTRACE(mat = mat_counts, ncores = ncores)

    base::save(results, file = results_file)
    base::message("Analysis saved to ", results_file)
  }

  # 5. Visualization
  # plotCytoTRACE automatically generates PDFs in the output directory
  base::message("Generating plots in ", output_figure_dir)

  CytoTRACE::plotCytoTRACE(
    results,
    outputDir = output_figure_dir,
    emb = scRNA@reductions[[emb]]@cell.embeddings,
    phenotype = phe
  )
}
