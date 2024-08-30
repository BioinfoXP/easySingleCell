# =========== 1. PrepareCpdb =============

#' @title Preprocess and Save Single-cell Data for CellPhoneDB Analysis
#' @description This function preprocesses single-cell data, normalizes it, and saves the results for CellPhoneDB analysis.
#' @param scRNA A Seurat object containing the single-cell RNA data.
#' @param output_dir A character string specifying the output directory where results will be saved. Default is './output_data/'.
#' @param celltype_column A character string specifying the column name in the metadata that contains cell type information. Default is 'celltype'.
#' @param normalization_threshold A numeric value specifying the threshold for filtering genes based on mean expression. Default is 0.001.
#' @param prefix A character string specifying the prefix to be added to the output file names. Default is "".
#' @return None. The function saves the normalized counts and metadata files in the specified output directory.
#' @export
#' @import Seurat
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' PrepareCpdb(
#'   scRNA = sce,
#'   output_dir = './output_data/',
#'   celltype_column = 'celltype',
#'   normalization_threshold = 0.001,
#'   prefix = ''
#' )
#' }

PrepareCpdb <- function(scRNA, output_dir = './output_data/', celltype_column = 'celltype', normalization_threshold = 0.001, prefix = "") {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(magrittr)

  # Create the output directory for CellPhoneDB if it doesn't exist
  cpdb_dir <- file.path(output_dir, 'cpdb')
  dir.create(cpdb_dir, showWarnings = FALSE, recursive = TRUE)

  # Normalize the data and filter genes with low expression
  Normalized_counts <- GetAssayData(scRNA, slot = 'data') %>% as.data.frame()
  Normalized_counts$mean_exp <- rowMeans(Normalized_counts)
  Normalized_counts <- Normalized_counts[Normalized_counts$mean_exp > normalization_threshold, ]
  Normalized_counts <- Normalized_counts[, -ncol(Normalized_counts)]
  Normalized_counts <- cbind(Gene = rownames(Normalized_counts), Normalized_counts)

  # Create metadata file
  metadata <- data.frame(Cell = rownames(scRNA@meta.data), celltype = scRNA@meta.data[, celltype_column])
  metadata$celltype <- gsub(' ', '_', metadata$celltype)

  # Save the normalized counts and metadata files
  write.table(Normalized_counts, file.path(cpdb_dir, paste0(prefix, 'Normalized_counts.txt')), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(metadata, file.path(cpdb_dir, paste0(prefix, 'cellphonedb_meta.txt')), row.names = FALSE, sep = '\t', quote = FALSE)
}

# Example usage
# PrepareCpdb(
#   scRNA = sce,
#'   output_dir = './output_data/',
#'   celltype_column = 'celltype',
#'   normalization_threshold = 0.001,
#'   prefix = ''
#' )


# ============ 2. PrepareCell2loc ==========

#' @title Prepare Data for Cell2loc
#' @description This function prepares data for Cell2loc analysis by exporting counts matrix, UMAP embeddings, and metadata from a Seurat object.
#' @param scRNA A Seurat object.
#' @param output_dir A character string specifying the directory to save the output files. Default is "output_data".
#' @param prefix A character string specifying the prefix to be added to the output file names. Default is "".
#' @return NULL
#' @export
#' @import Seurat
#' @import qs
#' @examples
#' \dontrun{
#' PrepareCell2loc(input_file = './output_data/Figure11/sce.tumor.所有细胞注释.qs', output_dir = 'output_data', prefix = 'Figure12_')
#' }

PrepareCell2loc <- function(scRNA, output_dir = "output_data", prefix = "") {
  # Load required libraries
  library(Seurat)
  library(qs)

  # Ensure the output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Export counts matrix
  counts_matrix <- t(as.matrix(scRNA@assays$RNA@counts))
  write.table(counts_matrix, file = file.path(output_dir, paste0(prefix, "scRNA_counts.csv")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Export UMAP embeddings
  embedding <- Embeddings(scRNA, "umap")
  write.table(embedding, file = file.path(output_dir, paste0(prefix, "scRNA_embedding.csv")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Export metadata
  metadata <- scRNA@meta.data
  write.table(metadata, file = file.path(output_dir, paste0(prefix, "metadata.csv")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  return(NULL)
}

# Example usage
# PrepareCell2loc(input_file = './output_data/Figure11/sce.tumor.所有细胞注释.qs', output_dir = 'output_data', prefix = 'Figure12_')



# =========== 3. PreparePyscenic =============
# makeMetaCells.R

#' Create Metacells
#'
#' This function creates metacells from a Seurat object.
#'
#' @param seu A Seurat object containing the single-cell RNA-seq data.
#' @param min.cells Minimum number of cells required to form a metacell. Default is 10.
#' @param reduction The reduction method to use. Default is "umap".
#' @param dims Dimensions to use for the reduction. Default is 1:2.
#' @param k.param Number of nearest neighbors to use in the clustering. Default is 10.
#' @param cores Number of cores to use for parallel computation. Default is 20.
#'
#' @return A list containing the metacell matrix and metadata.
#' @import Seurat
#' @import parallel
#' @import Matrix
#' @import dplyr
makeMetaCells <- function(seu, min.cells = 10, reduction = "umap", dims = 1:2, k.param = 10, cores = 20) {
  seu <- seu %>%
    FindNeighbors(reduction = reduction, dims = dims, k.param = k.param) %>%
    FindClusters(res = 50)
  metadata <- seu@meta.data
  metadata$METACELL_ID <- factor(metadata$seurat_clusters)
  dge_mat <- seu[["RNA"]]@counts

  dge_mat_mc <- parallel::mclapply(levels(metadata$METACELL_ID), function(xx) {
    cells <- rownames(subset(metadata, METACELL_ID == xx))
    Matrix::rowSums(dge_mat[, cells])
  }, mc.cores = cores)
  dge_mat_mc <- do.call(cbind, dge_mat_mc)

  metacell_metadata <- metadata[["METACELL_ID"]] %>% table() %>% as.data.frame()
  colnames(metacell_metadata) <- c("METACELL_ID", "CELL_COUNT")
  rownames(metacell_metadata) <- metacell_metadata[["METACELL_ID"]]

  kept.cells <- subset(metacell_metadata, CELL_COUNT >= min.cells)[["METACELL_ID"]]
  metacells <- list(
    mat = dge_mat_mc[, kept.cells],
    metadata = metacell_metadata[kept.cells, ]
  )
  colnames(metacells$mat) <- paste0(seu@project.name, ".METACELL_", kept.cells)
  rownames(metacells$metadata) <- colnames(metacells$mat)
  return(metacells)
}



#' Prepare Data for pySCENIC Analysis
#'
#' This function prepares the data for pySCENIC analysis by splitting the Seurat object by cell type,
#' creating meta cells, and generating the necessary files for pySCENIC.
#'
#' @param scRNA A Seurat object containing single-cell RNA data.
#' @param celltype A character string to split the Seurat object by cell type. Default is "celltype".
#' @param output_dir A character string specifying the output directory where results will be saved. Default is './output_data/'.
#' @param nCells An integer specifying the number of cells to downsample to if not using MetaCell. Default is 500.
#' @param use_MetaCell A logical value indicating whether to use MetaCell. If FALSE, the function will downsample the Seurat object using raw counts. Default is FALSE.
#' @param min_cells An integer specifying the minimum number of cells for meta cell creation. Only used if use_MetaCell is TRUE. Default is 10.
#' @param reduction A character string specifying the reduction method to use. Only used if use_MetaCell is TRUE. Default is 'umap'.
#' @param dims A numeric vector specifying the dimensions to use for the reduction. Only used if use_MetaCell is TRUE. Default is 1:2.
#' @param k_param An integer specifying the k parameter for meta cell creation. Only used if use_MetaCell is TRUE. Default is 10.
#' @param cores An integer specifying the number of cores to use. Only used if use_MetaCell is TRUE. Default is 10.
#' @return None. The function saves the meta cell matrix and other necessary files for pySCENIC analysis in the specified output directory.
#' @export
#' @import Seurat
#' @importFrom magrittr %>%
#' @import data.table
#' @import arrow
#' @import SCopeLoomR
#' @examples
#' \dontrun{
#' PreparePyscenic(
#'   scRNA = scRNA,
#'   output_dir = './output_data/',
#'   nCells = 500,
#'   use_MetaCell = FALSE,
#'   min_cells = 10,
#'   reduction = 'umap',
#'   dims = 1:2,
#'   k_param = 10,
#'   cores = 10
#' )
#' }
PreparePyscenic <- function(scRNA, celltype = 'celltype', output_dir = './output_data/',
                            nCells = 500, use_MetaCell = FALSE, min_cells = 10, reduction = 'umap',
                            dims = 1:2, k_param = 10, cores = 10) {
  # Load necessary libraries
  library(Seurat)
  library(magrittr)
  library(data.table)
  library(arrow)
  library(SCopeLoomR)

  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!use_MetaCell) {
    # Downsample Seurat object
    sce.sub <- subset(scRNA, downsample = nCells)
    mc.mat <- GetAssayData(sce.sub, slot = 'counts') %>%
      as.matrix()

    # Create a loom file for pySCENIC
    loom <- SCopeLoomR::build_loom(
      file.name = file.path(output_dir, "00-2.mc_mat_for_step1.loom"),
      dgem = mc.mat,
      default.embedding = NULL
    )
    loom$close()
  } else {
    # Load pre-prepared data
    load(system.file("data", "pyscenic_database.Rdata", package = "easySingleCell"))

    # Split the Seurat object by cell type
    seu.list <- SplitObject(scRNA, split.by = celltype)
    for (i in seq_along(seu.list)) {
      seu.list[[i]]@project.name <- names(seu.list)[i]
    }

    # Create meta cells for each cell type
    metacells.list <- lapply(seq_along(seu.list), function(ii) {
      makeMetaCells(
        seu = seu.list[[ii]],
        min.cells = min_cells,
        reduction = reduction,
        dims = dims,
        k.param = k_param,
        cores = cores
      )
    })

    # Combine meta cell matrices and metadata
    mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
    mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)

    # Save the meta cell matrix
    saveRDS(mc.mat, file.path(output_dir, "00-1.mc.mat.rds"))

    # Read the saved meta cell matrix
    mc.mat <- readRDS(file.path(output_dir, "00-1.mc.mat.rds"))

    # Filter low expression genes
    expr.in.cells <- rowSums(mc.mat > 0)
    mc.mat <- mc.mat[expr.in.cells >= 5, ]

    # Intersect genes with cisDB
    genes.use <- intersect(cisdb.genes, rownames(mc.mat))
    mc.mat <- mc.mat[genes.use, ]

    # Create a loom file for pySCENIC
    loom <- SCopeLoomR::build_loom(
      file.name = file.path(output_dir, "00-2.mc_mat_for_step1.loom"),
      dgem = mc.mat,
      default.embedding = NULL
    )
    loom$close()
  }
}


# =========== 4. ProcessPyscenic =============
#' Calculate gene module scores
#'
#' @param x gene expression matrix, rows are genes, columns are cells. Can be any format, UMI, CPM, TPM, etc.
#' @param ... Arguments passed to other methods.
#' @return A signature score matrix or Seurat object.
#'
#' @export
ComputeModuleScore <- function(x, ...) UseMethod('ComputeModuleScore')


#' @param gene.sets a list of gene sets in data.frame or named list.
#' @param min.size The minimal genes of the gene sets. The size of gene sets less than this value were ignored.  Default: 20
#' @param batch.size The number of cells were calculated for each batch. Default: 500
#' @param cores number of threads for parallel computing. Default: 1
#'
#' @rdname ComputeModuleScore
#' @export
ComputeModuleScore.default <- function(x, gene.sets, min.size=20, batch.size=500, cores=1, ...) {
  if (!is.list(gene.sets)) {
    stop("'gene.sets' should be a list or data.frame!")
  }
  gene.sets <- gene.sets[sapply(gene.sets, length) >= min.size]
  n.cells <- ncol(x)
  batches <- floor((1:n.cells-1) / batch.size)
  batch.levels <- unique(batches)
  aucell <- function(i) {
    dge.tmp <- x[, batches == i]
    cr <- AUCell::AUCell_buildRankings(dge.tmp, nCores=1, plotStats=F, verbose = F)
    auc <- AUCell::AUCell_calcAUC(gene.sets, cr, nCores=1, verbose = F)
    AUCell::getAUC(auc)
  }
  auc_scores <- parallel::mclapply(batch.levels, aucell, mc.cores = cores)
  do.call(cbind, auc_scores)
}


#' @rdname ComputeModuleScore
#' @param assay Name of the seurat object assay.
#' @export
ComputeModuleScore.Seurat <- function(x, gene.sets, min.size=20, batch.size=500, cores=1, assay = Seurat::DefaultAssay(x), ...) {
  dge <- x[[assay]]@counts
  ras_mat <- ComputeModuleScore.default(x = dge, gene.sets, min.size, batch.size, cores)
  x[["AUCell"]] <- Seurat::CreateAssayObject(data = ras_mat)
  return(x)
}


library(magrittr)

#' Find Differential Expressed Genes via Variance Deomposition using Mixed Linear Model
#' @param data A data.frame or matrix of gene expression matrix
#' @param meta.data A data.frame contains cell meta data
#' @param vd.vars variable in meta.data for variance decomposition.
#' @param genes Genes for variance decomposition. Default: "all"
#' @param cores The number of threads for calculation, -1 means all available
#' threads. Default: -1
#'
#' @return A data.frame contains the results of variance decomposition.
#' @export
#'
VarDecompose <- function(data, meta.data, vd.vars, genes = "all", cores = -1) {
  ## check params
  if (missing(data) || missing(meta.data) || missing(vd.vars)) {
    stop("Must provide 'data', 'meta.data', and 'vd.vars'.")
  }
  if (is.null(colnames(meta.data)) || is.null(rownames(meta.data)) ) {
    stop("The row and column names of 'meta.data' should be provided.")
  }
  if (is.null(colnames(data)) || is.null(rownames(data)) ) {
    stop("The row and column names of 'data' should be provided.")
  }
  if (!all(rownames(data) == rownames(meta.data)) ) {
    stop("The row names of 'data' and 'meta.data' should be matched.")
  }
  if (!all(vd.vars %in% colnames(meta.data))) {
    vd.vars.404 <- setdiff(vd.vars, colnames(meta.data))
    stop(paste("vd.vars:", vd.vars.404, "is(are) not found in 'meta.data'"))
  }
  if (length(genes) == 1) {
    if (genes == "all") {
      genes <- colnames(data)
    } else {
      stop("'genes' should be 'all' or a vector.")
    }
  } else {
    genes.404 <- setdiff(genes, colnames(data))
    if (length(genes.404) > 0) {
      warning(paste(length(genes.404), "gene(s) are not found in 'data'."))
      genes <- setdiff(genes, genes.404)
    }
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## prepare data for VD
  vd.vars.str <- sapply(vd.vars, function(xx) sprintf("(1|%s)", xx))
  modelFormulaStr <- paste("expression ~", paste(vd.vars.str, collapse = "+"))
  data.use <- cbind(data[, genes], meta.data)
  ## exe VD
  vd.res <- do.call(rbind, parallel::mclapply(genes, function(genename) {
    data.model <- data.use[, c(vd.vars, genename)]
    colnames(data.model) <- c(vd.vars, "expression")
    tryCatch({
      model <- suppressWarnings(lme4::lmer(stats::as.formula(modelFormulaStr), data = data.model, REML = TRUE, verbose = FALSE))
      results <- as.data.frame(lme4::VarCorr(model))
      rownames(results) <- results$grp
      results <- results[c(vd.vars, "Residual"), ]
      frac.var <- results$vcov / sum(results$vcov)

      res.tmp <- c("OK", frac.var)
    },
    error = function(e) {
      print(e)
      res.tmp <- c("FAIL", rep(-1, length(vd.vars)+1))
    })
    names(res.tmp) <- c("status", vd.vars, "residual")
    as.data.frame(as.list(res.tmp)) # return
  }, mc.cores = cores)
  )
  rownames(vd.res) <- genes
  vd.res %<>% as.data.frame()
  vd.res <- vd.res %>% dplyr::mutate(gene = rownames(vd.res), .before=1)
  # vd.res <- vd.res %>% as.data.frame() %>% dplyr::mutate(gene = rownames(.), .before=1)
  for (i in 3:ncol(vd.res)) {
    vd.res[[i]] %<>% as.numeric()
  }
  return(vd.res)
}




# process_scenic.R

#' Process SCENIC Output and Compute Module Scores
#'
#' This function processes the output from SCENIC, computes module scores, and performs dimensionality reduction.
#'
#' @param sce Seurat object containing the single-cell data.
#' @param inputdir Path to the input directory containing the SCENIC output files. Ignored if `regulons` is provided.
#' @param outputdir Path to the output directory where results will be saved. If not specified, defaults to the input directory.
#' @param gmtfile Name of the GMT file containing regulons. Ignored if `regulons` is provided.
#' @param regulons A list of regulons. If provided, `inputdir` and `gmtfile` are ignored.
#' @param min.size Minimum size of the regulons to be considered. Default is 10.
#' @param cores Number of cores to use for computation. Default is 15.
#' @param vd.vars Variable for variance decomposition. Default is 'celltype'.
#' @param prefix Prefix to add to the output file names. Default is "03-1".
#'
#' @return None. The function outputs files to the specified directory.
#' @import Seurat
#' @import AUCell
#' @import dplyr
#' @import future
#' @import furrr
#' @import readr
#' @import tidyr
#' @import clusterProfiler
#' @import qs
#' @export
ProcessPyscenic <- function(sce,
                            inputdir = NULL,
                            outputdir = inputdir,
                            gmtfile = NULL,
                            regulons = NULL,
                            min.size = 10,
                            cores = 15,
                            vd.vars = 'celltype',
                            prefix = "03-1") {
  # Check if regulons are provided
  if (is.null(regulons)) {
    if (is.null(inputdir) || is.null(gmtfile)) {
      stop("Either 'regulons' must be provided, or both 'inputdir' and 'gmtfile' must be specified.")
    }

    # Check if the input file exists
    gmt_path <- file.path(inputdir, gmtfile)
    if (!file.exists(gmt_path)) {
      stop("GMT file not found: ", gmt_path)
    }

    # Read GMT file
    regulons <- read.gmt(gmt_path)
  }

  # Create regulon list
  regulon.list <- split(regulons$gene, sub("[0-9]+g", "\\+", regulons$term))

  # Save regulon list
  regulon_rds_path <- file.path(outputdir, paste0(prefix, ".regulons.rds"))
  saveRDS(regulon.list, regulon_rds_path)

  # Set the default assay of the Seurat object
  DefaultAssay(sce) <- 'RNA'

  # Compute module scores
  sce <- ComputeModuleScore(sce, regulon.list, min.size = min.size, cores = cores)

  # Set the default assay of the Seurat object to AUCell
  DefaultAssay(sce) <- "AUCell"

  # Run UMAP
  sce <- RunUMAP(sce, features = rownames(sce), metric = "correlation", reduction.name = "umapRAS", reduction.key = "umapRAS_")

  # Compute PCA and visualize
  sce <- ScaleData(sce)
  sce <- RunPCA(sce, features = rownames(sce), reduction.name = "pcaRAS", reduction.key = "pcaRAS_")

  # Variance decomposition
  meta.data <- sce@meta.data[, vd.vars, drop = FALSE]
  meta.data[, vd.vars] <- factor(meta.data[, vd.vars])
  ras.data <- FetchData(sce, vars = rownames(sce))
  vd.res <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = cores)

  # Save variance decomposition results
  vd_res_path <- file.path(outputdir, paste0("04-1.VD_res.rds"))
  saveRDS(vd.res, vd_res_path)

  # Save Seurat object
  sce_rds_path <- file.path(outputdir, paste0(prefix, ".sce.qs"))
  qsave(sce, sce_rds_path)
}



# =========== 5. ProcessScFEA =============
#' Prepare data for scFEA
#'
#' This function prepares single-cell RNA sequencing data for scFEA analysis.
#' It exports the raw counts data to a specified output directory.
#'
#' @param scRNA A Seurat object containing single-cell RNA sequencing data.
#' @param output_dir A character string specifying the directory where the output files will be saved. Default is './output_data/'.
#' @param prefix A character string to be prefixed to the output file names. Default is an empty string.
#'
#' @return None. The function writes the counts data to a file in the specified output directory.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Assuming 'scRNA' is your Seurat object
#' PrepareScFEA(scRNA, output_dir = './output_data/', prefix = "")
#'
#' # Example to run full pipeline:
#' # ===> 1. Prepare environment
#' # git clone https://github.com/changwn/scFEA
#' # cd scFEA
#' # conda create -n scFEA -y \
#' # && conda activate scFEA \
#' # && conda install python=3.8 -y \
#' # && conda install --file requirements -y \
#' # && conda install pytorch torchvision -c pytorch -y \
#' # && pip install --user magic-impute \
#' # && conda install pandas=1.5.3 -y
#' # pip install ipykernel
#' # python -m ipykernel install --user --name scFEA --display-name "scFEA"
#'
#'
#'# pip install --user magic-impute
#' # Use "python src/scFEA.py --help" to check install.
#'
#'
#' # ===> 2. Prepare file
#' # PrepareScFEA(scRNA, output_dir = './output_data/', prefix = "")
#'
#' # ===> 3. Run scFEA analysis
#' # cd scFEA
#' # ===> Enter the scFEA folder and move your .csv file into "input" folder.
#' # python src/scFEA.py --data_dir data --input_dir input \
#' # --test_file scFEA_counts.txt \
#' # --moduleGene_file module_gene_m168.csv \
#' # --stoichiometry_matrix cmMat_c70_m168.csv \
#' # --output_flux_file output/flux.csv \
#' # --output_balance_file output/flux_balance.csv \
#' # --sc_imputation True
#' }
#'
#' @export
PrepareScFEA <- function(scRNA, output_dir = './output_data/', prefix = "") {
  # Ensure necessary libraries are loaded
  library(Seurat)

  # Create the output directory
  scFEA_dir <- file.path(output_dir, 'scFEA')
  dir.create(scFEA_dir, showWarnings = FALSE, recursive = TRUE)

  # Export the counts
  counts <- GetAssayData(scRNA, slot = 'counts')
  write.csv(counts, file = file.path(scFEA_dir, paste0(prefix, 'scFEA_counts.txt')), row.names = TRUE)
}
