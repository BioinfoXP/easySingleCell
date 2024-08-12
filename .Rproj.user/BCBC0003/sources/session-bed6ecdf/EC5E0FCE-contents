# Load necessary libraries
library(tidyverse)
library(UCSCXenaTools)
library(IOBR)

# Define cancer types
cancer_types <- c('PANCAN', 'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL',
                  'COAD', 'DLBC', 'ESCA', 'GBM', 'KIRP', 'LAML',
                  'HNSC', 'KICH', 'KIRC', 'LGG', 'LIHC',
                  'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD',
                  'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                  'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UVM')

# ======================== 1. Index object============================
setClass('Index',
         slots = list(save_path = "character",
                      data_type = "character",
                      cancer_type = "character"),
         prototype = list(save_path = 'output_xena',
                          data_type = 'counts',
                          cancer_type = 'NULL'))

# Validity check for Index class
setValidity("Index", function(object) {
  if (object@save_path == 'NULL') {
    warning("No path provided, it will create an 'output_xena' folder!", call. = FALSE)
  } else if (!dir.exists(object@save_path)) {
    cli::cli_abort(message = {
      cli::col_red("Error: Invalid path!")
    })
  } else {
    message('Your save path: ', object@save_path)
  }

  if (object@cancer_type  == 'NULL') {
    cli::cli_abort(message = {
      cli::col_red("Error: Please provide valid cancer type!")
    })
  } else if (!object@cancer_type %in% cancer_types) {
    cli::cli_abort(message = {
      cli::col_red("Error: Please provide valid cancer type!")
    })
  } else {
    message('We will download: ', object@cancer_type)
  }
})

# ======================== 2. PanCancer & TCGA object============================
setClass('PanCancer', contains = 'Index',
         slots = list(target_info = "character"),
         prototype = list(target_info = "tcga_target_no_normal_RSEM_hugo_norm_count",
                          cancer_type = 'PANCAN',
                          save_path = 'output_xena',
                          data_type = 'counts'))

setClass('TCGA', contains = 'Index',
         slots = list(target_info = "character"),
         prototype = list(target_info = 'NULL',
                          cancer_type = 'NULL',
                          save_path = 'output_xena',
                          data_type = 'counts'))

# ========================= 3. XenaTools_xp Interface============================
setGeneric("XenaTools_xp", function(object, ...) {
  standardGeneric("XenaTools_xp")
})

# ========================= 4. Methods============================
# ========================= PanCancer =========================
setMethod("XenaTools_xp", "PanCancer", function(object, ...) {
  {
    library(tidyverse)
    library(UCSCXenaTools)
    options(timeout = 10000)
    Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
    dir.create(object@save_path, showWarnings = FALSE)
  }

  {
    filterDataset1 <- 'TCGA_TARGET_phenotype'
    filterDataset2 <- 'TCGA_survival_data_2.txt'
  }

  # Retrieve data
  message("=> Starting query!\n")
  project_data <- XenaScan() %>%
    XenaGenerate()
  message("=> Starting download!\n")
  exp <- project_data %>%
    XenaFilter(filterDatasets = object@target_info) %>%
    XenaQuery() %>%
    XenaDownload() %>%
    XenaPrepare()

  # Download clinical data
  cli <- project_data %>%
    XenaFilter(filterDatasets = filterDataset1) %>%
    XenaQuery() %>%
    XenaDownload() %>%
    XenaPrepare()

  # Download survival data
  surv <- project_data %>%
    XenaFilter(filterDatasets = filterDataset2) %>%
    XenaQuery() %>%
    XenaDownload() %>%
    XenaPrepare()

  # Save data
  save(exp, cli, surv, file = paste0(object@save_path, '/', object@data_type, "_xena_", object@cancer_type, ".Rdata"))
  message("=> Success! Data saved.")

  # Clean
  message("=> Start! Data clean.")
  coding_gene <- system.file("data", "protein_coding.Rdata", package = "easySingleCell")
  load(coding_gene)

  coding_gene <- coding_gene %>%
    filter(gene_type == "protein_coding") %>%
    distinct(gene_name, .keep_all = TRUE)

  # Common genes
  com_gene <- intersect(exp$Ensembl_ID, coding_gene$gene_id)

  # Filter (keep only coding genes)
  coding_gene <- coding_gene[coding_gene$gene_id %in% com_gene, ]
  exp <- exp[exp$Ensembl_ID %in% com_gene, ]

  exp <- exp %>% tibble::column_to_rownames(var = "Ensembl_ID")
  exp <- exp[coding_gene$gene_id, ]
  identical(row.names(exp), coding_gene$gene_id)  # TRUE
  row.names(exp) <- NULL
  row.names(exp) <- coding_gene$gene_name

  # Extract coding genes from expression matrix
  com_gene <- intersect(colnames(exp), cli$submitter_id.samples)
  com_gene <- intersect(com_gene, surv$sample)

  exp <- exp[, com_gene]
  exp <- 2^exp - 1
  cli <- cli[match(com_gene, cli$submitter_id.samples), ]
  surv <- surv[match(com_gene, surv$sample), ]

  # Save cleaned data
  identical(colnames(exp), surv$sample)  # TRUE
  identical(colnames(exp), cli$submitter_id.samples)  # TRUE

  # Save data
  save(exp, cli, surv, file = paste0(object@save_path, '/', object@data_type, "_xena_clean_", object@cancer_type, ".Rdata"))
  message("=> Success! Data saved.")
})

# =====================TCGA=====================
setMethod("XenaTools_xp", "TCGA", function(object, ...) {
  {
    library(tidyverse)
    library(UCSCXenaTools)
    options(timeout = 10000)
    Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
    dir.create(object@save_path, showWarnings = FALSE)
  }

  # Retrieve data
  message("=> Starting query!\n")
  project <- object@cancer_type
  project_data <- XenaScan(pattern = project)
  # Only GDC
  cohort <- str_extract(project_data$XenaCohorts, pattern = "^GDC.*") %>%
    na.omit() %>% unique()

  # Download expression matrix
  dataset <- project_data %>%
    XenaGenerate(subset = XenaCohorts == cohort & DataSubtype == "gene expression RNAseq")

  message("=> Starting download!\n")
  exp <- dataset %>%
    XenaFilter(filterDatasets = dataset@datasets[str_detect(dataset@datasets, object@data_type)]) %>% # Filter datasets
    XenaQuery() %>% # Query
    XenaDownload() %>% # Download
    XenaPrepare() # Prepare

  # Download clinical data
  # Retrieve
  sur_pheno <- project_data %>% XenaGenerate(subset = XenaCohorts == cohort & DataSubtype == "phenotype")

  dataset2 <- sur_pheno@datasets # Index

  filterDataset1 <- dataset2[1]
  filterDataset2 <- dataset2[2]

  # Download survival and cli data
  data1 <- sur_pheno %>%
    XenaFilter(filterDatasets = filterDataset1) %>% # Filter datasets
    XenaQuery() %>% # Query
    XenaDownload() %>% # Download
    XenaPrepare() # Prepare

  data2 <- sur_pheno %>%
    XenaFilter(filterDatasets = filterDataset2) %>% # Filter datasets
    XenaQuery() %>% # Query
    XenaDownload() %>% # Download
    XenaPrepare() # Prepare

  # results judge
  # if cols > 4 is cli，else if surv
  cols <- c(ncol(data1), ncol(data2))
  result <- ifelse(cols > 4, "cli", "surv")
  assign(result[1], data1)
  assign(result[2], data2)

  # Save data
  save(exp, cli, surv, file = paste0(object@save_path, '/', object@data_type, "_xena_", object@cancer_type, ".Rdata"))
  message("=> Success! Data saved.")

  # Clean
  message("=> Start! Data clean.")
  coding_gene <- system.file("data", "protein_coding.Rdata", package = "OneGene")
  load(coding_gene)

  coding_gene <- coding_gene %>%
    filter(gene_type == "protein_coding") %>%
    distinct(gene_name, .keep_all = TRUE)

  # Common genes
  com_gene <- intersect(exp$Ensembl_ID, coding_gene$gene_id)

  # Filter (keep only coding genes)
  coding_gene <- coding_gene[coding_gene$gene_id %in% com_gene, ]
  exp <- exp[exp$Ensembl_ID %in% com_gene, ]

  exp <- exp %>% tibble::column_to_rownames(var = "Ensembl_ID")
  exp <- exp[coding_gene$gene_id, ]
  identical(row.names(exp), coding_gene$gene_id)  # TRUE
  row.names(exp) <- NULL
  row.names(exp) <- coding_gene$gene_name

  # Extract coding genes from expression matrix
  com_gene <- intersect(colnames(exp), cli$submitter_id.samples)
  com_gene <- intersect(com_gene, surv$sample)

  exp <- exp[, com_gene]
  exp <- 2^exp - 1
  cli <- cli[match(com_gene, cli$submitter_id.samples), ]
  surv <- surv[match(com_gene, surv$sample), ]

  # Save cleaned data
  identical(colnames(exp), surv$sample)  # TRUE
  identical(colnames(exp), cli$submitter_id.samples)  # TRUE

  # TPM conversion and log transformation
  library(IOBR)
  tpm_data <- count2tpm(exp, idType = 'Symbol')
  tpm_data <- log2(tpm_data + 1)

  # Save data
  save(tpm_data, exp, cli, surv, file = paste0(object@save_path, '/', object@data_type, "_xena_clean_", object@cancer_type, ".Rdata"))
  message("=> Success! Data saved.")
})

# =====================5. get_TCGA ==============
#' Create an Object for TCGA or PanCancer Project then download and clean the data
#'
#' This function creates an object for either the TCGA or PanCancer project.
#'
#' @param project_type The type of project ('TCGA' or 'PanCancer').
#' @param cancer_type The type of cancer to study.
#' @return An object of class 'TCGA' or 'PanCancer'.
#' @export
#' @examples
#'     get_TCGA(cancer_type = 'COAD', project_type = 'TCGA')
get_TCGA <- function(cancer_type = 'NULL',
                     project_type = 'TCGA') {
  save_path = 'output_xena'
  dir.create(save_path, showWarnings = FALSE)


  # Create the object
  obj <- if (!project_type %in% c('TCGA', 'PanCancer')) {
    cli::cli_abort(message = {
      cli::col_red("Error: Please provide valid project type!\n\t'TCGA' or 'PanCancer'!")
    })
  } else if (project_type == 'PanCancer') {
    new(Class = 'PanCancer', save_path = save_path)
  } else {
    new(Class = project_type, cancer_type = cancer_type, save_path = save_path)
  }

  # Run the XenaTools_xp function
  XenaTools_xp(obj)
}



#' @title TCGA Group Classification
#' @description This function classifies samples into 'tumor' and 'normal' based on TCGA sample barcodes.
#' @param expr A matrix or data frame containing expression data with genes as rows and samples as columns.
#' @return A factor vector indicating the group ('tumor' or 'normal') for each sample.
#' @export
#' @examples
#' \dontrun{
#' expr <- readRDS('./path_to_your_expression_data.rds')
#' group <- TCGA_group(expr)
#' }

# ======== 6. TCGA Group =============
TCGA_group <- function(expr) {
  group <- sapply(strsplit(colnames(expr), "\\-"), "[", 4)
  group <- sapply(strsplit(group, ""), "[", 1)
  group_list <- ifelse(group == "0", 'tumor', 'normal')
  group_list <- factor(group_list, levels = c('normal', 'tumor'))
  return(group_list)
}


# ======== 7. TCGA DEG limma=============
#' @title Differential Expression Analysis
#' @description This function performs differential expression analysis using limma.
#' @param expr.tpm A matrix or data frame containing TPM expression data with genes as rows and samples as columns.
#' @param group A factor vector indicating the group ('tumor' or 'normal') for each sample. If not provided, the function will generate it.
#' @param contrast1 A character string specifying the first group for comparison. Default is 'tumor'.
#' @param contrast2 A character string specifying the second group for comparison. Default is 'normal'.
#' @param output_dir A character string specifying the directory to save the results. Default is "./output_data".
#' @return A data frame containing the differential expression analysis results.
#' @export
#' @import dplyr
#' @import limma
#' @import IOBR
#' @examples
#' \dontrun{
#' expr <- readRDS('./path_to_your_expression_data.rds')
#' group <- TCGA_group(expr)
#' results <- TCGADegLimma(expr, group, output_dir = "./output_data")
#' }

TCGADegLimma <- function(expr.tpm, group = NULL, contrast1 = 'tumor', contrast2 = 'normal', output_dir = "./output_data") {
  # Ensure necessary libraries are loaded
  library(dplyr)
  library(limma)
  library(IOBR)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate group if not provided
  if (is.null(group)) {
    group <- TCGA_group(expr.tpm)
  }

  # Differential expression analysis
  group_info <- data.frame(id = colnames(expr.tpm), group)
  row.names(group_info) <- group_info$id
  design <- model.matrix(~0 + factor(group_info$group))
  colnames(design) <- levels(factor(group_info$group))
  contrast.matrix <- makeContrasts(contrasts = paste0(contrast1, "-", contrast2), levels = design)
  fit <- lmFit(expr.tpm, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  deg_list <- topTable(fit2, coef = 1, n = Inf) %>% na.omit()
  deg_list$gene <- row.names(deg_list)

  # Save results to output directory
  output_file <- file.path(output_dir, "differential_expression_results.csv")
  write.csv(deg_list, output_file, row.names = FALSE)

  return(deg_list)
}

# Example usage
# expr <- readRDS('./path_to_your_expression_data.rds')
# group <- TCGA_group(expr)
# results <- differential_expression_analysis(expr, group, output_dir = "./results")


# ============== 8. TCGA tumor extract =================
#' @title Data Preprocessing for Tumor Samples
#' @description This function preprocesses the expression, survival, and clinical data for tumor samples.
#' @param expr A matrix or data frame containing TPM expression data with genes as rows and samples as columns.
#' @param surv A data frame containing survival data with columns 'sample', 'OS', and 'OS.time'.
#' @param cli A data frame containing clinical data with a column 'submitter_id.samples'.
#' @param output_file A character string specifying the file path to save the preprocessed data as an RData file. Default is "./output_data/preprocessed_data.RData".
#' @return A list containing the preprocessed expression, survival, and clinical data for tumor samples.
#' @export
#' @import dplyr
#' @import IOBR
#' @examples
#' \dontrun{
#' expr <- readRDS('./path_to_your_expression_data.rds')
#' surv <- read.csv('./path_to_your_survival_data.csv')
#' cli <- read.csv('./path_to_your_clinical_data.csv')
#' TCGA_Tumor <- TCGAExtractTumor(expr, surv, cli, output_file = "./output_data/TCGAExtractTumor.RData")
#' }

TCGAExtractTumor <- function(expr, surv, cli, output_file = "./output_data/TCGAExtractTumor.RData") {
  # Ensure necessary libraries are loaded
  library(dplyr)
  library(IOBR)

  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate group if not provided
  group <- TCGA_group(expr)

  # Filter data for tumor samples
  tumorIndex <- group == 'tumor'
  exp <- expr[, tumorIndex]
  surv <- surv[tumorIndex, ]
  cli <- cli[tumorIndex, ]

  # Ensure sample IDs match across datasets
  if (!identical(colnames(exp), surv$sample) || !identical(colnames(exp), cli$submitter_id.samples)) {
    stop("Sample IDs do not match across expression, survival, and clinical data.")
  }

  # Save preprocessed data
  save(exp, surv, cli, file = output_file)

  return(list(expression = exp, survival = surv, clinical = cli))
}

# Example usage
# expr <- readRDS('./path_to_your_expression_data.rds')
# surv <- read.csv('./path_to_your_survival_data.csv')
# cli <- read.csv('./path_to_your_clinical_data.csv')
# preprocessed_data <- preprocess_tumor_data(expr, surv, cli, output_file = "./output_data/preprocessed_data.RData")

# ============== 9. Unicox Regression =================

#' @title Differential Expression and Survival Analysis
#' @description This function performs differential expression analysis and identifies prognostic genes using Cox proportional hazards model.
#' @param genes A vector of gene names to be included in the survival analysis.
#' @param expr.tpm A matrix or data frame containing only tumor samples' TPM expression data with genes as rows and samples as columns.
#' @param surv A data frame containing survival data with columns 'sample', 'OS', and 'OS.time'.
#' @param cli A data frame containing clinical data with a column 'submitter_id.samples'.
#' @param output_dir A character string specifying the directory to save the results. Default is "./output_data".
#' @return A list containing the differential expression analysis results and the survival analysis results.
#' @export
#' @import dplyr
#' @import limma
#' @import survival
#' @examples
#' \dontrun{
#' genes <- c("gene1", "gene2", "gene3")
#' results <- TCGAUniCox(genes, expr.tpm = expr, surv = surv, cli = cli, output_dir = "./output_data")
#' }

TCGAUniCox <- function(genes, expr.tpm, surv, cli, width = 8, height = 6,output_dir = "./output_data") {
  # Ensure necessary libraries are loaded
  library(dplyr)
  library(limma)
  library(survival)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Ensure sample IDs match across datasets
  if (!identical(colnames(expr.tpm), surv$sample) || !identical(colnames(expr.tpm), cli$submitter_id.samples)) {
    stop("Sample IDs do not match across expression, survival, and clinical data.")
  }

  # Prepare data for survival analysis
  rt <- data.frame('fustat' = surv$OS, 'futime' = surv$OS.time / 30, t(expr.tpm[genes, ]))

  # Perform Cox proportional hazards analysis
  outTab <- data.frame()
  sigGenes <- c()
  for(i in colnames(rt[,3:ncol(rt)])){
    #cox分析
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<0.05){
      sigGenes=c(sigGenes,i)
      outTab=rbind(outTab,
                   cbind(id=i,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }

  # Save results to output directory
  write.table(outTab, file = file.path(output_dir, "uniCox.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(rt, file = file.path(output_dir, "rt.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Plot forest
  source(system.file("data", "bioForest.R", package = "easySingleCell"))
  bioForest(coxFile = file.path(output_dir, "uniCox.txt"),
            forestFile = file.path(output_dir, "forest.pdf"))

  return(list(cox_results = outTab, survival_data = rt))
}

# Example usage
# expr <- readRDS('./path_to_your_expression_data.rds')
# surv <- read.csv('./path_to_your_survival_data.csv')
# cli <- read.csv('./path_to_your_clinical_data.csv')
# genes <- c("gene1", "gene2", "gene3")
# results <- TCGAUniCox(genes, expr.tpm = expr, surv = surv, cli = cli, output_dir = "./output_data")

# ============== 10. TNM Stage =================

#' @title Plot Gene Expression by Pathologic Stages
#' @description This function plots gene expression levels by pathologic stages (T, N, M) using boxplots.
#' @param genes A character vector specifying the genes of interest.
#' @param exp A data frame or matrix containing gene expression data.
#' @param cli A data frame containing clinical data with columns for pathologic_T, pathologic_N, and pathologic_M.
#' @param output_dir A string specifying the output directory. Default is "./output_figure/".
#' @param output_filename A string specifying the output filename for the plots. Default is "gene_expression_by_stage.pdf".
#' @param palette A character vector specifying the colors for the boxplots. Default is c("#F1788D", "#54990F", "#E6550D", "#843C39").
#' @param width Width of the PDF file. Default is 10.
#' @param height Height of the PDF file. Default is 6.
#' @param method A string specifying the statistical test method for `stat_compare_means`. Default is "kruskal.test".
#' @return None. The function saves the boxplots in a PDF file.
#' @export
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import stringr
#' @import ggpubr
#' @examples
#' \dontrun{
#' plot_gene_expression_by_stage(
#'   genes = c("Gene1", "Gene2"),
#'   exp = expression_data,
#'   cli = clinical_data,
#'   output_dir = "./output_figure/",
#'   output_filename = "gene_expression_by_stage.pdf",
#'   palette = c("#F1788D", "#54990F", "#E6550D", "#843C39"),
#'   method = "kruskal.test"
#' )
#' }

TCGAStagePlot <- function(genes, exp, cli,
                          output_dir = "./output_figure/",
                          output_filename = "gene_expression_by_stage.pdf",
                          palette = c("#F1788D", "#54990F", "#E6550D", "#843C39"),
                          width = 10, height = 6,
                          method = "kruskal.test") {

  # Ensure necessary libraries are loaded
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(ggpubr)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Prepare clinical data
  rt <- cli %>% dplyr::select(pathologic_T, pathologic_N, pathologic_M)
  rt <- data.frame(t(exp[genes, ]), rt)
  rt <- rt %>%
    drop_na() %>%
    filter(pathologic_M != "MX") %>%
    filter(pathologic_T != "Tis") %>%
    filter(pathologic_N != "NX")

  # Process pathologic_T
  rt$pathologic_T <- str_extract(rt$pathologic_T, "T[0-9]")
  rt <- rt %>% mutate("pathologic_T" =
                        case_when(
                          pathologic_T == 'T1' ~ 'T1',
                          pathologic_T == 'T2' ~ 'T2',
                          TRUE ~ 'T3/4'
                        ))
  rt$pathologic_T <- factor(rt$pathologic_T, levels = c("T1", "T2", "T3/4"))

  # Process pathologic_N
  rt$pathologic_N <- str_extract(rt$pathologic_N, "N[0-9]")
  rt$pathologic_N <- factor(rt$pathologic_N, levels = paste0('N', seq(length(unique(rt$pathologic_N)))-1))

  # Process pathologic_M
  rt$pathologic_M <- str_extract(rt$pathologic_M, "M[0-9]")
  rt$pathologic_M <- factor(rt$pathologic_M, levels = paste0('M', seq(length(unique(rt$pathologic_M)))-1))

  # Pivot longer for gene expression
  rt <- rt %>% pivot_longer(cols = all_of(genes),
                            values_to = 'gene_expression',
                            names_to = 'gene')

  # Plot function
  plot_boxplot <- function(data, stage, palette, method) {
    ggplot(data, aes(x = gene, y = gene_expression, fill = !!sym(stage))) +
      geom_boxplot(outlier.shape = 21, color = "black") +
      scale_fill_manual(values = palette) +
      theme_bw() +
      labs(x = NULL, y = "Gene expression") +
      theme(legend.position = "top") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(color = "black", size = 12)) +
      stat_compare_means(aes(group = !!sym(stage), label = ..p.signif..),
                         method = method)
  }

  # Generate plots
  pdf(file.path(output_dir, output_filename), width = width, height = height)

  # Plot for pathologic_T
  plot2 <- plot_boxplot(rt, "pathologic_T", palette, method)
  print(plot2)

  # Plot for pathologic_N
  plot1 <- plot_boxplot(rt, "pathologic_N", palette, method)
  print(plot1)

  # Plot for pathologic_M
  plot3 <- plot_boxplot(rt, "pathologic_M", palette, method)
  print(plot3)

  dev.off()
}

# ============== 11. TCGA survival =================
#' @title Perform Survival Analysis and Plot Survival Curves
#' @description This function performs survival analysis for a list of genes and plots the survival curves.
#' @param rt A data frame containing survival data and gene expression data. The first two columns should be 'futime' and 'fustat'.
#' @param genes A character vector specifying the genes for which to perform survival analysis.
#' @param palette A character vector specifying the colors for the survival plots.
#' @param output_dir A string specifying the output directory for the survival plots. Default is "./output_figure/".
#' @param width Width of the PDF file. Default is 5.5.
#' @param height Height of the PDF file. Default is 4.
#' @return None. The function saves the survival plots in PDF files.
#' @export
#' @import survminer
#' @import survival
#' @examples
#' \dontrun{
#' rt <- read.table('./output_data/rt.txt', header=T, sep="\t", check.names=F)
#' genes <- c("PLA2R1", "FAM20A")
#' palette <- c("#F1788D", "#54990F")
#' TCGASurvivalPlot(rt, genes, palette, output_dir = "./output_figure/")
#' }

TCGASurvivalPlot <- function(rt, genes, palette = c("#F1788D", "#54990F"),
                             output_dir = "./output_figure/", width = 5.5, height = 4) {
  # Ensure necessary libraries are loaded
  library(survminer)
  library(survival)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Determine the best cutoff values for gene expression
  cutoff <- surv_cutpoint(rt, time = "futime", event = "fustat", variables = genes)
  summary(cutoff) # Output the summary of cutoff

  # Categorize continuous variables based on the cutoff values
  groups <- surv_categorize(cutoff)
  str(groups)
  head(groups)

  # Loop through each gene in the genes vector
  for (gene in genes) {
    message("Processing gene: ", gene)

    # Assign gene to a global variable
    gene <<- gene # It is necessary, or the as.formula can not recongonize!!

    # Construct the survival analysis formula and fit the model
    fit <- survfit(as.formula(paste0('Surv(futime, fustat) ~ ', gene)), data = groups)

    # Plot the survival curve using ggsurvplot
    p <- ggsurvplot(fit,
                    data = groups,                   # Use the categorized data
                    pval = TRUE, # Show p-value
                    legend = c(0.75, 0.75),
                    xlab='Time(Months)',
                    ylab='OS(Overall Survival)',
                    ggtheme =  ggpubr::theme_pubr(base_size = 16),
                    pval.method = TRUE,              # Show p-value method
                    palette = palette,             # Use the specified palette
                    risk.table = TRUE,               # Show risk table
                    conf.int = FALSE)              # Show confidence interval (CI)

    # Save the plot to a PDF file
    ggsave(filename = paste0(output_dir, 'Surv_', gene, '.pdf'),
           plot = p$plot,                           # Use the plot object returned by ggsurvplot
           width = width,
           height = height)

    # Remove the global gene variable
    gene <<- NULL
  }

}

# Example usage
# rt <- read.table('./output_data/rt.txt', header=T, sep="\t", check.names=F)
# genes <- c("PLA2R1", "FAM20A")
# palette <- c("#F1788D", "#54990F")
# TCGASurvivalPlot(rt, genes, palette, output_dir = "./output_figure/")
