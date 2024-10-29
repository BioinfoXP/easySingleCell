# Load necessary libraries
library(DT)
library(tidyverse)
library(data.table)

# ========================= 1. Define process function============================
# ======== 定义下载函数
download_data <- function(file_link, output_dir, overwrite = FALSE) {
  path <- file.path(output_dir, basename(file_link))
  if (!file.exists(path) || overwrite) {
    download.file(url = file_link, destfile = path)
    message("Downloaded: ", path)
  } else {
    message("File already exists and will not be overwritten: ", path)
  }
}

# ======== 定义处理表达数据的函数
process_expression_data <- function(exp_data, coding_genes) {
  # 过滤和去重 protein_coding 数据
  coding_gene <- coding_genes %>%
    filter(gene_type == "protein_coding") %>%
    distinct(gene_name, .keep_all = TRUE)

  # 获取共同的基因
  com_gene <- intersect(exp_data$Ensembl_ID, coding_gene$gene_id)
  exp_data <- exp_data[exp_data$Ensembl_ID %in% com_gene, ]
  row.names(exp_data) <- NULL
  exp_data <- exp_data %>% tibble::column_to_rownames(var = "Ensembl_ID")
  exp_data <- exp_data[coding_gene$gene_id, ]

  if (identical(row.names(exp_data), coding_gene$gene_id)) {
    row.names(exp_data) <- coding_gene$gene_name
  } else {
    stop("Row names do not match gene IDs in coding_gene.")
  }

  return(exp_data)
}

# =====================2. get_TCGA ==============
#' Get TCGA Data for a Specific Cancer Type
#'
#' This function downloads and processes TCGA (The Cancer Genome Atlas) data for a given cancer type.
#' It retrieves expression, clinical, and survival data, processes the expression data to include only
#' protein-coding genes, and saves the processed data to an Rdata file.
#'
#' @param cancer_type A character string representing the cancer type abbreviation (e.g., "BRCA" for breast cancer).
#' @param overwrite A logical value indicating whether to overwrite existing files during the download process. Default is FALSE.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{exp_tpm}{The raw expression data in TPM format.}
#'   \item{exp_tpm_coding_only}{The processed expression data, including only protein-coding genes (TPM).}
#'   \item{exp_counts}{The raw expression data in counts format.}
#'   \item{exp_counts_coding_only}{The processed expression data, including only protein-coding genes (counts).}
#'   \item{clinical}{The clinical data for the cancer type.}
#'   \item{survival}{The survival data for the cancer type.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   result <- get_TCGA("BRCA", overwrite = TRUE)
#' }
#'
#' @export
get_TCGA <- function(cancer_type, overwrite = FALSE) {

  load(system.file("data", "protein_coding.Rdata", package = "easySingleCell"))
  load(system.file("data", "TCGA_Table.Rdata", package = "easySingleCell"))

  # 检查 cancer_type 是否在 tcga.table 中
  index <- which(tcga.table$Abbreviation == cancer_type)
  if (length(index) == 0) {
    message(paste("Cancer type", cancer_type, "not found in tcga.table. Please select a valid cancer type from the table below:"))
    datatable(tcga.table[c(1, 2)])
    return(NULL)
  }
  names(index) <- cancer_type

  # 提取 URL 信息
  url <- tcga.table[index, , drop = TRUE]

  # 创建输出目录
  output_dir <- paste0('./output_xena/', names(index))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 定义要下载的文件类型
  file_types <- c("star_tpm.Link", "clinical.Link", "survival.Link", "star_counts.Link")

  # 循环下载每种文件类型
  message("Downloading files...")
  for (file_type in file_types) {
    download_data(url[[file_type]], output_dir, overwrite = overwrite)
  }

  # 获取下载的文件列表
  files <- list.files(path = paste0('./output_xena/', cancer_type), full.names = TRUE)

  # 读取表达数据文件
  message("Reading expression data files...")
  exp.tpm <- fread(files[grep('star_tpm', files)], data.table = FALSE)
  exp.counts <- fread(files[grep('star_counts', files)], data.table = FALSE)

  # 使用函数处理表达数据
  message("Processing expression data...")
  exp.tpm.CodingOnly <- process_expression_data(exp.tpm, protein_coding)
  exp.counts.CodingOnly <- process_expression_data(exp.counts, protein_coding)

  # 读取临床和生存数据文件
  message("Reading clinical and survival data files...")
  cli <- fread(files[grep('clinical', files)], data.table = FALSE)
  surv <- fread(files[grep('survival', files)], data.table = FALSE)

  # 匹配临床和生存数据
  message("Matching clinical and survival data...")
  cli <- cli[match(colnames(exp.tpm.CodingOnly), cli$sample), ];row.names(cli)=NULL
  surv <- surv[match(colnames(exp.tpm.CodingOnly), surv$sample), ];row.names(surv)=NULL

  # 保存处理后的数据为 Rdata 文件
  output_file <- file.path(output_dir, paste0(cancer_type, "_processed.Rdata"))
  save(exp.tpm.CodingOnly, exp.counts.CodingOnly, cli, surv, file = output_file)
  message("Processed data saved to: ", output_file)

  # 返回结果
  list(
    exp_tpm = exp.tpm,
    exp_tpm_coding_only = exp.tpm.CodingOnly,
    exp_counts = exp.counts,
    exp_counts_coding_only = exp.counts.CodingOnly,
    clinical = cli,
    survival = surv
  )
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

# ======== 3. TCGA Group =============
TCGA_group <- function(expr) {
  group <- sapply(strsplit(colnames(expr), "\\-"), "[", 4)
  group <- sapply(strsplit(group, ""), "[", 1)
  group_list <- ifelse(group == "0", 'tumor', 'normal')
  group_list <- factor(group_list, levels = c('normal', 'tumor'))
  return(group_list)
}


# ======== 4. TCGA DEG limma=============
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


# ============== 5. TCGA tumor extract =================
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

# ============== 6. Unicox Regression =================

# bioForest.R

#' Create a Forest Plot
#'
#' This function creates a forest plot based on Cox proportional hazards model results.
#'
#' @param coxFile Path to the input file containing Cox model results. The file should be a tab-separated text file with columns: "HR", "HR.95L", "HR.95H", and "pvalue".
#' @param forestFile Path to the output PDF file where the forest plot will be saved.
#'
#' @return None. The function outputs a PDF file with the forest plot.
#' @importFrom graphics arrows axis box plot points text
#' @importFrom utils read.table
#' @export
bioForest <- function(coxFile = NULL, forestFile = NULL) {
  # 读取输入文件
  rt <- read.table(coxFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f", rt$HR)
  hrLow  <- sprintf("%.3f", rt$HR.95L)
  hrHigh <- sprintf("%.3f", rt$HR.95H)
  Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
  pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))

  # 输出图形
  pdf(file = forestFile, width = 6.6, height = 14)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))

  # 绘制森林图左边的临床信息
  xlim = c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
  text.cex = 0.8
  text(0, n:1, gene, adj = 0, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard.ratio, adj = 1, cex = text.cex)
  text(3.1, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)

  # 绘制右边森林图
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
  abline(v = 1, col = "black", lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, "#F1788D", "#54990F")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()
}















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

# ============== 7. TNM Stage =================

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

# ============== 8. TCGA survival =================
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
