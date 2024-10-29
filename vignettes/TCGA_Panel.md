# TCGAPanel: A Comprehensive Tool for TCGA Data Analysis

The `TCGAPanel` R package provides a comprehensive suite of tools for downloading, cleaning, and analyzing TCGA (The Cancer Genome Atlas) data. This package is designed to facilitate various tasks such as differential expression analysis, survival analysis, and visualization of gene expression data by clinical stages. Below is a detailed description of the functions available in the `TCGAPanel` package:

## 1. get_TCGA

### Description
This function creates an object for either the TCGA or PanCancer project, downloads the relevant data, and cleans it.

### Usage
```r
get_TCGA(cancer_type = 'NULL', project_type = 'TCGA')
```

### Arguments
- **cancer_type**: The type of cancer to study. Default is 'NULL'.
- **project_type**: The type of project ('TCGA' or 'PanCancer'). Default is 'TCGA'.

### Example
```r
get_TCGA(cancer_type = 'COAD', project_type = 'TCGA')
```

## 2. TCGA_group

### Description
This function classifies samples into 'tumor' and 'normal' based on TCGA sample barcodes.

### Usage
```r
TCGA_group(expr)
```

### Arguments
- **expr**: A matrix or data frame containing expression data with genes as rows and samples as columns.

### Example
```r
expr <- readRDS('./path_to_your_expression_data.rds')
group <- TCGA_group(expr)
```

## 3. TCGADegLimma

### Description
This function performs differential expression analysis using limma.

### Usage
```r
TCGADegLimma(expr.tpm, group = NULL, contrast1 = 'tumor', contrast2 = 'normal', output_dir = "./output_data")
```

### Arguments
- **expr.tpm**: A matrix or data frame containing TPM expression data with genes as rows and samples as columns.
- **group**: A factor vector indicating the group ('tumor' or 'normal') for each sample. If not provided, the function will generate it.
- **contrast1**: A character string specifying the first group for comparison. Default is 'tumor'.
- **contrast2**: A character string specifying the second group for comparison. Default is 'normal'.
- **output_dir**: A character string specifying the directory to save the results. Default is "./output_data".

### Example
```r
expr <- readRDS('./path_to_your_expression_data.rds')
group <- TCGA_group(expr)
results <- TCGADegLimma(expr, group, output_dir = "./output_data")
```

## 4. TCGAExtractTumor

### Description
This function preprocesses the expression, survival, and clinical data for tumor samples.

### Usage
```r
TCGAExtractTumor(expr, surv, cli, output_file = "./output_data/TCGAExtractTumor.RData")
```

### Arguments
- **expr**: A matrix or data frame containing TPM expression data with genes as rows and samples as columns.
- **surv**: A data frame containing survival data with columns 'sample', 'OS', and 'OS.time'.
- **cli**: A data frame containing clinical data with a column 'submitter_id.samples'.
- **output_file**: A character string specifying the file path to save the preprocessed data as an RData file. Default is "./output_data/TCGAExtractTumor.RData".

### Example
```r
expr <- readRDS('./path_to_your_expression_data.rds')
surv <- read.csv('./path_to_your_survival_data.csv')
cli <- read.csv('./path_to_your_clinical_data.csv')
TCGA_Tumor <- TCGAExtractTumor(expr, surv, cli, output_file = "./output_data/TCGAExtractTumor.RData")
```

## 5. TCGAUniCox

### Description
This function performs differential expression analysis and identifies prognostic genes using Cox proportional hazards model.

### Usage
```r
TCGAUniCox(genes, expr.tpm, surv, cli, width = 8, height = 6, output_dir = "./output_data")
```

### Arguments
- **genes**: A vector of gene names to be included in the survival analysis.
- **expr.tpm**: A matrix or data frame containing only tumor samples' TPM expression data with genes as rows and samples as columns.
- **surv**: A data frame containing survival data with columns 'sample', 'OS', and 'OS.time'.
- **cli**: A data frame containing clinical data with a column 'submitter_id.samples'.
- **width**: Width of the PDF file. Default is 8.
- **height**: Height of the PDF file. Default is 6.
- **output_dir**: A character string specifying the directory to save the results. Default is "./output_data".

### Example
```r
genes <- c("gene1", "gene2", "gene3")
results <- TCGAUniCox(genes, expr.tpm = expr, surv = surv, cli = cli, output_dir = "./output_data")
```

## 6. TCGAStagePlot

### Description
This function plots gene expression levels by pathologic stages (T, N, M) using boxplots.

### Usage
```r
TCGAStagePlot(genes, exp, cli, output_dir = "./output_figure/", output_filename = "gene_expression_by_stage.pdf", palette = c("#F1788D", "#54990F", "#E6550D", "#843C39"), width = 10, height = 6, method = "kruskal.test")
```

### Arguments
- **genes**: A character vector specifying the genes of interest.
- **exp**: A data frame or matrix containing gene expression data.
- **cli**: A data frame containing clinical data with columns for pathologic_T, pathologic_N, and pathologic_M.
- **output_dir**: A string specifying the output directory. Default is "./output_figure/".
- **output_filename**: A string specifying the output filename for the plots. Default is "gene_expression_by_stage.pdf".
- **palette**: A character vector specifying the colors for the boxplots. Default is c("#F1788D", "#54990F", "#E6550D", "#843C39").
- **width**: Width of the PDF file. Default is 10.
- **height**: Height of the PDF file. Default is 6.
- **method**: A string specifying the statistical test method for `stat_compare_means`. Default is "kruskal.test".

### Example
```r
plot_gene_expression_by_stage(
  genes = c("Gene1", "Gene2"),
  exp = expression_data,
  cli = clinical_data,
  output_dir = "./output_figure/",
  output_filename = "gene_expression_by_stage.pdf",
  palette = c("#F1788D", "#54990F", "#E6550D", "#843C39"),
  method = "kruskal.test"
)
```

## 7. TCGASurvivalPlot

### Description
This function performs survival analysis for a list of genes and plots the survival curves.

### Usage
```r
TCGASurvivalPlot(rt, genes, palette = c("#F1788D", "#54990F"), output_dir = "./output_figure/", width = 6, height = 4)
```

### Arguments
- **rt**: A data frame containing survival data and gene expression data. The first two columns should be 'futime' and 'fustat'.
- **genes**: A character vector specifying the genes for which to perform survival analysis.
- **palette**: A character vector specifying the colors for the survival plots.
- **output_dir**: A string specifying the output directory for the survival plots. Default is "./output_figure/".
- **width**: Width of the PDF file. Default is 6.
- **height**: Height of the PDF file. Default is 4.

### Example
```r
rt <- read.table('./output_data/rt.txt', header=T, sep="\t", check.names=F)
genes <- c("PLA2R1", "FAM20A")
palette <- c("#F1788D", "#54990F")
TCGASurvivalPlot(rt, genes, palette, output_dir = "./output_figure/")
```

This documentation provides a comprehensive guide to using the `TCGAPanel` package for TCGA data analysis. Each function is described in detail with its usage, arguments, and example code. This should help users to effectively utilize the package for their research needs.
