# easySingleCell

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

`easySingleCell` is an R package designed to streamline bioinformatics analysis for single-cell RNA sequencing (scRNA-seq) data. It provides a suite of functions to perform various analyses, including colocalization analysis using MistyR and Seurat, normalization, clustering, trajectory analysis, and more. The package aims to simplify the handling of complex single-cell data, making it accessible for researchers working with data from sources such as The Cancer Genome Atlas (TCGA) and other single-cell datasets.

## Features

- **Colocalization Analysis**: Utilize MistyR and Seurat for colocalization analysis, facilitating the study of spatial gene expression patterns.
- **Normalization and Clustering**: Functions for normalizing scRNA-seq data and performing clustering to identify distinct cell populations.
- **Trajectory Analysis**: Tools for trajectory inference to study cell differentiation processes.
- **Metabolic Pathway Analysis**: Analyze metabolic pathways and generate visualizations to understand cellular metabolism.
- **Gene Set Enrichment Analysis (GSEA)**: Perform GSEA to identify enriched pathways and gene sets.
- **Doublet Detection**: Identify and mark doublets in scRNA-seq data.
- **Integration with TCGA Data**: Handle and analyze TCGA data based on single-gene expressions.

## Environment prepare

To install the development version of `easySingleCell`, use the following commands:

Tips: For R version 4.2.2 (2022-10-31), it is recommed.
- Seurat 4.4.0
- SeuratObject 4.1.4
- Matrix 1.5-1
- Matrix 1.6-1
> Seurat: https://github.com/satijalab/seurat/releases
> 
> SeuratObject: https://github.com/satijalab/seurat-object/releases
> 
> Matrix: https://cran.r-project.org/src/contrib/Archive/Matrix/


## Installation

### Cran R packages

```r
# Install CRAN packages one by one, checking if they are already installed
if (!requireNamespace('tibble', quietly = TRUE)) install.packages('tibble')
if (!requireNamespace('survival', quietly = TRUE)) install.packages('survival')
if (!requireNamespace('survminer', quietly = TRUE)) install.packages('survminer')
if (!requireNamespace('limma', quietly = TRUE)) install.packages('limma')
if (!requireNamespace('DESeq2', quietly = TRUE)) install.packages('DESeq2')
if (!requireNamespace('limSolve', quietly = TRUE)) install.packages('limSolve')
if (!requireNamespace('GSVA', quietly = TRUE)) install.packages('GSVA')
if (!requireNamespace('e1071', quietly = TRUE)) install.packages('e1071')
if (!requireNamespace('preprocessCore', quietly = TRUE)) install.packages('preprocessCore')
if (!requireNamespace('tidyHeatmap', quietly = TRUE)) install.packages('tidyHeatmap')
if (!requireNamespace('caret', quietly = TRUE)) install.packages('caret')
if (!requireNamespace('glmnet', quietly = TRUE)) install.packages('glmnet')
if (!requireNamespace('ppcor', quietly = TRUE)) install.packages('ppcor')
if (!requireNamespace('timeROC', quietly = TRUE)) install.packages('timeROC')
if (!requireNamespace('pracma', quietly = TRUE)) install.packages('pracma')
if (!requireNamespace('factoextra', quietly = TRUE)) install.packages('factoextra')
if (!requireNamespace('FactoMineR', quietly = TRUE)) install.packages('FactoMineR')
if (!requireNamespace('patchwork', quietly = TRUE)) install.packages('patchwork')
if (!requireNamespace('ggplot2', quietly = TRUE)) install.packages('ggplot2')
if (!requireNamespace('biomaRt', quietly = TRUE)) install.packages('biomaRt')
if (!requireNamespace('ggpubr', quietly = TRUE)) install.packages('ggpubr')
if (!requireNamespace('ComplexHeatmap', quietly = TRUE)) install.packages('ComplexHeatmap')
if (!requireNamespace('export', quietly = TRUE)) install.packages('export')
if (!requireNamespace('harmony', quietly = TRUE)) install.packages('harmony')
if (!requireNamespace('vioplot', quietly = TRUE)) install.packages('vioplot')
```
### Biocmanager R packages

```r
# Install Bioconductor manager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install Bioconductor packages one by one, checking if they are already installed
if (!requireNamespace('mistyR', quietly = TRUE)) BiocManager::install('mistyR')
if (!requireNamespace('CoGAPS', quietly = TRUE)) BiocManager::install('CoGAPS')
if (!requireNamespace('Mfuzz', quietly = TRUE)) BiocManager::install('Mfuzz')
if (!requireNamespace('UCSCXenaTools', quietly = TRUE)) BiocManager::install('UCSCXenaTools')
if (!requireNamespace('WGCNA', quietly = TRUE)) BiocManager::install('WGCNA')
if (!requireNamespace('igraph', quietly = TRUE)) BiocManager::install('igraph')
if (!requireNamespace('GeneOverlap', quietly = TRUE)) BiocManager::install('GeneOverlap')
if (!requireNamespace('ggrepel', quietly = TRUE)) BiocManager::install('ggrepel')
if (!requireNamespace('UCell', quietly = TRUE)) BiocManager::install('UCell')
```
### Github R packages


```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install GitHub packages one by one using devtools, checking if they are already installed
if (!requireNamespace('ROGUE', quietly = TRUE)) devtools::install_github('PaulingLiu/ROGUE')
if (!requireNamespace('ClusterGVis', quietly = TRUE)) devtools::install_github('junjunlab/ClusterGVis')
if (!requireNamespace('ComplexHeatmap', quietly = TRUE)) devtools::install_github('jokergoo/ComplexHeatmap')
if (!requireNamespace('Scillus', quietly = TRUE)) devtools::install_github('xmc811/Scillus@development')
if (!requireNamespace('CoGAPS', quietly = TRUE)) devtools::install_github('FertigLab/CoGAPS')
if (!requireNamespace('easySingleCell', quietly = TRUE)) devtools::install_github('BioinfoXP/easySingleCell')
if (!requireNamespace('ggunchull', quietly = TRUE)) devtools::install_github('sajuukLyu/ggunchull')
if (!requireNamespace('scRNAtoolVis', quietly = TRUE)) devtools::install_github('junjunlab/scRNAtoolVis')
if (!requireNamespace('IOBR', quietly = TRUE)) devtools::install_github('IOBR/IOBR')
if (!requireNamespace('CellChat', quietly = TRUE)) devtools::install_github('jinworks/CellChat')
if (!requireNamespace('SCopeLoomR', quietly = TRUE)) devtools::install_github('aertslab/SCopeLoomR')
if (!requireNamespace('hdWGCNA', quietly = TRUE)) devtools::install_github('smorabit/hdWGCNA')
if (!requireNamespace('ktplots', quietly = TRUE)) devtools::install_github('zktuong/ktplots')
if (!requireNamespace('scMetabolism', quietly = TRUE)) devtools::install_github('wu-yc/scMetabolism')
if (!requireNamespace('VISION', quietly = TRUE)) devtools::install_github('YosefLab/VISION@v2.1.0')
if (!requireNamespace('copykat', quietly = TRUE))
devtools::install_github("navinlabcode/copykat")
```

### Install easySingleCell

```r
# Finally, install easySingleCell
remotes::install_github("BioinfoXP/easySingleCell", upgrade = F, dependencies = F)
```

## Vignettes

For more details, please read the files below!

### [TCGA_Panel](vignettes/TCGA_Panel.md)

The `TCGA_Panel` vignette provides detailed instructions and examples on how to use `easySingleCell` for analyzing TCGA data. It includes steps for data preprocessing, normalization, and various downstream analyses such as differential expression and pathway enrichment.

### [SingleCell_Panel](vignettes/SingleCell_Panel.md)

The `SingleCell_Panel` vignette focuses on single-cell RNA-seq data analysis. It covers the entire workflow from raw data import, quality control, normalization, clustering, trajectory analysis, and visualization. This vignette is designed to help users get started with single-cell data analysis using the `easySingleCell` package, providing practical examples and best practices.

By following these vignettes, users can effectively utilize the `easySingleCell` package to perform comprehensive analyses on both bulk and single-cell RNA-seq data, gaining deeper insights into their biological questions.
