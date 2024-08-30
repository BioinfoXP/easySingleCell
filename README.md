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

## Installation

To install the development version of `easySingleCell`, use the following commands:

Tips: For R version 4.2.2 (2022-10-31), it is recommed.
- Seurat 4.4.0
- SeuratObject 4.1.4
- Matrix 1.5-1
- Matrix 1.6-1
> Seurat: https://github.com/satijalab/seurat/releases
> SeuratObject: https://github.com/satijalab/seurat-object/releases
> Matrix: https://cran.r-project.org/src/contrib/Archive/Matrix/

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install dependencies from Bioconductor
BiocManager::install("mistyR")
BiocManager::install("CoGAPS")
BiocManager::install("Mfuzz")

# Install additional dependencies from GitHub
devtools::install_github("PaulingLiu/ROGUE")
devtools::install_github("junjunlab/ClusterGVis")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("xmc811/Scillus", ref = "development")
remotes::install_github("FertigLab/CoGAPS")

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
