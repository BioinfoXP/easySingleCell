# easySingleCell

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

`easySingleCell` is an R package designed to facilitate bioinformatics analysis for single-cell data. It provides functions to run colocalization analysis using MistyR and Seurat, making it easier to handle TCGA data based on a single gene.

## Installation

To install the development version of `easySingleCell`, use the following commands:

```r
# Install devtools if you haven't already
install.packages("devtools")

BiocManager::install("mistyR")
BiocManager::install("CoGAPS")
BiocManager::install("Mfuzz")

devtools::install_github("PaulingLiu/ROGUE")
devtools::install_github("junjunlab/ClusterGVis")
devtools::install_github("jokergoo/ComplexHeatmap")
remotes::install_github("BioinfoXP/easySingleCell",upgrade = F,dependencies = F)
remotes::install_github("FertigLab/CoGAPS")
```


# [TCGA_Panel](vignettes/TCGA_Panel.md)
