# Single cell pipeline Panel

This document provides detailed descriptions and usage examples for various functions designed to facilitate single-cell RNA-seq data analysis using the Seurat package and other related tools.

## Table of Contents

1. [Filter and Plot scRNA-seq Data](#filter-and-plot-scrna-seq-data)
2. [Run Normalization on a Seurat Object](#run-normalization-on-a-seurat-object)
3. [Generate and Save Clustree Plot](#generate-and-save-clustree-plot)
4. [Mark Doublets in Seurat Object](#mark-doublets-in-seurat-object)
5. [Perform Metabolism Analysis and Generate DotPlot](#perform-metabolism-analysis-and-generate-dotplot)
6. [Run GSEA Analysis For ScRNA](#run-gsea-analysis-for-scrna)
7. [Run CytoTRACE Analysis](#run-cytotrace-analysis)
8. [Run CoGAPS Analysis on Single-cell Data](#run-cogaps-analysis-on-single-cell-data)
9. [Process scRNA Data with Monocle](#process-scrna-data-with-monocle)
10. [Run CellChat Analysis](#run-cellchat-analysis)
11. [Perform Initial HdWGCNA Analysis to Find Soft Threshold](#perform-initial-hdwgcna-analysis-to-find-soft-threshold)
12. [Perform HdWGCNA Analysis Step 2](#perform-hdwgcna-analysis-step-2)
13. [Run HdWGCNA Step 3](#run-hdwgcna-step-3)
14. [Run MistyR Analysis](#run-mistyr-analysis)
15. [Run InferCNV Analysis](#run-infercnv-analysis)
16. [Gene Naming Conversion Functions](#gene-naming-conversion-functions)
17. [Import SCENIC Loom Files into Seurat](#import-scenic-loom-files-into-seurat)

## Filter and Plot scRNA-seq Data

### Description

This function filters scRNA-seq data based on specified criteria, generates violin plots for QC metrics, and saves the plots as PDF files.

### Usage

```r
runScRNAQC(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000,
           species = 'human', pal = c("#F1788D", "#54990F","#E6550D","#843C39", 
           "#3182BD","#8C6D31", "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", 
           "#6BAED6", "#9ECAE1", "#AD494A", "#A1D99B", "#C7E9C0", "#99600F",
           "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B", "#66A61E",
           "#E6550D", "#E7969C",'#53A85F'), output_dir = "./output_figure/", 
           width = 6, height = 4)
```

### Example

```r
sce <- CreateSeuratObject(counts = your_data)
pal <- c("#F1788D", "#54990F")
filtered_sce <- runScRNAQC(sce, minGene = 200, maxGene = 5000, pctMT = 15, maxCounts = 20000, species = 'human', pal = pal)
```

## Run Normalization on a Seurat Object

### Description

This function normalizes a Seurat object and performs dimensionality reduction.

### Usage

```r
run_normalize(seurat_obj, dims = 1:30, batch_var = "orig.ident")
```

### Example

```r
seurat_obj <- CreateSeuratObject(counts = your_data)
normalized_sce <- run_normalize(seurat_obj = seurat_obj, dims = 1:30, batch_var = "batch")
```

## Generate and Save Clustree Plot

### Description

This function generates and saves a clustree plot for a Seurat object across specified clustering resolutions.

### Usage

```r
run_clustree(sce, resolutions = seq(0.1, 1.4, 0.2), prefix = "RNA_snn_res.", output_dir = "./output_figure/", output_filename = "clustree_plot.pdf")
```

### Example

```r
sce <- CreateSeuratObject(counts = your_data)
clustree_plot <- run_clustree(sce, resolutions = seq(0.1, 1.4, 0.2), prefix = "RNA_snn_res.")
```

## Mark Doublets in Seurat Object

### Description

This function wraps the MarkDoublets function from the DoubletFinder package to identify doublets in a Seurat object.

### Usage

```r
MarkDoubletsWrapper(seu, PCs = 1:10, split.by = "orig.ident")
```

### Example

```r
seu <- CreateSeuratObject(counts = your_data)
seu <- MarkDoubletsWrapper(seu = seu, PCs = 1:10, split.by = "orig.ident")
```

## Perform Metabolism Analysis and Generate DotPlot

### Description

This function performs metabolism analysis on Seurat object data and generates a DotPlot of metabolic pathways.

### Usage

```r
runScMetabolismAnalysis(scRNA, output_file = './output_data/scMetabolism.Rdata', pdf_file = './output_figure/scMetabolism.pdf', npathways = 20, cell_type_column = "celltype", width = 6, height = 10, method = "VISION")
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
scRNA <- runScMetabolismAnalysis(scRNA = scRNA, npathways = 20, cell_type_column = "celltype")
```

## Run GSEA Analysis For ScRNA

### Description

This function performs Gene Set Enrichment Analysis (GSEA) on Seurat object data.

### Usage

```r
runGSEA(scRNA, reference_group, comparison_group, output_file = './output_data/GSEA_results.Rdata', pdf_file = './output_figure/GSEA.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runGSEA(scRNA = scRNA, reference_group = "group1", comparison_group = "group2")
```

## Run CytoTRACE Analysis

### Description

This function runs CytoTRACE analysis on single-cell RNA-seq data.

### Usage

```r
runCytoTRACE(scRNA, output_file = './output_data/CytoTRACE_results.Rdata', pdf_file = './output_figure/CytoTRACE.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runCytoTRACE(scRNA = scRNA)
```

## Run CoGAPS Analysis on Single-cell Data

### Description

This function runs CoGAPS analysis on single-cell RNA-seq data.

### Usage

```r
runCoGAPS(scRNA, output_file = './output_data/CoGAPS_results.Rdata', pdf_file = './output_figure/CoGAPS.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runCoGAPS(scRNA = scRNA)
```

## Process scRNA Data with Monocle

### Description

This function processes single-cell RNA-seq data using Monocle for trajectory analysis.

### Usage

```r
runMonocle(scRNA, output_file = './output_data/Monocle_results.Rdata', pdf_file = './output_figure/Monocle.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runMonocle(scRNA = scRNA)
```

## Run CellChat Analysis

### Description

This function runs CellChat analysis on single-cell RNA-seq data to analyze cell-cell communication.

### Usage

```r
runCellChat(scRNA, output_file = './output_data/CellChat_results.Rdata', pdf_file = './output_figure/CellChat.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runCellChat(scRNA = scRNA)
```

## Perform Initial HdWGCNA Analysis to Find Soft Threshold

### Description

This function performs the initial step of HdWGCNA analysis to find the soft threshold.

### Usage

```r
runHdWGCNA_Step1(scRNA, output_file = './output_data/HdWGCNA_Step1_results.Rdata', pdf_file = './output_figure/HdWGCNA_Step1.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runHdWGCNA_Step1(scRNA = scRNA)
```

## Perform HdWGCNA Analysis Step 2

### Description

This function performs the second step of HdWGCNA analysis.

### Usage

```r
runHdWGCNA_Step2(scRNA, output_file = './output_data/HdWGCNA_Step2_results.Rdata', pdf_file = './output_figure/HdWGCNA_Step2.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runHdWGCNA_Step2(scRNA = scRNA)
```

## Run HdWGCNA Step 3

### Description

This function performs the third step of HdWGCNA analysis.

### Usage

```r
runHdWGCNA_Step3(scRNA, output_file = './output_data/HdWGCNA_Step3_results.Rdata', pdf_file = './output_figure/HdWGCNA_Step3.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runHdWGCNA_Step3(scRNA = scRNA)
```

## Run MistyR Analysis

### Description

This function runs MistyR analysis on single-cell RNA-seq data.

### Usage

```r
runMistyR(scRNA, output_file = './output_data/MistyR_results.Rdata', pdf_file = './output_figure/MistyR.pdf', width = 6, height = 10)
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
runMistyR(scRNA = scRNA)
```

## Run InferCNV Analysis

### Description

This function runs InferCNV analysis on single-cell RNA-seq data to infer copy number variations.

### Example

```r
runInferCNVPipeline(sce_epi = sce.hepa,
                    sce_refer = sce.refer,
                    celltype = 'celltype',
                    infercnv_path = './output_data/Figure6/inferCNV/',
                    name = 'infer_run',
                    ref_group_names = c("Tcell", "Bcell"),
                    ref_cell_name = c("Tcell", "Bcell"),
                    obs_cell_name = "Hepatocyte",
                    ref_group = "immune",
                    obs_group = "Hepatocyte",
                    k_clusters = 5,
                    heatmap_colors = c("#2166ac", "white", "#b2182b"),
                    cluster_colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"),
                    num_threads = 16
                    )
```

## Gene Naming Conversion Functions

### Description

This function converts gene names between different formats (e.g., Ensembl IDs to gene symbols).

### Usage

```r
convertGeneNames(gene_list, from_format = "ENSEMBL", to_format = "SYMBOL")
```

### Example

```r
gene_list <- c("ENSG00000141510", "ENSG00000171862")
converted_genes <- convertGeneNames(gene_list, from_format = "ENSEMBL", to_format = "SYMBOL")
```

## Import SCENIC Loom Files into Seurat

### Description

This function imports SCENIC loom files into a Seurat object for further analysis.

### Usage

```r
importSCENIC(scRNA, loom_file, output_file = './output_data/SCENIC_results.Rdata')
```

### Example

```r
scRNA <- CreateSeuratObject(counts = your_data)
importSCENIC(scRNA = scRNA, loom_file = "path_to_loom_file.loom")
```

By following these guidelines, you should be able to effectively use the provided functions for single-cell RNA-seq data analysis using Seurat and other related tools.
