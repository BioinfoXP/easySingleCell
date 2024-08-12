# ============ 1. scVisFeaturePlot =============

#' @title Feature Plot for Single-cell Data
#' @description This function creates a feature plot for single-cell data using a specified reduction method.
#' @param scRNA A Seurat object containing single-cell data.
#' @param feature A character string specifying the feature to be plotted.
#' @param reduction A character string specifying the reduction method to be used (default is "umap").
#' @param pt.size A numeric value specifying the size of the points (default is 0.0001).
#' @param max.cutoff A numeric value specifying the maximum cutoff for the feature values (default is 1.5).
#' @param cols A vector of colors to be used for the plot.
#' @param title A character string specifying the title of the plot. If NULL, the feature name will be used as the title.
#' @return A ggplot object representing the feature plot.
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(ggplot2)
#'
#' # Example Seurat object
#' seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
#' seurat_obj <- NormalizeData(seurat_obj)
#' seurat_obj <- FindVariableFeatures(seurat_obj)
#' seurat_obj <- ScaleData(seurat_obj)
#' seurat_obj <- RunPCA(seurat_obj)
#' seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#'
#' # Example feature plot
#' plot <- scVisFeaturePlot(
#'   scRNA = seurat_obj,
#'   feature = "CD3E",
#'   reduction = "umap",
#'   pt.size = 0.5,
#'   max.cutoff = 1.5,
#'   cols =  c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),
#'   title = NULL
#' )
#' print(plot)
#' }

scVisFeaturePlot <- function(scRNA, feature, reduction = "umap", pt.size = 0.0001, max.cutoff = 1.5, cols =  c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"), title = NULL) {
  # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Set title to feature name if title is NULL
  if (is.null(title)) {
    title <- feature
  }

  # Create the feature plot
  plot <- FeaturePlot(
    object = scRNA,
    features = feature,
    reduction = reduction,
    pt.size = pt.size,
    max.cutoff = max.cutoff,
    cols = cols
  ) +
    scale_x_continuous("") +
    scale_y_continuous("") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    ggtitle(title)

  return(plot)
}

# Example usage
# library(Seurat)
# library(ggplot2)
#
# # Example Seurat object
# seurat_obj <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 100, ncol = 10))
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#
# # Example feature plot
# plot <- scVisFeaturePlot(
#   scRNA = seurat_obj,
#   feature = "PC_1",
#   reduction = "umap",
#   pt.size = 0.5,
#   max.cutoff = 1.5,
#   cols = c("blue", "red"),
#   title = NULL
# )
# print(plot)



# ============ 2. scVisDimPlot ==========
#' @title Plot Dimensionality Reduction Visualization
#' @description This function plots a dimensionality reduction visualization (e.g., UMAP, TSNE) of a Seurat object.
#' @param scRNA A Seurat object containing scRNA-seq data.
#' @param reduction A character string specifying the reduction method to use. Default is "umap".
#' @param group.by A character string specifying the metadata column to group by. Default is 'celltype'.
#' @param split.by A character string specifying the metadata column to split the plot by. Default is NULL.
#' @param label A logical value indicating whether to label the clusters. Default is TRUE.
#' @param repel A logical value indicating whether to use label repelling. Default is TRUE.
#' @param pt.size A numeric value specifying the size of points. Default is 0.8.
#' @param cols A vector of colors to use for the plot. Default is a predefined set of colors.
#' @return A ggplot2 object representing the dimensionality reduction plot (UMAP/TSNE).
#' @export
#' @import Seurat
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Assuming `seurat_obj` is a pre-existing Seurat object
#' scVisDimPlot(seurat_obj, group.by = "celltype")
#' }

scVisDimPlot <- function(scRNA,
                         reduction = "umap",
                         group.by = 'celltype',
                         split.by = NULL,
                         label = TRUE,
                         repel = TRUE,
                         pt.size = 0.8,
                         cols = c("#F1788D", "#54990F", "#E6550D", "#843C39", "#3182BD", "#8C6D31",
                                           "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", "#6BAED6", "#9ECAE1",
                                           "#AD494A", "#A1D99B", "#C7E9C0", "#99600F", "#C3BC3F", "#D6616B",
                                           "#FF7F00", "#1B9E77", "#FDAE6B", "#66A61E", "#E6550D", "#E7969C",
                                           '#53A85F', '#EFD4B4', '#9BCFB9', '#F0C0BD', '#FB6264', '#4E9C99',
                                           '#70C055', '#E98A27', '#FEBC28', "#F1788D", "#54990F", "#E6550D",
                                           "#843C39", "#66A61E", "#E6550D", "#E7969C")) {
                                           # Ensure necessary libraries are loaded
  library(Seurat)
  library(ggplot2)

  # Plot UMAP/TSNE
  plot <- Seurat::DimPlot(scRNA,
                          reduction = reduction,
                          group.by = group.by,
                          split.by = split.by,
                          label = label,
                          repel = repel,
                          pt.size = pt.size,
                          cols = cols)

  return(plot)
}

# Example usage
# Assuming `seurat_obj` is a pre-existing Seurat object
# scVisDimPlot(seurat_obj, group.by = "celltype")





# ============ 3. scVisTissueOR ==========
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data
  library(data.table)
  dir.create(dirname(out.prefix),F,T)

  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)

  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}

  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)

  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }

  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


# 加载必要的包
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

# 定义主分析函数
analyze_tissue_dist <- function(meta_data, output_prefix, pdf_width = 8, pdf_height = 4, verbose = 1) {

  # 调用 do.tissueDist 函数进行主要分析
  OR_immune_list <- do.tissueDist(cellInfo.tb = meta_data,
                                  out.prefix = sprintf("%s.Immune_cell", output_prefix),
                                  pdf.width = pdf_width, pdf.height = pdf_height, verbose = verbose)

  # 返回分析结果
  return(OR_immune_list)
}

# 定义绘图函数
plot_heatmap <- function(OR_list) {
  # 提取 OR 值结果
  a <- OR_list[["OR.dist.tb"]] %>%
    as.data.frame() %>%
    column_to_rownames(var = "rid") %>%
    na.omit()

  # 提取 P 值结果
  b <- OR_list$count.dist.melt.ext.tb[, c(1, 2, 6)] %>%
    spread(key = "cid", value = "adj.p.value") %>%
    column_to_rownames(var = "rid")

  # 只选择在a中的行
  b <- b[rownames(a),]

  # 调整 P 值符号表示
  col <- viridis(11, option = "D")
  b <- ifelse(b >= 0.05 & (a > 1.5 | a < 0.5), "",
              ifelse(b < 0.0001 & (a > 1.5 | a < 0.5), "****",
                     ifelse(b < 0.001 & (a > 1.5 | a < 0.5), "***",
                            ifelse(b < 0.01 & (a > 1.5 | a < 0.5), "**",
                                   ifelse(b < 0.05 & (a > 1.5 | a < 0.5), "*", "")))))

  bk <- c(seq(0, 0.99, by = 0.01), seq(1, 2, by = 0.01))

  # 绘制热图
  pheatmap(a, border_color = NA, fontsize = 9, cellheight = 12, cellwidth = 20,
           clustering_distance_rows = "correlation", display_numbers = b,
           number_color = "black", fontsize_number = 10, cluster_col = FALSE,
           cluster_rows = TRUE, breaks = bk, treeheight_row = 20, treeheight_col = 20,
           color = c(colorRampPalette(colors = col[1:6])(length(bk) / 2),
                     colorRampPalette(colors = col[6:11])(length(bk) / 2)))
}






#' @title Visualize Tissue Odds Ratio (OR) Analysis
#' @description This function performs tissue odds ratio analysis and plots the results.
#' @param scRNA A Seurat object containing scRNA-seq data.
#' @param group A character string specifying the grouping variable in the metadata. Default is 'orig.ident'.
#' @param celltype A character string specifying the cell type variable in the metadata. Default is 'celltype'.
#' @param output_prefix A character string specifying the prefix for output files. Default is "./output_figure".
#' @param output_file A character string specifying the output file name. Default is "/tissue_OR.pdf".
#' @param width Numeric value specifying the width of the output plot. Default is 5.
#' @param height Numeric value specifying the height of the output plot. Default is 4.
#' @return A list with tissue OR analysis results and a heatmap plot.
#' @export
#' @import Seurat
#' @import ggplot2
#' @import export
#' @examples
#' \dontrun{
#' # Assuming `sce` is a pre-existing Seurat object
#' scVisTissueOR(sce, group = "group", celltype = "celltype", output_prefix = "./output_figure")
#' }

scVisTissueOR <- function(scRNA, group = 'orig.ident', celltype = 'celltype', output_prefix = "./output_figure", output_file = "/tissue_OR.pdf", width = 5, height = 4) {

  # Extract metadata
  meta <- scRNA@meta.data

  # Assign group and celltype to new columns
  meta$loc <- meta[, group]
  meta$meta.cluster <- meta[, celltype]

  # Perform tissue OR analysis
  OR_immune_list <- analyze_tissue_dist(meta_data = meta, output_prefix = output_prefix)

  # Plot heatmap
  p <- plot_heatmap(OR_immune_list)

  # Ensure the output directory exists
  dir.create(output_prefix, showWarnings = FALSE, recursive = TRUE)

  # Export the plot to a PDF file
  export::graph2pdf(p, file = paste0(output_prefix, output_file), width = width, height = height)

  return(list(OR_immune_list = OR_immune_list, heatmap_plot = p))
}

# Example usage
# Assuming `sce` is a pre-existing Seurat object
# scVisTissueOR(sce, group = "group", celltype = "celltype", output_prefix = "./output_figure", output_file = "/tissue_OR.pdf", width = 5, height = 4)




# ============ 4. CPDB Visualization =============


#' @title CPDB Visualization
#' @description This function visualizes the results from CellPhoneDB analysis.
#' @param sce A Seurat object containing the single-cell RNA data.
#' @param pvals_path A string specifying the path to the p-values file.
#' @param means_path A string specifying the path to the means file.
#' @param decon_path A string specifying the path to the deconvolution file.
#' @param celltype_key A string specifying the key for cell type in the metadata. Default is "celltype".
#' @param sender A character vector specifying the sender cell types.
#' @param receiver A character vector specifying the receiver cell types.
#' @param output_dir A string specifying the output directory. Default is "./output/".
#' @param output_filename A string specifying the output filename. Default is "cpdb.pdf".
#' @param palette1 A character vector specifying the color palette for dot plot. Default is c("darkblue", "yellow", "red").
#' @param palette2 A character vector specifying the color palette for heatmap. Default is c("navy", "white", "firebrick3").
#' @return None. The function saves the visualizations in a PDF file.
#' @export
#' @import Seurat
#' @import pheatmap
#' @import ktplots
#' @import ggplot2
#' @examples
#' \dontrun{
#' cpdb_visualization(
#'   sce = seurat_obj,
#'   pvals_path = "path/to/pvals.txt",
#'   means_path = "path/to/means.txt",
#'   decon_path = "path/to/decon.txt",
#'   sender = c("CellTypeA"),
#'   receiver = c("CellTypeB"),
#'   output_dir = "./output_data/",
#'   output_filename = "cpdb.pdf"
#' )
#' }

cpdb_visualization <- function(sce, pvals_path, means_path, decon_path,
                               celltype_key = "celltype",
                               sender, receiver,
                               output_dir = "./output_data/",
                               output_filename = "cpdb.pdf",
                               palette1 = c("darkblue", "yellow", "red"),
                               palette2 = c("navy", "white", "firebrick3")) {

  # Ensure necessary libraries are loaded
  library(pheatmap)
  library(ktplots)
  library(Seurat)
  library(ggplot2)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set up the PDF output
  pdf(file.path(output_dir, output_filename))

  # Define color palettes
  my_palette1 <- colorRampPalette(palette1)(n = 100)
  my_palette2 <- colorRampPalette(palette2)(100)

  # Load data
  pvals <- read.delim(pvals_path, check.names = FALSE)
  means <- read.delim(means_path, check.names = FALSE)
  decon <- read.delim(decon_path, check.names = FALSE)

  # Prepare cell types for plotting
  cell_type1 <- paste(sender, collapse = '|')
  cell_type2 <- paste(receiver, collapse = '|')

  # Generate interaction table
  table <- plot_cpdb(scdata = sce, keep_significant_only = TRUE,
                     return_table = TRUE,
                     cell_type1 = cell_type1,
                     cell_type2 = cell_type2,
                     celltype_key = celltype_key,
                     means = means, pvals = pvals)

  # Prepare data for dot plot
  plot_data1 <- table[c(2:5)]
  colnames(plot_data1)[1:4] <- c('interacting_pair', 'cell_pairs', 'means', 'pvals')
  plot_data1$interacting_pair <- str_split(plot_data1$interacting_pair, ">@<>@<>@<", simplify = TRUE)[, 2]

  # Data preprocessing
  plot_data1$means[is.na(plot_data1$means)] <- 0
  plot_data1$pvals[is.na(plot_data1$pvals)] <- 1

  # Dot plot
  plot1 <- ggplot(plot_data1, aes(x = cell_pairs, y = interacting_pair)) +
    geom_point(aes(size = -log10(pvals), color = means)) +
    scale_size_continuous(range = c(0, 3), breaks = c(0, 1.0, 2.0)) +
    scale_color_gradientn('Mean expression', colors = my_palette1, limits = c(0, 1)) +
    theme_bw() +
    theme(axis.text = element_text(size = 8, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 7),
          axis.text.y = element_text(size = 7, colour = "black"),
          axis.title = element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  print(plot1)

  # Heatmap
  plot_data2 <- plot_cpdb_heatmap(pvals = pvals, return_tables = TRUE,
                                  cellheight = 20, cellwidth = 20,
                                  cluster_cols = FALSE, cluster_rows = FALSE,
                                  cell_types = unique(sce@meta.data[,celltype_key]))$count_network
  pheatmap(plot_data2, color = my_palette2,
           cluster_rows = FALSE, cluster_cols = FALSE,
           border_color = "white",
           cellwidth = 30, cellheight = 30)

  # Chord diagram
  plot_cpdb3(scdata = as.SingleCellExperiment(sce),
             cell_type1 = cell_type1, cell_type2 = cell_type2,
             celltype_key = celltype_key, means = means, pvals = pvals,
             keep_significant_only = TRUE, deconvoluted = decon)

  dev.off()
}
