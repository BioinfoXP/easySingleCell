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
#' @import patchwork
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

  if (length(feature)>=2) {
    plot <- lapply(feature, function(x){
      FeaturePlot(
        object = scRNA,
        features = x,
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
        ggtitle(x)
    })

    return(patchwork::wrap_plots(plot))
  } else {
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

# =============== 2.scVisFeaturePlot3  ===============
#' @title Simultaneous Visualization of Three Features in a Single Plot
#' @description This function visualizes three distinct features on a single dimension reduction plot using a color blending system. It allows for the quantitative display of gene expressions or other continuous variables by mixing colors according to the RYB or RGB color models, providing a unique perspective on feature interactions and expression levels within individual cells.
#' @param seu A Seurat object that contains the data for plotting. This object should have precomputed dimensionality reduction coordinates.
#' @param feature.1 The name of the first feature (gene or other variable) to be plotted. Default: NA.
#' @param feature.2 The name of the second feature. Default: NA.
#' @param feature.3 The name of the third feature. Default: NA.
#' @param color The color model used to blend the expression data of the three features. Options include "ryb" (red-yellow-blue) and "rgb" (red-green-blue), affecting how expression intensities are represented through color. Default: c("ryb", "rgb").
#' @param color.range The range of expression intensity that is represented by the color spectrum in the plot, helping to enhance visibility of lower expressions and prevent oversaturation at high expression levels. Default: c(0.1, 0.9).
#' @param reduction The type of dimension reduction used to display the data, such as 'umap' or 'tsne'. This choice determines the underlying plot layout. Default: 'umap'.
#' @param order A logical value indicating whether to plot cells with higher expressions on top of those with lower expressions, which can help prevent significant data points from being obscured in dense areas of the plot. Default: TRUE.
#' @param pt.size Point size for plotting individual cells in the grid. Smaller values are typically used for large datasets or dense plots, whereas larger values enhance visibility for plots with fewer cells or less overlap. Default: 0.1.
#' @return A ggplot object that represents a dimension reduction plot incorporating three features with color blending, showing how each feature contributes to the overall expression patterns observed.
#' @details `scVisFeaturePlot3` is designed for detailed exploratory analysis where understanding the interplay between multiple variables is crucial. This function is particularly useful for researchers looking to explore gene expressions in complex datasets, such as those involving interactions between different cell types or conditions.
#' @examples
#' library(Seurat)
#' library(easySingleCell)
#'
#' scVisFeaturePlot3(
#'   pbmc,
#'   feature.1 = "CD3D",
#'   feature.2 = "CD14",
#'   feature.3 = "CD79A",
#'   color = "ryb"
#' )
#'
#' @rdname scVisFeaturePlot3
#' @export

scVisFeaturePlot3 <- function(
    seu,
    feature.1 = NA,
    feature.2 = NA,
    feature.3 = NA,
    color = c("ryb","rgb"),
    color.range = c(0.1,0.9),
    reduction = "umap",
    order = T,
    pt.size = 0.1
){
  library(Seurat)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  features <- c(feature.1, feature.2, feature.3)
  color <- color[1]
  if(all(is.na(features))) stop("No feature to plot.")
  tp <- matrix(rep(color.range[1], times = (ncol(seu) * 3)), ncol = 3) %>%
    as.data.frame()
  l <- color.range[1]
  h <- color.range[2]
  colors <- switch (color,
                    "ryb" = c(
                      ryb2rgb(c(h,l,l)),
                      ryb2rgb(c(l,h,l)),
                      ryb2rgb(c(l,l,h)),
                      ryb2rgb(c(l,l,l))
                    ),
                    "rgb" = c(
                      rgb(h,l,l),
                      rgb(l,h,l),
                      rgb(l,l,h),
                      rgb(l,l,l)
                    )
  )
  lgd <- function(title, col) {
    value = tp[,title]
    df <- data.frame(
      value = c(min(value),max(value))
    )
    if(max(value) == min(value)) col[2] <- col[1]
    p.tmp <-
      ggplot(df, aes(x = 1, y = 1, color = value)) +
      geom_point() +
      scale_color_gradient(low = col[1], high = col[2]) +
      labs(color = title) +
      theme(legend.justification = c(0,1))
    p.tmp <- get_legend(p.tmp)
    return(p.tmp)
  }
  p.leg <- list()
  for (i in 1:3) {
    if(!is.na(features[i])) {
      tp[,i] <- FetchData(seu, features[i])
      colnames(tp)[i] <- features[i]
      p.leg[[features[i]]] <- lgd(features[i], col = c(colors[4], colors[i]))
    }
  }
  tp.c <- as.data.frame(apply(tp, 2, function(x){
    if(max(x) == min(x)){
      y <- rep(color.range[1], length(x))
    } else {
      y <- sca(x, color.range)
    }
  }))
  tp.c$color <- apply(tp.c, 1, function(x){
    col <- switch(
      color,
      "ryb" = ryb2rgb(x),
      "rgb" = rgb(x[1],x[2],x[3]))
    return(col)
  })
  tp.c <- cbind(Embeddings(seu, reduction = reduction), tp.c)
  if(order) tp.c <- tp.c[order(rowMeans(tp.c[3:5])), ]
  p <-
    ggplot(tp.c, aes_string(
      x = colnames(tp.c)[1],
      y = colnames(tp.c)[2])) +
    geom_point(color = tp.c$color, size = pt.size) +
    theme_classic()
  p <- ggarrange(
    p, ggarrange(plotlist = p.leg, ncol = 1),
    widths = c(8,1)
  )
  return(p)
}

#' @title Simultaneous Visualization of Three Features on a Grid of Plots
#' @description This function allows for the simultaneous visualization of three features across multiple plots, utilizing a grid layout. It supports two color blending systems (RYB or RGB) to represent the intensity of each feature within a single plot, providing a comprehensive overview of expression patterns across a dataset.
#' @param seu A Seurat object containing single-cell RNA sequencing data.
#' @param features A vector of feature names (genes or other continuous variables) to be displayed. This vector should be divisible by three, with each triplet of features being plotted in a separate subplot within the grid. If a triplet includes `NA`, that position will not display a feature, allowing for flexibility in visualization.
#' @param color Specifies the color blending system used to display the features. The available options are "ryb" for red-yellow-blue and "rgb" for red-green-blue. This parameter controls how the three features are visually represented based on their expression levels.
#' @param color.range Adjusts the luminance range used for feature visualization to ensure that low expressions are visible and not obscured by the background color. Default: c(0.1, 0.95), where 0.1 prevents the lowest expressions from being pure white and 0.95 keeps the highest expressions from saturating to pure color.
#' @param reduction The dimensionality reduction technique to use for the plot layout. Typically 'umap' or 'tsne' are used, with 'umap' being the default.
#' @param order Controls whether cells with higher feature expressions are plotted above those with lower expressions. This is useful for ensuring that cells with significant expression levels are visible and not obscured by those with lower levels. Default: TRUE.
#' @param pt.size Point size for plotting individual cells in the grid. Smaller values are typically used for large datasets or dense plots, whereas larger values enhance visibility for plots with fewer cells or less overlap. Default: 0.1.
#' @param legend Determines whether to display a legend describing the features and color scales. Default: FALSE.
#' @return A ggplot object displaying a grid of dimension reduction plots, each illustrating the expression patterns of three features using the specified color blending system.
#' @details The `scVisFeaturePlot3.grid` function is particularly useful for exploratory data analysis where visualization of multiple gene interactions or expression patterns across different cell populations is required. It effectively combines data from multiple features into a single coherent visual representation.
#' @examples
#' library(Seurat)
#' library(easySingleCell)
#'
#' scVisFeaturePlot3.grid(
#'   pbmc,
#'   features = c("CD3D", "CD14", "CD79A", "FCGR3A", NA, "LYZ"),
#'   color = "ryb",
#'   pt.size = 0.5
#' )
#'
#' @rdname scVisFeaturePlot3-grid
#' @export

scVisFeaturePlot3.grid <- function(
    seu,
    features,
    color = c("ryb","rgb"),
    color.range = c(0.1,0.95),
    reduction = "umap",
    order = T,
    pt.size = 0.1,
    legend = F
){
  library(Seurat)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  library(magrittr)
  library(rlist)
  import("ggtext")

  # features
  add3 <- function(x) {
    y <- switch(length(x) %% 3 + 1, 0, 2, 1)
    x <- c(x, rep(NA, times = y))
    return(x)
  }
  trim3 <- function(x) {
    if(length(x) > 3) {
      warning("Vector length > 3. Only use the first 3 elements")
      x <- x[1:3]
    }
    if(length(x) < 3) x <- add3(x)
    return(x)
  }
  if(is.data.frame(features)){
    ft <- features[1:3,] %>% set_rownames(NULL)
  }else if(is.list(features)){
    ft <- lapply(features, trim3)
    ft <- as.data.frame(ft)
  }else if(is.vector(features)){
    ft <- as.data.frame(matrix(add3(features), nrow = 3))
  }else stop("'features' should be class 'vector', 'data.frame' or 'list'.")

  # colors
  color <- color[1]
  l <- color.range[1]
  h <- color.range[2]
  colors <- switch (
    color,
    "ryb" = c(
      ryb2rgb(c(h,l,l)),
      ryb2rgb(c(l,h,l)),
      ryb2rgb(c(l,l,h)),
      ryb2rgb(c(l,l,l))
    ),
    "rgb" = c(
      rgb(h,l,l),
      rgb(l,h,l),
      rgb(l,l,h),
      rgb(l,l,l)
    )
  )

  # data for plot
  tp <- list()
  tp.c <- list()
  for (i in 1:ncol(ft)) {
    tp[[i]] <-
      matrix(rep(color.range[1], times = (ncol(seu) * 3)), ncol = 3) %>%
      as.data.frame()
    for (j in 1:3) {
      f <- ft[j,i]
      if(!is.na(f)) {
        tp[[i]][,j] <- FetchData(seu, f)
        colnames(tp[[i]])[j] <- f
      }
    }
    tp.c[[i]] <- as.data.frame(
      apply(tp[[i]], 2, function(x){
        if(max(x) == min(x)){
          y <- rep(color.range[1], length(x))
        } else {
          y <- sca(x, color.range)
        }
      })
    )
    tp.c[[i]]$color <- apply(
      tp.c[[i]], 1,
      function(x){
        col <- switch(
          color,
          "ryb" = ryb2rgb(x),
          "rgb" = rgb(x[1],x[2],x[3]))
        return(col)
      }
    )
    tp.c[[i]] <- cbind(Embeddings(seu, reduction = reduction), tp.c[[i]])
    if(order) tp.c[[i]] <- tp.c[[i]][order(rowMeans(tp.c[[i]][,3:5])), ]
    tp.c[[i]]$index <- i
    tp.c[[i]] <- tp.c[[i]][,-c(3:5)]
  }
  tp.c <- list.rbind(tp.c)

  # strip title
  ft.col.str <- function(f, col){
    if(is.na(f)){
      s <- NULL
    }else{
      s <- paste0(
        "<span style='color:",
        col, ";'>",
        f, "</span>")
    }
    return(s)
  }
  lab.str <-
    apply(ft, 2, function(x){
      tmp <- c(
        ft.col.str(x[1], colors[1]),
        ft.col.str(x[2], colors[2]),
        ft.col.str(x[3], colors[3])
      )
      tmp <- paste(tmp, collapse = " ")
      return(tmp)
    })
  names(lab.str) <- 1:ncol(ft)

  # main plot
  p <-
    ggplot(tp.c, aes_string(
      x = colnames(tp.c)[1],
      y = colnames(tp.c)[2])) +
    facet_wrap(~index, labeller = labeller(index = lab.str)) +
    geom_point(color = tp.c$color, size = pt.size) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_markdown(face = "bold", size = 11),
          strip.background = element_rect(fill = NA))

  # figure legend
  if(legend){
    p.leg <- list()
    for (i in 1:3) {
      if(!all(is.na(ft[i,]))){
        p.tmp <-
          ggplot(data.frame(value = 0:1), aes(x = 1, y = 1, color = value)) +
          geom_point() +
          scale_color_gradient(
            low = colors[4], high = colors[i],
            breaks = c(0,1), labels = c("min", "max")) +
          labs(color = NULL) +
          theme(legend.justification = c(0,1))
        p.leg[[as.character(i)]] <- get_legend(p.tmp)
      }
    }
    p <- ggarrange(
      p, ggarrange(plotlist = p.leg, ncol = 1),
      widths = c(8,1)
    )
  }
  return(p)
}
















# ============ 3. scVisDimPlot ==========
#' @title Create an Enhanced Dimensional Reduction Plot
#' @description This function creates a dimension reduction plot that can handle both discrete and continuous variables seamlessly. It incorporates additional customization options for visual representation and automatically recognizes input variable types to optimize visualization.
#' @param seu Seurat object containing single-cell data for visualization.
#' @param features Variables to be visualized in the plot, accepting both discrete and continuous variables. Default: NULL.
#' @param group.by Alias for `features`. Default: NULL.
#' @param split.by A metadata column name by which to split the plot, creating separate plots for each unique value. This can be useful for visualizing differences across conditions or experiments.
#'   Default: NULL.
#' @param cells A vector specifying a subset of cells to include in the plot.
#'   Default: all cells are included.
#' @param slot Which data slot to use for pulling expression data. Accepts 'data', 'scale.data', or 'counts'. Default: 'data'.
#' @param assay Specify the assay from which to retrieve data.
#'   Default: NULL, which will use the default assay.
#' @param dims A two-length numeric vector specifying which dimensions to use for the x and y axes, typically from a PCA, tSNE, or UMAP reduction.
#'   Default: c(1, 2).
#' @param reduction Which dimensionality reduction to use. If not specified, will search in order of 'umap', 'tsne', then 'pca'.
#'   Default: NULL.
#' @param priority Specifies which to prioritize when metadata column names conflict with gene names: 'expr' for expression, 'none' for metadata.
#'   Default: c("expr", "none").
#' @param ncol Number of columns to display when combining multiple plots into a single patchworked ggplot object.
#'   Default: NULL.
#' @param nrow Number of rows to display when combining multiple plots into a single patchworked ggplot object.
#'   Default: NULL.
#' @param nrow.each Specifies the number of rows each split plot should have when using the split.by parameter.
#'   Default: NULL.
#' @param ncol.legend Integer specifying the number of columns in the plot legend of discrete variables. Default: NULL.
#' @param cols Flexible color settings for the plot, accepting a variety of inputs:
#'
#'   - A vector specifying a global color setting similar to Seurat's `DimPlot`/`FeaturePlot`.
#'
#'   - A list specifying colors for each variable type (discrete/continuous) or for each individual variable. For example, `list(discrete = "auto", continuous = "A")` applies automatic styling from `color_pro()` for discrete variables and `viridis` "A" style for continuous variables. More detailed setups can include `list("cluster" = "pro_blue", "CD14" = c("#EEEEEE", "black"))`.
#'
#'   For continuous variables:
#'
#'     - Predefined color schemes from the `viridis` package ("A", "B", "C", "D", "E").
#'
#'     - Named vector with keys "low", "mid", and "high" for three-point gradients. Example: `c(low = "blue", mid = "white", high = "red")`.
#'
#'     - Two-point gradient with keys "low" and "high". Example: `c(low = "blue", high = "red")`.
#'
#'     - Custom color gradient using a vector of colors.
#'
#'   For discrete variables:
#'
#'     - Seven color_pro styles: "default", "light", "pro_red", "pro_yellow", "pro_green", "pro_blue", "pro_purple".
#'
#'     - Five color_iwh styles: "iwh_default", "iwh_intense", "iwh_pastel", "iwh_all", "iwh_all_hard".
#'
#'     - Brewer color scales as specified by `brewer.pal.info`.
#'
#'     - Any manually specified colors.
#'
#'   Default: list(discrete = "auto", continuous = "A").
#' @param load.cols When TRUE, automatically loads pre-stored color information for variables from `seu@misc[["var_colors"]]`.
#'   Default: TRUE.#' @param pt.size Point size for plotting, adjusts the size of each cell in the plot.
#'   Default: NULL.
#' @param shape.by Metadata column or expression data used to specify different shapes for cells in the plot, allowing for additional visual distinctions.
#'   Default: NULL.
#' @param alpha.by Transparency of points in the plot, which can be helpful in densely plotted areas.
#'   Default: NULL.
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if points representing higher expression levels are being buried.
#'   Default: c(discrete = FALSE, continuous = TRUE).
#' @param shuffle Randomly shuffles the order of points to prevent burial under other points.
#'   Default: c(discrete = TRUE, continuous = FALSE).
#' @param label Whether to label clusters or other features in the plot.
#'   Default: FALSE.
#' @param label.color Customize the color of labels; defaults to the same as cluster colors unless specified, such as "black".
#'   Default: NULL.
#' @param box Whether to draw a box around labels to enhance visibility.
#'   Default: FALSE.
#' @param index.title Specify a prefix for cluster indices when labels are replaced by numerical indices to simplify the plot.
#'   Default: NULL.
#' @param repel Whether to use a repelling algorithm to avoid overlapping text labels.
#'   Default: FALSE.
#' @param label.size Size of the text labels used for clusters or features.
#'   Default: 4.
#' @param theme Allows customization of ggplot themes, for example, to remove axes or adjust text.
#'   Default: NULL.
#' @param cells.highlight A vector of cell names to highlight; simpler input than Seurat's approach, focusing on ease of use.
#'   Default: NULL.
#' @param cols.highlight A color or vector of colors to use for highlighting specified cells; will repeat to match the number of groups in cells.highlight.
#'   Default: '#DE2D26'.
#' @param sizes.highlight Size of highlighted points, providing emphasis where needed.
#'   Default: 1.
#' @param na.value Color to use for NA points when using a custom color scale.
#'   Default: 'grey80'.
#' @param raster Whether to convert the plot points to a raster format, which can help with performance on large datasets.
#'   Default: NULL.
#' @param raster.dpi The resolution for rasterized plots, useful for maintaining detail in dense plots.
#'   Default: NULL.
#' @param combine Whether to combine multiple plots into a single ggplot object using patchwork.
#'   Default: TRUE.
#' @param align Specifies how plots should be aligned if combined, accepting 'h' for horizontal, 'v' for vertical, or 'hv' for both.
#'   Default: 'hv'.
#' @return A ggplot object if `combine` is TRUE; otherwise, a list of ggplot objects, allowing for flexible plot arrangements or combined visualizations.
#' @details `scVisDimPlot` extends the functionality of Seurat's visualization tools by combining the features of `DimPlot` and `FeaturePlot` into a single, more versatile function. It automatically recognizes whether the input features are discrete or continuous, adjusting the visualization accordingly. This makes `scVisDimPlot` ideal for exploring complex scRNA-seq data without the need to switch between different plotting functions based on variable types. The function also offers advanced customization options for colors, themes, and labeling, making it highly adaptable to various data visualization needs.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Create a basic dimensional reduction plot with default settings
#' scVisDimPlot(pbmc)
#'
#' # Visualize different variables, including both discrete and continuous types
#' scVisDimPlot(pbmc, features = c("cluster", "orig.ident", "CD14", "CD3D"))
#'
#' # Split the visualization by a specific variable for comparative analysis
#' scVisDimPlot(pbmc, features = c("cluster", "CD14"), split.by = "orig.ident", ncol = 1)
#'
#' # Highlight specific cells, such as a particular cluster
#' b_cells <- colnames(pbmc)[pbmc$cluster == "B cell"]
#' scVisDimPlot(pbmc, cells.highlight = b_cells)
#'
#' # Apply advanced customization for colors and themes
#' scVisDimPlot(
#'   pbmc,
#'   features = c("cluster", "orig.ident", "CD14", "CD3D"),
#'   cols = list(
#'     "cluster" = "pro_blue",
#'     "CD14" = "D",
#'     "CD3D" = c("#EEEEEE", "black")
#'   ),
#'   theme = NoAxes())
#'
#' # Enhance the plot with labels and bounding boxes
#' scVisDimPlot(pbmc, label = TRUE, box = TRUE, label.color = "black", repel = TRUE, theme = NoLegend())
#'
#' # Use indices instead of long cluster names to simplify labels in the plot
#' scVisDimPlot(pbmc, index.title = "C", box = TRUE, label.color = "black")
#' @rdname scVisDimPlot
#' @export

scVisDimPlot <- function(
    seu,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    dims = c(1, 2),
    reduction = NULL,
    priority = c("expr", "none"),
    ncol = NULL,
    nrow = NULL,
    nrow.each = NULL,
    ncol.legend = NULL,
    cols = list(discrete = "auto", continuous = "A"),
    load.cols = TRUE,
    pt.size = NULL,
    shape.by = NULL,
    alpha.by = NULL,
    order = c(discrete = FALSE, continuous = TRUE),
    shuffle = c(discrete = TRUE, continuous = FALSE),
    label = FALSE,
    label.color = NULL,
    box = FALSE,
    index.title = NULL,
    repel = FALSE,
    label.size = 4,
    theme = NULL,
    cells.highlight = NULL,
    cols.highlight = "#DE2D26",
    sizes.highlight = 1,
    na.value = "grey80",
    raster = NULL,
    raster.dpi = NULL,
    combine = TRUE,
    align = "hv"
) {
  plot_data <- scVisDimPlot_GetData(
    seu = seu,
    features = features,
    group.by = group.by,
    split.by = split.by,
    cells = cells,
    slot = slot,
    assay = assay,
    dims = dims,
    reduction = reduction,
    shape.by = shape.by,
    alpha.by = alpha.by,
    cells.highlight = cells.highlight
  )

  data.dim <- plot_data$data.dim
  data.var <- plot_data$data.var
  vars <- colnames(data.var)

  get_disc_cont_par <- function(par, type, default) {
    if (is.vector(par) && is.logical(par) && type %in% names(par)) {
      par.d <- par[type]
    } else if (is.logical(par) && length(par) == 1) {
      par.d <- par
    } else {
      par.d <- default
    }
    return(par.d)
  }

  library(ggplot2)
  p <- list()
  for (i in vars) {
    data.single <- data.dim
    data.single$var <- data.var[[i]]
    if(is_continuous(data.single$var)) {
      scale_color <- scVisDimPlot_SelColCont(
        seu = seu,
        var = i,
        cols = cols,
        load.cols = load.cols
      )
      order.c <- get_disc_cont_par(order, "continuous", TRUE)
      shuffle.c <- get_disc_cont_par(shuffle, "continuous", FALSE)

      p[[i]] <- scVisDimPlot_PlotSingle(
        data.plot = data.single,
        title = i,
        cols = scale_color,
        pt.size = pt.size,
        order = order.c,
        shuffle = shuffle.c,
        label = FALSE,
        label.color = NULL,
        repel = FALSE,
        box = FALSE,
        label.size = 4,
        cols.highlight = "#DE2D26",
        sizes.highlight = 1,
        na.value = na.value,
        nrow.each = nrow.each,
        theme = theme,
        raster = raster,
        raster.dpi = raster.dpi)
    } else {
      if(i != "Selected_cells") data.single$var <- factor(data.single$var)
      n <- nlevels(data.single$var)
      if(!is.null(index.title)) {
        data.single$var_orig <- factor(data.single$var)
        index <- paste0(index.title, seq(nlevels(data.single$var_orig)))
        index_cluster <- paste(index, levels(data.single$var_orig))
        data.single$var <- factor(index[data.single$var_orig], levels = index)
        data.single$index_cluster <- factor(index_cluster[data.single$var_orig], levels = index_cluster)
        label <- TRUE
        labels <- index_cluster
      } else {labels <- waiver()}
      if(n > 100) stop("Discrete variable '",i,"' has > 100 values. Unable to plot.")
      if(i != "Selected_cells") {
        scale_color <- scVisDimPlot_SelColDisc(
          seu = seu,
          n = n,
          var = i,
          cols = cols,
          load.cols = load.cols,
          label = label,
          labels = labels,
          box = box
        )
      } else scale_color <- NULL

      order.d <- get_disc_cont_par(order, "discrete", FALSE)
      shuffle.d <- get_disc_cont_par(shuffle, "discrete", TRUE)

      p[[i]] <- scVisDimPlot_PlotSingle(
        data.plot = data.single,
        title = i,
        cols = scale_color,
        pt.size = pt.size,
        order = order.d,
        shuffle = shuffle.d,
        label = label,
        label.color = label.color,
        repel = repel,
        box = box,
        label.size = label.size,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        nrow.each = nrow.each,
        theme = theme,
        raster = raster,
        raster.dpi = raster.dpi,
        ncol.legend = ncol.legend)
    }
  }
  if(!combine) {
    return(p)
  } else {
    import("cowplot")
    p <- plot_grid(plotlist = p, ncol = ncol, nrow = nrow, align = align)
    return(p)
  }
}

is_continuous <- function(vec) {
  return(is.numeric(vec) & !is.factor(vec))
}

scVisDimPlot_GetData <- function(
    seu,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    dims = c(1, 2),
    reduction = NULL,
    shape.by = NULL,
    alpha.by = NULL,
    cells.highlight = NULL
) {
  if(!require(SeuratObject)) library(Seurat)

  if(!is.null(assay)) DefaultAssay(seu) <- assay

  # cells
  cells <- cells %||% colnames(seu)
  if(is.logical(cells)) {
    if(length(cells) != ncol(seu)) {
      stop("Logical value of 'cells' should be the same length as cells in ",
           "Seurat object")
    } else {
      cells.l <- cells
      cells <- colnames(seu)[cells]
    }
  } else {
    if(all(!cells %in% colnames(seu))) {
      stop("'cells' not found in Seurat object")
    }else if(any(!cells %in% colnames(seu))) {
      cells.out <- setdiff(cells, colnames(seu))
      stop(length(cells.out), " cell(s) not found in Seurat object: '",
           cells.out[1], "'...")
    }
    cells.l <- colnames(seu) %in% cells
  }

  # reduction
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = seu)
  data.dim <- Embeddings(object = seu[[reduction]])[cells, dims]
  data.dim <- as.data.frame(x = data.dim)
  dims <- colnames(data.dim)

  # vars
  if(!is.null(cells.highlight)) {
    seu@meta.data[["Selected_cells"]] <- factor(
      ifelse(
        colnames(seu) %in% cells.highlight,
        "Selected", "Unselected"),
      levels = c("Unselected","Selected")
    )
  }
  if(is.null(features) & is.null(group.by) & is.null(cells.highlight)) {
    data.var <- data.frame(Idents = factor(Idents(seu)[cells]))
  } else {
    vars <- c(features, group.by)
    if(!is.null(cells.highlight)) vars = c(vars, "Selected_cells")
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      data.var <- FetchData(object = seu, vars = vars, cells = cells, layer = slot, clean = "none")
    } else {
      data.var <- FetchData(object = seu, vars = vars, cells = cells, slot = slot)
    }
  }

  # split.by, shape.by, alpha.by
  plot_vars <- list(
    split.by = split.by,
    shape.by = shape.by,
    alpha.by = alpha.by)
  for (i in names(plot_vars)) {
    if(!is.null(plot_vars[[i]])) {
      data.dim[[i]] <- scVisDimPlot_PlotVars(
        seu = seu,
        var = plot_vars[[i]],
        var_name = i,
        cells = cells,
        cells.l = cells.l)
    }
  }

  return(list(
    data.dim = data.dim,
    data.var = data.var
  ))
}

scVisDimPlot_PlotVars <- function(
    seu,
    var,
    var_name,
    cells,
    cells.l
) {
  if(is.null(var)) {
    f2 <- NULL
  }else if(length(var) == 1) {
    if(var %in% colnames(seu@meta.data)){
      f2 <- factor(seu[[var]][cells,])
      names(f2) <- cells
    }else{
      stop("Cannot find '", var, "' in meta.data")
    }
  }else if(length(var) == ncol(seu)) {
    f2 <- factor(var[cells.l])
    names(f2) <- cells
  }else if(length(var) == length(cells)) {
    f2 <- factor(var)
    names(f2) <- cells
  }else{
    stop("'",var_name,"' should be variable name in 'meta.data' or ",
         "string with the same length of cells")
  }
  return(f2)
}

scVisDimPlot_PlotSingle <- function (
    data.plot,
    title = NULL,
    cols = NULL,
    pt.size = NULL,
    order = NULL,
    shuffle = NULL,
    label = FALSE,
    label.color = NULL,
    repel = FALSE,
    box = FALSE,
    label.size = 4,
    cols.highlight = "#DE2D26",
    sizes.highlight = 1,
    na.value = "grey80",
    nrow.each = NULL,
    theme = NULL,
    raster = NULL,
    raster.dpi = NULL,
    ncol.legend = NULL)
{
  library(ggplot2)
  pt.size <- pt.size %||% scVisDimPlot_AutoPointSize(data = data.plot, raster = raster)
  dims <- colnames(data.plot)[1:2]
  if ((nrow(x = data.plot) > 1e+05) & !isFALSE(raster)) {
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data.plot) > 1e+05)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  } else raster.dpi <- c(512, 512)

  if (isTRUE(x = shuffle)) {
    data.plot <- data.plot[sample(x = 1:nrow(x = data.plot)), ]
  }
  istrue.cell.highlight <- (title == "Selected_cells" & identical(levels(data.plot$var), c("Unselected","Selected")))
  if (isTRUE(x = order | istrue.cell.highlight)) {
    data.plot <- data.plot[order(data.plot$var), ]
  }

  if(!"shape.by" %in% colnames(data.plot)) shape.by <- NULL
  if(!"alpha.by" %in% colnames(data.plot)) alpha.by <- NULL

  plot <- ggplot(data = data.plot, mapping = aes(x = .data[[dims[1]]], y = .data[[dims[2]]], color = var, shape = shape.by, alpha = alpha.by))
  if(istrue.cell.highlight) pt.size <- sizes.highlight
  plot <-
    if (isTRUE(x = raster)) {
      import("scattermore")
      plot + geom_scattermore(pointsize = pt.size, pixels = raster.dpi)
    } else {
      plot + geom_point(size = pt.size)
    }
  if(!is_continuous(data.plot$var)) {
    plot <- plot +
      guides(color = guide_legend(override.aes = list(size = 3), ncol = ncol.legend))
  }
  plot <- plot + labs(color = NULL, title = title)
  if (label & title != "Selected_cells") {
    plot <- scVisDimPlot_LabelClusters(
      plot = plot,
      id = "var",
      repel = repel,
      label.color = label.color,
      box = box,
      size = label.size)
  }
  if("split.by" %in% colnames(data.plot)) {
    plot <- plot + facet_wrap(vars(split.by), nrow = nrow.each)
  }
  import("cowplot")
  plot <- plot + theme_cowplot() + CenterTitle()
  if(istrue.cell.highlight) {
    plot <- plot +
      scale_color_manual(values = c("#C3C3C3", cols.highlight)) +
      theme
  } else {
    plot <- plot + cols + theme
  }
  return(plot)
}

scVisDimPlot_AutoPointSize <- function (data, raster = NULL) {
  return(ifelse(test = isTRUE(x = raster), yes = 1, no = min(1583/nrow(x = data), 1)))
}

scVisDimPlot_LabelClusters <- function(
    plot,
    id,
    clusters = NULL,
    labels = NULL,
    label.color = NULL,
    split.by = NULL,
    repel = TRUE,
    box = FALSE,
    geom = 'GeomPoint',
    position = "median",
    ...
) {
  if(repel) import("ggrepel")
  xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == 'GeomSpatial') {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] = max(data[, xynames["y"]]) - data[, xynames["y"]] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      data.medians$color <- data.use$color[1]
      return(data.medians)
    }
  )
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) == as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1,2)]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[, id]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel, no = geom_label)
    if(is.null(label.color)) {
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
        show.legend = FALSE,
        ...
      )
    } else {
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
        show.legend = FALSE,
        color = label.color,
        ...
      )
    }
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    plot <- plot + geom.use(
      data = labels.loc,
      mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
      show.legend = FALSE,
      color = "black",
      ...
    )
  }
  # restore old axis ranges
  if (geom == 'GeomSpatial') {
    plot <- suppressMessages(expr = plot + coord_fixed(xlim = xrange.save, ylim = yrange.save))
  }
  return(plot)
}

GetXYAesthetics <- function (plot, geom = "GeomPoint", plot.first = TRUE) {
  geoms <- sapply(X = plot$layers, FUN = function(layer) {
    return(class(x = layer$geom)[1])
  })
  if (geom == "GeomPoint" && "GeomScattermore" %in% geoms) {
    geom <- "GeomScattermore"
  }
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as_label(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
    y <- as_label(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
  } else {
    x <- as_label(x = plot$layers[[geoms]]$mapping$x %||%
                    plot$mapping$x)
    y <- as_label(x = plot$layers[[geoms]]$mapping$y %||%
                    plot$mapping$y)
  }
  return(list(x = x, y = y))
}

scVisDimPlot_SelColCont <- function(
    seu,
    var,
    cols,
    load.cols
) {
  list_l <- is.list(cols)
  cols_var_l <- var %in% names(cols)
  load_var <- seu@misc[["var_colors"]][[var]]
  cols_cont_l <- "continuous" %in% names(cols)
  load_cont <- seu@misc[["var_colors"]][["continuous"]]
  if(list_l & cols_var_l) {
    cols <- cols[[var]]
  } else if(load.cols & !is.null(load_var)) {
    cols <- load_var
  } else if(list_l & cols_cont_l) {
    cols <- cols[["continuous"]]
  } else if(load.cols & !is.null(load_cont)) {
    cols <- load_cont
  } else if (list_l) {
    cols <- "A"
  }
  scale_color <- scale_color_cont_auto(cols)
  return(scale_color)
}

scVisDimPlot_SelColDisc <- function(
    seu,
    var,
    n,
    cols,
    load.cols,
    label,
    labels = waiver(),
    box
) {
  list_l <- is.list(cols)
  cols_var_l <- var %in% names(cols)
  load_var <- seu@misc[["var_colors"]][[var]]
  cols_disc_l <- "discrete" %in% names(cols)
  load_disc <- seu@misc[["var_colors"]][["discrete"]]
  label_l <- label & !box
  if(list_l & cols_var_l) {
    cols <- cols[[var]]
  } else if(load.cols & !is.null(load_var)) {
    cols <- load_var
  } else if(list_l & cols_disc_l) {
    cols <- cols[["discrete"]]
  } else if(load.cols & !is.null(load_disc)) {
    cols <- load_disc
  } else if(list_l) {
    cols <- "pro_default"
  }
  if(is.null(cols)) return(NULL)
  if(cols[1] == "auto") cols <- ifelse(label_l, "pro_light","pro_default")
  scale_color <- scale_color_disc_auto(cols, n = n, labels = labels)
  return(scale_color)
}

CenterTitle <- function () {
  return(theme(plot.title = element_text(hjust = 0.5), validate = TRUE))
}

# =========== 4.scVisVlnPlot =============
#' @include generics.R
#'
NULL

#' @param seu A Seurat object. Only applicable for the Seurat method.
#' @param group.by A variable from `meta.data` for grouping or a character vector of equal length as the number of cells. Only applicable for the Seurat method.
#' @param split.by A variable from `meta.data` to bifurcate the violin plots. Only applicable for the Seurat method.
#' @param cells Cell identifiers for use. Defaults to all cells. Only applicable for the Seurat method.
#' @param slot Slot to retrieve feature data from. Only applicable for the Seurat method.
#' @param assay Name of the assay to employ. Defaults to the active assay. Only applicable for the Seurat method.
#' @param priority If set to "expr", extracts data from the expression matrix over `meta.data`. Only applicable for the Seurat method.
#' @param load.cols When TRUE, automatically loads pre-stored color information for variables from `seu@misc[["var_colors"]]`.
#' @rdname scVisVlnPlot
#' @export

scVisVlnPlot.Seurat <- function(
    seu,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    priority = c("expr","none"),
    cols = "auto",
    load.cols = TRUE,
    ncol = NULL,
    lab_fill = "group",
    scales = "free_y",
    violin = T,
    box = T,
    width = 0.9,
    pt = T,
    hide.outlier = F,
    pt.style = c("jitter","quasirandom"),
    pt.size = 0.2,
    pt.alpha = 1,
    strip.position = "top",
    stat.method = c("none", "wilcox.test", "t.test"),
    p.adjust.method = "holm",
    label = c("p.signif","p","p.adj","p.format"),
    comparisons = NULL,
    hide.ns = TRUE,
    step.increase = 0.12,
    tip.length = 0.03,
    ...
) {
  Std.matr <- Seu2Matr(
    seu = seu,
    features = features,
    group.by = group.by,
    split.by = split.by,
    cells = cells,
    slot = slot,
    assay = assay,
    priority = priority
  )

  cols <- scVisVlnPlot_SelColDisc(
    seu = seu,
    group.by = group.by,
    split.by = split.by,
    cols = cols,
    load.cols = load.cols
  )

  p <- scVisVlnPlot.default(
    matr = Std.matr$matr,
    f = Std.matr$f,
    f2 = Std.matr$f2,
    t = T,
    cols = cols,
    ncol = ncol,
    lab_fill = lab_fill,
    scales = scales,
    violin = violin,
    box = box,
    width = width,
    pt = pt,
    hide.outlier = hide.outlier,
    pt.style = pt.style,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    strip.position = strip.position,
    stat.method = stat.method,
    p.adjust.method = p.adjust.method,
    label = label,
    comparisons = comparisons,
    hide.ns = hide.ns,
    step.increase = step.increase,
    tip.length = tip.length,
    ...)
  return(p)
}

#' @param matr A matrix or data frame with rows as features and columns as cells.
#' @param f A factor or vector indicating the identity of each cell. Should match the column length of `matr`.
#' @param f2 A factor or vector akin to `f` for splitting the violin plots. Default: NULL.
#' @param features Features to depict, such as gene expression, metrics, PC scores, or any data obtainable via `FetchData()`. Default: NULL (all features in matrix).
#' @param t If the matrix has features in columns and cells in rows, transpose the matrix first. Default: FALSE.
#' @param cols Flexible color settings for the plot, accepting a variety of inputs:
#'
#'     - Seven color_pro styles: "default", "light", "pro_red", "pro_yellow", "pro_green", "pro_blue", "pro_purple".
#'
#'     - Five color_iwh styles: "iwh_default", "iwh_intense", "iwh_pastel", "iwh_all", "iwh_all_hard".
#'
#'     - Brewer color scales as specified by `brewer.pal.info`.
#'
#'     - Any manually specified colors.
#' @param ncol Specifies the number of columns for display if multiple plots are shown. Default: NULL.
#' @param lab_fill Label for the figure legend. Default: 'group'.
#' @param scales Scales parameter passed to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}. Default: 'free_y'.
#' @param violin Indicates whether to generate a violin plot. Default: TRUE.
#' @param box Indicates whether to depict a box plot. Default: TRUE.
#' @param width Width of the box plot. Default: 0.9.
#' @param pt Indicates if points should be plotted. Default: TRUE.
#' @param hide.outlier Conceals outlier points from the box plot. Default: FALSE.
#' @param pt.style Position adjustment. Default choices: "jitter", "quasirandom".
#' @param pt.size Point size setting. Default: 0.2.
#' @param pt.alpha Adjusts the transparency of points. Default: 1.
#' @param strip.position Positions the strip ("top" (default), "bottom", "left", or "right"). Only used when `f2 = NULL`.
#' @param stat.method Determines if pairwise statistics are added to the plot. Either "wilcox.test" or "t.test". Default: "none".
#' @param p.adjust.method Method for adjusting p-values, especially when conducting multiple pairwise tests or dealing with multiple grouping variables. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none". Note: Adjustments are independently conducted for each variable in formulas containing multiple variables. Default: 'holm'.
#' @param label Specifies label type. Options include "p.signif" (showing significance levels), "p.format" (formatted p value), or "p", "p.adj". Default: "p.signif".
#' @param comparisons List of length-2 vectors, each containing either names of two x-axis values or two integers pointing to groups of interest for comparison. Default: all groups.
#' @param hide.ns If TRUE, the 'ns' symbol is concealed when displaying significance levels. Default: TRUE.
#' @param step.increase Numeric vector denoting the increase in fraction of total height for each additional comparison, minimizing overlap. Default: 0.12.
#' @param tip.length Numeric vector indicating the fraction of total height the bar descends to specify the exact column. For a line display instead of a bracket, set to 0. Default: 0.03.
#' @param ... Further arguments passed to the \code{\link[ggpubr:stat_pvalue_manual]{ggpubr::stat_pvalue_manual()}}.
#' @rdname scVisVlnPlot
#' @export

scVisVlnPlot.default <- function(
    matr, f, f2 = NULL,
    features = NULL,
    t = F,
    cols = "pro_default",
    ncol = NULL,
    lab_fill = "group",
    scales = "free_y",
    violin = T,
    box = T,
    width = 0.9,
    pt = T,
    hide.outlier = F,
    pt.style = c("jitter","quasirandom"),
    pt.size = 0.2,
    pt.alpha = 1,
    strip.position = "top",
    stat.method = c("none", "wilcox.test", "t.test"),
    p.adjust.method = "holm",
    label = c("p.signif","p","p.adj","p.format"),
    comparisons = NULL,
    hide.ns = TRUE,
    step.increase = 0.12,
    tip.length = 0.03,
    ...
) {

  scores <- scVisVlnPlot_Calc(
    matr = matr,
    f = f,
    f2 = f2,
    features = features,
    t = t
  )

  p <- scVisVlnPlot_Plot(
    scores = scores,
    cols = cols,
    ncol = ncol,
    lab_fill = lab_fill,
    scales = scales,
    violin = violin,
    box = box,
    width = width,
    pt = pt,
    hide.outlier = hide.outlier,
    pt.style = pt.style,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    strip.position = strip.position
  )

  p <- scVisVlnPlot_Stat(
    p = p,
    stat.method = stat.method,
    p.adjust.method = p.adjust.method,
    label = label,
    comparisons = comparisons,
    hide.ns = hide.ns,
    step.increase = step.increase,
    tip.length = tip.length,
    ...
  )

  return(p)
}

#' @title StackedViolin
#' @description Alias of \code{\link[SeuratExtend:scVisVlnPlot]{scVisVlnPlot()}}
#' @seealso \code{\link[SeuratExtend:scVisVlnPlot]{scVisVlnPlot()}}
#' @rdname StackedViolin
#' @export

StackedViolin <- scVisVlnPlot.default


#' @rdname StackedViolin
#' @export

StackedViolin_v3 <- scVisVlnPlot.Seurat

# Internal ----------------------------------------------------------------

scVisVlnPlot_Calc <- function(
    matr,
    f,
    f2,
    features,
    t
) {
  library(reshape2)
  if(!t) matr <- t(matr)
  features <- features %||% colnames(matr)
  if(length(setdiff(features, colnames(matr))) > 0){
    message(paste0(setdiff(features, colnames(matr)), collapse = ", "), " not found")
    features <- intersect(features, colnames(matr))
  }
  f <- factor(f)
  f2 <- f2 %||% data.frame(row.names = rownames(matr))
  scores <- cbind(f, f2, as.data.frame(matr[,features,drop = F]))
  scores <- melt(scores, measure.vars = features, variable.name = "feature")
  return(scores)
}

scVisVlnPlot_Plot <- function(
    scores,
    cols,
    ncol,
    lab_fill,
    scales,
    violin,
    box,
    width,
    pt,
    hide.outlier,
    pt.style,
    pt.size,
    pt.alpha,
    strip.position
) {
  library(ggplot2)
  x <- ifelse(!"f2" %in% colnames(scores), "f", "f2")
  p <- ggplot(scores, aes(x = .data[[x]], y = value))
  n <- nlevels(factor(scores[[x]]))

  if(violin) {
    p <- p + geom_violin(mapping = aes(fill = .data[[x]]), scale = "width", width = width)
  }
  if(box & !violin) {
    if(pt | hide.outlier) {
      p <- p + geom_boxplot(mapping = aes(fill = .data[[x]]), outlier.shape = NA, width = width)
    } else {
      p <- p + geom_boxplot(mapping = aes(fill = .data[[x]]), outlier.size = pt.size, width = width)
    }
  }
  if(pt) {
    pt.style <- pt.style[1]
    if(!pt.style %in% c("quasirandom", "jitter")) stop('"pt.style" shoule be "quasirandom" or "jitter"')
    if(pt.style == "jitter") p <- p + geom_jitter(width = width/2.2, size = pt.size, alpha= pt.alpha)
    if(pt.style == "quasirandom") {
      import("ggbeeswarm")
      p <- p + geom_quasirandom(size = pt.size, width = width/2, alpha= pt.alpha)
    }
  }
  if(box & violin) {
    if(pt | hide.outlier) {
      p <- p + geom_boxplot(outlier.shape = NA, width = 0.12, fill = "white")
    } else {
      p <- p + geom_boxplot(fill = "white", outlier.size = pt.size, width = 0.12, outlier.alpha = pt.alpha)
    }
  }

  p <- p + scale_fill_disc_auto(color_scheme = cols, n = n)

  if(x == "f"){
    p <- p +
      facet_wrap(vars(feature), ncol = ncol, strip.position=strip.position, scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none",
            axis.text.x = element_text(angle = 45,hjust = 1),
            strip.text = element_text(face = "bold", size = 10)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.08)))
  }else{
    p <- p +
      facet_grid(vars(feature), vars(f), switch = c("both"), scales = scales)+
      ylab(NULL) +
      xlab(NULL) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            axis.text.x = element_blank(),
            strip.text.y = element_text(face = "bold", size = 10)) +
      labs(fill = lab_fill) +
      scale_y_continuous(expand = expansion(mult = c(0,0.08)))
  }
  return(p)
}

scVisVlnPlot_Stat <- function(
    p,
    stat.method = c("none", "wilcox.test", "t.test"),
    p.adjust.method = "holm",
    label = c("p.signif","p","p.adj","p.format"),
    comparisons = NULL,
    hide.ns = TRUE,
    step.increase = 0.12,
    tip.length = 0.03,
    ...
) {
  stat.method <- stat.method[1]
  if(stat.method %in% c("wilcox.test", "t.test")) {
    library(ggpubr)
    library(dplyr)

    scores <- p$data

    if (!"f2" %in% colnames(scores)) {
      formula <- value ~ f
      group_by_arg <- "feature"
    } else {
      formula <- value ~ f2
      group_by_arg <- c("feature", "f")
    }
    stat.test <- compare_means(
      formula,
      data = scores,
      method = stat.method,
      group.by = group_by_arg,
      p.adjust.method = p.adjust.method
    )

    if(hide.ns == TRUE) {
      stat.test <- filter(stat.test, p.signif != "ns")
      if(nrow(stat.test) == 0) {
        message("No statistical significance.")
        return(p)
      }
    }
    stat.test$groups <- apply(stat.test, 1, function(x) c(x[["group1"]], x[["group2"]]), simplify = F)
    if(!is.null(comparisons)) {
      level_comparisons <- lapply(comparisons, function(pair) {
        if(is.numeric(pair)) {
          if (!"f2" %in% colnames(scores)) {
            levels(scores$f)[pair]
          } else {
            levels(scores$f2)[pair]
          }
        } else if(is.character(pair)) {
          pair
        } else {
          stop("Invalid comparison type")
        }
      })
      stat.test <- filter(stat.test, groups %in% level_comparisons)
    }
    stat.test <- scVisVlnPlot_Stat_add_y(stat.test, scores = scores, step.increase = step.increase)
    p <- p +
      stat_pvalue_manual(
        stat.test,
        label = label[1],
        tip.length = tip.length,
        ...)
  }
  return(p)
}

scVisVlnPlot_Stat_add_y <- function(stat.test, scores, step.increase) {
  summary_data <- scores %>%
    group_by(feature) %>%
    dplyr::summarize(
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE)
    )
  summary_data <- summary_data %>%
    dplyr::mutate(
      step = (max_value - min_value) * step.increase,
      start = (max_value - min_value) * 0.1 + max_value
    )
  grouping_cols <- if ("f2" %in% colnames(scores)) c("feature", "f") else "feature"
  stat.test <- stat.test %>%
    left_join(summary_data, by = "feature") %>%
    group_by(across(all_of(grouping_cols))) %>%
    dplyr::mutate(y.position = start + step * (row_number() - 1)) %>%
    select(-min_value, -max_value, -step, -start) %>%
    ungroup()
  return(stat.test)
}

scVisVlnPlot_SelColDisc <- function(
    seu,
    group.by,
    split.by,
    cols,
    load.cols
) {
  if(is.null(cols)) return(NULL)
  if(cols[1] != "auto") return(cols)
  if(is.null(group.by)) group.by <- "ident"
  if(!is.null(split.by)) {
    var <- split.by
  } else {
    var <- group.by
  }
  if(length(var) == 1) {
    load_var <- seu@misc[["var_colors"]][[var]]
    if(!is.null(load_var)) {
      cols <- load_var
    } else {
      cols <- "pro_default"
    }
  }
  return(cols)
}



# ============ 5. scVisTissueOR ==========
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
  OR_list <- do.tissueDist(cellInfo.tb = meta_data,
                                  out.prefix = sprintf("%s.Immune_cell", output_prefix),
                                  pdf.width = pdf_width, pdf.height = pdf_height, verbose = verbose)

  # 返回分析结果
  return(OR_list)
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
#' @param output_dir A character string specifying the prefix for output files. Default is "./output_figure".
#' @param output_file A character string specifying the output file name. Default is "tissue_OR.pdf".
#' @param width Numeric value specifying the width of the output plot. Default is 5.
#' @param height Numeric value specifying the height of the output plot. Default is 4.
#' @return A list with tissue OR analysis results and a heatmap plot.
#' @export
#' @import Seurat
#' @import ggplot2
#' @import export
#' @import plyr
#' @examples
#' \dontrun{
#' # Assuming `sce` is a pre-existing Seurat object
#' scVisTissueOR(sce, group = "group", celltype = "celltype", output_dir = "./output_figure", output_file = "tissue_OR.pdf", width = 5, height = 4)
#' }

scVisTissueOR <- function(scRNA, group = 'orig.ident',
                          celltype = 'celltype',
                          output_dir = "./output_figure",
                          output_file = "tissue_OR.pdf", width = 5, height = 4) {

  # Load necessary libraries
  library(Seurat)
  library(ggplot2)
  library(export)
  library(plyr)

  # Ensure the output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Extract metadata
  meta <- scRNA@meta.data

  # Check if group and celltype columns exist in metadata
  if (!(group %in% colnames(meta))) {
    stop(paste("Grouping variable", group, "not found in metadata"))
  }
  if (!(celltype %in% colnames(meta))) {
    stop(paste("Cell type variable", celltype, "not found in metadata"))
  }

  # Assign group and celltype to new columns
  meta$loc <- meta[, group]
  meta$meta.cluster <- meta[, celltype]

  # Perform tissue OR analysis
  OR_list <- analyze_tissue_dist(meta_data = meta, output_prefix = output_dir)

  # Plot heatmap and save to file
  pdf(file = file.path(output_dir, output_file), width = width, height = height)
  plot_heatmap(OR_list)
  dev.off()

  return(OR_list)
}

# Example usage
# Assuming `sce` is a pre-existing Seurat object
# scVisTissueOR(sce, group = "group", celltype = "celltype", output_dir = "./output_figure", output_file = "tissue_OR.pdf", width = 5, height = 4)

# ============= 6.scVisHeatmap =============
# R/scVisHeatmap.R
#' Plot the heatmap of single cell dataset
#'
#' This function generates a heatmap for a single cell dataset using specified markers and annotation variables.
#'
#' @param dataset A Seurat object.
#' @param markers A dataframe generated by the FindAllMarkers function in Seurat, or a character vector specifying the names of genes.
#' @param sort_var A character vector specifying the variables to be sorted in the heatmap. These variables should exist in the metadata of the dataset. Default is "seurat_clusters".
#' @param n An integer specifying the number of genes to be plotted for each seurat_cluster. This parameter will not be used if the markers are directly specified. Default is 8.
#' @param anno_var A character vector specifying the variables in the annotation. These variables should exist in the metadata of the dataset.
#' @param anno_colors A list specifying the color specification of annotation bars. The length should be the same as anno_var, with each list element corresponding to each variable. If not provided, a default color set will be used.
#' @param hm_limit A numeric vector of three values dictating the lowest limit, midpoint, and highest limit of the color gradient. Default is c(-2, 0, 2).
#' @param hm_colors A character vector of length 3 specifying the colors corresponding to the lowest limit, midpoint, and highest limit specified in hm_limit. Default is c("#4575b4", "white", "#d73027").
#' @param row_font_size An integer specifying the font size of row names. Default is 12.
#'
#' @return A heatmap plot.
#' @examples
#' \dontrun{
#' scVisHeatmap(dataset = sce2,
#'              markers = top5$gene,
#'              sort_var = c("celltype","sample"),
#'              anno_var = c("celltype","sample","percent.mt","AUCell"),
#'              anno_colors = list("Set2", my36colors, "Reds", "Greens"))
#' }
#' @export
scVisHeatmap <- function(dataset, markers, sort_var = c("seurat_clusters"), n = 8,
                         anno_var, anno_colors = NULL, hm_limit = c(-2, 0, 2),
                         hm_colors = c("#4575b4", "white", "#d73027"), row_font_size = 12) {

  library(Scillus)

  # default color
  default_colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                               '#968175')

  # if not provide anno_colors, it will use the default
  if (is.null(anno_colors)) {
    anno_colors <- list(default_colors)
  }

  plot_heatmap(dataset = dataset,
               markers = markers,
               sort_var = sort_var,
               n = n,
               anno_var = anno_var,
               anno_colors = anno_colors,
               hm_limit = hm_limit,
               hm_colors = hm_colors,
               row_font_size = row_font_size)
}



# ============ 7. CPDB Visualization =============


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
