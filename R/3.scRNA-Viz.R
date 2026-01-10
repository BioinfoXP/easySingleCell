# =============== FeaturePlot  ================
# =============== 1.scVisFeaturePlot  ================
#' @title Feature Plot for Single-cell Data
#' @description This function creates customizable feature plots for single-cell data.
#' It allows for multi-color gradients, clean themes, and flexible layouts.
#'
#' @param scRNA A Seurat object.
#' @param features Vector of features to plot (e.g., gene names).
#' @param reduction Reduction method (default "umap").
#' @param pt.size Point size (default 0.1).
#' @param max.cutoff Value for quantile clipping (default "q95").
#' @param cols A vector of colors for the gradient.
#' @param ncol Number of columns for layout.
#' @param nrow Number of rows for layout.
#' @param show.legend Logical, whether to show the legend (default FALSE).
#' @param plot.title Optional global title for the combined plot.
#' @param ... Additional arguments passed to \code{Seurat::FeaturePlot}.
#'
#' @return A patchwork object.
#' @export
#'
#' @importFrom Seurat FeaturePlot
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_blank element_text ggtitle scale_x_continuous scale_y_continuous
#' @importFrom patchwork wrap_plots plot_annotation
#'
scVisFeaturePlot <- function(scRNA,
                             features,
                             reduction = "umap",
                             pt.size = 0.1,
                             max.cutoff = "q95", # 推荐用分位数，比固定数值1.5更健壮
                             cols = c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),
                             ncol = NULL,
                             nrow = NULL,
                             show.legend = FALSE,
                             plot.title = NULL,
                             ...) {

  # 内部辅助函数：处理单个 Feature
  create_single_plot <- function(feature) {

    # 检查 Feature 是否存在 (避免报错中断)
    if (!feature %in% rownames(scRNA) && !feature %in% colnames(scRNA@meta.data)) {
      warning(paste("Feature", feature, "not found in object."))
      return(NULL)
    }

    # 基础绘图 (注意：这里把 cols留空，我们用 ggplot 图层手动加颜色)
    p <- Seurat::FeaturePlot(
      object = scRNA,
      features = feature,
      reduction = reduction,
      pt.size = pt.size,
      max.cutoff = max.cutoff,
      ...
    )

    # 应用自定义颜色渐变和主题
    p <- p +
      ggplot2::scale_color_gradientn(colors = cols) + # 关键：正确应用多色渐变
      ggplot2::scale_x_continuous("") +
      ggplot2::scale_y_continuous("") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = if(show.legend) "right" else "none" # 灵活控制图例
      ) +
      ggplot2::ggtitle(feature) # 子图标题始终为基因名

    return(p)
  }

  # 循环生成图形列表
  plot_list <- lapply(features, create_single_plot)

  # 移除 NULL (即没找到的基因)
  plot_list <- plot_list[!sapply(plot_list, is.null)]

  if (length(plot_list) == 0) {
    stop("No valid features found to plot.")
  }

  # 使用 patchwork 拼图
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol, nrow = nrow)

  # 如果有总标题，添加总标题
  if (!is.null(plot.title)) {
    combined_plot <- combined_plot + patchwork::plot_annotation(title = plot.title)
  }

  return(combined_plot)
}


# =============== DimPlot  ================
# =============== 2.scVisDimPlot  ================
#' @title Cloud/Misty Style DimPlot (Zemin Zhang Style)
#' @description Creates a DimPlot with a "cloud/misty" aesthetic commonly seen in Zemin Zhang's lab publications.
#' Key features include small point sizes, zero stroke, and arrowed axes.
#'
#' @param scRNA A Seurat object.
#' @param reduction Reduction method (default "umap").
#' @param group.by Vector of meta.data columns to group by (color).
#' @param split.by Vector of meta.data columns to split by (facet).
#' @param colors Vector of colors. If NULL, uses a default "Zhang Lab" style palette.
#' @param pt.size Point size. Default is 0.05 (crucial for the "misty" look).
#' @param stroke Point stroke thickness. Default is 0 (crucial for the "misty" look).
#' @param alpha Point transparency. Default is 1.
#' @param show.arrow Logical, whether to show arrows on axes. Default TRUE.
#' @param strip.color Background color for facet strips. Default "#e6bac5".
#' @param legend.position Position of legend ("right", "bottom", "none", etc.). Default "right".
#' @param ... Additional arguments passed to ggplot2 theme.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic theme element_blank element_line element_rect element_text arrow unit facet_wrap facet_grid labs guide_legend guides
scVisDimPlot <- function(scRNA,
                         reduction = "umap",
                         group.by = NULL,
                         split.by = NULL,
                         colors = NULL,
                         pt.size = 0.05,
                         stroke = 0,
                         alpha = 1,
                         show.arrow = TRUE,
                         strip.color = "#e6bac5",
                         legend.position = "right",
                         aspect.ratio = 1,
                         ...) {

  # 1. 检查 Reduction 并获取坐标轴名称
  if (!reduction %in% names(scRNA)) {
    stop(paste("Reduction", reduction, "not found in object."))
  }

  # 获取坐标列名 (例如 UMAP_1, UMAP_2)
  emb <- scRNA[[reduction]]
  dims <- colnames(emb)[1:2]
  key <- emb@key

  # 2. 准备绘图数据
  # 如果没有指定 group.by，默认使用当前的 Idents
  if (is.null(group.by)) {
    group.by <- "ident"
    plot_data <- Seurat::FetchData(scRNA, vars = c(dims, split.by))
    plot_data$ident <- Seurat::Idents(scRNA)
  } else {
    plot_data <- Seurat::FetchData(scRNA, vars = c(dims, group.by, split.by))
  }

  # 确保分组变量是因子
  group_col <- if (is.null(group.by)) "ident" else group.by
  plot_data[[group_col]] <- as.factor(plot_data[[group_col]])

  # 3. 设置默认颜色 (参考原文配色)
  if (is.null(colors)) {
    # 原文中的低饱和度高级灰配色
    colors <- c("#919ac2","#ffac98","#70a4c8","#a5a9af","#63917d",
                "#dbd1b4","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6",
                "#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")

    # 如果类别数超过颜色数，扩展颜色
    n_groups <- length(levels(plot_data[[group_col]]))
    if (n_groups > length(colors)) {
      colors <- scales::hue_pal()(n_groups)
    } else {
      colors <- colors[1:n_groups]
    }
  }

  # 4. 构建 ggplot 基础图层
  # 注意：核心在于 shape=16, stroke=0
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = dims[1], y = dims[2], color = group_col)) +
    ggplot2::geom_point(size = pt.size, shape = 16, stroke = stroke, alpha = alpha) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = dims[1], y = dims[2])

  # 5. 应用“张泽民团队”风格主题

  # 定义坐标轴样式 (带箭头)
  axis_line_setting <- if (show.arrow) {
    ggplot2::element_line(colour = "black", size = 0.3,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"), type = "closed"))
  } else {
    ggplot2::element_line(colour = "black", size = 0.3)
  }

  p <- p + ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),

    # 坐标轴设置
    axis.title = ggplot2::element_blank(), # 原文通常不显示 Label，只显示箭头
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = axis_line_setting,

    # 图例和比例
    legend.position = legend.position,
    aspect.ratio = aspect.ratio,

    # 分面标题背景 (Strip)
    strip.background = ggplot2::element_rect(fill = strip.color, color = NA),
    strip.text = ggplot2::element_text(size = 8, face = "bold"),
    strip.placement = "outside"
  )

  # 6. 处理分面 (Split.by)
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", split.by)))
  }

  # 7. 优化图例点的显示 (防止图例点太小看不清)
  p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  return(p)
}
