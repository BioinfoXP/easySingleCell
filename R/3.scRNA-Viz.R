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
#' @param show.axis.title Logical, whether to show axis titles (e.g. UMAP_1). Default FALSE.
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
                         pt.size = 0.15,
                         stroke = 0,
                         alpha = 1,
                         show.arrow = TRUE,
                         show.axis.title = FALSE, # 新增参数：控制轴标题显示
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
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = dims[1], y = dims[2], color = group_col)) +
    ggplot2::geom_point(size = pt.size, shape = 16, stroke = stroke, alpha = alpha) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = dims[1], y = dims[2]) # 这里已经设置了标题内容

  # 5. 应用“张泽民团队”风格主题

  # 定义坐标轴线条样式 (带箭头)
  axis_line_setting <- if (show.arrow) {
    ggplot2::element_line(colour = "black", size = 0.3,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"), type = "closed"))
  } else {
    ggplot2::element_line(colour = "black", size = 0.3)
  }

  # 定义坐标轴标题样式 (新增逻辑)
  axis_title_setting <- if (show.axis.title) {
    ggplot2::element_text(size = 10, face = "plain", color = "black") # 如果True，显示文字
  } else {
    ggplot2::element_blank() # 如果False，隐藏文字
  }

  p <- p + ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),

    # 坐标轴设置
    axis.title = axis_title_setting,   # 应用刚才定义的标题样式
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

  # 7. 优化图例点的显示
  p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  return(p)
}



# =============== CellRatioPlot  ================
# =============== 3.scVisCellRatioPlot  ================
#' @title Visualize Cell Type Proportion Differences (Fixed Stats)
#' @description Calculates cell type proportions and performs statistical comparisons.
#' Fixed the "missing value" error by correctly handling grouped statistical testing.
#'
#' @param sce A Seurat object.
#' @param group.by Column name for conditions (e.g., "Treatment").
#' @param cell.type Column name for cell types.
#' @param sample.by Column name for samples.
#' @param test.method Default "wilcox.test".
#' @param label.type Default "p.signif" (*).
#' @param cols Colors.
#' @param pt.size Point size.
#' @param width Box width.
#'
#' @export
scVisCellRatioPlot <- function(sce,
                               group.by,
                               cell.type = "celltype",
                               sample.by = "orig.ident",
                               test.method = "wilcox.test",
                               label.type = "p.signif",
                               cols = NULL,
                               pt.size = 1.5,
                               width = 0.6,
                               ...) {

  # 1. Check input
  if (!all(c(group.by, cell.type, sample.by) %in% colnames(sce@meta.data))) {
    stop("Column not found in meta.data")
  }

  # 2. Calculate Proportions
  message("Calculating cell proportions per sample...")
  meta_df <- Seurat::FetchData(sce, vars = c(group.by, cell.type, sample.by))
  colnames(meta_df) <- c("Group", "CellType", "Sample")

  ratio_data <- meta_df %>%
    dplyr::group_by(Sample, Group, CellType) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    tidyr::complete(Sample, CellType, fill = list(n = 0)) %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Group = unique(stats::na.omit(Group))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Ratio = n / sum(n)) %>%
    dplyr::ungroup()

  # 3. Colors
  if (is.null(cols)) {
    n_groups <- length(unique(ratio_data$Group))
    if (n_groups <= 10) {
      cols <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF")
    } else {
      cols <- scales::hue_pal()(n_groups)
    }
  }

  # 4. Plotting
  p <- ggplot2::ggplot(ratio_data, ggplot2::aes(x = CellType, y = Ratio, fill = Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8, width = width,
                          position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_jitter(ggplot2::aes(color = Group),
                         position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                         size = pt.size, alpha = 0.8, show.legend = FALSE) +

    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "", y = "Cell Proportion", fill = group.by) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dashed"),
      ...
    )

  # 5. Statistical Tests (关键修复部分)
  # 在分组箱线图中，不要传递 comparisons 参数，直接让 ggpubr 识别 group
  tryCatch({
    p <- p + ggpubr::stat_compare_means(
      aes(group = Group), # 显式告诉 ggplot 按 Group 分组进行比较
      method = test.method,
      label = label.type,
      hide.ns = TRUE,    # 是否隐藏无意义的结果
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                         symbols = c("****", "***", "**", "*", "ns"))
    )
  }, error = function(e) {
    message("Statistical annotation failed. Returning plot without stats.")
    print(e)
  })

  return(p)
}
