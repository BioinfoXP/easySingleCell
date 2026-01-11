# =============== BoxPlot  ================
# =============== 1.VizBox  ================

#' @title Generate Publication-Ready Boxplots (Smart Statistics)
#' @description A flexible wrapper for creating Nature-style boxplots.
#' Features smart positioning of significance bars based on data range, 
#' custom aesthetics (jitter + clean theme), and automatic pairwise comparisons.
#'
#' @param data A data frame.
#' @param x_col String. Column name for the X-axis grouping variable.
#' @param y_col String. Column name for the Y-axis numeric variable.
#' @param fill_col String. Column name for fill color. Default is same as x_col.
#' @param comparisons A character vector (e.g., \code{c("A", "B", "C")}) or list of pairs.
#' If a vector is provided, pairwise comparisons are generated automatically.
#' @param test_method String. "wilcox.test" (default) or "t.test".
#' @param add_dots Logical. Whether to show jitter points. Default TRUE.
#' @param title String. Plot title.
#' @param xlab String. X-axis label.
#' @param ylab String. Y-axis label.
#' @param pal Character vector. Custom color palette. 
#' Default is \code{c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")}.
#' @param jitter.size Numeric. Size of the jitter points. Default 1.
#' @param jitter.width Numeric. Width of the jitter points. Default 0.1.
#' @param font.size Numeric. Font size for axis text. Default 15.
#' @param step.percent Numeric. Percentage of Y-range to step up for each significance bar. Default 0.12 (12%).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter ylim scale_fill_manual labs theme element_text unit
#' @importFrom cowplot theme_cowplot
#' @importFrom ggsignif geom_signif
#' @importFrom rlang .data
#' @importFrom utils combn
#'
#' @examples
#' \dontrun{
#'   # 1. 默认显示散点
#'   VizBox(df, x_col = "Group", y_col = "Value")
#'   
#'   # 2. 隐藏散点 (适用于大数据量)
#'   VizBox(df, x_col = "Group", y_col = "Value", add_dots = FALSE)
#' }
VizBox <- function(data, 
                   x_col, 
                   y_col, 
                   fill_col = NULL, 
                   comparisons = NULL,
                   test_method = "wilcox.test",
                   add_dots = TRUE,          # 新增控制参数
                   title = NULL, 
                   xlab = "", 
                   ylab = y_col,
                   pal = NULL,
                   jitter.size = 1,
                   jitter.width = 0.1,
                   font.size = 15,
                   step.percent = 0.12) {
  
  # 1. 参数默认值处理
  if (is.null(fill_col)) fill_col <- x_col
  
  # 2. 默认配色 (更新为最新提供的配色)
  if (is.null(pal)) {
    pal <- c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")
  }
  
  # 3. 智能生成比较列表
  final_comparisons <- NULL
  if (!is.null(comparisons)) {
    if (is.list(comparisons)) {
      final_comparisons <- comparisons
    } else if (is.vector(comparisons) && length(comparisons) >= 2) {
      final_comparisons <- utils::combn(comparisons, 2, simplify = FALSE)
    }
  }
  
  # 4. === 核心智能算法：计算 Y 轴位置 ===
  y_vals <- data[[y_col]]
  y_max_data <- max(y_vals, na.rm = TRUE)
  y_min_data <- min(y_vals, na.rm = TRUE)
  y_rng <- y_max_data - y_min_data
  
  if (y_rng == 0) y_rng <- abs(y_max_data) * 0.5
  if (y_rng == 0) y_rng <- 1 
  
  step_height <- y_rng * step.percent
  
  y_pos_vec <- NULL
  if (!is.null(final_comparisons)) {
    n_comp <- length(final_comparisons)
    y_pos_vec <- y_max_data + seq(from = step_height * 0.5, by = step_height, length.out = n_comp)
  }
  
  top_line_height <- if(is.null(y_pos_vec)) y_max_data else max(y_pos_vec)
  ylim_custom <- c(y_min_data - y_rng * 0.1, top_line_height + step_height * 0.5)
  
  # 5. 绘图
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    
    # 基础箱线图 (无离群点)
    ggplot2::geom_boxplot(outlier.shape = NA) +
    
    # 智能 Y 轴范围
    ggplot2::ylim(ylim_custom) +
    
    # 配色与标签
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    
    # 主题美化
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(size = font.size),
      axis.text = ggplot2::element_text(size = font.size),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0),
      plot.margin = ggplot2::unit(c(0.5, 1, 0, 1), "cm")
    )
  
  # 6. 可选：添加抖动点
  if (add_dots) {
    p <- p + ggplot2::geom_jitter(width = jitter.width, color = "grey20", size = jitter.size)
  }
  
  # 7. 添加显著性检验
  if (!is.null(final_comparisons)) {
    p <- p + ggsignif::geom_signif(
      comparisons = final_comparisons,
      test = test_method,
      y_position = y_pos_vec,
      map_signif_level = TRUE,
      textsize = 5,
      tip_length = 0,
      vjust = 0.4
    )
  }
  
  return(p)
}


# =============== BarPlot  ================
# =============== 2.VizBar  ================

#' @title Generate Publication-Ready Bar Plots (Mean + Error Bar)
#' @description A flexible wrapper for creating Nature-style bar plots.
#' Features smart positioning of significance bars based on error bar heights,
#' automatic pairwise comparisons, and options for SD/SE error bars.
#'
#' @param data A data frame.
#' @param x_col String. Column name for the X-axis grouping variable.
#' @param y_col String. Column name for the Y-axis numeric variable.
#' @param fill_col String. Column name for fill color. Default is same as x_col.
#' @param comparisons A character vector (e.g., \code{c("A", "B", "C")}) or list of pairs.
#' If a vector is provided, pairwise comparisons are generated automatically.
#' @param test_method String. "t.test" (default) or "wilcox.test".
#' @param error_type String. "sd" (Mean +/- SD) or "se" (Mean +/- SE). Default "sd".
#' @param add_dots Logical. Whether to add jitter points over the bars. Default TRUE.
#' @param title String. Plot title.
#' @param xlab String. X-axis label.
#' @param ylab String. Y-axis label.
#' @param pal Character vector. Custom color palette. 
#' Default is \code{c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")}.
#' @param jitter.size Numeric. Size of the jitter points. Default 1.
#' @param jitter.width Numeric. Width of the jitter points. Default 0.2.
#' @param bar.width Numeric. Width of the bars. Default 0.7.
#' @param font.size Numeric. Font size for axis text. Default 15.
#' @param step.percent Numeric. Percentage of Y-range to step up for each significance bar. Default 0.12 (12%).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes stat_summary geom_jitter ylim scale_fill_manual labs theme element_text unit position_dodge geom_errorbar
#' @importFrom cowplot theme_cowplot
#' @importFrom ggsignif geom_signif
#' @importFrom rlang .data
#' @importFrom utils combn
#' @importFrom stats aggregate sd
#'
#' @examples
#' \dontrun{
#'   # 1. 默认显示散点
#'   VizBar(df, x_col = "Group", y_col = "Value")
#'   
#'   # 2. 隐藏散点 (add_dots = FALSE)
#'   VizBar(df, x_col = "Group", y_col = "Value", add_dots = FALSE)
#' }
VizBar <- function(data, 
                   x_col, 
                   y_col, 
                   fill_col = NULL, 
                   comparisons = NULL,
                   test_method = "t.test", 
                   error_type = "sd", 
                   add_dots = TRUE,        # 控制参数
                   title = NULL, 
                   xlab = "", 
                   ylab = y_col,
                   pal = NULL,
                   jitter.size = 1,
                   jitter.width = 0.2,
                   bar.width = 0.7,
                   font.size = 15,
                   step.percent = 0.12) {
  
  # 1. 参数默认值处理
  if (is.null(fill_col)) fill_col <- x_col
  
  # 2. 默认配色 (与 VizBox 保持一致)
  if (is.null(pal)) {
    pal <- c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")
  }
  
  # 3. 智能生成比较列表
  final_comparisons <- NULL
  if (!is.null(comparisons)) {
    if (is.list(comparisons)) {
      final_comparisons <- comparisons
    } else if (is.vector(comparisons) && length(comparisons) >= 2) {
      final_comparisons <- utils::combn(comparisons, 2, simplify = FALSE)
    }
  }
  
  # 4. === 核心智能算法 ===
  y_vals <- data[[y_col]]
  x_vals <- data[[x_col]]
  
  # 4.1 计算统计量
  agg_mean <- stats::aggregate(y_vals ~ x_vals, FUN = mean)
  colnames(agg_mean) <- c("group", "mean")
  
  agg_sd <- stats::aggregate(y_vals ~ x_vals, FUN = stats::sd)
  colnames(agg_sd) <- c("group", "sd")
  
  # 4.2 计算误差
  if (error_type == "se") {
    agg_n <- stats::aggregate(y_vals ~ x_vals, FUN = length)
    error_val <- agg_sd$sd / sqrt(agg_n[,2])
  } else {
    error_val <- agg_sd$sd
  }
  
  # 4.3 找到视觉最高点 (BarTop vs MaxDataPoint)
  max_bar_top <- max(agg_mean$mean + error_val, na.rm = TRUE)
  min_data_val <- min(y_vals, na.rm = TRUE)
  
  # 如果有点，最高点取 max(bar_top, max_data)
  # 如果没点，最高点就是 max(bar_top)
  if (add_dots) {
    max_visual_point <- max(max_bar_top, max(y_vals, na.rm=TRUE))
  } else {
    max_visual_point <- max_bar_top
  }
  
  y_rng <- max_visual_point - min_data_val
  if (y_rng == 0) y_rng <- 1
  
  step_height <- y_rng * step.percent
  
  # 4.4 显著性高度
  y_pos_vec <- NULL
  if (!is.null(final_comparisons)) {
    n_comp <- length(final_comparisons)
    y_pos_vec <- max_visual_point + seq(from = step_height * 0.5, by = step_height, length.out = n_comp)
  }
  
  top_line_height <- if(is.null(y_pos_vec)) max_visual_point else max(y_pos_vec)
  ylim_custom <- c(min(0, min_data_val * 1.1), top_line_height + step_height * 0.5)
  
  # 6. 定义统计函数
  mean_sdl_custom <- function(x) {
    m <- mean(x)
    s <- sd(x)
    return(c(y = m, ymin = m - s, ymax = m + s))
  }
  
  mean_se_custom <- function(x) {
    m <- mean(x)
    se <- sd(x) / sqrt(length(x))
    return(c(y = m, ymin = m - se, ymax = m + se))
  }
  
  chosen_fun_data <- if (error_type == "se") mean_se_custom else mean_sdl_custom
  
  # 7. 绘图
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    
    # 柱子
    ggplot2::stat_summary(
      fun = mean, 
      geom = "bar", 
      width = bar.width, 
      color = "black", 
      size = 0.5
    ) +
    
    # 误差棒
    ggplot2::stat_summary(
      fun.data = chosen_fun_data, 
      geom = "errorbar", 
      width = 0.2, 
      size = 0.5, 
      color = "black"
    ) +
    
    # Y轴范围
    ggplot2::ylim(ylim_custom) +
    
    # 样式
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(size = font.size),
      axis.text = ggplot2::element_text(size = font.size),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0),
      plot.margin = ggplot2::unit(c(0.5, 1, 0, 1), "cm")
    )
  
  # 可选: 散点
  if (add_dots) {
    p <- p + ggplot2::geom_jitter(
      width = jitter.width, 
      color = "black", 
      size = jitter.size, 
      alpha = 0.6
    )
  }
  
  # 显著性
  if (!is.null(final_comparisons)) {
    p <- p + ggsignif::geom_signif(
      comparisons = final_comparisons,
      test = test_method,
      y_position = y_pos_vec,
      map_signif_level = TRUE,
      textsize = 5,
      tip_length = 0,
      vjust = 0.4
    )
  }
  
  return(p)
}