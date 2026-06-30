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
#'   # 1. Default plot with jitter points.
#'   VizBox(df, x_col = "Group", y_col = "Value")
#'
#'   # 2. Hide jitter points for large datasets.
#'   VizBox(df, x_col = "Group", y_col = "Value", add_dots = FALSE)
#' }
VizBox <- function(data,
                   x_col,
                   y_col,
                   fill_col = NULL,
                   comparisons = NULL,
                   test_method = "wilcox.test",
                   add_dots = TRUE,          # \u65B0\u589E\u63A7\u5236\u53C2\u6570
                   title = NULL,
                   xlab = "",
                   ylab = y_col,
                   pal = NULL,
                   jitter.size = 1,
                   jitter.width = 0.1,
                   font.size = 15,
                   step.percent = 0.12) {

  # 1. \u53C2\u6570\u9ED8\u8BA4\u503C\u5904\u7406
  if (is.null(fill_col)) fill_col <- x_col

  # 2. \u9ED8\u8BA4\u914D\u8272 (\u66F4\u65B0\u4E3A\u6700\u65B0\u63D0\u4F9B\u7684\u914D\u8272)
  if (is.null(pal)) {
    pal <- c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")
  }

  # 3. \u667A\u80FD\u751F\u6210\u6BD4\u8F83\u5217\u8868
  final_comparisons <- NULL
  if (!is.null(comparisons)) {
    if (is.list(comparisons)) {
      final_comparisons <- comparisons
    } else if (is.vector(comparisons) && length(comparisons) >= 2) {
      final_comparisons <- utils::combn(comparisons, 2, simplify = FALSE)
    }
  }

  # 4. === \u6838\u5FC3\u667A\u80FD\u7B97\u6CD5\uFF1A\u8BA1\u7B97 Y \u8F74\u4F4D\u7F6E ===
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

  # 5. \u7ED8\u56FE
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +

    # \u57FA\u7840\u7BB1\u7EBF\u56FE (\u65E0\u79BB\u7FA4\u70B9)
    ggplot2::geom_boxplot(outlier.shape = NA) +

    # \u667A\u80FD Y \u8F74\u8303\u56F4
    ggplot2::ylim(ylim_custom) +

    # \u914D\u8272\u4E0E\u6807\u7B7E
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +

    # \u4E3B\u9898\u7F8E\u5316
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(size = font.size),
      axis.text = ggplot2::element_text(size = font.size),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0),
      plot.margin = ggplot2::unit(c(0.5, 1, 0, 1), "cm")
    )

  # 6. \u53EF\u9009\uFF1A\u6DFB\u52A0\u6296\u52A8\u70B9
  if (add_dots) {
    p <- p + ggplot2::geom_jitter(width = jitter.width, color = "grey20", size = jitter.size)
  }

  # 7. \u6DFB\u52A0\u663E\u8457\u6027\u68C0\u9A8C
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
#'   # 1. Default plot with jitter points.
#'   VizBar(df, x_col = "Group", y_col = "Value")
#'
#'   # 2. Hide jitter points.
#'   VizBar(df, x_col = "Group", y_col = "Value", add_dots = FALSE)
#' }
VizBar <- function(data,
                   x_col,
                   y_col,
                   fill_col = NULL,
                   comparisons = NULL,
                   test_method = "t.test",
                   error_type = "sd",
                   add_dots = TRUE,        # \u63A7\u5236\u53C2\u6570
                   title = NULL,
                   xlab = "",
                   ylab = y_col,
                   pal = NULL,
                   jitter.size = 1,
                   jitter.width = 0.2,
                   bar.width = 0.7,
                   font.size = 15,
                   step.percent = 0.12) {

  # 1. \u53C2\u6570\u9ED8\u8BA4\u503C\u5904\u7406
  if (is.null(fill_col)) fill_col <- x_col

  # 2. \u9ED8\u8BA4\u914D\u8272 (\u4E0E VizBox \u4FDD\u6301\u4E00\u81F4)
  if (is.null(pal)) {
    pal <- c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")
  }

  # 3. \u667A\u80FD\u751F\u6210\u6BD4\u8F83\u5217\u8868
  final_comparisons <- NULL
  if (!is.null(comparisons)) {
    if (is.list(comparisons)) {
      final_comparisons <- comparisons
    } else if (is.vector(comparisons) && length(comparisons) >= 2) {
      final_comparisons <- utils::combn(comparisons, 2, simplify = FALSE)
    }
  }

  # 4. === \u6838\u5FC3\u667A\u80FD\u7B97\u6CD5 ===
  y_vals <- data[[y_col]]
  x_vals <- data[[x_col]]

  # 4.1 \u8BA1\u7B97\u7EDF\u8BA1\u91CF
  agg_mean <- stats::aggregate(y_vals ~ x_vals, FUN = mean)
  colnames(agg_mean) <- c("group", "mean")

  agg_sd <- stats::aggregate(y_vals ~ x_vals, FUN = stats::sd)
  colnames(agg_sd) <- c("group", "sd")

  # 4.2 \u8BA1\u7B97\u8BEF\u5DEE
  if (error_type == "se") {
    agg_n <- stats::aggregate(y_vals ~ x_vals, FUN = length)
    error_val <- agg_sd$sd / sqrt(agg_n[,2])
  } else {
    error_val <- agg_sd$sd
  }

  # 4.3 \u627E\u5230\u89C6\u89C9\u6700\u9AD8\u70B9 (BarTop vs MaxDataPoint)
  max_bar_top <- max(agg_mean$mean + error_val, na.rm = TRUE)
  min_data_val <- min(y_vals, na.rm = TRUE)

  # \u5982\u679C\u6709\u70B9\uFF0C\u6700\u9AD8\u70B9\u53D6 max(bar_top, max_data)
  # \u5982\u679C\u6CA1\u70B9\uFF0C\u6700\u9AD8\u70B9\u5C31\u662F max(bar_top)
  if (add_dots) {
    max_visual_point <- max(max_bar_top, max(y_vals, na.rm=TRUE))
  } else {
    max_visual_point <- max_bar_top
  }

  y_rng <- max_visual_point - min_data_val
  if (y_rng == 0) y_rng <- 1

  step_height <- y_rng * step.percent

  # 4.4 \u663E\u8457\u6027\u9AD8\u5EA6
  y_pos_vec <- NULL
  if (!is.null(final_comparisons)) {
    n_comp <- length(final_comparisons)
    y_pos_vec <- max_visual_point + seq(from = step_height * 0.5, by = step_height, length.out = n_comp)
  }

  top_line_height <- if(is.null(y_pos_vec)) max_visual_point else max(y_pos_vec)
  ylim_custom <- c(min(0, min_data_val * 1.1), top_line_height + step_height * 0.5)

  # 6. \u5B9A\u4E49\u7EDF\u8BA1\u51FD\u6570
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

  # 7. \u7ED8\u56FE
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +

    # \u67F1\u5B50
    ggplot2::stat_summary(
      fun = mean,
      geom = "bar",
      width = bar.width,
      color = "black",
      size = 0.5
    ) +

    # \u8BEF\u5DEE\u68D2
    ggplot2::stat_summary(
      fun.data = chosen_fun_data,
      geom = "errorbar",
      width = 0.2,
      size = 0.5,
      color = "black"
    ) +

    # Y\u8F74\u8303\u56F4
    ggplot2::ylim(ylim_custom) +

    # \u6837\u5F0F
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

  # \u53EF\u9009: \u6563\u70B9
  if (add_dots) {
    p <- p + ggplot2::geom_jitter(
      width = jitter.width,
      color = "black",
      size = jitter.size,
      alpha = 0.6
    )
  }

  # \u663E\u8457\u6027
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


# =============== Grouped BoxPlot  ================
# =============== 3.VizGroupedBox  ================

#' @title Visualize Grouped Boxplots with Adjusted P-values (Native Implementation)
#' @description A wrapper to create grouped boxplots (e.g., Genes on X, Conditions as Groups)
#' with manually calculated statistical tests. Uses native ggplot2 layers for significance bars
#' to avoid compatibility warnings and ensure precise alignment.
#'
#' @param data A data frame.
#' @param x_col String. Column name for the X-axis (e.g., "Gene").
#' @param y_col String. Column name for the Y-axis value (e.g., "Expression").
#' @param group_col String. Column name for the grouping variable (e.g., "Condition").
#' @param test_method String. "t.test" (default) or "wilcox.test".
#' @param p_adjust_method String. Correction method, passed to \code{p.adjust}. Default "bonferroni" or "BH".
#' @param symnum.args List. Arguments for significance symbols. Default list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")).
#' @param pal Character vector. Custom color palette. Default is c("#B2DF8A", "#33A02C").
#' @param add_dots Logical. Whether to show jitter points. Default TRUE.
#' @param jitter.size Numeric. Size of jitter points. Default 1.5.
#' @param jitter.alpha Numeric. Transparency of jitter points. Default 0.6.
#' @param box.width Numeric. Width of boxplots. Default 0.6.
#' @param legend.position String or Vector. Position of the legend. Default "top". Examples: "right", "none", c(0.8, 0.9).
#' @param title String. Plot title.
#' @param xlab String. X-axis label.
#' @param ylab String. Y-axis label.
#' @param font.size Numeric. Font size for axis text. Default 15.
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter position_jitterdodge scale_fill_manual labs theme element_text unit theme_classic element_blank geom_segment geom_text
#' @importFrom cowplot theme_cowplot
#' @importFrom stats wilcox.test t.test p.adjust aggregate symnum
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   # 1. Default legend at the top.
#'   VizGroupedBox(df, x_col = "gene", y_col = "value", group_col = "group")
#'
#'   # 2. Put legend on the right.
#'   VizGroupedBox(df, x_col = "gene", y_col = "value", group_col = "group",
#'                 legend.position = "right")
#'
#'   # 3. Hide the legend.
#'   VizGroupedBox(df, x_col = "gene", y_col = "value", group_col = "group",
#'                 legend.position = "none")
#' }
VizGroupedBox <- function(data,
                          x_col,
                          y_col,
                          group_col,
                          test_method = "t.test",
                          p_adjust_method = "bonferroni",
                          symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                             symbols = c("****", "***", "**", "*", "ns")),
                          pal = NULL,
                          add_dots = TRUE,
                          jitter.size = 1.5,
                          jitter.alpha = 0.6,
                          box.width = 0.6,
                          legend.position = "top", # \u65B0\u589E\u53C2\u6570
                          title = NULL,
                          xlab = "",
                          ylab = "Expression Level",
                          font.size = 15) {

  # 1. \u9ED8\u8BA4\u914D\u8272
  if (is.null(pal)) {
    pal <- c("#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C")
  }

  # 2. \u9501\u5B9A\u56E0\u5B50\u987A\u5E8F
  if (!is.factor(data[[x_col]])) {
    data[[x_col]] <- factor(data[[x_col]], levels = unique(data[[x_col]]))
  }
  x_levels <- levels(data[[x_col]])

  if (!is.factor(data[[group_col]])) {
    data[[group_col]] <- factor(data[[group_col]])
  }
  group_levels <- levels(data[[group_col]])

  if (length(group_levels) != 2) {
    stop("VizGroupedBox currently supports exactly 2 groups for pairwise comparison.")
  }

  # 3. === \u6838\u5FC3\u7EDF\u8BA1\u8BA1\u7B97 ===
  res_list <- lapply(x_levels, function(curr_x) {
    sub_df <- data[data[[x_col]] == curr_x, ]
    vals_g1 <- sub_df[[y_col]][sub_df[[group_col]] == group_levels[1]]
    vals_g2 <- sub_df[[y_col]][sub_df[[group_col]] == group_levels[2]]

    pval <- NA
    if (length(vals_g1) > 0 && length(vals_g2) > 0) {
      if (test_method == "wilcox.test") {
        try({pval <- stats::wilcox.test(vals_g1, vals_g2)$p.value}, silent=TRUE)
      } else {
        try({pval <- stats::t.test(vals_g1, vals_g2)$p.value}, silent=TRUE)
      }
    }

    max_val <- max(sub_df[[y_col]], na.rm = TRUE)
    return(data.frame(x = curr_x, p = pval, max_val = max_val, stringsAsFactors = FALSE))
  })

  stat_df <- do.call(rbind, res_list)

  # P\u503C\u6821\u6B63
  stat_df$p.adj <- stats::p.adjust(stat_df$p, method = p_adjust_method)
  stat_df$p.signif <- stats::symnum(stat_df$p.adj, corr = FALSE, na = FALSE,
                                    cutpoints = symnum.args$cutpoints,
                                    symbols = symnum.args$symbols)
  stat_df$p.signif <- as.character(stat_df$p.signif)

  # \u5750\u6807\u8BA1\u7B97
  stat_df$x_numeric <- as.numeric(factor(stat_df$x, levels = x_levels))
  stat_df$xmin <- stat_df$x_numeric - 0.2
  stat_df$xmax <- stat_df$x_numeric + 0.2

  global_range <- max(data[[y_col]], na.rm = TRUE) - min(data[[y_col]], na.rm = TRUE)
  if (global_range == 0) global_range <- 1
  stat_df$y.position <- stat_df$max_val + global_range * 0.1

  # 4. \u7ED8\u56FE
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[group_col]])) +

    # 4.1 \u7BB1\u7EBF\u56FE
    ggplot2::geom_boxplot(
      width = box.width,
      outlier.shape = NA,
      position = ggplot2::position_dodge(width = 0.8),
      size = 0.5
    ) +

    # 4.2 \u6296\u52A8\u70B9
    {if(add_dots)
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
        size = jitter.size,
        alpha = jitter.alpha,
        color = "black"
      )
    } +

    # 4.3 \u6807\u7B7E\u4E0E\u914D\u8272
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +

    # 4.4 \u4E3B\u9898 (\u4F7F\u7528 legend.position \u53C2\u6570)
    cowplot::theme_cowplot() +
    ggplot2::theme(
      # \u81EA\u5B9A\u4E49\u56FE\u4F8B\u4F4D\u7F6E
      legend.position = legend.position,
      legend.justification = "center",
      legend.title = ggplot2::element_blank(),

      # \u6807\u9898\u5C45\u4E2D
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = font.size + 2),

      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black", size = font.size),
      axis.text.y = ggplot2::element_text(color = "black", size = font.size),
      axis.title = ggplot2::element_text(size = font.size),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # 5. \u539F\u751F\u7ED8\u5236\u663E\u8457\u6027
  if (nrow(stat_df) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = stat_df,
        ggplot2::aes(x = xmin, xend = xmax, y = y.position, yend = y.position),
        inherit.aes = FALSE,
        size = 0.6,
        color = "black"
      ) +
      ggplot2::geom_text(
        data = stat_df,
        ggplot2::aes(x = (xmin + xmax) / 2, y = y.position, label = p.signif),
        inherit.aes = FALSE,
        size = 5,
        vjust = -0.2
      )
  }

  return(p)
}
