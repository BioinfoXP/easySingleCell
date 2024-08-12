#  =============== 1. PlotViolin ==================
#' PlotViolin
#'
#' This function creates a violin plot with optional significance testing.
#'
#' @param df Data frame containing the data to be plotted.
#' @param x Variable name for grouping (column name).
#' @param y Numeric variable name (column name).
#' @param comparisons List of groups to compare. If not provided, no significance testing is performed.
#' @param fill.col Vector of colors for filling (default is colorRampPalette(brewer.pal(9, "Set1"))(6)).
#' @param color Vector of background colors (default is colorRampPalette(brewer.pal(11, "BrBG"))(30)).
#' @param title Title of the plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param angle_x_text Angle of the x-axis text (default is 45).
#' @param legend_position Position of the legend (default is "none").
#' @param signif_test Method for significance testing (default is "t.test").
#' @param signif_map Whether to use asterisks to show significance (default is TRUE).
#' @param signif_tip_length Length of the significance markers (default is c(0.01)).
#' @param x_limits Limits for the x-axis (default is NULL, which means ggplot2 will determine the limits), e.g., c(0, 20).
#' @param y_limits Limits for the y-axis (default is NULL, which means ggplot2 will determine the limits).
#' @param x_breaks Interval for the x-axis breaks (default is 1).
#' @param y_breaks Interval for the y-axis breaks (default is 1).
#' @return A ggplot object.
#' @examples
#' df <- data.frame(samples = rep(c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2"), each = 10),
#'                  values = rnorm(60))
#' comparisons <- list(c("A_1", "A_2"), c("B_1", "B_2"), c("C_1", "C_2"))
#' PlotViolin(df, x = "samples", y = "values", comparisons = comparisons)
#' @export
#' @import ggplot2
#' @import ggpubr
#' @import ggsignif
#' @import tidyverse
#' @import ggprism
#' @import vioplot
#' @import RColorBrewer
#' @import grid
#' @import scales
PlotViolin <- function(df, x, y, comparisons = NULL,
                       fill.col = colorRampPalette(brewer.pal(9, "Set1"))(6),
                       color = colorRampPalette(brewer.pal(11, "BrBG"))(30),
                       title = NULL, xlab = NULL, ylab = NULL,
                       angle_x_text = 45, legend_position = "none",
                       signif_test = "t.test", signif_map = TRUE,
                       signif_tip_length = c(0.01),
                       x_limits = NULL, y_limits = NULL,
                       x_breaks = 1, y_breaks = 1) {

  # Create the plot
  p <- ggplot(df, aes_string(x = x, y = y, fill = x)) +
    geom_violin(trim = TRUE, position = position_dodge(width = 0.1), scale = 'width') +
    geom_boxplot(alpha = 1, outlier.size = 0, size = 0.3, width = 0.2, fill = "white") +
    stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, fill = "blue") +
    labs(x = xlab, y = ylab, title = title) +
    theme_prism() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          legend.position = legend_position,
          axis.text = element_text(color = 'black', size = 12),
          legend.text = element_text(color = 'black', size = 12),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = angle_x_text, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = fill.col)

  # Set x and y axis limits if provided
  if (!is.null(x_limits)) {
    p <- p + scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], x_breaks))
  }
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], y_breaks))
  }

  # Add significance layer if comparisons are provided
  if (!is.null(comparisons)) {
    p <- p + geom_signif(comparisons = comparisons,
                         map_signif_level = signif_map,
                         test = signif_test,
                         tip_length = signif_tip_length,
                         size = 0.8, color = "black")
  }

  return(p)
}


