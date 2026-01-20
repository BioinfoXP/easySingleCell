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
#' @description Creates a DimPlot with a "cloud/misty" aesthetic.
#'
#' @param scRNA A Seurat object.
#' @param reduction Reduction method (default "umap").
#' @param group.by Vector of meta.data columns to group by (color).
#' @param split.by Vector of meta.data columns to split by (facet).
#' @param colors Vector of colors.
#' @param pt.size Point size. Default is 0.05.
#' @param stroke Point stroke thickness. Default is 0.
#' @param alpha Point transparency. Default is 1.
#' @param show.arrow Logical, whether to show arrows on axes. Default TRUE.
#' @param show.axis.title Logical, whether to show axis titles. Default FALSE.
#' @param label Logical, whether to label the clusters. Default FALSE.
#' @param label.size Numeric, size of the label text. Default 4.
#' @param label.color String, color of the label text. Default "black".
#' @param repel Logical, whether to use ggrepel to avoid overlap. Default FALSE.
#' @param label.box Logical, whether to draw a box around the label. Default FALSE.
#' @param strip.color Background color for facet strips.
#' @param legend.position Position of legend.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData Idents
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic theme labs element_blank element_line element_rect element_text arrow unit facet_wrap guides guide_legend geom_text geom_label
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom stats aggregate median
scVisDimPlot <- function(scRNA,
                         reduction = "umap",
                         group.by = NULL,
                         split.by = NULL,
                         colors = NULL,
                         pt.size = 0.15,
                         stroke = 0,
                         alpha = 1,
                         show.arrow = TRUE,
                         show.axis.title = FALSE,
                         label = FALSE,        # 新增：是否标记
                         label.size = 4,       # 新增：字体大小
                         label.color = "black",# 新增：字体颜色
                         repel = FALSE,        # 新增：是否防重叠
                         label.box = FALSE,    # 新增：是否加背景框
                         strip.color = "#e6bac5",
                         legend.position = "right",
                         aspect.ratio = 1,
                         ...) {

  # 1. 检查 Reduction 并获取坐标轴名称
  if (!reduction %in% names(scRNA)) {
    stop(paste("Reduction", reduction, "not found in object."))
  }

  emb <- scRNA[[reduction]]
  dims <- colnames(emb)[1:2]
  key <- emb@key

  # 2. 准备绘图数据
  if (is.null(group.by)) {
    group.by <- "ident"
    plot_data <- Seurat::FetchData(scRNA, vars = c(dims, split.by))
    plot_data$ident <- Seurat::Idents(scRNA)
  } else {
    plot_data <- Seurat::FetchData(scRNA, vars = c(dims, group.by, split.by))
  }

  group_col <- if (is.null(group.by)) "ident" else group.by
  plot_data[[group_col]] <- as.factor(plot_data[[group_col]])

  # 3. 设置默认颜色
  if (is.null(colors)) {
    colors <- c("#919ac2","#ffac98","#70a4c8","#a5a9af","#63917d",
                "#dbd1b4","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6",
                "#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")
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
    ggplot2::labs(x = dims[1], y = dims[2])

  # ================== 新增 Label 处理逻辑 ==================
  if (label) {
    # 计算每个群体的中心点 (Median)
    # 如果有 split.by，需要按 group + split 分组计算
    if (!is.null(split.by)) {
      group_vars <- c(group_col, split.by)
    } else {
      group_vars <- c(group_col)
    }

    # 使用 base R aggregate 计算中位数，避免增加 dplyr 依赖
    label_data <- stats::aggregate(
      list(x = plot_data[[dims[1]]], y = plot_data[[dims[2]]]),
      by = plot_data[group_vars],
      FUN = stats::median
    )

    # 定义绘图参数
    geom_fun <- if (label.box) {
      if (repel) ggrepel::geom_label_repel else ggplot2::geom_label
    } else {
      if (repel) ggrepel::geom_text_repel else ggplot2::geom_text
    }

    # 添加标签图层
    p <- p + geom_fun(
      data = label_data,
      ggplot2::aes_string(x = "x", y = "y", label = group_col),
      color = label.color,
      size = label.size,
      show.legend = FALSE,
      inherit.aes = FALSE # 关键：不继承主图的aes，防止报错
    )
  }
  # =======================================================

  # 5. 应用“张泽民团队”风格主题
  axis_line_setting <- if (show.arrow) {
    ggplot2::element_line(colour = "black", size = 0.3,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"), type = "closed"))
  } else {
    ggplot2::element_line(colour = "black", size = 0.3)
  }

  axis_title_setting <- if (show.axis.title) {
    ggplot2::element_text(size = 10, face = "plain", color = "black")
  } else {
    ggplot2::element_blank()
  }

  p <- p + ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = axis_title_setting,
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = axis_line_setting,
    legend.position = legend.position,
    aspect.ratio = aspect.ratio,
    strip.background = ggplot2::element_rect(fill = strip.color, color = NA),
    strip.text = ggplot2::element_text(size = 8, face = "bold"),
    strip.placement = "outside"
  )

  # 6. 处理分面
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", split.by)))
  }

  # 7. 优化图例点的显示
  p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  return(p)
}



# =============== DotPlot  ================
# =============== 3.scVisDotPlot  ================

#' @title Visualize Marker Expression using a Styled DotPlot
#' @description A wrapper around \code{Seurat::DotPlot} with enhanced aesthetics.
#' Features customizable color palettes, refined legends (hollow circles), and consistent axis styling.
#'
#' @param object A Seurat object.
#' @param features A vector of features (genes) to plot.
#' @param group.by String. Name of the meta.data column to group the cells by. Default is NULL (uses active idents).
#' @param pal Character vector. Custom color palette for the gradient.
#' If NULL, defaults to \code{rev(RColorBrewer::brewer.pal(7, "RdBu"))}.
#' @param rot_angle Numeric. Angle of x-axis text labels. Default is 45.
#' @param font.size Numeric. Font size for axis labels. Default is 12.
#' @param ... \strong{Additional arguments passed to Seurat::DotPlot}.
#' Examples: \code{split.by}, \code{assay}, \code{cols} (will be overwritten by pal), \code{scale}.
#'
#' @return A ggplot object.
#' @export
#' @importFrom Seurat DotPlot
#' @importFrom ggplot2 scale_color_gradientn scale_size_continuous guide_colorbar guide_legend guides theme_bw theme element_blank element_rect element_line element_text margin unit
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#'   markers <- c("CD3D", "CD79A", "MS4A1", "CD14", "FCGR3A")
#'
#'   # 1. 默认用法 (45度旋转，RdBu配色，字体大小12)
#'   scVisDotPlot(pbmc, features = markers)
#'
#'   # 2. 如果想改回90度
#'   scVisDotPlot(pbmc, features = markers, rot_angle = 90)
#'
#'   # 3. 自定义配色
#'   my_pal <- c("#2166ac", "#f7f7f7", "#b2182b")
#'   scVisDotPlot(pbmc, features = markers, pal = my_pal)
#' }
scVisDotPlot <- function(object,
                         features,
                         group.by = NULL,
                         pal = NULL,
                         rot_angle = 45,  # 默认修改为 45 度
                         font.size = 12,
                         ...) {

  # 1. 确定配色方案
  if (is.null(pal)) {
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
    } else {
      pal <- c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B")
    }
  }

  # 2. 调用 Seurat::DotPlot
  p <- Seurat::DotPlot(
    object = object,
    features = features,
    group.by = group.by,
    ...
  )

  # 3. 移除 Seurat 默认的颜色标尺
  p$scales$scales <- list()

  # 4. 应用美化图层
  p <- p +
    # --- 颜色设置 ---
    ggplot2::scale_color_gradientn(
      colors = pal,
      guide = ggplot2::guide_colorbar(
        title = "Mean expression",
        barwidth = 0.8,
        barheight = 4,
        ticks = TRUE
      )
    ) +

    # --- 大小设置 ---
    ggplot2::scale_size_continuous(
      range = c(0, 5),
      breaks = c(25, 50, 75, 100),
      labels = c("25", "50", "75", "100")
    ) +

    # --- 图例样式 ---
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title = "Fraction of cells (%)",
        override.aes = list(
          shape = 21,
          colour = "black",
          fill = NA,
          stroke = 0.5
        )
      )
    ) +

    # --- 主题设置 ---
    ggplot2::theme_bw() +
    ggplot2::theme(
      # 去除分面背景
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),

      # 面板间距与边框
      panel.spacing = ggplot2::unit(0.1, "lines"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),

      # 网格线
      panel.grid.major = ggplot2::element_line(color = "grey95", linewidth = 0.2),

      # --- 字体美化核心区域 ---

      # X轴: 基因名 (斜体), 黑色, 45度角优化
      axis.text.x = ggplot2::element_text(
        size = font.size,
        color = "black",
        face = "italic",
        angle = rot_angle,
        hjust = 1, # 45度时，hjust=1 保证文字末端对齐刻度
        vjust = 1  # 45度时，vjust=1 保证文字不与轴线重叠
      ),

      # Y轴: 细胞类型 (常规), 黑色
      axis.text.y = ggplot2::element_text(
        size = font.size,
        color = "black",
        face = "plain"
      ),

      # 去除轴标题
      axis.title = ggplot2::element_blank(),

      # 页边距
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  return(p)
}

# =============== scVisRatioBox  ================
# =============== 4.scVisRatioBox  ================
#' @title scVis: Cell Proportion Boxplot with Stats
#' @description Calculates cell type proportions per sample and performs statistical comparisons between groups using Boxplots.
#' Automatically handles missing cell types (fills 0) and maps samples back to their groups correctly.
#' Uses base R pipe (R >= 4.1.0).
#'
#' @param sce A Seurat object.
#' @param group_col String. Column name for conditions (e.g., "Treatment").
#' @param celltype_col String. Column name for cell types.
#' @param sample_col String. Column name for biological samples (e.g., "orig.ident").
#' @param comparisons List of vectors. Pairs to compare (e.g., \code{list(c("A", "B"))}). If NULL, performs global test.
#' @param sign_method String. Test method: "wilcox.test" (default), "t.test", "kruskal.test", or "anova".
#' @param sign_label String. Label type: "p.signif" (*) or "p.format" (numbers).
#' @param palette Vector. Custom colors. If NULL, uses default vibrant palette.
#' @param pt_size Numeric. Size of jitter points. Default 1.5.
#' @param box_width Numeric. Width of boxplots. Default 0.6.
#' @param base_size Numeric. Base font size. Default 14.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData
#' @importFrom dplyr group_by summarise mutate ungroup n select distinct left_join
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter scale_fill_manual scale_color_manual scale_y_continuous labs theme_classic theme element_text element_line position_dodge position_jitterdodge
#' @importFrom ggpubr stat_compare_means
#' @importFrom scales percent hue_pal
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   library(Seurat)
#'   library(ggplot2)
#'
#'   # --- 1. Mock Data ---
#'   set.seed(123)
#'   n_cells <- 1000
#'   counts <- matrix(rpois(n_cells * 10, 1), nrow = 10, ncol = n_cells)
#'   rownames(counts) <- paste0("Gene_", 1:10)
#'   colnames(counts) <- paste0("Cell_", 1:n_cells)
#'   sce <- CreateSeuratObject(counts = counts)
#'
#'   # Simulate Metadata
#'   sce$Group <- sample(c("Ctrl", "Treat"), n_cells, replace = TRUE)
#'   sce$Sample <- paste0(sce$Group, "_Rep", sample(1:3, n_cells, replace = TRUE))
#'   sce$CellType <- sample(c("T_cell", "B_cell", "Macro"), n_cells, replace = TRUE)
#'
#'   # --- 2. Plot ---
#'   p <- scVisRatioBox(
#'     sce = sce,
#'     group_col = "Group",
#'     celltype_col = "CellType",
#'     sample_col = "Sample",
#'     comparisons = list(c("Treat", "Ctrl")), # Pairwise comparison
#'     sign_label = "p.signif"
#'   )
#'   print(p)
#' }
scVisRatioBox <- function(sce,
                          group_col,
                          celltype_col,
                          sample_col,
                          comparisons = NULL,
                          sign_method = "wilcox.test",
                          sign_label = "p.signif",
                          palette = NULL,
                          pt_size = 1.5,
                          box_width = 0.6,
                          base_size = 14) {

  # 1. Input Validation
  if (!inherits(sce, "Seurat")) stop("Input 'sce' must be a Seurat object.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required.")

  # Check columns
  req_cols <- c(group_col, celltype_col, sample_col)
  if (!all(req_cols %in% colnames(sce@meta.data))) {
    missing <- req_cols[!req_cols %in% colnames(sce@meta.data)]
    stop(paste("❌ Columns not found in meta.data:", paste(missing, collapse = ", ")))
  }

  # 2. Data Preparation
  # Extract raw data
  meta_df <- Seurat::FetchData(sce, vars = c(group_col, celltype_col, sample_col))
  colnames(meta_df) <- c("Group", "CellType", "Sample")

  # Create a lookup table for Sample -> Group mapping
  # This is SAFER than na.omit() after complete(), which can lose group info
  sample_info <- meta_df |>
    dplyr::select("Sample", "Group") |>
    dplyr::distinct()

  # Calculate Proportions
  ratio_data <- meta_df |>
    dplyr::group_by(.data$Sample, .data$CellType) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    # Fill missing combinations with 0 counts
    tidyr::complete(.data$Sample, .data$CellType, fill = list(n = 0)) |>
    # Re-attach Group information
    dplyr::left_join(sample_info, by = "Sample") |>
    # Calculate Ratio
    dplyr::group_by(.data$Sample) |>
    dplyr::mutate(Ratio = .data$n / sum(.data$n)) |>
    dplyr::ungroup()

  # 3. Color Palette Logic (Standard 15 colors)
  n_groups <- length(unique(ratio_data$Group))

  if (is.null(palette)) {
    default_pal <- c(
      "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
      "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
      "#f47c7c", "#aadea7", "#f7f494", "#72ccff", "#c999ff"
    )
    if (n_groups > length(default_pal)) {
      palette <- scales::hue_pal()(n_groups)
    } else {
      palette <- default_pal[1:n_groups]
    }
  } else {
    # Handle user palette recycling
    if (length(palette) < n_groups) palette <- rep(palette, length.out = n_groups)
  }

  # 4. Plotting
  p <- ggplot2::ggplot(ratio_data, ggplot2::aes(x = .data$CellType, y = .data$Ratio, fill = .data$Group)) +

    # Boxplot
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      width = box_width,
      position = ggplot2::position_dodge(width = 0.8)
    ) +

    # Jitter Points
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data$Group),
      position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = pt_size,
      alpha = 0.8,
      show.legend = FALSE
    ) +

    # Aesthetics
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(x = NULL, y = "Cell Proportion", fill = group_col) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dashed")
    )

  # 5. Statistical Tests
  tryCatch({
    if (is.null(comparisons)) {
      # Global Test (e.g., ANOVA/Kruskal if >2 groups, or basic comparison)
      p <- p + ggpubr::stat_compare_means(
        ggplot2::aes(group = .data$Group),
        method = sign_method,
        label = sign_label,
        hide.ns = TRUE,
        label.y.npc = "top" # Auto position at top
      )
    } else {
      # Pairwise Comparisons (Specified by user)
      p <- p + ggpubr::stat_compare_means(
        mapping = ggplot2::aes(group = .data$Group),
        comparisons = comparisons,
        method = sign_method,
        label = sign_label,
        hide.ns = TRUE,
        symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                           symbols = c("****", "***", "**", "*", "ns"))
      )
    }
  }, error = function(e) {
    message("⚠️ Statistical annotation failed. Returning plot without stats.")
    message(e$message)
  })

  return(p)
}


# =============== Ro/e 分布偏好  ================
# =============== 5.scVisRoePlot  ================
# https://mp.weixin.qq.com/s/Fhz72Cjbd3zzCi1tB9xP5Q
#' @title Visualize Ro/e (Tissue Enrichment) Heatmap
#' @description Calculates and visualizes the Ratio of observed to expected (Ro/e) cell numbers,
#' a metric popularized by Zemin Zhang's lab to quantify tissue enrichment/depletion preference.
#'
#' @param sce A Seurat object.
#' @param group.by Column name in meta.data representing the tissue/group (e.g., "Tissue", "Group").
#' @param cell.type Column name in meta.data representing cell types (e.g., "celltype").
#' @param sample.by Column name in meta.data representing biological samples/patients (e.g., "orig.ident", "Patient").
#' @param title The title of the heatmap. Default is "Ro/e Tissue Enrichment".
#' @param method Statistical method for calculation. Default "chisq".
#' @param min.rowSum Minimum row sum to keep a cell type. Default 0.
#' @param display.mode Mode for cell annotation ("symbol", "numeric", "none").
#' @param cluster_rows Logical, whether to cluster rows. Default TRUE.
#' @param cluster_cols Logical, whether to cluster columns. Default FALSE.
#' @param font.size Font size for numbers/symbols. Default 10.
#' @param ... Additional arguments passed to \code{pheatmap}.
#'
#' @return A pheatmap object (invisibly returns the Ro/e matrix).
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#'
scVisRoePlot <- function(sce,
                         group.by,
                         cell.type = "celltype",
                         sample.by = "orig.ident",
                         title = "Ro/e Tissue Enrichment", # 新增标题参数
                         method = "chisq",
                         min.rowSum = 0,
                         display.mode = c("symbol", "numeric", "none"),
                         cluster_rows = TRUE,
                         cluster_cols = FALSE,
                         font.size = 10,
                         ...) {

  # 1. 检查依赖包 Startrac
  if (!requireNamespace("Startrac", quietly = TRUE)) {
    stop("Package 'Startrac' is required for Ro/e calculation.\nPlease install it using: devtools::install_github('Japrin/STARTRAC')")
  }

  # 2. 检查输入列
  if (!all(c(group.by, cell.type, sample.by) %in% colnames(sce@meta.data))) {
    stop("One or more specified columns not found in meta.data.")
  }

  display.mode <- match.arg(display.mode)

  # 3. 计算 Ro/e 矩阵
  message("Calculating Ro/e matrix using Startrac...")
  meta_data <- sce@meta.data

  # 调用 Startrac (注意：直接传 meta_data 作为第一个参数)
  roe_mat <- Startrac::calTissueDist(meta_data,
                                     byPatient = FALSE,
                                     colname.cluster = cell.type,
                                     colname.patient = sample.by,
                                     colname.tissue = group.by,
                                     method = method,
                                     min.rowSum = min.rowSum)

  # 4. 设置视觉风格
  my_palette <- grDevices::colorRampPalette(c("#483D8B", "#00FFFF", "#F8F8FF", "#FF69B4", "#8B008B"))(100)

  # 设置对称断点
  max_abs_deviation <- max(abs(roe_mat - 1), na.rm = TRUE)
  if (max_abs_deviation == 0) max_abs_deviation <- 0.1

  my_breaks <- seq(1 - max_abs_deviation,
                   1 + max_abs_deviation,
                   length.out = 101)

  # 5. 处理标注
  if (display.mode == "symbol") {
    display_mat <- ifelse(roe_mat > 1, "+++",
                          ifelse(roe_mat > 0.8 & roe_mat <= 1, "++",
                                 ifelse(roe_mat > 0.2 & roe_mat <= 0.8, "+",
                                        ifelse(roe_mat > 0 & roe_mat <= 0.2, "+/-", "-"))))
  } else if (display.mode == "numeric") {
    display_mat <- TRUE
  } else {
    display_mat <- FALSE
  }

  # 6. 绘图 (将 title 传给 main 参数)
  p <- pheatmap::pheatmap(roe_mat,
                          color = my_palette,
                          breaks = my_breaks,
                          display_numbers = display_mat,
                          fontsize_number = font.size,
                          cluster_rows = cluster_rows,
                          cluster_cols = cluster_cols,
                          border_color = "grey90",
                          scale = "none",
                          main = title,  # 这里设置标题
                          ...)

  invisible(roe_mat)
}





# =============== scVisCellFC  ================
# =============== 6.scVisCellFC  ================
#' @title scVis: Cell Type Log2FC Plot
#' @description Calculates and visualizes the Log2 Fold Change of cell type proportions between groups.
#' Optimized for flexibility, supporting both single vector comparisons and list of comparisons.
#'
#' @param sce A Seurat object.
#' @param group_col String. Column name in \code{meta.data} for grouping (e.g., "group").
#' @param celltype_col String. Column name in \code{meta.data} for cell types (e.g., "cell_type").
#' @param comparisons A single vector (e.g., \code{c("Treat", "Ctrl")}) OR a list of vectors (e.g., \code{list(c("A", "B"), c("C", "D"))}).
#' The first element is the Case (Numerator), the second is the Control (Denominator).
#' @param palette Vector. Custom colors. Can be a named vector (names = cell types) or a simple vector of hex codes.
#' @param show_labels Logical. Whether to show text labels on bars. Default TRUE.
#' @param show_legend Logical. Whether to show the legend. Default FALSE.
#' @param legend_position String. Position of legend ("right", "top", "bottom", "left"). Default "right".
#' @param label_size Numeric. Font size for labels. Default 4.
#' @param base_size Numeric. Base font size for the plot theme. Default 14.
#' @param bar_width Numeric. Width of the bars (0-1). Default 0.7.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom dplyr group_by summarise mutate ungroup select filter rename full_join bind_rows
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_col geom_hline facet_wrap scale_fill_manual theme_classic labs theme element_blank element_text element_line geom_text scale_y_continuous expansion
#' @importFrom stats reorder
#' @importFrom Seurat FetchData
#' @importFrom scales hue_pal
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'    # 1. Simple Vector Input (Case vs Control)
#'    scVisCellFC(sce, "Group", "CellType", comparisons = c("Treat", "Ctrl"))
#'
#'    # 2. List Input (Multiple Comparisons)
#'    scVisCellFC(sce, "Group", "CellType",
#'                comparisons = list(c("TreatA", "Ctrl"), c("TreatB", "Ctrl")))
#' }
scVisCellFC <- function(sce,
                        group_col,
                        celltype_col,
                        comparisons,
                        palette = NULL,
                        show_labels = TRUE,
                        show_legend = FALSE,
                        legend_position = "right",
                        label_size = 4,
                        base_size = 14,
                        bar_width = 0.7) {

  # --- 1. Input Validation & Standardization ---

  if (!inherits(sce, "Seurat")) stop("Input 'sce' must be a Seurat object.")

  # Check if columns exist
  if (!all(c(group_col, celltype_col) %in% colnames(sce@meta.data))) {
    stop(paste0("❌ Columns '", group_col, "' or '", celltype_col, "' not found in meta.data."))
  }

  # Smart Comparisons Handling: Normalize to a list
  # Logic: If it's not a list, wrap it. If it is a list, check contents.
  if (!is.list(comparisons)) {
    if (length(comparisons) != 2) {
      stop("❌ For a single comparison, 'comparisons' must be a vector of length 2: c('Case', 'Control').")
    }
    comparisons <- list(comparisons)
  } else {
    # Check if any element in the list is not length 2
    if (any(sapply(comparisons, length) != 2)) {
      stop("❌ Every element in the 'comparisons' list must be a vector of length 2.")
    }
  }

  # Validate Groups Exist in Data
  available_groups <- unique(sce@meta.data[[group_col]])
  all_requested_groups <- unique(unlist(comparisons))
  missing_groups <- setdiff(all_requested_groups, available_groups)

  if (length(missing_groups) > 0) {
    stop(paste0("❌ The following groups defined in 'comparisons' are missing from meta.data: ",
                paste(missing_groups, collapse = ", ")))
  }

  # --- 2. Data Preparation ---

  # Fetch Data safely
  meta <- Seurat::FetchData(sce, vars = c(group_col, celltype_col))
  colnames(meta) <- c("group", "celltype")
  meta$celltype <- as.character(meta$celltype)

  # Calculate Proportions (with 0-filling)
  props <- meta |>
    dplyr::group_by(.data$group, .data$celltype) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(freq = .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    tidyr::complete(.data$group, .data$celltype, fill = list(n = 0, freq = 0))

  # --- 3. Calculate Log2FC ---

  res_list <- list()
  pseudo_count <- 1e-4 # Avoid division by zero

  for (pair in comparisons) {
    g1 <- pair[1] # Case
    g2 <- pair[2] # Control

    d1 <- props |> dplyr::filter(.data$group == g1) |> dplyr::select("celltype", freq1 = "freq")
    d2 <- props |> dplyr::filter(.data$group == g2) |> dplyr::select("celltype", freq2 = "freq")

    merged <- dplyr::full_join(d1, d2, by = "celltype") |>
      dplyr::mutate(
        freq1 = ifelse(is.na(.data$freq1), 0, .data$freq1),
        freq2 = ifelse(is.na(.data$freq2), 0, .data$freq2),
        # Log2FC Calculation
        log2fc = log2((.data$freq1 + pseudo_count) / (.data$freq2 + pseudo_count)),
        comp_name = paste0(g1, " vs. ", g2)
      )

    res_list[[paste(g1, g2)]] <- merged
  }

  plot_data <- dplyr::bind_rows(res_list)

  # --- 4. Color Palette Handling ---

  unique_cells <- sort(unique(plot_data$celltype))
  n_colors <- length(unique_cells)
  default_palette <- c("#00F672","#C8A4F9","#E64B35","#FF00DB","#4DBBD5", "#00A087",
                       "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148")

  if (is.null(palette)) {
    # Use default extended palette or hue_pal if too many cells
    if (n_colors > length(default_palette)) {
      my_colors <- scales::hue_pal()(n_colors)
    } else {
      my_colors <- default_palette[1:n_colors]
    }
    names(my_colors) <- unique_cells
  } else {
    # If user provides palette
    if (!is.null(names(palette))) {
      # Named vector: strictly match
      my_colors <- palette
    } else {
      # Unnamed vector: interpolate or subset
      if (length(palette) < n_colors) {
        warning("Provided palette has fewer colors than cell types. Recycling colors.")
        palette <- rep(palette, length.out = n_colors)
      }
      my_colors <- palette[1:n_colors]
      names(my_colors) <- unique_cells
    }
  }

  # --- 5. Plotting ---

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = stats::reorder(.data$celltype, -.data$log2fc),
                                    y = .data$log2fc,
                                    fill = .data$celltype)) +

    # Zero Line
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +

    # Bars
    ggplot2::geom_col(color = "black", width = bar_width, linewidth = 0.4) +

    # Faceting (Only if multiple comparisons, but safe to always use)
    ggplot2::facet_wrap(~comp_name, scales = "free", ncol = 2) +

    # Colors
    ggplot2::scale_fill_manual(values = my_colors) +

    # Axis Expansion (15% padding for labels)
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +

    # Theme
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(y = "Log2 Fold Change", x = NULL) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.text.y = ggplot2::element_text(color = "black", size = base_size),
      axis.title.y = ggplot2::element_text(color = "black", size = base_size + 1),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size * 1.1, face = "bold", color = "black"),
      legend.position = if (show_legend) legend_position else "none"
    )

  # --- 6. Add Labels ---
  if (show_labels) {
    # Calculate smart offset relative to bar height direction
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        label = .data$celltype,
        vjust = ifelse(.data$log2fc >= 0, -0.5, 1.5) # Push up for pos, down for neg
      ),
      size = label_size,
      color = "black"
    )
  }

  return(p)
}


# =============== scVisDimSplit ================
# =============== 7.scVisDimSplit (Dimension Reduction Version) ================
#' @title scVis: Split Dimension Reduction Plot (Smart Labels & Global Contour)
#' @description A comprehensive Dimension Reduction (DimPlot) visualizer matching high-impact paper styles.
#' Features: Global dashed contour (Ghost), Filtered highlights, Smart label repulsion (ggrepel).
#' Applicable to UMAP, t-SNE, or PCA coordinates.
#' Defaults to showing sample size (n=xxx) at the TOP to avoid overlap.
#'
#' @param sce A Seurat object.
#' @param group_col String. Grouping column (e.g., "group").
#' @param celltype_col String. Cell type column (coloring/labeling).
#' @param groups Vector. Specific groups to plot.
#' @param sub_celltypes Vector. Specific cell types to highlight.
#' @param reduction String. Reduction to use (e.g., "umap", "tsne", "pca"). Default "umap".
#' @param palette Vector. Custom colors.
#' @param pt_size Numeric. Point size. Default 0.1.
#' @param global_contour Logical. Show dashed outline of ALL cells? Default TRUE.
#' @param contour_color String. Contour color. Default "grey70".
#' @param contour_bins Numeric. Contour density. Default 5.
#' @param label_cells Logical. Label cell types? Default TRUE.
#' @param use_repel Logical. Use ggrepel to avoid label overlap? Default TRUE.
#' @param label_size Numeric. Label font size.
#' @param n_pos String. Position of 'n=xxx'. Defaults to "top". Options: "top", "bottom".
#' @param n_shift_y Numeric. Manual vertical offset for 'n' label. Default 0.
#' @param base_size Numeric. Base theme size.
#' @param xlims,ylims Vector. Fixed axis limits.
#' @param ncol Numeric. Grid columns.
#'
#' @return A cowplot grid object.
#' @export
#'
#' @importFrom Seurat Embeddings FetchData
#' @importFrom dplyr group_by summarise filter n mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_text theme_void theme element_rect element_text annotate ggtitle coord_cartesian margin scale_color_manual stat_density_2d after_stat
#' @importFrom cowplot plot_grid
#' @importFrom scales hue_pal comma
#' @importFrom stats median
#' @importFrom rlang .data
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#'   # Works for UMAP
#'   scVisDimSplit(sce, "group", "celltype", reduction = "umap")
#'
#'   # Works for t-SNE
#'   scVisDimSplit(sce, "group", "celltype", reduction = "tsne")
#' }
scVisDimSplit <- function(sce,
                          group_col,
                          celltype_col,
                          groups = NULL,
                          sub_celltypes = NULL,
                          reduction = "umap", # 通用参数
                          palette = NULL,
                          pt_size = 0.1,
                          global_contour = TRUE,
                          contour_color = "grey70",
                          contour_bins = 5,
                          label_cells = TRUE,
                          use_repel = TRUE,
                          label_size = 4,
                          n_pos = "top",
                          n_shift_y = 0,
                          base_size = 14,
                          xlims = c(-15, 15),
                          ylims = c(-15, 15),
                          ncol = NULL) {

  # 1. Validation & Data Prep
  if (!inherits(sce, "Seurat")) stop("Input must be a Seurat object.")
  if (!reduction %in% names(sce@reductions)) stop(paste0("Reduction '", reduction, "' not found."))

  # Extract Global Data for Contour (Generic Dim Reduction)
  coords <- as.data.frame(Seurat::Embeddings(sce, reduction = reduction))[, 1:2]
  colnames(coords) <- c("dim_1", "dim_2")

  meta_cols <- c(celltype_col, group_col)
  meta <- tryCatch({
    Seurat::FetchData(sce, vars = meta_cols)
  }, error = function(e) stop(paste0("Columns not found: ", paste(meta_cols, collapse=", "))))

  global_data <- cbind(coords, meta)
  global_data$celltype_internal <- as.character(global_data[[celltype_col]])
  global_data$group_internal <- as.character(global_data[[group_col]])

  # Filter Data for Dots
  plot_data <- global_data
  if (!is.null(sub_celltypes)) {
    plot_data <- plot_data |> dplyr::filter(.data$celltype_internal %in% sub_celltypes)
    if (nrow(plot_data) == 0) stop("No cells remaining after filtering.")
  }

  target_groups <- if (is.null(groups)) sort(unique(plot_data$group_internal)) else groups

  # Color Logic
  unique_cells <- sort(unique(plot_data$celltype_internal))
  plot_data$celltype_internal <- factor(plot_data$celltype_internal, levels = unique_cells)
  n_types <- length(unique_cells)

  if (is.null(palette)) {
    default_pal <- c("#fb8072", "#8dd3c7", "#ffffb3", "#bebada", "#80b1d3", "#fdb462",
                     "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
    my_colors <- if (n_types > length(default_pal)) scales::hue_pal()(n_types) else default_pal[1:n_types]
  } else {
    my_colors <- if (is.null(names(palette))) rep(palette, length.out = n_types) else palette
  }
  if (is.null(names(my_colors))) names(my_colors) <- unique_cells

  # 2. Plotting Loop
  plot_list <- list()

  # Calculate N Label Position
  if (n_pos == "top") {
    n_label_y_base <- ylims[1] + (diff(ylims) * 0.92)
  } else {
    n_label_y_base <- ylims[1] + (diff(ylims) * 0.08)
  }
  n_label_y <- n_label_y_base + n_shift_y

  for (grp in target_groups) {
    if (!grp %in% unique(plot_data$group_internal)) next

    sub_data <- plot_data |> dplyr::filter(.data$group_internal == grp)
    n_text <- paste0("n = ", scales::comma(nrow(sub_data)))

    # Centroids
    centroids <- NULL
    if (label_cells) {
      centroids <- sub_data |> dplyr::group_by(.data$celltype_internal) |>
        dplyr::summarise(x = stats::median(.data$dim_1), y = stats::median(.data$dim_2), .groups = "drop")
    }

    # Base Plot
    p <- ggplot2::ggplot()

    # Layer 1: Global Contour (The "Map")
    if (global_contour) {
      p <- p + ggplot2::stat_density_2d(
        data = global_data,
        ggplot2::aes(x = .data$dim_1, y = .data$dim_2),
        color = contour_color, linetype = "dashed", bins = contour_bins, size = 0.3, show.legend = FALSE
      )
    }

    # Layer 2: Highlights (The "Dots")
    p <- p +
      ggplot2::geom_point(data = sub_data,
                          ggplot2::aes(x = .data$dim_1, y = .data$dim_2, color = .data$celltype_internal),
                          size = pt_size, show.legend = FALSE) +
      scale_color_manual(values = my_colors) +
      coord_cartesian(xlim = xlims, ylim = ylims) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = base_size, face = "plain", margin = margin(b = 5)),
        panel.border = element_rect(colour = "grey60", fill = NA, linetype = "dashed", linewidth = 0.8),
        legend.position = "none"
      ) +
      ggtitle(grp) +
      # Annotation
      annotate("text", x = mean(xlims), y = n_label_y, label = n_text,
               size = base_size * 0.3, fontface = "italic")

    # Layer 3: Smart Labels
    if (label_cells && !is.null(centroids)) {
      if (use_repel) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = centroids,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$celltype_internal),
            color = "black",
            size = label_size,
            fontface = "bold",
            bg.color = "white",
            bg.r = 0.15,
            seed = 42,
            inherit.aes = FALSE
          )
        } else {
          warning("ggrepel not found. Using standard text.")
          p <- p + geom_text(data = centroids, aes(x=x, y=y, label=celltype_internal),
                             color="black", size=label_size, fontface="bold")
        }
      } else {
        p <- p + geom_text(data = centroids, aes(x=x, y=y, label=celltype_internal),
                           color="black", size=label_size, fontface="bold")
      }
    }
    plot_list[[grp]] <- p
  }

  if (length(plot_list) == 0) stop("❌ No plots generated.")

  final_col <- if (is.null(ncol)) ceiling(sqrt(length(plot_list))) else ncol
  return(cowplot::plot_grid(plotlist = plot_list, nrow = NULL, ncol = final_col))
}



# =============== scVisCellRatio ================
# =============== 8.scVisCellRatio ================
#' @title scVis: Cell Ratio Alluvial Plot
#' @description Visualizes cell type proportion dynamics across groups using an alluvial (Sankey-like) plot.
#' Features a vibrant, high-contrast 15-color palette and adjustable layout options.
#' Uses base R pipe (R >= 4.1.0).
#'
#' @param sce A Seurat object.
#' @param group_col String. Column name in \code{meta.data} for X-axis grouping (e.g., "Condition", "Timepoint").
#' @param celltype_col String. Column name in \code{meta.data} for cell types (used for stacking and coloring).
#' @param group_order Vector. Explicit order for the X-axis groups. Crucial for determining the flow direction.
#' @param palette Vector. Custom colors. If NULL, uses the extended 15-color vibrant palette.
#' @param legend_position String. Position of the legend ("right", "bottom", "top", "left", "none"). Default "right".
#' @param alpha_flow Numeric. Transparency of the connecting flow bands (0-1). Default 0.5.
#' @param bar_width Numeric. Width of the vertical bars (0-1). Default 0.5.
#' @param base_size Numeric. Base font size for the theme. Default 14.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom dplyr group_by summarise mutate ungroup n
#' @importFrom ggplot2 ggplot aes scale_y_continuous labs scale_fill_manual theme_classic theme element_text element_blank element_line
#' @importFrom ggalluvial geom_flow geom_stratum
#' @importFrom Seurat FetchData
#' @importFrom scales percent hue_pal
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   library(Seurat)
#'   library(ggplot2)
#'
#'   # --- 1. Generate Mock Data ---
#'   set.seed(123)
#'   n_cells <- 1000
#'   counts <- matrix(rpois(n_cells * 10, lambda = 1), nrow = 10, ncol = n_cells)
#'   colnames(counts) <- paste0("Cell_", 1:n_cells)
#'   rownames(counts) <- paste0("Gene_", 1:10)
#'
#'   sce <- CreateSeuratObject(counts = counts)
#'   sce$Stage <- sample(c("NC", "MASH", "HCC-N", "HCC-T"), n_cells, replace = TRUE)
#'   sce$CellType <- sample(paste0("Cluster_", 1:8), n_cells, replace = TRUE)
#'
#'   # --- 2. Run Function ---
#'   # Example: Legend at the bottom
#'   p <- scVisCellRatio(
#'     sce = sce,
#'     group_col = "Stage",
#'     celltype_col = "CellType",
#'     group_order = c("NC", "MASH", "HCC-N", "HCC-T"),
#'     legend_position = "bottom"  # <--- Change legend position here
#'   )
#'
#'   print(p)
#' }
scVisCellRatio <- function(sce,
                           group_col,
                           celltype_col,
                           group_order = NULL,
                           palette = NULL,
                           legend_position = "right", # New parameter
                           alpha_flow = 0.5,
                           bar_width = 0.5,
                           base_size = 14) {

  # 1. Input Validation
  if (!inherits(sce, "Seurat")) stop("Input 'sce' must be a Seurat object.")

  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required. Please install it: install.packages('ggalluvial')")
  }

  # Secure Data Extraction
  meta <- tryCatch({
    Seurat::FetchData(sce, vars = c(group_col, celltype_col))
  }, error = function(e) {
    stop(paste0("❌ Columns not found: '", group_col, "' or '", celltype_col, "'. Please check meta.data."))
  })

  colnames(meta) <- c("group", "celltype")

  # 2. Data Processing
  plot_data <- meta |>
    dplyr::group_by(.data$group, .data$celltype) |>
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(Frequency = .data$Count / sum(.data$Count)) |>
    dplyr::ungroup()

  # Handle Group Order
  if (!is.null(group_order)) {
    missing <- setdiff(group_order, unique(plot_data$group))
    if(length(missing) > 0) warning("⚠️ Groups in 'group_order' missing from data.")
    plot_data$group <- factor(plot_data$group, levels = group_order)
  }

  # 3. Color Palette Logic (15 Colors)
  unique_types <- sort(unique(plot_data$celltype))
  n_types <- length(unique_types)

  if (is.null(palette)) {
    # Custom Vibrant Palette (15 colors)
    ref_pal <- c(
      "#32CD32", "#C77CFF", "#FFA500", "#007AFF", "#FF3333",
      "#00C19F", "#FF61CC", "#FFD700", "#00FFFF", "#8A2BE2",
      "#FF4500", "#2E8B57", "#DC143C", "#4169E1", "#8B4513"
    )

    if (n_types > length(ref_pal)) {
      warning("⚠️ More cell types than default 15 colors. Extending palette with hue_pal.")
      my_colors <- scales::hue_pal()(n_types)
    } else {
      my_colors <- ref_pal[1:n_types]
    }
  } else {
    if (is.null(names(palette))) {
      if (length(palette) < n_types) {
        warning("⚠️ Not enough colors provided. Recycling palette.")
        palette <- rep(palette, length.out = n_types)
      }
      my_colors <- palette[1:n_types]
    } else {
      my_colors <- palette
    }
  }

  if (is.null(names(my_colors))) names(my_colors) <- unique_types

  # 4. Plotting
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$group,
                                    y = .data$Frequency,
                                    fill = .data$celltype,
                                    stratum = .data$celltype,
                                    alluvium = .data$celltype)) +

    # Layer 1: Flow
    ggalluvial::geom_flow(alpha = alpha_flow, width = bar_width, curve_type = "linear") +

    # Layer 2: Bar
    ggalluvial::geom_stratum(width = bar_width, color = NA) +

    # Styling
    ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    ggplot2::scale_fill_manual(values = my_colors) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(x = NULL, y = "Proportion") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black", size = base_size, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(color = "black", size = base_size),
      legend.position = legend_position, # <--- Used here
      legend.title = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5)
    )

  return(p)
}




# ========= scVisVlnPlot =========
# ========= 9. scVisVlnPlot =========
#' @title scVis: Intelligent VlnPlot (Split-Stats & Legend Support)
#' @description A robust violin plot builder that handles Seurat V4/V5 data.
#' Features "Smart Stats": Automatically detects if you want to compare X-axis groups
#' OR split-by groups (e.g., Control vs Disease within cell types).
#'
#' @param sce A Seurat object.
#' @param features Vector of feature names.
#' @param group_col String. X-axis grouping (e.g., cell type).
#' @param split_by String. Column to split/fill violins by (e.g., "Condition").
#' @param comparisons Vector/List. Pairs to compare. Can be groups (X-axis) OR conditions (Split-by).
#' @param palette Vector. Custom colors.
#' @param pt_size Numeric. Size of jitter points.
#' @param legend_position String. "right", "top", "bottom", "none", etc.
#' @param ncol Numeric. Number of columns.
#' @param sign_method String. "wilcox.test" or "t.test".
#' @param sign_label String. "p.signif" (stars) or "p.format" (numbers).
#' @param hide_ns Logical. Hide non-significant results.
#' @param ... Additional arguments (ignored safely).
#'
#' @return A patchwork ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData Idents
#' @importFrom ggplot2 ggplot aes geom_violin geom_point stat_summary theme_classic labs scale_fill_manual scale_y_continuous expansion position_jitterdodge position_dodge theme element_text
#' @importFrom ggpubr stat_compare_means compare_means
#' @importFrom patchwork wrap_plots
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   # 1. Split Comparison (Control vs PCOS within each cell type)
#'   scVisVlnPlot(sce, "Gdf15", group_col = "celltype", split_by = "group",
#'                comparisons = c("Control", "PCOS"), legend_position = "top")
#' }
scVisVlnPlot <- function(sce,
                         features,
                         group_col = NULL,
                         split_by = NULL,
                         comparisons = NULL,
                         palette = NULL,
                         pt_size = 0.1,
                         legend_position = "right",
                         ncol = NULL,
                         sign_method = "wilcox.test",
                         sign_label = "p.signif",
                         hide_ns = FALSE,
                         ...) {

  # --- 1. Catch split.by in ... ---
  args <- list(...)
  if (is.null(split_by) && "split.by" %in% names(args)) split_by <- args$split.by

  # --- 2. Validation ---
  if (missing(sce) || is.null(sce)) stop("❌ Input 'sce' is missing.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' required.")

  # Set Group Column
  if (is.null(group_col)) {
    sce$scVis_ident <- Seurat::Idents(sce)
    group_col <- "scVis_ident"
  } else if (!group_col %in% colnames(sce@meta.data)) {
    stop(paste0("❌ Group column '", group_col, "' not found."))
  }

  # Set Split Column & Colors
  color_by <- if (!is.null(split_by)) split_by else group_col
  if (!is.null(split_by) && !split_by %in% colnames(sce@meta.data)) stop(paste0("❌ Split column '", split_by, "' not found."))

  # Palette Logic
  unique_colors <- unique(sce@meta.data[[color_by]])
  n_colors <- length(unique_colors)
  if (is.null(palette)) {
    default_pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2")
    palette <- if (n_colors > length(default_pal)) scales::hue_pal()(n_colors) else default_pal[1:n_colors]
  }
  if (is.null(names(palette))) names(palette) <- unique_colors[1:length(palette)]

  # --- 3. Comparisons Pre-processing ---
  # Normalize comparisons to list
  if (!is.null(comparisons) && is.atomic(comparisons) && !is.list(comparisons)) {
    comparisons <- list(comparisons)
  }

  # --- 4. Plot Loop ---
  plot_list <- list()

  for (gene in features) {
    # Fetch Data (Crash-Proof)
    fetch_vars <- c(gene, group_col, split_by)
    plot_data <- tryCatch({ Seurat::FetchData(sce, vars = fetch_vars) }, error = function(e) NULL)

    if (is.null(plot_data)) { message(paste("⚠️ Skip:", gene)); next }

    # Normalize Columns
    colnames(plot_data)[1] <- "expression"
    colnames(plot_data)[2] <- "group"
    if (!is.null(split_by)) colnames(plot_data)[3] <- "split"

    plot_data <- plot_data[is.finite(plot_data$expression), ]
    plot_data$group <- as.factor(plot_data$group)
    if (!is.null(split_by)) plot_data$split <- as.factor(plot_data$split)

    # --- Build Plot ---
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$group, y = .data$expression))

    if (!is.null(split_by)) {
      p <- p + ggplot2::aes(fill = .data$split)
      # Split Mode Jitter
      if (pt_size > 0) {
        p <- p + ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
                                     size = pt_size, alpha = 0.6, show.legend = FALSE)
      }
    } else {
      p <- p + ggplot2::aes(fill = .data$group)
      # Normal Mode Jitter
      if (pt_size > 0) {
        p <- p + ggplot2::geom_jitter(width = 0.2, size = pt_size, alpha = 0.6, show.legend = FALSE)
      }
    }

    p <- p +
      ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.8,
                           position = ggplot2::position_dodge(width = 0.9),
                           size = 0.4, color = "black") +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white", color = "black",
                            position = ggplot2::position_dodge(width = 0.9), show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = gene, x = NULL, y = "Expression", fill = if(!is.null(split_by)) split_by else group_col) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold.italic"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = legend_position # <--- Applies user setting
      )

    # --- Smart Statistics Logic ---
    if (!is.null(comparisons)) {

      # Determine Compare Mode: "Split" vs "Group"
      mode <- "group" # Default
      flat_comps <- unique(unlist(comparisons))

      if (!is.null(split_by)) {
        # Check if requested comparisons exist in the split column (e.g. "Control", "PCOS")
        if (all(flat_comps %in% levels(plot_data$split))) {
          mode <- "split"
        }
      }

      if (mode == "split") {
        # --- Mode A: Compare Splits (Control vs PCOS) ---
        # Note: stat_compare_means with 'group' aesthetic works best with simple labels, not brackets, for dodged plots
        p <- p + ggpubr::stat_compare_means(
          ggplot2::aes(group = .data$split),
          method = sign_method,
          label = sign_label, # e.g. "****"
          label.y.npc = "top", # Put labels at top
          hide.ns = hide_ns,
          vjust = -0.5
        )
        # Increase Y limit slightly to fit labels
        p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15)))

      } else {
        # --- Mode B: Compare Groups (Endothelial vs Immune) ---
        # Existing logic for X-axis comparisons

        valid_comps <- list()
        counts <- table(plot_data$group)
        for (pair in comparisons) {
          # Ensure groups exist and have cells
          if (all(pair %in% names(counts)) && all(counts[pair] >= 3)) {
            valid_comps <- c(valid_comps, list(pair))
          }
        }

        # Filter Significant
        final_comps <- valid_comps
        if (hide_ns && length(valid_comps) > 0) {
          tryCatch({
            # Simple pre-test
            res <- ggpubr::compare_means(expression ~ group, data = plot_data, method = sign_method)
            sig_res <- res[res$p < 0.05, ]
            keep_idx <- which(sapply(valid_comps, function(pr) {
              any((sig_res$group1 == pr[1] & sig_res$group2 == pr[2]) |
                    (sig_res$group1 == pr[2] & sig_res$group2 == pr[1]))
            }))
            final_comps <- valid_comps[keep_idx]
          }, error = function(e) NULL)
        }

        if (length(final_comps) > 0) {
          p <- p + ggpubr::stat_compare_means(
            comparisons = final_comps,
            method = sign_method,
            label = sign_label,
            vjust = 0.5
          ) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1 + length(final_comps)*0.1)))
        }
      }
    }
    plot_list[[gene]] <- p
  }

  if (is.null(ncol)) ncol <- min(length(plot_list), 3)
  patchwork::wrap_plots(plot_list, ncol = ncol)
}
