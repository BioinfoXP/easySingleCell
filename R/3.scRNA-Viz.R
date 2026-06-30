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
                             max.cutoff = "q95", # \u63A8\u8350\u7528\u5206\u4F4D\u6570\uFF0C\u6BD4\u56FA\u5B9A\u6570\u503C1.5\u66F4\u5065\u58EE
                             cols = c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),
                             ncol = NULL,
                             nrow = NULL,
                             show.legend = FALSE,
                             plot.title = NULL,
                             ...) {

  # \u5185\u90E8\u8F85\u52A9\u51FD\u6570\uFF1A\u5904\u7406\u5355\u4E2A Feature
  create_single_plot <- function(feature) {

    # \u68C0\u67E5 Feature \u662F\u5426\u5B58\u5728 (\u907F\u514D\u62A5\u9519\u4E2D\u65AD)
    if (!feature %in% rownames(scRNA) && !feature %in% colnames(scRNA@meta.data)) {
      warning(paste("Feature", feature, "not found in object."))
      return(NULL)
    }

    # \u57FA\u7840\u7ED8\u56FE (\u6CE8\u610F\uFF1A\u8FD9\u91CC\u628A cols\u7559\u7A7A\uFF0C\u6211\u4EEC\u7528 ggplot \u56FE\u5C42\u624B\u52A8\u52A0\u989C\u8272)
    p <- Seurat::FeaturePlot(
      object = scRNA,
      features = feature,
      reduction = reduction,
      pt.size = pt.size,
      max.cutoff = max.cutoff,
      ...
    )

    # \u5E94\u7528\u81EA\u5B9A\u4E49\u989C\u8272\u6E10\u53D8\u548C\u4E3B\u9898
    p <- p +
      ggplot2::scale_color_gradientn(colors = cols) + # \u5173\u952E\uFF1A\u6B63\u786E\u5E94\u7528\u591A\u8272\u6E10\u53D8
      ggplot2::scale_x_continuous("") +
      ggplot2::scale_y_continuous("") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = if(show.legend) "right" else "none" # \u7075\u6D3B\u63A7\u5236\u56FE\u4F8B
      ) +
      ggplot2::ggtitle(feature) # \u5B50\u56FE\u6807\u9898\u59CB\u7EC8\u4E3A\u57FA\u56E0\u540D

    return(p)
  }

  # \u5FAA\u73AF\u751F\u6210\u56FE\u5F62\u5217\u8868
  plot_list <- lapply(features, create_single_plot)

  # \u79FB\u9664 NULL (\u5373\u6CA1\u627E\u5230\u7684\u57FA\u56E0)
  plot_list <- plot_list[!sapply(plot_list, is.null)]

  if (length(plot_list) == 0) {
    stop("No valid features found to plot.")
  }

  # \u4F7F\u7528 patchwork \u62FC\u56FE
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol, nrow = nrow)

  # \u5982\u679C\u6709\u603B\u6807\u9898\uFF0C\u6DFB\u52A0\u603B\u6807\u9898
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
#' @param aspect.ratio Plot aspect ratio. Default `1`.
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
                         pt.size = 1,
                         stroke = 0,
                         alpha = 1,
                         show.arrow = TRUE,
                         show.axis.title = FALSE,
                         label = FALSE,        # \u65B0\u589E\uFF1A\u662F\u5426\u6807\u8BB0
                         label.size = 4,       # \u65B0\u589E\uFF1A\u5B57\u4F53\u5927\u5C0F
                         label.color = "black",# \u65B0\u589E\uFF1A\u5B57\u4F53\u989C\u8272
                         repel = FALSE,        # \u65B0\u589E\uFF1A\u662F\u5426\u9632\u91CD\u53E0
                         label.box = FALSE,    # \u65B0\u589E\uFF1A\u662F\u5426\u52A0\u80CC\u666F\u6846
                         strip.color = "#e6bac5",
                         legend.position = "right",
                         aspect.ratio = 1,
                         ...) {

  # 1. \u68C0\u67E5 Reduction \u5E76\u83B7\u53D6\u5750\u6807\u8F74\u540D\u79F0
  if (!reduction %in% names(scRNA)) {
    stop(paste("Reduction", reduction, "not found in object."))
  }

  emb <- scRNA[[reduction]]
  dims <- colnames(emb)[1:2]
  key <- emb@key

  # 2. \u51C6\u5907\u7ED8\u56FE\u6570\u636E
  if (is.null(group.by)) {
    group.by <- "ident"
    plot_data <- Seurat::FetchData(scRNA, vars = c(dims, split.by))
    plot_data$ident <- Seurat::Idents(scRNA)
  } else {
    plot_data <- Seurat::FetchData(scRNA, vars = c(dims, group.by, split.by))
  }

  group_col <- if (is.null(group.by)) "ident" else group.by
  plot_data[[group_col]] <- as.factor(plot_data[[group_col]])

  # 3. \u8BBE\u7F6E\u9ED8\u8BA4\u989C\u8272
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

  # 4. \u6784\u5EFA ggplot \u57FA\u7840\u56FE\u5C42
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = dims[1], y = dims[2], color = group_col)) +
    ggplot2::geom_point(size = pt.size, shape = 16, stroke = stroke, alpha = alpha) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = dims[1], y = dims[2])

  # ================== \u65B0\u589E Label \u5904\u7406\u903B\u8F91 ==================
  if (label) {
    # \u8BA1\u7B97\u6BCF\u4E2A\u7FA4\u4F53\u7684\u4E2D\u5FC3\u70B9 (Median)
    # \u5982\u679C\u6709 split.by\uFF0C\u9700\u8981\u6309 group + split \u5206\u7EC4\u8BA1\u7B97
    if (!is.null(split.by)) {
      group_vars <- c(group_col, split.by)
    } else {
      group_vars <- c(group_col)
    }

    # \u4F7F\u7528 base R aggregate \u8BA1\u7B97\u4E2D\u4F4D\u6570\uFF0C\u907F\u514D\u589E\u52A0 dplyr \u4F9D\u8D56
    label_data <- stats::aggregate(
      list(x = plot_data[[dims[1]]], y = plot_data[[dims[2]]]),
      by = plot_data[group_vars],
      FUN = stats::median
    )

    # \u5B9A\u4E49\u7ED8\u56FE\u53C2\u6570
    geom_fun <- if (label.box) {
      if (repel) ggrepel::geom_label_repel else ggplot2::geom_label
    } else {
      if (repel) ggrepel::geom_text_repel else ggplot2::geom_text
    }

    # \u6DFB\u52A0\u6807\u7B7E\u56FE\u5C42
    p <- p + geom_fun(
      data = label_data,
      ggplot2::aes_string(x = "x", y = "y", label = group_col),
      color = label.color,
      size = label.size,
      show.legend = FALSE,
      inherit.aes = FALSE # \u5173\u952E\uFF1A\u4E0D\u7EE7\u627F\u4E3B\u56FE\u7684aes\uFF0C\u9632\u6B62\u62A5\u9519
    )
  }
  # =======================================================

  # 5. \u5E94\u7528\u201C\u5F20\u6CFD\u6C11\u56E2\u961F\u201D\u98CE\u683C\u4E3B\u9898
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

  # 6. \u5904\u7406\u5206\u9762
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", split.by)))
  }

  # 7. \u4F18\u5316\u56FE\u4F8B\u70B9\u7684\u663E\u793A
  p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  return(p)
}



# =============== DotPlot  ================
# =============== 3.scVisDotPlot  ================
#' @title scVis: Universal DotPlot (Vector & List Support)
#' @description A pure ggplot2 implementation of DotPlot that automatically adapts to input type.
#' - If 'features' is a VECTOR: Plots a standard DotPlot (no facets).
#' - If 'features' is a LIST: Plots a grouped DotPlot with facets (supports duplicate genes).
#'
#' @param object A Seurat object.
#' @param features A vector of genes OR a named list of vectors.
#' @param group.by String. Column to group cells by (y-axis). Default: active Idents.
#' @param pal Character vector. Gradient colors.
#' @param rot_angle Numeric. Angle of x-axis text. Default 45 for vector, 90 for list usually looks best.
#' @param font.size Numeric. Base font size.
#' @param scale.min Numeric. Minimum limit for scaled expression (default -2.5).
#' @param scale.max Numeric. Maximum limit for scaled expression (default 2.5).
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A ggplot object.
#' @export
#' @importFrom Seurat FetchData Idents
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn scale_size facet_grid theme_bw element_text element_blank unit element_rect
#' @importFrom dplyr group_by summarise mutate ungroup bind_rows
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' \dontrun{
#'   # Mode 1: Standard Vector
#'   scVisDotPlot(sce, features = c("Ptprc", "Lyz2", "CD14"))
#'
#'   # Mode 2: Grouped List (with facets)
#'   markers <- list("Immune" = c("Ptprc", "Lyz2"), "Myeloid" = c("CD14", "Lyz2"))
#'   scVisDotPlot(sce, features = markers)
#' }
scVisDotPlot <- function(object,
                         features,
                         group.by = NULL,
                         pal = NULL,
                         rot_angle = 45,
                         font.size = 12,
                         scale.min = NA,
                         scale.max = NA,
                         ...) {

  # --- 1. \u5224\u65AD\u6A21\u5F0F (\u5411\u91CF vs \u5217\u8868) ---
  is_grouped <- is.list(features)

  if (is_grouped && is.null(names(features))) {
    stop("\u274C If 'features' is a list, it must be named (e.g., list('GroupA' = c(...))).")
  }

  # --- 2. \u51C6\u5907\u57FA\u56E0\u5217\u8868 ---
  if (is_grouped) {
    all_genes_input <- unique(unlist(features))
  } else {
    all_genes_input <- unique(features)
  }

  # \u68C0\u67E5\u57FA\u56E0\u662F\u5426\u5B58\u5728
  avail_genes <- rownames(object)
  valid_genes <- all_genes_input[all_genes_input %in% avail_genes]
  missing_genes <- setdiff(all_genes_input, avail_genes)

  if (length(missing_genes) > 0) {
    warning(paste("\u26A0\uFE0F Skipped missing genes:", paste(missing_genes, collapse = ", ")), call. = FALSE)
  }
  if (length(valid_genes) == 0) stop("\u274C No valid genes found in the object.")

  # --- 3. \u8BBE\u7F6E\u5206\u7EC4 ---
  if (is.null(group.by)) {
    object$scVis_ident <- Seurat::Idents(object)
    group.by <- "scVis_ident"
  } else if (!group.by %in% colnames(object@meta.data)) {
    stop(paste0("\u274C Group column '", group.by, "' not found."))
  }

  # --- 4. \u8BA1\u7B97\u57FA\u7840\u7EDF\u8BA1\u91CF (Mean & Pct) ---
  # \u8FD9\u4E00\u6B65\u5BF9\u4E24\u79CD\u6A21\u5F0F\u662F\u901A\u7528\u7684
  dat <- Seurat::FetchData(object, vars = c(valid_genes, group.by))

  base_data <- dat |>
    tidyr::pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expr") |>
    dplyr::group_by(CellType = .data[[group.by]], Gene) |>
    dplyr::summarise(
      AvgExp = mean(expm1(.data$Expr)), # \u975ELog\u7A7A\u95F4\u6C42\u5747\u503C
      PctExp = sum(.data$Expr > 0) / dplyr::n() * 100,
      .groups = "drop"
    ) |>
    dplyr::mutate(AvgExp = log1p(.data$AvgExp)) # Log\u56DE\u53BB

  # Z-Score Scaling
  base_data <- base_data |>
    dplyr::group_by(Gene) |>
    dplyr::mutate(AvgExpScaled = as.numeric(scale(.data$AvgExp))) |>
    dplyr::ungroup()

  # \u622A\u65AD\u503C\u5904\u7406
  if (!is.na(scale.min)) base_data$AvgExpScaled[base_data$AvgExpScaled < scale.min] <- scale.min
  if (!is.na(scale.max)) base_data$AvgExpScaled[base_data$AvgExpScaled > scale.max] <- scale.max

  # --- 5. \u6784\u5EFA\u7ED8\u56FE\u6570\u636E (\u5206\u652F\u903B\u8F91) ---

  if (is_grouped) {
    # === \u6A21\u5F0F A: \u5217\u8868 (\u5E26\u5206\u9762\uFF0C\u652F\u6301\u91CD\u590D) ===
    final_df_list <- list()
    for (grp_name in names(features)) {
      grp_genes <- features[[grp_name]]
      grp_genes <- grp_genes[grp_genes %in% valid_genes]
      if (length(grp_genes) == 0) next

      sub_df <- base_data[base_data$Gene %in% grp_genes, ]
      sub_df$FeatureGroup <- grp_name
      # \u5F3A\u5236\u56E0\u5B50\u6C34\u5E73\u4EE5\u4FDD\u6301\u5217\u8868\u987A\u5E8F
      sub_df$Gene <- factor(sub_df$Gene, levels = grp_genes)
      final_df_list[[grp_name]] <- sub_df
    }
    final_plot_data <- dplyr::bind_rows(final_df_list)
    final_plot_data$FeatureGroup <- factor(final_plot_data$FeatureGroup, levels = names(features))

  } else {
    # === \u6A21\u5F0F B: \u5411\u91CF (\u65E0\u5206\u9762\uFF0C\u6807\u51C6\u5C55\u793A) ===
    final_plot_data <- base_data
    # \u5F3A\u5236\u56E0\u5B50\u6C34\u5E73\u4EE5\u4FDD\u6301\u8F93\u5165\u5411\u91CF\u987A\u5E8F
    final_plot_data$Gene <- factor(final_plot_data$Gene, levels = valid_genes)
  }

  # --- 6. \u914D\u8272\u65B9\u6848 ---
  if (is.null(pal)) {
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
    } else {
      pal <- c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B")
    }
  }

  # --- 7. \u7ED8\u56FE (ggplot2) ---
  p <- ggplot2::ggplot(final_plot_data, ggplot2::aes(x = .data$Gene, y = .data$CellType)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$PctExp, color = .data$AvgExpScaled)) +

    # \u6807\u5C3A
    ggplot2::scale_color_gradientn(
      colors = pal,
      guide = ggplot2::guide_colorbar(title = "Mean Exp\n(Scaled)", order = 1)
    ) +
    ggplot2::scale_size_continuous(
      range = c(0, 6),
      labels = c("20", "40", "60", "80", "100"),
      guide = ggplot2::guide_legend(title = "Fraction (%)", order = 2, override.aes = list(shape=21, fill="grey", color="black"))
    ) +

    # \u4E3B\u9898
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = rot_angle, hjust = 1, vjust = if(rot_angle==45) 1 else 0.5,
                                          color = "black", size = font.size),
      axis.text.y = ggplot2::element_text(color = "black", size = font.size),
      axis.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      panel.spacing = ggplot2::unit(0.1, "lines")
    )

  # --- 8. \u4EC5\u5728\u5217\u8868\u6A21\u5F0F\u4E0B\u6DFB\u52A0\u5206\u9762 ---
  if (is_grouped) {
    p <- p +
      ggplot2::facet_grid(cols = ggplot2::vars(.data$FeatureGroup), scales = "free_x", space = "free_x") +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
        strip.text = ggplot2::element_text(face = "bold", size = font.size)
      )
  }

  return(p)
}
# =============== scVisRatioBox  ================
# =============== 4.scVisRatioBox  ================
#' @title scVis: Publication-Ready Cell Proportion Boxplot
#' @description A beautified boxplot for cell type proportions with flexible legend controls and robust statistics.
#' Supports global tests (Kruskal/ANOVA), pairwise tests (Wilcox/T-test), and One-vs-All comparisons.
#'
#' @param sce A Seurat object.
#' @param group_col String. Grouping column (e.g., "Group").
#' @param celltype_col String. Cell type column (e.g., "cluster").
#' @param sample_col String. Sample column (e.g., "sample").
#' @param ref_group String. Optional. Reference group for One-vs-All stats (e.g., "Control").
#' @param comparisons Ignored. Auto-calculated based on groups.
#' @param sign_method String. "wilcox.test" (default) or "t.test".
#' @param sign_label String. "p.signif" (*) or "p.format" (value).
#' @param hide_ns Logical. Whether to hide non-significant results. Default FALSE.
#' @param legend_position String or Vector. e.g., "top", "right", "none", or c(0.8, 0.8).
#' @param palette Vector. Custom colors.
#' @param pt_size Numeric. Jitter point size. Default 1.2.
#' @param box_width Numeric. Box width. Default 0.5.
#' @param base_size Numeric. Font size. Default 14.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData
#' @importFrom dplyr group_by summarise mutate ungroup distinct left_join select
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter scale_fill_manual scale_color_manual scale_y_continuous labs theme_classic theme element_text element_line position_dodge position_jitterdodge expansion element_rect
#' @importFrom ggpubr stat_compare_means
#' @importFrom scales percent hue_pal
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   # 1. Polished style with legend on the right.
#'   scVisRatioBox(sce, "Group", "celltype", "sample", legend_position = "right")
#'
#'   # 2. Specify a reference group, hide NS labels, and put legend on top.
#'   scVisRatioBox(sce, "Group", "celltype", "sample",
#'                 ref_group = "Control", hide_ns = TRUE, legend_position = "top")
#' }
scVisRatioBox <- function(sce,
                          group_col,
                          celltype_col,
                          sample_col,
                          ref_group = NULL,
                          comparisons = NULL,
                          sign_method = "wilcox.test",
                          sign_label = "p.signif",
                          hide_ns = FALSE,
                          legend_position = "top", # <--- \u65B0\u589E\u53C2\u6570
                          palette = NULL,
                          pt_size = 1.2,
                          box_width = 0.5,
                          base_size = 14) {

  # --- 1. \u8F93\u5165\u68C0\u67E5 ---
  if (!inherits(sce, "Seurat")) stop("Input 'sce' must be a Seurat object.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required.")
  if (!all(c(group_col, celltype_col, sample_col) %in% colnames(sce@meta.data))) {
    stop("\u274C Columns not found in meta.data.")
  }

  # --- 2. \u6570\u636E\u51C6\u5907 ---
  meta_df <- Seurat::FetchData(sce, vars = c(group_col, celltype_col, sample_col))
  colnames(meta_df) <- c("Group", "CellType", "Sample")
  meta_df$Group <- as.factor(meta_df$Group)

  if (!is.null(ref_group) && !ref_group %in% levels(meta_df$Group)) {
    stop(paste0("\u274C ref_group '", ref_group, "' not found."))
  }

  sample_info <- meta_df |> dplyr::select("Sample", "Group") |> dplyr::distinct()

  ratio_data <- meta_df |>
    dplyr::group_by(.data$Sample, .data$CellType) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    tidyr::complete(.data$Sample, .data$CellType, fill = list(n = 0)) |>
    dplyr::left_join(sample_info, by = "Sample") |>
    dplyr::group_by(.data$Sample) |>
    dplyr::mutate(Ratio = .data$n / sum(.data$n)) |>
    dplyr::ungroup()

  # --- 3. \u667A\u80FD\u914D\u8272 ---
  groups <- unique(ratio_data$Group)
  n_groups <- length(groups)
  if (is.null(palette)) {
    # \u66F4\u52A0\u67D4\u548C\u4E14\u533A\u5206\u5EA6\u9AD8\u7684\u51FA\u7248\u7EA7\u914D\u8272
    default_pal <- c("#EE6677", "#4477AA", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
    palette <- if (n_groups > length(default_pal)) scales::hue_pal()(n_groups) else default_pal[1:n_groups]
  }

  # --- 4. \u7ED8\u56FE (\u7F8E\u5316\u6838\u5FC3) ---
  p <- ggplot2::ggplot(ratio_data, ggplot2::aes(x = .data$CellType, y = .data$Ratio, fill = .data$Group)) +

    # \u7BB1\u7EBF\u56FE (\u9ED1\u8272\u8FB9\u6846\uFF0C\u589E\u5F3A\u5BF9\u6BD4\u5EA6)
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      width = box_width,
      color = "black", # \u9ED1\u8272\u8FB9\u6846\u66F4\u4E13\u4E1A
      size = 0.4,      # \u8FB9\u6846\u7EBF\u6761\u7C97\u7EC6
      position = ggplot2::position_dodge(width = 0.8)
    ) +

    # \u6563\u70B9 (\u5E26\u9ED1\u8272\u63CF\u8FB9\uFF0C\u589E\u52A0\u8D28\u611F)
    ggplot2::geom_jitter(
      ggplot2::aes(fill = .data$Group), # \u7528 fill \u800C\u4E0D\u662F color\uFF0C\u914D\u5408 shape=21 \u5B9E\u73B0\u63CF\u8FB9\u70B9
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size = pt_size,
      shape = 21,      # \u7A7A\u5FC3\u5706\uFF0C\u53EF\u586B\u5145\u989C\u8272
      color = "black", # \u70B9\u7684\u8FB9\u6846\u989C\u8272
      stroke = 0.3,    # \u70B9\u7684\u8FB9\u6846\u7C97\u7EC6
      alpha = 0.9,
      show.legend = FALSE
    ) +

    # \u6807\u5C3A\u4E0E\u8F74
    ggplot2::scale_fill_manual(values = palette) +
    # \u589E\u52A0\u9876\u90E8 20% \u7559\u767D\uFF0C\u9632\u6B62\u661F\u661F\u663E\u793A\u4E0D\u5168
    ggplot2::scale_y_continuous(labels = scales::percent, expand = ggplot2::expansion(mult = c(0.05, 0.2))) +

    # \u4E3B\u9898\u7F8E\u5316
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(x = NULL, y = "Cell Proportion", fill = group_col) +
    ggplot2::theme(
      # \u5750\u6807\u8F74\u6587\u5B57\u7EAF\u9ED1
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black", face = "plain"),
      axis.text.y = ggplot2::element_text(color = "black", face = "plain"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black"),

      # \u56FE\u4F8B\u63A7\u5236 (\u54CD\u5E94\u53C2\u6570)
      legend.position = legend_position,
      legend.title = ggplot2::element_text(face = "bold"),
      legend.background = ggplot2::element_blank(),

      # \u7F51\u683C\u7EBF (\u6A2A\u5411\u865A\u7EBF\u8F85\u52A9\u9605\u8BFB)
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linetype = "dashed")
    )

  # --- 5. \u7EDF\u8BA1\u68C0\u9A8C ---
  tryCatch({
    if (!is.null(ref_group)) {
      # \u6A21\u5F0F A: \u53C2\u8003\u7EC4\u5BF9\u6BD4 (\u6253\u661F\u53F7)
      p <- p + ggpubr::stat_compare_means(
        ggplot2::aes(group = .data$Group),
        method = sign_method,
        label = sign_label,
        ref.group = ref_group,
        hide.ns = hide_ns,
        label.y.npc = "top",
        vjust = 0.5
      )
    } else {
      # \u6A21\u5F0F B: \u81EA\u52A8\u5224\u65AD (2\u7EC4=Wilcox, >2\u7EC4=Kruskal)
      p <- p + ggpubr::stat_compare_means(
        ggplot2::aes(group = .data$Group),
        method = sign_method,
        label = if (n_groups > 2) "p.format" else sign_label,
        hide.ns = hide_ns,
        label.y.npc = "top",
        vjust = -0.5
      )
    }
  }, error = function(e) message("\u26A0\uFE0F Stats failed: ", e$message))

  return(p)
}

# =============== Ro/e \u5206\u5E03\u504F\u597D  ================
# =============== 5.scVisRoePlot  ================
# https://mp.weixin.qq.com/s/Fhz72Cjbd3zzCi1tB9xP5Q

.scVisRoeLabels <- function(roe_mat, display.mode = c("symbol", "numeric", "none")) {
  display.mode <- match.arg(display.mode)

  if (display.mode == "none") {
    return(FALSE)
  }

  if (display.mode == "numeric") {
    display_mat <- matrix(
      formatC(roe_mat, format = "f", digits = 2),
      nrow = nrow(roe_mat),
      ncol = ncol(roe_mat),
      dimnames = dimnames(roe_mat)
    )
    return(display_mat)
  }

  display_mat <- ifelse(roe_mat >= 2, "+++",
                        ifelse(roe_mat >= 1.5, "++",
                               ifelse(roe_mat > 1.1, "+",
                                      ifelse(roe_mat >= 0.9 & roe_mat <= 1.1, "",
                                             ifelse(roe_mat >= 0.67, "-",
                                                    ifelse(roe_mat >= 0.5, "--", "---"))))))

  matrix(display_mat,
         nrow = nrow(roe_mat),
         ncol = ncol(roe_mat),
         dimnames = dimnames(roe_mat))
}

.scVisRoeHeatmap <- function(roe_mat,
                             title = "Ro/e Tissue Enrichment",
                             display.mode = c("symbol", "numeric", "none"),
                             cluster_rows = TRUE,
                             cluster_cols = FALSE,
                             font.size = 10,
                             font.size.row = font.size,
                             font.size.col = font.size,
                             ...) {
  display.mode <- match.arg(display.mode)

  my_palette <- grDevices::colorRampPalette(c("#483D8B", "#00FFFF", "#F8F8FF", "#FF69B4", "#8B008B"))(100)

  max_abs_deviation <- max(abs(roe_mat - 1), na.rm = TRUE)
  if (max_abs_deviation == 0) max_abs_deviation <- 0.1

  my_breaks <- seq(1 - max_abs_deviation,
                   1 + max_abs_deviation,
                   length.out = 101)

  display_mat <- .scVisRoeLabels(roe_mat, display.mode)

  p <- pheatmap::pheatmap(roe_mat,
                          color = my_palette,
                          breaks = my_breaks,
                          display_numbers = display_mat,
                          fontsize_number = font.size,
                          fontsize_row = font.size.row,
                          fontsize_col = font.size.col,
                          cluster_rows = cluster_rows,
                          cluster_cols = cluster_cols,
                          border_color = "grey90",
                          scale = "none",
                          main = title,
                          ...)
  attr(p, "roe_mat") <- roe_mat
  p
}

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
#'   In "symbol" mode, labels are centered on Ro/e = 1:
#'   "+++"/"++"/"+" indicate increasing enrichment, empty text indicates
#'   approximately expected abundance, and "-"/"--"/"---" indicate depletion.
#'   In "numeric" mode, Ro/e values are printed with two decimal places.
#' @param cluster_rows Logical, whether to cluster rows. Default TRUE.
#' @param cluster_cols Logical, whether to cluster columns. Default FALSE.
#' @param font.size Font size for numbers/symbols. Default 10.
#' @param font.size.row Font size for row names. Defaults to `font.size`.
#' @param font.size.col Font size for column names. Defaults to `font.size`.
#' @param ... Additional arguments passed to \code{pheatmap}.
#'
#' @return A pheatmap object with the Ro/e matrix attached as `attr(x, "roe_mat")`.
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#'
scVisRoePlot <- function(sce,
                         group.by,
                         cell.type = "celltype",
                         sample.by = "orig.ident",
                         title = "Ro/e Tissue Enrichment", # \u65B0\u589E\u6807\u9898\u53C2\u6570
                         method = "chisq",
                         min.rowSum = 0,
                         display.mode = c("symbol", "numeric", "none"),
                         cluster_rows = TRUE,
                         cluster_cols = FALSE,
                         font.size = 10,
                         font.size.row = font.size,
                         font.size.col = font.size,
                         ...) {

  # 1. \u68C0\u67E5\u4F9D\u8D56\u5305 Startrac
  if (!requireNamespace("Startrac", quietly = TRUE)) {
    stop("Package 'Startrac' is required for Ro/e calculation.\nPlease install it using: devtools::install_github('Japrin/STARTRAC')")
  }

  # 2. \u68C0\u67E5\u8F93\u5165\u5217
  if (!all(c(group.by, cell.type, sample.by) %in% colnames(sce@meta.data))) {
    stop("One or more specified columns not found in meta.data.")
  }

  display.mode <- match.arg(display.mode)

  # 3. \u8BA1\u7B97 Ro/e \u77E9\u9635
  message("Calculating Ro/e matrix using Startrac...")
  meta_data <- sce@meta.data

  # \u8C03\u7528 Startrac (\u6CE8\u610F\uFF1A\u76F4\u63A5\u4F20 meta_data \u4F5C\u4E3A\u7B2C\u4E00\u4E2A\u53C2\u6570)
  roe_mat <- Startrac::calTissueDist(meta_data,
                                     byPatient = FALSE,
                                     colname.cluster = cell.type,
                                     colname.patient = sample.by,
                                     colname.tissue = group.by,
                                     method = method,
                                     min.rowSum = min.rowSum)

  .scVisRoeHeatmap(roe_mat,
                   title = title,
                   display.mode = display.mode,
                   cluster_rows = cluster_rows,
                   cluster_cols = cluster_cols,
                   font.size = font.size,
                   font.size.row = font.size.row,
                   font.size.col = font.size.col,
                   ...)
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
    stop(paste0("\u274C Columns '", group_col, "' or '", celltype_col, "' not found in meta.data."))
  }

  # Smart Comparisons Handling: Normalize to a list
  # Logic: If it's not a list, wrap it. If it is a list, check contents.
  if (!is.list(comparisons)) {
    if (length(comparisons) != 2) {
      stop("\u274C For a single comparison, 'comparisons' must be a vector of length 2: c('Case', 'Control').")
    }
    comparisons <- list(comparisons)
  } else {
    # Check if any element in the list is not length 2
    if (any(sapply(comparisons, length) != 2)) {
      stop("\u274C Every element in the 'comparisons' list must be a vector of length 2.")
    }
  }

  # Validate Groups Exist in Data
  available_groups <- unique(sce@meta.data[[group_col]])
  all_requested_groups <- unique(unlist(comparisons))
  missing_groups <- setdiff(all_requested_groups, available_groups)

  if (length(missing_groups) > 0) {
    stop(paste0("\u274C The following groups defined in 'comparisons' are missing from meta.data: ",
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
                          reduction = "umap", # \u901A\u7528\u53C2\u6570
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

  if (length(plot_list) == 0) stop("\u274C No plots generated.")

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
    stop(paste0("\u274C Columns not found: '", group_col, "' or '", celltype_col, "'. Please check meta.data."))
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
    if(length(missing) > 0) warning("\u26A0\uFE0F Groups in 'group_order' missing from data.")
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
      warning("\u26A0\uFE0F More cell types than default 15 colors. Extending palette with hue_pal.")
      my_colors <- scales::hue_pal()(n_types)
    } else {
      my_colors <- ref_pal[1:n_types]
    }
  } else {
    if (is.null(names(palette))) {
      if (length(palette) < n_types) {
        warning("\u26A0\uFE0F Not enough colors provided. Recycling palette.")
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

#' @title scVis: Intelligent VlnPlot (Fully Automated & SCI Styled)
#' @description A high-end, minimalist violin plot builder designed for Seurat V4/V5 data.
#' Features automatic pairwise comparisons, flat significance brackets, and publication-ready
#' SCI styling (Arial font, bold borders, no inner points). Parameters are perfectly aligned
#' with Seurat's native ecosystem.
#'
#' @param sce A Seurat object.
#' @param features Vector of feature names to plot.
#' @param group.by String. Name of the metadata column for X-axis grouping. Defaults to Idents.
#' @param split.by String. Name of the metadata column to split/dodge violins.
#' @param comparisons Vector or List. Specific pairs to compare. If NULL and auto_compare is TRUE, pairs are auto-generated.
#' @param auto_compare Logical. If TRUE, automatically computes all pairwise combinations for the active groups.
#' @param palette Vector. Custom color palette. Defaults to a standard SCI publication palette.
#' @param pt_size Numeric. Size of jitter points. Default is 0 (hidden) for a clean visual.
#' @param legend_position String. Position of the legend ("none", "right", "bottom", "top"). Default is "none".
#' @param ncol Numeric. Number of columns for the patchwork layout.
#' @param sign_method String. Statistical test to use ("wilcox.test" or "t.test"). Default is "wilcox.test".
#' @param sign_label String. Format of significance labels ("p.signif" for stars, "p.format" for p-values).
#' @param hide_ns Logical. If TRUE, hides non-significant comparison brackets.
#' @param assay String. Specific assay to pull data from (e.g., "RNA", "SCT").
#' @param layer String. Specific layer to pull data from in Seurat V5 (e.g., "data", "counts").
#' @param base_size Numeric. Base font size for ggplot theme. Default is 15.
#' @param ... Additional arguments passed to specific internal functions.
#'
#' @return A combined patchwork ggplot object containing the styled violin plots.
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. Basic auto pairwise comparison
#' p <- scVisVlnPlot(sce.endo, features = c("PLVAP", "CLDN5"), group.by = "Group", auto_compare = TRUE)
#'
#' # 2. Split by condition
#' p <- scVisVlnPlot(sce.endo, features = "ANGPT2", group.by = "celltype", split.by = "Condition")
#' }
scVisVlnPlot <- function(sce,
                         features,
                         group.by = NULL,
                         split.by = NULL,
                         comparisons = NULL,
                         auto_compare = FALSE,
                         palette = NULL,
                         pt_size = 0,
                         legend_position = "none",
                         ncol = NULL,
                         sign_method = "wilcox.test",
                         sign_label = "p.signif",
                         hide_ns = FALSE,
                         assay = NULL,
                         layer = NULL,
                         base_size = 15,
                         ...) {

  # 1. \u53C2\u6570\u5BF9\u9F50\u4E0E\u57FA\u7840\u8BBE\u7F6E
  args <- list(...)
  if (is.null(group.by)) group.by <- "scVis_ident"
  if (is.null(sce@meta.data[[group.by]])) sce[[group.by]] <- Seurat::Idents(sce)

  if (is.null(split.by) && "split_by" %in% names(args)) split.by <- args$split_by
  if ("pt.size" %in% names(args)) pt_size <- args$pt.size

  # 2. \u8C03\u8272\u76D8\u903B\u8F91
  color_by <- if (!is.null(split.by)) split.by else group.by
  lvls <- levels(as.factor(sce@meta.data[[color_by]]))
  if (is.null(palette)) {
    palette <- c("Adjacent" = "#D5DBDB", "Tumor-MVI-" = "#5DADE2", "Tumor-MVI+" = "#E74C3C")
    if (!all(names(palette) %in% lvls)) palette <- grDevices::hcl.colors(length(lvls), "Dark 3")
  }

  # 3. \u81EA\u52A8\u5316\u6BD4\u8F83
  if (is.null(comparisons) && auto_compare) {
    comp_lvls <- if (is.null(split.by)) lvls else levels(as.factor(sce@meta.data[[split.by]]))
    if (length(comp_lvls) >= 2) comparisons <- utils::combn(comp_lvls, 2, simplify = FALSE)
  }
  if (!is.null(comparisons) && !is.list(comparisons)) comparisons <- list(comparisons)

  # 4. \u7ED8\u56FE\u5FAA\u73AF
  plot_list <- list()
  for (feature in features) {
    dat <- Seurat::FetchData(sce, vars = c(feature, group.by, split.by), assay = assay, layer = layer)
    colnames(dat)[1:2] <- c("expression", "group")
    if (!is.null(split.by)) colnames(dat)[3] <- "split"

    dat <- dat[is.finite(dat$expression), ]
    dat$group <- factor(dat$group)
    if (!is.null(split.by)) dat$split <- factor(dat$split)

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = group, y = expression))

    if (is.null(split.by)) {
      p <- p + ggplot2::aes(fill = group)
      if (pt_size > 0) p <- p + ggplot2::geom_jitter(width = 0.18, size = pt_size, alpha = 0.45, show.legend = FALSE)
    } else {
      p <- p + ggplot2::aes(fill = split)
      if (pt_size > 0) p <- p + ggplot2::geom_point(position = ggplot2::position_jitterdodge(0.15, dodge.width = 0.85), size = pt_size, alpha = 0.45, show.legend = FALSE)
    }

    p <- p +
      ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.9, linewidth = 0.6, position = ggplot2::position_dodge(0.85)) +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Arial", color = "black"),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(linewidth = 1.2, fill = NA),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = base_size - 1),
        axis.text.y = ggplot2::element_text(face = "bold", size = base_size - 1),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(linewidth = 1),
        axis.ticks.length = ggplot2::unit(0.2, "cm"),
        plot.title = ggplot2::element_text(hjust = 0, size = base_size + 3, margin = ggplot2::margin(b = 10)),
        legend.position = legend_position
      ) +
      ggplot2::labs(title = feature)

    # 5. \u7EDF\u8BA1
    if (!is.null(comparisons)) {
      v_just <- ifelse(sign_label %in% c("p.signif", "p.adj.signif"), 0.4, 0.55)

      if (is.null(split.by)) {
        p <- p + ggpubr::stat_compare_means(comparisons = comparisons, method = sign_method, label = sign_label,
                                            hide.ns = hide_ns, tip.length = 0, step.increase = 0.12, vjust = v_just, size = 4) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1 + length(comparisons) * 0.1)))
      } else {
        p <- p + ggpubr::stat_compare_means(ggplot2::aes(group = split), method = sign_method, label = sign_label,
                                            hide.ns = hide_ns, label.y.npc = "top", vjust = -0.3, size = 4) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25)))
      }
    }
    plot_list[[feature]] <- p
  }

  if (is.null(ncol)) ncol <- min(length(plot_list), 4)
  patchwork::wrap_plots(plot_list, ncol = ncol)
}
