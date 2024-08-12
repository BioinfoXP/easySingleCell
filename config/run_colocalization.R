run_colocalization <- function(slide, assay, useful_features, out_label, misty_out_alias) {
  # 定义每个视图的assay
  view_assays <- list(
    "main" = assay,
    "juxta" = assay,
    "para" = assay
  )

  # 定义每个视图的特征
  view_features <- list(
    "main" = useful_features,
    "juxta" = useful_features,
    "para" = useful_features
  )

  # 定义每个视图的空间上下文
  view_types <- list(
    "main" = "intra",
    "juxta" = "juxta",
    "para" = "para"
  )

  # 定义额外的参数
  view_params <- list(
    "main" = NULL,
    "juxta" = 2,
    "para" = 5
  )

  # 输出路径
  misty_out <- paste0(misty_out_alias, out_label, "_", assay)

  # 运行Misty Seurat分析
  run_misty_seurat(
    visium.slide = slide,
    view.assays = view_assays,
    view.features = view_features,
    view.types = view_types,
    view.params = view_params,
    spot.ids = NULL,
    out.alias = misty_out
  )

  return(misty_out)
}
