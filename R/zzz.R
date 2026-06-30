#' easySingleCell: Seurat-style single-cell, bulk, GEO, DepMap, and model-assisted helpers
#'
#' @description
#' Streamlines bioinformatics analysis for bulk RNA-seq, scRNA-seq, GEO
#' metadata screening, DepMap visualization, gene ID conversion, and
#' function-level model-assisted annotation or metadata screening workflows.
#'
#' @keywords internal
#' @importFrom methods new
#' @importFrom stats as.formula model.matrix na.omit
#' @importFrom utils assignInNamespace capture.output head read.table tail write.csv
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "AI_SeqType", "AvgExp", "AvgExpScaled", "CellType", "Count", "Expr",
    "FeatureGroup", "Frequency", "Gene", "ModelID", "PctExp", "Ratio",
    "Sample", "Score", "YVal", "avg_log2FC", "avg_logFC", "cell_id",
    "celltype", "celltype_internal", "cleaned_id", "cluster", "copykat.pred",
    "dispersion_empirical", "dispersion_fit", "feature", "freq", "freq1",
    "freq2", "gene", "group", "group_internal", "mean_expression", "n",
    "p.signif", "p_val_adj", "raw_id", "sample_id", "tmp_sort_N", "x",
    "xmax", "xmin", "y", "y.position"
  ))
}
