#' Generate Pairwise Comparisons from Group Vector
#'
#' This function takes a vector of group labels and generates all pairwise comparisons
#' between the unique groups. It is useful for statistical testing or plotting purposes
#' where pairwise comparisons are required.
#'
#' @param group_vector A vector of group labels. The function will identify unique groups
#' and generate all possible pairwise comparisons.
#'
#' @return A list of pairwise comparisons, where each element is a vector of length two
#' representing a comparison between two groups.
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   group_vector <- c("A", "B", "C", "A", "B", "C")
#'   comparisons <- generate_comparisons(group_vector)
#'   print(comparisons)
#' }
#'
#' @export
generate_comparisons <- function(group_vector) {
  # 将 group_vector 转换为字符向量并提取唯一的组
  unique_groups <- unique(as.character(group_vector))

  # 生成两两组合
  comparisons <- combn(unique_groups, 2, simplify = FALSE)

  return(comparisons)
}
