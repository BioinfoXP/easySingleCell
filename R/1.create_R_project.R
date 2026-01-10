#' Create Project Directory Structure
#'
#' This function initializes a standardized directory structure for a new bioinformatics project,
#' including subdirectories for scripts, data, and results. It also generates an RStudio
#' project (.Rproj) file to facilitate structured analysis.
#'
#' @param ProjectName A character string specifying the name of the project.
#' @param BasePath A character string specifying the root directory where the project
#' will be created. Defaults to the current working directory (".").
#'
#' @return Invisible NULL. The function is called for its side effects (creating directories and files).
#' @export
#'
#' @examples
#' \dontrun{
#' creat_project(ProjectName = "scRNA_Study")
#' }
creat_project <- function(ProjectName, BasePath = ".") {
  # 1. 严格参数检查 (借鉴 Internal.R 的校验风格)
  if (missing(ProjectName) || !is.character(ProjectName) || nchar(trimws(ProjectName)) == 0) {
    stop("Error: 'ProjectName' must be a non-empty character string.")
  }

  # 清理项目名称（移除首尾空格）
  ProjectName <- trimws(ProjectName)
  full_project_path <- file.path(BasePath, ProjectName)

  # 2. 防止覆盖现有数据
  if (dir.exists(full_project_path)) {
    stop(sprintf("Creation aborted: Directory '%s' already exists.", full_project_path))
  }

  # 3. 定义子目录 (修正了原代码中的 'resourse' 拼写错误)
  # 增加 'raw_data' 和 'temp' 目录是生信分析的推荐做法
  subdirs <- c("script", "input_data", "output_figure", "output_data", "other_resource", "temp")

  # 4. 递归创建目录 (recursive = TRUE 确保父目录不存在时也能成功)
  message(sprintf("Initializing project at: %s", full_project_path))
  dir.create(full_project_path, recursive = TRUE, showWarnings = FALSE)

  # 批量创建子目录
  for (subdir in subdirs) {
    target_dir <- file.path(full_project_path, subdir)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 5. 定义并写入 R project 文件内容
  rproj_content <- c(
    "Version: 1.0",
    "RestoreWorkspace: Default",
    "SaveWorkspace: Default",
    "AlwaysSaveHistory: Default",
    "EnableCodeIndexing: Yes",
    "UseSpacesForTab: Yes",
    "NumSpacesForTab: 2",
    "Encoding: UTF-8",
    "RnwWeave: Sweave",
    "LaTeX: pdfLaTeX"
  )

  rproj_file_path <- file.path(full_project_path, paste0(basename(ProjectName), ".Rproj"))
  writeLines(rproj_content, con = rproj_file_path)

  # 6. 可选：初始化一个简单的 .gitignore (如果用户使用 Git)
  gitignore_content <- c(".Rhistory", ".RData", ".Rproj.user/", "temp/", "input_data/large_files/")
  writeLines(gitignore_content, con = file.path(full_project_path, ".gitignore"))

  # 使用规范的 message 输出
  message(sprintf("Project '%s' created successfully with %d subdirectories.", ProjectName, length(subdirs)))

  return(invisible(NULL))
}
