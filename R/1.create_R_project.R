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
  # 1. \u4E25\u683C\u53C2\u6570\u68C0\u67E5 (\u501F\u9274 Internal.R \u7684\u6821\u9A8C\u98CE\u683C)
  if (missing(ProjectName) || !is.character(ProjectName) || nchar(trimws(ProjectName)) == 0) {
    stop("Error: 'ProjectName' must be a non-empty character string.")
  }

  # \u6E05\u7406\u9879\u76EE\u540D\u79F0\uFF08\u79FB\u9664\u9996\u5C3E\u7A7A\u683C\uFF09
  ProjectName <- trimws(ProjectName)
  full_project_path <- file.path(BasePath, ProjectName)

  # 2. \u9632\u6B62\u8986\u76D6\u73B0\u6709\u6570\u636E
  if (dir.exists(full_project_path)) {
    stop(sprintf("Creation aborted: Directory '%s' already exists.", full_project_path))
  }

  # 3. \u5B9A\u4E49\u5B50\u76EE\u5F55 (\u4FEE\u6B63\u4E86\u539F\u4EE3\u7801\u4E2D\u7684 'resourse' \u62FC\u5199\u9519\u8BEF)
  # \u589E\u52A0 'raw_data' \u548C 'temp' \u76EE\u5F55\u662F\u751F\u4FE1\u5206\u6790\u7684\u63A8\u8350\u505A\u6CD5
  subdirs <- c("script", "input_data", "output_figure", "output_data", "other_resource", "temp")

  # 4. \u9012\u5F52\u521B\u5EFA\u76EE\u5F55 (recursive = TRUE \u786E\u4FDD\u7236\u76EE\u5F55\u4E0D\u5B58\u5728\u65F6\u4E5F\u80FD\u6210\u529F)
  message(sprintf("Initializing project at: %s", full_project_path))
  dir.create(full_project_path, recursive = TRUE, showWarnings = FALSE)

  # \u6279\u91CF\u521B\u5EFA\u5B50\u76EE\u5F55
  for (subdir in subdirs) {
    target_dir <- file.path(full_project_path, subdir)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 5. \u5B9A\u4E49\u5E76\u5199\u5165 R project \u6587\u4EF6\u5185\u5BB9
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

  # 6. \u53EF\u9009\uFF1A\u521D\u59CB\u5316\u4E00\u4E2A\u7B80\u5355\u7684 .gitignore (\u5982\u679C\u7528\u6237\u4F7F\u7528 Git)
  gitignore_content <- c(".Rhistory", ".RData", ".Rproj.user/", "temp/", "input_data/large_files/")
  writeLines(gitignore_content, con = file.path(full_project_path, ".gitignore"))

  # \u4F7F\u7528\u89C4\u8303\u7684 message \u8F93\u51FA
  message(sprintf("Project '%s' created successfully with %d subdirectories.", ProjectName, length(subdirs)))

  return(invisible(NULL))
}
