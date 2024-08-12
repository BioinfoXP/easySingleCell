#' Create Project Directory Structure
#'
#' This function creates a directory structure for a new project, including subdirectories for scripts, input data, output figures, output data, and other resources. It also generates a basic R project file.
#'
#' @param ProjectName A character string specifying the name of the project.
#'
#' @return NULL. The function is called for its side effects.
#' @export
#'
#' @examples
#' creat_project(ProjectName = "ATAC")
creat_project <- function(ProjectName) {
  if (missing(ProjectName) || !is.character(ProjectName) || nchar(ProjectName) == 0) {
    stop("Please provide a valid project name as a non-empty character string.")
  }

  # Define the subdirectories to create
  subdirs <- c("script", "input_data", "output_figure", "output_data", "other_resourse")

  # Create the main project directory
  dir.create(ProjectName, showWarnings = FALSE)

  # Create subdirectories
  for (subdir in subdirs) {
    dir.create(file.path(ProjectName, subdir), showWarnings = FALSE)
  }

  # Define the content of the R project file
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

  # Write the R project file
  writeLines(rproj_content, con = file.path(ProjectName, paste0(ProjectName, ".Rproj")))

  message("Project '", ProjectName, "' created successfully with the necessary subdirectories and R project file.")
}

# Example usage
# creat_project(ProjectName = "ATAC")
