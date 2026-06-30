with_no_openai_key <- function(expr) {
  old_key <- Sys.getenv("OPENAI_API_KEY", unset = NA)
  Sys.unsetenv("OPENAI_API_KEY")
  on.exit({
    if (is.na(old_key)) {
      Sys.unsetenv("OPENAI_API_KEY")
    } else {
      Sys.setenv(OPENAI_API_KEY = old_key)
    }
  }, add = TRUE)
  eval.parent(substitute(expr))
}
