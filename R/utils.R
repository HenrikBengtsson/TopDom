#' @importFrom utils capture.output
mcat <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "", file = stderr())
}
