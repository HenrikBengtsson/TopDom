#' @importFrom utils str
#' @export
print.TopDomData <- function(x, ...) {
  cat(sprintf("%s:\n", class(x)))
  cat("bins:\n")
  str(x$bins)
  cat("counts:\n")
  str(x$counts)
}
