#' @importFrom utils str
#' @export
print.TopDomData <- function(x, ...) {
  cat(sprintf("%s:\n", class(x)))
  cat("bins:\n")
  str(x$bins)
  cat("counts:\n")
  str(x$counts)
}

#' @export
dim.TopDomData <- function(x) {
  dim(x$counts)
}

#' @export
`[.TopDomData` <- function(x, i, ...) {
  structure(list(
    bins = x$bins[i, , drop = FALSE],
    counts = x$counts[i, i, drop = FALSE]
  ), class = "TopDomData")
}

#' @importFrom graphics image
#' @export
image.TopDomData <- function(x, transform = log2, ...) {
  image(transform(x$counts), ...)
}
