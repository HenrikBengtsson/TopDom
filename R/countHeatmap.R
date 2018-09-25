#' Produce a HiC Count Heatmap
#'
#' @param data A TopDomData object.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @seealso See [TopDom] for an example.
#'
#' @export
countHeatmap <- function(data) UseMethod("countHeatmap")

#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_raster coord_fixed theme_void scale_fill_gradient2
#' @export
countHeatmap.TopDomData <- function(data) {
  ## To please R CMD check
  x <- y <- counts <- NULL
  
  upper_tri <- function(x) {
    x[lower.tri(x)] <- NA
    x
  }
  
  dd <- upper_tri(data$counts)
  dd <- melt(dd, varnames = c("x", "y"), na.rm = TRUE, value.name = "counts")

  ## Genomic coordinate system
  mid <- with(data$bins, (from.coord + to.coord) / 2)
  dd$x <- mid[dd$x]
  dd$y <- mid[dd$y]

  gg <- ggplot(dd)
  gg <- gg + aes(x, y, fill = log2(counts))
  gg <- gg + geom_raster(show.legend = FALSE)
  gg <- gg + coord_fixed()
  gg <- gg + theme_void()
  gg <- gg + scale_fill_gradient2(limit = c(0, NA),
                                  mid = "blue", high = "yellow")
  gg
}
