#' Produce a Count Heatmap
#'
#' @param data A TopDomData object.
#'
#' @param transform A function applied to the counts prior to generating
#'        heatmap colors.
#'
#' @param colors A named list to control to color scale.
#'
#' @param \ldots Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @author Henrik Bengtsson.
#'
#' @seealso See [TopDom] for an example.
#'
#' @export
ggCountHeatmap <- function(data, transform, colors, ...) UseMethod("ggCountHeatmap")

#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_raster geom_segment coord_fixed theme_void scale_fill_gradient2
#' @export
ggCountHeatmap.TopDomData <- function(data, transform = function(x) log2(x + 1), colors = c(na = "white", mid = "blue", high = "yellow"), ...) {
  ## To please R CMD check
  x <- y <- counts <- NULL
  
  keepUpperTriangle <- function(x) {
    x[lower.tri(x)] <- NA
    x
  }

  cols <- eval(formals(ggCountHeatmap.TopDomData)$colors)
  keys <- names(cols)
  keys <- keys[is.na(colors[keys])]
  colors[keys] <- cols[keys]
    
  dd <- keepUpperTriangle(data$counts)
  dd <- melt(dd, varnames = c("x", "y"), na.rm = TRUE, value.name = "counts")

  ## Genomic coordinate system
  mid <- with(data$bins, (from.coord + to.coord) / 2)
  dd$x <- mid[dd$x]
  dd$y <- mid[dd$y]

  gg <- ggplot(dd)
  gg <- gg + aes(x, y, fill = transform(counts))
  gg <- gg + geom_raster(show.legend = FALSE)
  gg <- gg + coord_fixed()
  gg <- gg + theme_void()
  gg <- gg + scale_fill_gradient2(limit = c(0, NA), na.value = colors["na"],
                                  mid = colors["mid"], high = colors["high"])
  gg
}


#' Add a Topological Domain to a Count Heatmap
#'
#' @param td A single-row TopDomData object.
#'
#' @param delta Relative distance to heatmap.
#'
#' @param size,color The thickness and color of the domain line.
#'
#' @return A [ggplot2::geom_segment] object to be added to the count heatmap.
#'
#' @export
ggDomain <- function(td, delta = 0.04, size = 2.0, color = "#666666") {
  x0 <- td$from.coord
  x1 <- td$to.coord
  dx <- delta * (x1 - x0)
  geom_segment(aes(x = x0+dx, y = x0-dx, xend = x1+dx, yend = x1-dx),
               color = color, size = size)
}
