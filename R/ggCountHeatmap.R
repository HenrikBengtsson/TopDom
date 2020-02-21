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
#' @importFrom ggplot2 ggplot aes geom_raster coord_fixed theme_void scale_fill_gradient2
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
  gg <- gg + scale_fill_gradient2(limits = c(0, NA), na.value = colors["na"],
                                  mid = colors["mid"], high = colors["high"])
  gg
}


#' Add a Topological Domain to a Count Heatmap
#'
#' @param td A single-row data.frame.
#'
#' @param dx,delta,vline Absolute distance to heatmap.
#'   If `dx = NULL` (default), then `dx = delta * w + vline` where `w` is
#'   the width of the domain.
#'
#' @param size,color The thickness and color of the domain line.
#'
#' @return A [ggplot2::geom_segment] object to be added to the count heatmap.
#'
#' @importFrom ggplot2 geom_segment
#' @export
ggDomain <- function(td, dx = NULL, delta = 0.04, vline = 0, size = 2.0, color = "#666666") {
  x0 <- td$from.coord
  x1 <- td$to.coord
  if (is.null(dx)) dx <- delta * (x1 - x0) + vline
  gg <- geom_segment(aes(x = x0+dx, y = x0-dx, xend = x1+dx, yend = x1-dx),
                     color = color, size = size)
  attr(gg, "gg_params") <- list(x0 = x0, x1 = x1, width = x1 - x0, dx = dx, delta = delta, vline = vline)
  gg		     
}


#' Add a Topological Domain Label to a Count Heatmap
#'
#' @param td A single-row data.frame.
#'
#' @param fmt The [base::sprintf]-format string taking (chromosome, start, stop) as
#'        (string, numeric, numeric) input.
#'
#' @param rot The amount of rotation in \[0,360\] of label.
#'
#' @param dx,vjust The vertical adjustment of the label (relative to rotation)
#'
#' @param cex The scale factor of the label.
#'
#' @return A [ggplot2::ggproto] object to be added to the count heatmap.
#'
#' @importFrom ggplot2 annotation_custom
#' @importFrom grid gpar textGrob
#' @export
ggDomainLabel <- function(td, fmt = "%s: %.2f - %.2f Mbp", rot = 45, dx = 0, vjust = 2.5, cex = 1.5) {
  chr <- td$chr
  x0 <- td$from.coord + dx
  x1 <- td$to.coord + dx
  label <- sprintf(fmt, chr, x0/1e6, x1/1e6)
  grob <- textGrob(label = label, rot = rot, hjust = 0.5, vjust = vjust, gp = gpar(cex = cex))
  annotation_custom(grob = grob, ymin = x0, ymax = x1, xmin = x0, xmax = x1)
}
