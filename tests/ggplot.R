library("TopDom")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
  chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")

  pathname <- file.path(path, sprintf("nij.%s.gz", chr))
  data <- readHiC(pathname, chr = chr, binSize = 40e3)
  data <- data[1:150]

  fit <- TopDom(data, window.size = 5L)
  domains <- subset(fit$domain, tag == "domain")
  td <- domains[which.max(domains$size), ]

  ## ggCountHeatmap() - default arguments
  gg <- ggCountHeatmap(data)
  stopifnot(inherits(gg, "ggplot"))
  b <- ggplot2::ggplot_build(gg)
  stopifnot(is.list(b$data), nrow(b$data[[1]]) > 0L)

  ## ggCountHeatmap() - custom transform and partial 'colors' (fills in defaults)
  gg2 <- ggCountHeatmap(data, transform = function(x) x, colors = c(high = "red"))
  stopifnot(inherits(gg2, "ggplot"))
  invisible(ggplot2::ggplot_build(gg2))

  ## ggDomain() - default 'dx' derived from 'delta'/'vline'
  layer <- ggDomain(td)
  stopifnot(inherits(layer, "Layer") || inherits(layer, "LayerInstance"))
  params <- attr(layer, "gg_params")
  stopifnot(
    is.list(params),
    params$width == td$to.coord - td$from.coord,
    params$dx == params$delta * params$width + params$vline
  )
  gg3 <- gg + layer
  stopifnot(inherits(gg3, "ggplot"))

  ## ggDomain() - explicit 'dx'
  layer2 <- ggDomain(td, dx = 1000, color = "black", size = 1)
  stopifnot(attr(layer2, "gg_params")$dx == 1000)

  ## ggDomainLabel() - returns something addable to a ggplot
  lab <- ggDomainLabel(td)
  gg4 <- gg + lab
  stopifnot(inherits(gg4, "ggplot"))
  invisible(ggplot2::ggplot_build(gg4))

  lab2 <- ggDomainLabel(td, fmt = "%s [%.1f-%.1f]", rot = 0, dx = 10, cex = 1)
  stopifnot(inherits(gg + lab2, "ggplot"))
}
