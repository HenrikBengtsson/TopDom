if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## Original count data
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  print(data)  ## a TopDomData object

  ## Find topological domains using the TopDom method
  fit <- TopDom(data, window.size = 5L)
  print(fit)  ## a TopDom object

  ## Display the largest domain
  td <- subset(subset(fit$domain, tag == "domain"), size == max(size))
  print(td) ## a data.frame

  data_s <- subsetByRegion(data, region = td, margin = 0.9999)
  print(data_s)  ## a TopDomData object
  
  vp <- grid::viewport(angle = -45, width = 0.7, y = 0.1)
  gg <- ggCountHeatmap(data_s)
  gg <- gg + ggDomain(td, color = "#cccc00") + ggDomainLabel(td)
  print(gg, newpage = TRUE, vp = vp)

  gg <- ggCountHeatmap(data_s, colors = list(mid = "white", high = "black"))
  gg <- gg + ggDomain(td) + ggDomainLabel(td)
  gg <- gg + ggDomain(td, vline = 2)
  print(gg, newpage = TRUE, vp = vp)

  fit_s <- subsetByRegion(fit, region = td, margin = 0.9999)
  dx <- 0.04 * (td$to.coord - td$from.coord)
  print(fit_s)  ## a TopDom object
  for (kk in seq_len(nrow(fit_s$domain))) {
    gg <- gg + ggDomain(fit_s$domain[kk, ], dx = dx * (2 + kk %% 2), color = "red", size = 1)
  }

  print(gg, newpage = TRUE, vp = vp)
}
