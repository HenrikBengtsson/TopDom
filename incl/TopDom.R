if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## Original count data
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  print(data)

  ## Find topological domains using TopDom method
  fit <- TopDom(data, window.size = 5L)
  str(fit$domain)

  ## Display the largest domain
  td <- subset(fit$domain, tag == "domain" & size == max(size))
  print(td)

  data_s <- subsetByRegion(data, region = td, margin = 0.9999)

  vp <- grid::viewport(angle = -45, width = 0.7, y = 0.1)
  gg <- ggCountHeatmap(data_s)
  gg <- gg + ggDomain(td, color = "#cccc00") + ggDomainLabel(td)
  print(gg, newpage = TRUE, vp = vp)

  gg <- ggCountHeatmap(data_s, colors = list(mid = "white", high = "black"))
  gg <- gg + ggDomain(td) + ggDomainLabel(td)
  print(gg, newpage = TRUE, vp = vp)
}
