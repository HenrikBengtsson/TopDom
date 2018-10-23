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

  ## Subset TopDomData object
  data_s <- subsetByRegion(data, region = td, margin = 0.9999)
  print(data_s)  ## a TopDomData object
  
  vp <- grid::viewport(angle = -45, width = 0.7, y = 0.3)
  gg <- ggCountHeatmap(data_s)
  gg <- gg + ggDomain(td, color = "#cccc00") + ggDomainLabel(td)
  print(gg, newpage = TRUE, vp = vp)

  gg <- ggCountHeatmap(data_s, colors = list(mid = "white", high = "black"))
  gg_td <- ggDomain(td, delta = 0.08)
  dx <- attr(gg_td, "gg_params")$dx
  gg <- gg + gg_td + ggDomainLabel(td, vjust = 2.5)
  print(gg, newpage = TRUE, vp = vp)

  ## Subset TopDom object
  fit_s <- subsetByRegion(fit, region = td, margin = 0.9999)
  print(fit_s)  ## a TopDom object
  for (kk in seq_len(nrow(fit_s$domain))) {
    gg <- gg + ggDomain(fit_s$domain[kk, ], dx = dx * (4 + kk %% 2), color = "red", size = 1)
  }

  print(gg, newpage = TRUE, vp = vp)


  gg <- ggCountHeatmap(data_s)
  gg_td <- ggDomain(td, delta = 0.08)
  dx <- attr(gg_td, "gg_params")$dx
  gg <- gg + gg_td + ggDomainLabel(td, vjust = 2.5)
  fit_s <- subsetByRegion(fit, region = td, margin = 0.9999)
  for (kk in seq_len(nrow(fit_s$domain))) {
    gg <- gg + ggDomain(fit_s$domain[kk, ], dx = dx * (4 + kk %% 2), color = "blue", size = 1)
  }

  print(gg, newpage = TRUE, vp = vp)
}
