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
  
  data_s <- subsetByRegion(data, region = td, margin = 1/2)
  image(data_s)
}
