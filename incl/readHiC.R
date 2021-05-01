if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## Original count data
  chr <- "chr10"
  pathname <- file.path(path, sprintf("nij.%s.gz", chr))
  data <- readHiC(pathname, chr = chr, binSize = 40e3)
  print(data)
  str(data)
}
