if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## Original count data
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  print(data)

  ## Find topological domains using TopDom method for two window sizes
  tds_5 <- TopDom(data, window.size = 5L)
  tds_6 <- TopDom(data, window.size = 6L)

  ## Overlap scores (in both directions)
  overlap_56 <- overlapScores(tds_6, reference = tds_5)
  print(overlap_56)
  
  overlap_65 <- overlapScores(tds_5, reference = tds_6)
  print(overlap_65)
}
