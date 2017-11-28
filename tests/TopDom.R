library("TopDom")
if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## From Supplementary Materials of TopDom article
  pathname <- file.path(path, "mESC_5w_chr10.nij.HindIII.comb.40kb.domain")
  truth <- read.table(pathname, sep = "\t", header = TRUE,
                      colClasses = c("factor", "integer", "numeric", "integer",
                                     "numeric", "factor", "numeric"))
  str(truth)
  
  ## Original count data
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  str(data)

  ## Find topological domains using TopDom method
  fit <- TopDom(data, window.size = 5L)
  str(fit$domain)
}
