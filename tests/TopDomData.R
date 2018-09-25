library("TopDom")
if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## From Supplementary Materials of TopDom article
  ## Source: http://zhoulab.usc.edu/TopDom/topdom_sup.htm
  message("Loading truth ...")
  pathname <- file.path(path, "mESC_5w_chr10.nij.HindIII.comb.40kb.domain")
  truth <- read.table(pathname, sep = "\t", header = TRUE,
                      colClasses = c("factor", "integer", "numeric", "integer",
                                     "numeric", "factor", "numeric"))
  str(truth)
  
  ## Original count data
  message("Loading count data ...")
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  str(data)

  ## print()
  print(data)
  
  ## dim()
  print(dim(data))

  ## Subsetting via [()
  data_s <- data[101:200]
  str(data_s)
  stopifnot(identical(dim(data_s), c(100L, 100L)))

  ## image()
  image(data_s)
}
