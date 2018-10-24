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

  ## Find topological domains using TopDom method
  message("TopDom() ...")
  fit <- TopDom(data, window.size = 5L)
  print(fit)
  str(fit$domain)

  if (requireNamespace("diffobj", quietly = TRUE)) {
    message("TopDom versus published results ...")
    diff <- diffobj::diffPrint(fit$domain, truth,
                               extra = list(row.names = FALSE))
    print(diff)
  }

  ## The largest domain found
  td <- subset(subset(fit$domain, tag == "domain"), size == max(size))
  stopifnot(nrow(td) == 1L)

  data_s <- subsetByRegion(data, region = td, margin = 1/2)
  print(data_s)
  image(data_s)
}
