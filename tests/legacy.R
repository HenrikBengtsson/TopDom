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
  str(fit$domain)

  message("TopDom() v0.0.1 ...")
  fit1 <- TopDom::legacy("0.0.1")$TopDom(data, window.size = 5L)
  str(fit1$domain)

  message("TopDom() v0.0.2 ...")
  fit2 <- TopDom::legacy("0.0.2")$TopDom(data, window.size = 5L)
  str(fit2$domain)

  if (requireNamespace("diffobj", quietly = TRUE)) {
    message("TopDom v0.0.1 versus published results ...")
    diff <- diffobj::diffPrint(fit1$domain, truth,
                               extra = list(row.names = FALSE))
    print(diff)

    message("TopDom v0.0.2 versus published results ...")
    diff <- diffobj::diffPrint(fit2$domain, truth,
                               extra = list(row.names = FALSE))
    print(diff)

    message("TopDom v0.0.2 versus TopDom v0.0.1 results ...")
    diff <- diffobj::diffPrint(fit2$domain, fit1$domain,
                               extra = list(row.names = FALSE))
    print(diff)

    message("TopDom (package version) versus TopDom v0.0.2 results ...")
    diff <- diffobj::diffPrint(fit$domain, fit2$domain,
                               extra = list(row.names = FALSE))
    print(diff)
  }
}
