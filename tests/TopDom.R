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

  ## Compare to TopDom v0.0.1.
  fit1 <- TopDom:::TopDom_0.0.1(data, window.size = 5L)
  str(fit1$domain)
  
  ## Comparing fit with published results
  ## WORKAROUND: https://github.com/brodieG/diffobj/issues/111
  rownames(truth) <- NULL
  rownames(fit1$domain) <- NULL
  ## diff <- diffPrint(fit1$domain, truth)
  ## print(diff)
  ## < fit$domain                                                   
  ## > truth                                                        
  ## @@ 241,6 / 241,5 @@                                            
  ## ~       chr from.id from.coord to.id  to.coord      tag    size
  ##   240 chr10    3188  127480000  3194 127760000   domain  280000
  ##   241 chr10    3195  127760000  3200 128000000   domain  240000
  ## < 242 chr10    3201  128000000  3203 128120000 boundary  120000
  ## < 243 chr10    3204  128120000  3211 128440000 boundary  320000
  ## > 242 chr10    3201  128000000  3211 128440000 boundary  440000
  ##   244 chr10    3212  128440000  3241 129640000   domain 1200000
  ##   245 chr10    3242  129640000  3250 130000000   domain  360000

  ## Note that TD chr10:128000000-128440000 in Suppl. Mat. is split
  ## up in two in the fitted results.
}
