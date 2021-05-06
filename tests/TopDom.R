library("TopDom")

path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")

## From Supplementary Materials of TopDom article
## Source: http://zhoulab.usc.edu/TopDom/topdom_sup.htm
message("Loading truth ...")
pathname <- file.path(path, sprintf("mESC_5w_%s.nij.HindIII.comb.40kb.domain", chr))
truth <- read.table(pathname, sep = "\t", header = TRUE,
                    colClasses = c("factor", "integer", "numeric", "integer",
                                   "numeric", "factor", "numeric"))
str(truth)

## Original count data
message("Loading count data ...")
pathname <- file.path(path, sprintf("nij.%s.gz", chr))
data <- readHiC(pathname, chr = chr, binSize = 40e3)
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

