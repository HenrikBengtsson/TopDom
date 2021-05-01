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

