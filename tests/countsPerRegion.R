library("TopDom")

path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")

## Original count data (subset to keep the test fast)
pathname <- file.path(path, sprintf("nij.%s.gz", chr))
data <- readHiC(pathname, chr = chr, binSize = 40e3)
data <- data[1:200]
str(data)

## Find topological domains
fit <- TopDom(data, window.size = 5L)

## Counts per region for the identified domains
domains <- subset(fit$domain, tag == "domain")
counts <- countsPerRegion(data, regions = domains)
str(counts)
stopifnot(
  is.numeric(counts),
  length(counts) == nrow(domains),
  all(counts >= 0),
  !anyNA(counts)
)

## A single region works too
counts1 <- countsPerRegion(data, regions = domains[1, ])
stopifnot(length(counts1) == 1L, counts1 == counts[1])

## Regions spanning more than one chromosome are not allowed
domains2 <- domains[1:2, ]
domains2$chr <- c(chr, "chrX")
res <- tryCatch(countsPerRegion(data, regions = domains2), error = identity)
stopifnot(inherits(res, "error"))

## 'regions' must be a data.frame
res <- tryCatch(countsPerRegion(data, regions = "not-a-data-frame"), error = identity)
stopifnot(inherits(res, "error"))
