library("TopDom")

path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")

pathname <- file.path(path, sprintf("nij.%s.gz", chr))
data <- readHiC(pathname, chr = chr, binSize = 40e3)
data <- data[1:200]

fit <- TopDom(data, window.size = 5L)
domains <- subset(fit$domain, tag == "domain")
region <- domains[which.max(domains$size), ]

## (a) Subsetting a TopDomData object -----------------------------------------

## Relative margin (margin < 1)
sub_rel <- subsetByRegion(data, region = region, margin = 1/2)
print(sub_rel)
stopifnot(inherits(sub_rel, "TopDomData"))

## Absolute margin (margin >= 1)
sub_abs <- subsetByRegion(data, region = region, margin = 80e3)
stopifnot(inherits(sub_abs, "TopDomData"))

## Zero margin gives exactly the region
sub_zero <- subsetByRegion(data, region = region, margin = 0)
stopifnot(
  inherits(sub_zero, "TopDomData"),
  all(sub_zero$bins$from.coord >= region$from.coord),
  all(sub_zero$bins$to.coord <= region$to.coord),
  nrow(sub_zero$bins) <= nrow(sub_rel$bins)
)

## Negative margins are rejected
res <- tryCatch(subsetByRegion(data, region = region, margin = -1), error = identity)
stopifnot(inherits(res, "error"))

## 'region' must be a data.frame
res <- tryCatch(subsetByRegion(data, region = 42), error = identity)
stopifnot(inherits(res, "error"))

## (b) Subsetting a TopDom object -------------------------------------------

sub_td <- subsetByRegion(fit, region = region, margin = 1/2)
str(sub_td$domain)
stopifnot(
  inherits(sub_td, "TopDom"),
  is.data.frame(sub_td$domain),
  nrow(sub_td$domain) >= 1L
)

sub_td0 <- subsetByRegion(fit, region = region, margin = 0)
stopifnot(inherits(sub_td0, "TopDom"), nrow(sub_td0$domain) >= 1L)

sub_td_abs <- subsetByRegion(fit, region = region, margin = 200e3)
stopifnot(nrow(sub_td_abs$domain) >= nrow(sub_td0$domain))

res <- tryCatch(subsetByRegion(fit, region = "x"), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(subsetByRegion(fit, region = region, margin = -0.1), error = identity)
stopifnot(inherits(res, "error"))
