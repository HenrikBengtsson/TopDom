library("TopDom")

path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")
pathname <- file.path(path, sprintf("nij.%s.gz", chr))

data <- readHiC(pathname, chr = chr, binSize = 40e3)
data <- data[1:200]

## Reference fit -----------------------------------------------------------
fit <- TopDom(data, window.size = 5L)
stopifnot(inherits(fit, "TopDom"))
params <- attr(fit, "parameters")
stopifnot(params$window.size == 5L, isTRUE(params$statFilter))

## print.TopDom() --------------------------------------------------------
out <- capture.output(print(fit))
stopifnot(length(out) > 0L)

## dim.TopDom() and [.TopDom -------------------------------------------
stopifnot(identical(dim(fit), dim(fit$domain)))
fit_sub <- fit[1:3, ]
stopifnot(inherits(fit_sub, "TopDom"), nrow(fit_sub$domain) == 3L)

## print.TopDom() when no parameters are attached
fit_noparam <- fit
attr(fit_noparam, "parameters") <- NULL
out <- capture.output(print(fit_noparam))
stopifnot(any(grepl("N/A", out)))

## A larger window size yields 'boundary' domains, exercising the
## boundary-collapsing code path in Convert.Bin.To.Domain.TMP()
fit_w7 <- TopDom(data, window.size = 7L)
stopifnot("boundary" %in% fit_w7$domain$tag)
fit_w7_nf <- TopDom(data, window.size = 7L, statFilter = FALSE)
stopifnot(inherits(fit_w7_nf, "TopDom"))

## statFilter = FALSE branch -------------------------------------------
fit_nf <- TopDom(data, window.size = 5L, statFilter = FALSE)
stopifnot(
  inherits(fit_nf, "TopDom"),
  isFALSE(attr(fit_nf, "parameters")$statFilter),
  all(fit_nf$binSignal$pvalue == 0)
)

## Passing a file path as 'data' (dispatches to readHiC() internally) ---
fit_chr <- TopDom(pathname, window.size = 5L, chr = chr, binSize = 40e3,
                  bins = 1:200)
stopifnot(inherits(fit_chr, "TopDom"))
stopifnot(identical(fit_chr$domain, fit$domain))

## debug = TRUE exercises the messaging code paths --------------------
msg <- capture.output(
  fit_dbg <- TopDom(data, window.size = 5L, debug = TRUE),
  type = "message"
)
stopifnot(inherits(fit_dbg, "TopDom"), any(grepl("Step 1", msg)))
stopifnot(identical(fit_dbg$domain, fit$domain))

## outFile = <prefix> writes .binSignal, .domain and .bed files -------
tmp <- tempfile("TopDom-")
fit_out <- TopDom(data, window.size = 5L, outFile = tmp)
for (ext in c(".binSignal", ".domain", ".bed")) {
  f <- paste0(tmp, ext)
  stopifnot(file_test("-f", f), file.size(f) > 0L)
}
bed <- read.table(paste0(tmp, ".bed"), sep = "\t", header = FALSE)
stopifnot(nrow(bed) == nrow(fit_out$bed))
unlink(paste0(tmp, c(".binSignal", ".domain", ".bed")))

## Invalid arguments -----------------------------------------------------
res <- tryCatch(TopDom(data, window.size = -1L), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(TopDom(data, window.size = 5L, debug = NA), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(TopDom(42, window.size = 5L), error = identity)
stopifnot(inherits(res, "error"))
