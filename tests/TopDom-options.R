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
out <- utils::capture.output(print(fit))
stopifnot(length(out) > 0L)

## dim.TopDom() and [.TopDom -------------------------------------------
stopifnot(identical(dim(fit), dim(fit$domain)))
fit_sub <- fit[1:3, ]
stopifnot(inherits(fit_sub, "TopDom"), nrow(fit_sub$domain) == 3L)

## print.TopDom() when no parameters are attached
fit_noparam <- fit
attr(fit_noparam, "parameters") <- NULL
out <- utils::capture.output(print(fit_noparam))
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
msg <- utils::capture.output(
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
  stopifnot(utils::file_test("-f", f), file.size(f) > 0L)
}
bed <- utils::read.table(paste0(tmp, ".bed"), sep = "\t", header = FALSE)
stopifnot(nrow(bed) == nrow(fit_out$bed))
unlink(paste0(tmp, c(".binSignal", ".domain", ".bed")))

## Corner case: input with no processable region -----------------------
## When every non-gap run is shorter than the minimum region size, the
## internal 'idxs' loop variable is never created.  TopDom() used to warn
## "object 'idxs' not found" from rm(list = "idxs"); make sure it does not.
n <- 10L
z <- matrix(0, nrow = n, ncol = n)
diag(z)[c(1L, n)] <- 5  ## isolated signal, no run of >= 3 bins
edge <- structure(list(
  bins = data.frame(id = seq_len(n), chr = "chr1",
                    from.coord = (0:(n - 1)) * 40000L,
                    to.coord = (1:n) * 40000L,
                    stringsAsFactors = FALSE),
  counts = z
), class = "TopDomData")
warns <- NULL
fit_edge <- withCallingHandlers(
  TopDom(edge, window.size = 3L, statFilter = FALSE),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(inherits(fit_edge, "TopDom"))
stopifnot(!any(grepl("object 'idxs' not found", warns)))

## Corner case: extremely sparse data and statFilter = TRUE -----------
## Standardizing an all-zero near-diagonal band yields non-finite values
## that used to surface as an obscure wilcox.test() error
## <https://github.com/HenrikBengtsson/TopDom/issues/17>.  TopDom() should
## now fail early with an informative message that mentions statFilter.
m <- 40L
sparse <- structure(list(
  bins = data.frame(id = seq_len(m), chr = "22",
                    from.coord = (seq_len(m) - 1L) * 500000L,
                    to.coord = seq_len(m) * 500000L,
                    stringsAsFactors = FALSE),
  counts = { z <- matrix(0, nrow = m, ncol = m); diag(z) <- 5; z }
), class = "TopDomData")
res <- tryCatch(TopDom(sparse, window.size = 5L), error = identity)
stopifnot(
  inherits(res, "error"),
  grepl("statFilter", conditionMessage(res)),
  grepl("sparse", conditionMessage(res))
)

## Invalid arguments -----------------------------------------------------
res <- tryCatch(TopDom(data, window.size = -1L), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(TopDom(data, window.size = 5L, debug = NA), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(TopDom(42, window.size = 5L), error = identity)
stopifnot(inherits(res, "error"))
