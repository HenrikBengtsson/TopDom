library("TopDom")

options(TopDom.debug = FALSE)

## A tiny synthetic contact matrix -------------------------------------------
n <- 6L
set.seed(42L)
counts <- matrix(as.double(sample(0:9, size = n * n, replace = TRUE)), nrow = n)
counts[lower.tri(counts)] <- t(counts)[lower.tri(counts)]  ## symmetric
from <- seq(from = 0L, by = 40000L, length.out = n)
to <- from + 40000L

tmpdir <- tempfile("readHiC-")
dir.create(tmpdir)

write_ws <- function(df, file) {
  write.table(df, file = file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  file
}

## (1) N-by-(N+3): chr, from.coord, to.coord, <counts> ----------------------
f3 <- write_ws(data.frame(chr = "chr1", from = from, to = to, counts),
               file.path(tmpdir, "mat3.txt"))
d3 <- readHiC(f3)
str(d3)
stopifnot(
  inherits(d3, "TopDomData"),
  identical(dim(d3), c(n, n)),
  identical(d3$bins$id, seq_len(n)),
  all(d3$bins$chr == "chr1"),
  is.numeric(d3$counts),
  all(d3$counts == counts)
)

## (2) N-by-(N+4): id, chr, from.coord, to.coord, <counts> ----------------
f4 <- write_ws(data.frame(id = seq_len(n), chr = "chr1", from = from, to = to, counts),
               file.path(tmpdir, "mat4.txt"))
d4 <- readHiC(f4)
stopifnot(inherits(d4, "TopDomData"), identical(dim(d4), c(n, n)))
stopifnot(all(d4$counts == counts))

## (3) Unknown format -> error -------------------------------------------
f5 <- write_ws(data.frame(a = 1:n, b = 1:n, chr = "chr1", from = from, to = to, counts),
               file.path(tmpdir, "mat5.txt"))
res <- tryCatch(readHiC(f5), error = identity)
stopifnot(inherits(res, "error"), grepl("Unknown format", conditionMessage(res)))

## (4) Non-existing file -> error --------------------------------------
res <- tryCatch(readHiC(file.path(tmpdir, "does-not-exist.txt")), error = identity)
stopifnot(inherits(res, "error"))

## (5) Pure N-by-N count matrix + (chr, binSize) ------------------------
fm <- write_ws(as.data.frame(counts), file.path(tmpdir, "matN.txt"))
dm <- readHiC(fm, chr = "chrZ", binSize = 40e3)
stopifnot(
  inherits(dm, "TopDomData"),
  identical(dim(dm), c(n, n)),
  all(dm$bins$chr == "chrZ"),
  all(dm$bins$from.coord == from),
  all(dm$bins$to.coord == to)
)

## ... with debug output enabled
dm_dbg <- readHiC(fm, chr = "chrZ", binSize = 40e3,
                  debug = TRUE)
stopifnot(identical(dim(dm_dbg), c(n, n)))

## (6) 'bins' subsetting of a pure count matrix -----------------------
keep <- c(1L, 3L, 5L)
db <- readHiC(fm, chr = "chrZ", binSize = 40e3, bins = keep)
stopifnot(
  identical(dim(db), c(length(keep), length(keep))),
  identical(db$bins$id, keep),
  all(db$counts == counts[keep, keep])
)

## Out-of-range and non-positive 'bins' -> error
res <- tryCatch(readHiC(fm, chr = "chrZ", binSize = 40e3, bins = c(1L, 999L)),
                error = identity)
stopifnot(inherits(res, "error"), grepl("out of range", conditionMessage(res)))

res <- tryCatch(readHiC(fm, chr = "chrZ", binSize = 40e3, bins = c(0L, 2L)),
                error = identity)
stopifnot(inherits(res, "error"), grepl("non-positive", conditionMessage(res)))

## Invalid 'debug' argument -> error
res <- tryCatch(readHiC(f3, debug = "yes"), error = identity)
stopifnot(inherits(res, "error"))

unlink(tmpdir, recursive = TRUE)
