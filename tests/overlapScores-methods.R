library("TopDom")

path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")
pathname <- file.path(path, sprintf("nij.%s.gz", chr))

data <- readHiC(pathname, chr = chr, binSize = 40e3)
data <- data[1:200]

tds_5 <- TopDom(data, window.size = 5L)
tds_6 <- TopDom(data, window.size = 6L)

## overlapScores() with debug output ------------------------------------
msg <- capture.output(
  ov <- overlapScores(tds_6, reference = tds_5, debug = TRUE),
  type = "message"
)
stopifnot(inherits(ov, "TopDomOverlapScores"), any(grepl("Chromosome", msg)))

## Self-overlap should be perfect (score 1) for every reference domain --
ov_self <- overlapScores(tds_5, reference = tds_5)
scores <- ov_self[[chr]]$best_score
stopifnot(all(scores[!is.na(scores)] == 1))

## print.TopDomOverlapScores() ----------------------------------------
out <- capture.output(print(ov_self))
stopifnot(any(grepl("Chromosomes:", out)),
          any(grepl("best scores", out)))

## as_tibble.TopDomOverlapScores() -----------------------------------
tbl <- tibble::as_tibble(ov_self)
stopifnot(
  inherits(tbl, "data.frame"),
  all(c("chromosome", "best_score", "best_length", "best_set") %in% names(tbl)),
  nrow(tbl) == nrow(tds_5$domain)
)

## Non-TopDom inputs are rejected ----------------------------------
res <- tryCatch(overlapScores(tds_5, reference = list()), error = identity)
stopifnot(inherits(res, "error"))
res <- tryCatch(overlapScores("x", reference = tds_5), error = identity)
stopifnot(inherits(res, "error"))
