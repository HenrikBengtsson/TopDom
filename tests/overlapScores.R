library("TopDom")

path <- system.file("exdata", package = "TopDom", mustWork = TRUE)
chr <- Sys.getenv("R_TOPDOM_TESTS_CHROMOSOME", "chr19")

## Original count data
pathname <- file.path(path, sprintf("nij.%s.gz", chr))
data <- readHiC(pathname, chr = chr, binSize = 40e3)
str(data)

## Find topological domains using TopDom method for two window sizes
tds_5 <- TopDom(data, window.size = 5L)
tds_6 <- TopDom(data, window.size = 6L)

## Overlap scores (in both directions)
overlap_56 <- overlapScores(tds_6, reference = tds_5)
print(overlap_56)
overlap_65 <- overlapScores(tds_5, reference = tds_6)
print(overlap_65)

if (FALSE) {  
  ## Find topological domains as a function of windows size
  window.size <- c(5L, 6L, 8L, 10L, 12L, 14L, 16L, 18L, 20L)
  tds <- lapply(window.size, FUN = function(w) TopDom(data, window.size = w))
  names(tds) <- sprintf("window.size=%d", window.size)

  ## Overlap scores relative to the first window.size
  overlaps <- lapply(tds, FUN = overlapScores, reference = tds[[1]])
  print(overlaps)

  scores <- lapply(overlaps, FUN = function(overlaps) {
    unlist(lapply(overlaps, FUN = `[[`, "best_score"), use.names = FALSE)
  })

  avg_scores <- sapply(scores, FUN = mean, na.rm = TRUE)
  print(avg_scores)

  plot(window.size, avg_scores)
  lines(window.size, avg_scores, lwd = 2)
}

