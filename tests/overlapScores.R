library("TopDom")
if (require("TopDomData")) {
  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## Original count data
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  str(data)

  ## Find topological domains using TopDom method for two window sizes
  tds_5 <- TopDom(data, window.size = 5L)
  tds_6 <- TopDom(data, window.size = 6L)

  ## Overlap scores (in both directions)
  overlap_56 <- overlapScores(tds_5, tds_6)
  print(overlap_56)
  overlap_65 <- overlapScores(tds_6, tds_5)
  print(overlap_65)

  if (FALSE) {  
    ## Find topological domains as a function of windows size
    window.size <- c(5L, 6L, 8L, 10L, 12L, 14L, 16L, 18L, 20L)
    tds <- lapply(window.size, FUN = function(w) TopDom(data, window.size = w))
    names(tds) <- sprintf("window.size=%d", window.size)
  
    ## Overlap scores relative to the first window.size
    overlaps <- lapply(tds, FUN = overlapScores, tds[[1]])
    print(overlaps)
  
    scores <- lapply(overlaps, FUN = function(overlaps) {
      unlist(lapply(overlaps, FUN = `[[`, "best_scores"), use.names = FALSE)
    })
  
    avg_scores <- sapply(scores, FUN = mean, na.rm = TRUE)
    print(avg_scores)
  
    plot(window.size, avg_scores)
    lines(window.size, avg_scores, lwd = 2)
  }
}