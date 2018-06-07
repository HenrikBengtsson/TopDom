library("TopDom")
if (require("TopDomData")) {
  ## Randomize a bit
  resample <- function(data, z = 0.01, zr = 0:9, scale = 4.0) {
    counts <- data$counts
    is_zero <- (counts == 0L)
    ### (a) tweak some zeros
    idxs <- which(is_zero)
    idxs <- idxs[sample(length(idxs), size = z * length(idxs))]
    counts[idxs] <- sample(zr, size = length(idxs), replace = TRUE)
    ### (a) tweak non-zeros
    idxs <- which(!is_zero)
    counts[idxs] <- as.integer(counts[idxs] * runif(length(idxs), min = 0.0, max = scale))
    data$counts <- counts
    data
  }

  path <- system.file("exdata", package = "TopDomData", mustWork = TRUE)

  ## Original count data
  pathname <- file.path(path, "nij.chr10.gz")
  data <- readHiC(pathname, chr = "chr10", binSize = 40e3)
  str(data)

  ## Find topological domains using TopDom method
  tds_A <- TopDom(data, window.size = 5L)

  ## Slightly different
  dataB <- resample(data)
  tds_B <- TopDom(dataB, window.size = 5L)

  ab_scores <- overlapScores(tds_A, tds_B)
  print(ab_scores)

  ba_scores <- overlapScores(tds_B, tds_A)
  print(ba_scores)
}
