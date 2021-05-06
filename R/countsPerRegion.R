#' Calculates Counts per Region in a TopDomData Object
#'
#' @param data A TopDomData object.
#'
#' @param regions TopDom regions (a data.frame), e.g. domains.
#'
#' @return A numeric vector of length `nrow(regions)`.
#'
#' @author Henrik Bengtsson.
#'
#' @export
countsPerRegion <- function(data, regions) {
  UseMethod("countsPerRegion")
}


#' @importFrom matrixStats colSums2
#' @export
countsPerRegion.TopDomData <- function(data, regions) {
  chr <- NULL  ## To please R CMD check
  
  stopifnot(is.data.frame(regions))
  uchr <- unique(regions$chr)
  stopifnot(length(uchr) == 1L)

  bins <- subset(data$bins, chr == uchr)
  counts <- data$counts
  
  nregions <- nrow(regions)
  total <- vector(storage.mode(counts), length = nregions)

  for (kk in seq_len(nregions)) {
    region <- regions[kk, ]
    idxs_kk <- with(bins, which(from.coord >= region$from.coord & to.coord <= region$to.coord))
    total[kk] <- sum(colSums2(counts, rows = idxs_kk, cols = idxs_kk, na.rm = TRUE), na.rm = TRUE)
  }

  total
}
