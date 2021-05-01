#' Subset a TopDomData Object by Region
#'
#' @param data A TopDomData object.
#'
#' @param region A TopDom domain (a data.frame).
#'
#' @param margin An non-negative numeric specifying the additional margin
#'               extracted around the domain.
#'               If `margin < 1`, then the size of the margin is relative
#'               to the size of the domain.
#'
#' @return A TopDomData object.
#'
#' @author Henrik Bengtsson.
#'
#' @export
subsetByRegion <- function(data, region, margin = 1/2) {
  UseMethod("subsetByRegion")
}


#' @export
subsetByRegion.TopDomData <- function(data, region, margin = 1/2) {
  stopifnot(is.data.frame(region))
  stopifnot(margin >= 0)
  
  if (margin < 1) {
    margin <- margin * (region$to.coord - region$from.coord)
  }
  
  idxs <- with(data$bins, which(chr == region$chr & from.coord >= region$from.coord - margin & to.coord <= region$to.coord + margin))
  
  data[idxs]
}


#' @export
subsetByRegion.TopDom <- function(data, region, margin = 1/2) {
  stopifnot(is.data.frame(region))
  stopifnot(margin >= 0)
  
  if (margin < 1) {
    margin <- margin * (region$to.coord - region$from.coord)
  }

  idxs <- with(data$domain, which(chr == region$chr & from.coord >= region$from.coord - margin & to.coord <= region$to.coord + margin))
  
  data[idxs, ]
}
