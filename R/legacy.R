#' Easy Access to the Original TopDom 0.0.1 and 0.0.2 Implementations
#'
#' @param version A version string.
#'
#' @return An environment containing the legacy TopDom API.
#'
#' @examples
#' TopDom::legacy("0.0.2")$TopDom
#' TopDom::legacy("0.0.1")$Detect.Local.Extreme
#' 
#' @export
legacy <- function(version = c("0.0.1", "0.0.2")) {
  if (version == "0.0.1") {
    api <- environment(TopDom_0.0.1)
  } else if (version == "0.0.2") {
    api <- environment(TopDom_0.0.2)
  } else {
    stop("Unknown TopDom legacy version: ", sQuote(version))
  }
  api
}
