#' @importFrom utils capture.output
mcat <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "", file = stderr())
}


#' Asserts the Truth of R Expressions
#'
#' @param \dots Zero or more \R expressions to be asserted to be TRUE.
#'
#' @return Nothing.
#'
#' @details
#' A bare bone, faster version of [base::stopifnot].
#'
#' @keywords internal
stop_if_not <- function(...) {
  res <- list(...)
  for (ii in seq_along(res)) {
    res_ii <- .subset2(res, ii)
    if (length(res_ii) != 1L || is.na(res_ii) || !res_ii) {
        mc <- match.call()
        call <- deparse(mc[[ii + 1]], width.cutoff = 60L)
        if (length(call) > 1L) call <- paste(call[1L], "....")
        stop(sprintf("%s is not TRUE", sQuote(call)),
             call. = FALSE, domain = NA)
    }
  }
  invisible()
}
