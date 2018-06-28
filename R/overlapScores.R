#' Calculates Overlap Scores between Two Sets of Topological Domains
#' 
#' @param a,b Topological domain (TD) sets \eqn{A} and \eqn{B}
#' as returned by [TopDom()].
#'
#' @param \ldots ...
#'
#' @param debug If `TRUE`, debug output is produced.
#' 
#' @return
#' Returns a named list of class `TopDomOverlapScores`, where the names
#' correspond to the chromosomes in domain set \eqn{A}.
#' Each of these chromosome elements contains a named list of elements
#' `best_score` (\eqn{D_A_c} numerics in \eqn{[0,1]}) and
#' `best_sets' (\eqn{D_A_c} list of index vectors), where
#' \eqn{D_A_c} is the number of TDs in chromosome \eqn{c} in set \eqn{A}.
#'
#' @details
#' The _overlap score_, \eqn{overlap(a_i, B')}, represents how well topological
#' domain (TD) \eqn{a_i} in set \eqn{A} overlap with a _consecutive_ subset
#' \eqn{B'} of TDs in \eqn{B}.
#' For each TD \eqn{a_i}, the _best match_ \eqn{B'_max} is identified, that
#' is, the \eqn{B'} subset that maximize \eqn{overlap(a_i, B')}.
#' For exact definitions, see Page 8 in Shin et al. (2016).
#'
#' Note that the overlap score is an asymmetric score, which means that
#' `overlapScores(a, b) != overlapScores(b, a)`.
#' 
#' @example incl/overlapScores.R
#' 
#' @references
#' * Shin et al.,
#'   TopDom: an efficient and deterministic method for identifying
#'   topological domains in genomes,
#'   _Nucleic Acids Research_, 44(7): e70, April 2016.
#'   doi: 10.1093/nar/gkv1505,
#'   PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/),
#'   PMID: [26704975](https://www.ncbi.nlm.nih.gov/pubmed/26704975).
#'
#' @author Henrik Bengtsson - based on the description in Shin et al. (2016).
#' 
#' @seealso [TopDom].
#'
#' @export
overlapScores <- function(a, b, ..., debug = getOption("TopDom.debug", FALSE)) {
  stopifnot(inherits(a, "TopDom"), inherits(b, "TopDom"))
  stopifnot(is.logical(debug), length(debug) == 1L, !is.na(debug))

  domains_A <- a$domain
  domains_B <- b$domain
  
  chromosomes <- unique(domains_A$chr)
  scores <- vector("list", length = length(chromosomes))
  names(scores) <- chromosomes
  for (chr in chromosomes) {
    if (debug) message(sprintf("Chromosome %s ...", chr))
    doms_A <- domains_A[domains_A$chr == chr, ]
    doms_B <- domains_B[domains_B$chr == chr, ]
    scores[[chr]] <- overlapScoresOneChromosome(doms_A, doms_B, ...,
                                                debug = debug)
    if (debug) message(sprintf("Chromosome %s ... done", chr))
  }

  structure(scores, class = "TopDomOverlapScores")
}


overlapScoresOneChromosome <- function(doms_A, doms_B, ..., debug = getOption("TopDom.debug", FALSE)) {
  stopifnot(is.logical(debug), length(debug) == 1L, !is.na(debug))
  ## FIXME: Assert that A and B are sorted by (chr, pos)

  dtags <- diff(c(0L, as.integer(doms_B$tag == "domain"), 0L))
  sets <- data.frame(
    first = which(dtags == +1L),
    last  = which(dtags == -1L) - 1L
  )
  sets$from.coord <- doms_B$from.coord[sets$first]
  sets$to.coord <- doms_B$to.coord[sets$last]

  doms_A$length <- doms_A$to.coord - doms_A$from.coord
  doms_B$length <- doms_B$to.coord - doms_B$from.coord

  best_scores <- rep(NA_real_, length = nrow(doms_A))
  best_sets <- vector("list", length = nrow(doms_A))
  idxs_td <- which(doms_A$tag == "domain")
  for (ii in seq_along(idxs_td)) {
    idx_td <- idxs_td[ii]
    if (debug) message(sprintf("TD #%d of %d ...", ii, length(idxs_td)))
    td <- doms_A[idx_td, ]
    
    ## Identify sets to consider
    ## Q. Now many sets can match this? [0,1], [0,2], [0,3], or even more?
    sets_t <- sets[(sets$to.coord >= td$from.coord &
                    sets$from.coord <= td$to.coord), ]

    best_score <- 0.0
    best_set <- integer(0L)
    
    ## For each possible B' set ...
    for (kk in seq_len(nrow(sets_t))) {
      set <- sets_t[kk,]
      doms <- doms_B[set$first:set$last,]
      doms$cap <- doms$length
      
      ## TD in B' that is overlapping part of the beginning of A 
      before <- which(doms$from.coord <= td$from.coord)
      if (length(before) > 0L) {
        before <- before[length(before)]
        doms$cap[before] <- doms$to.coord[before] - td$from.coord
      }
      
      ## TD in B' that is overlapping part of the end of A 
      after <- which(doms$to.coord >= td$to.coord)
      if (length(after) > 0L) {
        after <- after[1L]
        doms$cap[after] <- td$to.coord - doms$from.coord[after]
      }
      
      ## TDs in B' that are strictly overlapping with A
      is_inside <- (doms$from.coord >= td$from.coord &
                    doms$to.coord <= td$to.coord)
      idxs_t <- c(before, which(is_inside), after)
      caps <- doms$cap[idxs_t]
      cups <- doms$length[idxs_t]
      n <- length(idxs_t)
      if (n == 1L) {
        idxs_u <- list(1L)
        max_score <- caps / cups
      } else if (n == 2L) {
        idxs_u <- list(1L, 1:2, 2L)
      } else if (n >= 3L) {
        idxs_u <- list(1L, 1L:(n-1L), 2L:(n-1L), 2L:n, n)
      }
      scores <- sapply(idxs_u, FUN = function(idxs) {
        sum(caps[idxs]) / sum(cups[idxs])
      })
      max_idx <- which.max(scores)

      max_score <- scores[max_idx]
      if (max_score > best_score) {
        best_score <- max_score
        best_set <- idxs_u[[max_idx]]
      }
    }

    best_scores[ii] <- best_score
    best_sets[[ii]] <- best_set
    
    if (debug) message(sprintf("TD #%d of %d ... done", ii, length(idxs_td)))
  } ## for (ii ...)

  list(best_scores = best_scores, best_sets = best_sets)
} ## overlapScores()


#' @export
print.TopDomOverlapScores <- function(x, ...) {
  cat(sprintf("%s:\n", class(x)))
  cat(sprintf("Chromosomes: [n = %d] %s\n",
              length(x), paste(sQuote(names(x)), collapse = ", ")))
  scores <- lapply(x, FUN = `[[`, "best_scores")
  scores[["whole genome"]] <- unlist(scores, use.names = FALSE)
  cat("Summary of best scores:\n")
  t <- t(sapply(scores, FUN = function(x) {
    c(summary(x), count = length(x))
  }))
  print(t)
}
