#' Calculates Overlap Scores Between Two Sets of Topological Domains
#' 
#' @param a,reference Topological domain (TD) set \eqn{A} and TD reference
#' set \eqn{R} both in a format as returned by [TopDom()].
#'
#' @param debug If `TRUE`, debug output is produced.
#' 
#' @return
#' Returns a named list of class `TopDomOverlapScores`, where the names
#' correspond to the chromosomes in domain reference set \eqn{R}.
#' Each of these chromosome elements contains a data.frame with fields:
#'
#' * `chromosome`  - \eqn{D_{R,c}} character strings
#' * `best_score`  - \eqn{D_{R,c}} numerics in \eqn{[0,1]}
#' * `best_length` - \eqn{D_{R,c}} positive integers
#' * `best_set`    - list of \eqn{D_{R,c}} index vectors
#'
#' where \eqn{D_{R,c}} is the number of TDs in reference set \eqn{R} on
#' chromosome \eqn{c}.  If a TD in reference \eqn{R} is not a `"domain"`,
#' then the corresponding `best_score` and `best_length` values are
#' `NA_real_` and `NA_integer_`, respectively, while `best_set` is an empty
#' list.
#'
#' @details
#' The _overlap score_, \eqn{overlap(A', r_i)}, represents how well a
#' _consecutive_ subset \eqn{A'} of topological domains (TDs) in \eqn{A}
#' overlap with topological domain \eqn{r_i} in reference set \eqn{R}.
#' For each reference TD \eqn{r_i}, the _best match_ \eqn{A'_max} is
#' identified, that is, the \eqn{A'} subset that maximize
#' \eqn{overlap(A', r_i)}.
#' For exact definitions, see Page 8 in Shin et al. (2016).
#'
#' Note that the overlap score is an asymmetric score, which means that
#' `overlapScores(a, b) != overlapScores(b, a)`.
#'
#' @section Warning - This might differ not be the correct implementation:
#' The original TopDom scripts do not provide an implementation for
#' calculating overlap scores.  Instead, the implementation of
#' `TopDom::overlapScores()` is based on the textual description of
#' overlap scores provided in Shin et al. (2016).  It is not known if this
#' is the exact same algorithm and implementation as the authors of the
#' TopDom article used.
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
#'   PMID: [26704975](https://pubmed.ncbi.nlm.nih.gov/26704975/)
#'
#' @author Henrik Bengtsson - based on the description in Shin et al. (2016).
#' 
#' @seealso [TopDom].
#'
#' @importFrom tibble as_tibble tibble
#' @export
overlapScores <- function(a, reference, debug = getOption("TopDom.debug", FALSE)) {
  stopifnot(inherits(reference, "TopDom"), inherits(a, "TopDom"))
  stopifnot(is.logical(debug), length(debug) == 1L, !is.na(debug))

  ## Extract the 'domain' fields from the 'reference' and the 'a' TopDom objects
  domains_R <- reference$domain
  domains_A <- a$domain
  
  chromosomes <- unique(domains_R$chr)
  scores <- vector("list", length = length(chromosomes))
  names(scores) <- chromosomes
  for (chr in chromosomes) {
    if (debug) message(sprintf("Chromosome %s ...", chr))
    doms_R <- domains_R[domains_R$chr == chr, ]
    doms_A <- domains_A[domains_A$chr == chr, ]
    scores_chr <- overlapScoresOneChromosome(doms_A, doms_R = doms_R, debug = debug)
    scores_chr <- as_tibble(cbind(tibble(chromosome = chr), scores_chr))
    scores[[chr]] <- scores_chr
    if (debug) message(sprintf("Chromosome %s ... done", chr))
  }

  class(scores) <- c("TopDomOverlapScores", class(scores))
  scores
}


#' @importFrom tibble as_tibble
overlapScoresOneChromosome <- function(doms_A, doms_R, debug = getOption("TopDom.debug", FALSE)) {
  stopifnot(is.logical(debug), length(debug) == 1L, !is.na(debug))
  ## FIXME: Assert that A and reference R are sorted by (chr, pos)

  dtags <- diff(c(0L, as.integer(doms_A$tag == "domain"), 0L))
  sets <- data.frame(
    first = which(dtags == +1L),
    last  = which(dtags == -1L) - 1L,
    stringsAsFactors = FALSE
  )
  sets$from.coord <- doms_A$from.coord[sets$first]
  sets$to.coord <- doms_A$to.coord[sets$last]

  doms_R$length <- as.integer(doms_R$to.coord - doms_R$from.coord)
  doms_A$length <- as.integer(doms_A$to.coord - doms_A$from.coord)

  best_scores <- rep(NA_real_, times = nrow(doms_R))
  best_lengths <- rep(NA_integer_, times = nrow(doms_R))
  best_sets <- vector("list", length = nrow(doms_R))
  idxs_td <- which(doms_R$tag == "domain")
  for (ii in seq_along(idxs_td)) {
    idx_td <- idxs_td[ii]
    if (debug) message(sprintf("TD \"domain\" #%d of %d ...", ii, length(idxs_td)))
    td_R <- doms_R[idx_td, ]
    
    ## Identify sets to consider
    ## Q. How many sets can match this? [0,1], [0,2], [0,3], or even more?
    sets_t <- sets[(sets$to.coord >= td_R$from.coord &
                    sets$from.coord <= td_R$to.coord), ]

    best_score <- 0.0
    best_set <- integer(0L)
    
    ## For each possible A' set ...
    for (kk in seq_len(nrow(sets_t))) {
      set <- sets_t[kk,]
      doms <- doms_A[set$first:set$last,]
      doms$cap <- doms$length
      
      ## TD in A' that is overlapping part of the beginning of reference R 
      before <- which(doms$from.coord <= td_R$from.coord)
      if (length(before) > 0L) {
        before <- before[length(before)]
        doms$cap[before] <- doms$to.coord[before] - td_R$from.coord
      }
      
      ## TD in A' that is overlapping part of the end of reference R
      after <- which(doms$to.coord >= td_R$to.coord)
      if (length(after) > 0L) {
        after <- after[1L]
        doms$cap[after] <- td_R$to.coord - doms$from.coord[after]
      }
      
      ## TDs in A' that are strictly overlapping with reference R
      is_inside <- (doms$from.coord >= td_R$from.coord &
                    doms$to.coord <= td_R$to.coord)
      idxs_t <- c(before, which(is_inside), after)
      n <- length(idxs_t)
      stop_if_not(n > 0L)
      caps <- doms$cap[idxs_t]
      stop_if_not(length(caps) > 0L, any(is.finite(caps)))
      stop_if_not(!anyNA(caps))
      cups <- doms$length[idxs_t]
      stop_if_not(length(cups) > 0L, any(is.finite(cups)))
      if (n == 1L) {
        idxs_u <- list(1L)
        max_score <- caps / cups
      } else if (n == 2L) {
        idxs_u <- list(1L, 1:2, 2L)
      } else if (n >= 3L) {
        idxs_u <- list(1L, 1L:(n-1L), 2L:(n-1L), 2L:n, n)
      }
      stop_if_not(!anyNA(idxs_u))
      scores <- sapply(idxs_u, FUN = function(idxs) {
        sum(caps[idxs]) / sum(cups[idxs])
      })
      stop_if_not(any(is.finite(scores)))
      max_idx <- which.max(scores)

      max_score <- scores[max_idx]
      if (max_score > best_score) {
        best_score <- max_score
        best_set <- idxs_u[[max_idx]]
      }
    } ## for (kk ...)

    ## Sanity checks
    stop_if_not(length(best_score) == 1L, length(td_R$length) == 1L)

    best_scores[ii] <- best_score
    best_lengths[ii] <- td_R$length
    best_sets[[ii]] <- best_set
    
    if (debug) message(sprintf("TD \"domain\" #%d of %d ... done", ii, length(idxs_td)))
  } ## for (ii ...)

  ## Sanity checks
  stop_if_not(length(best_scores) == nrow(doms_R),
              length(best_lengths) == nrow(doms_R),
              length(best_sets) == nrow(doms_R))

  res <- data.frame(best_score = best_scores, best_length = best_lengths, stringsAsFactors = FALSE)
  res$best_set <- best_sets

  res
} ## overlapScoresOneChromosome()


#' @importFrom tibble as_tibble
#' @export
as_tibble.TopDomOverlapScores <- function(x, ...) {
  do.call(rbind, args = x)
}


#' @export
print.TopDomOverlapScores <- function(x, ...) {
  cat(sprintf("%s:\n", paste(class(x), collapse = ", ")))

  cat(sprintf("Chromosomes: [n = %d] %s\n",
              length(x), paste(sQuote(names(x)), collapse = ", ")))

  lengths <- lapply(x, FUN = `[[`, "best_length")
  lengths[["whole genome"]] <- unlist(lengths, use.names = FALSE)
  cat("Summary of reference domain lengths:\n")
  t <- t(sapply(lengths, FUN = function(x) {
    c(summary(x), count = length(x))
  }))
  print(t)
  
  scores <- lapply(x, FUN = `[[`, "best_score")
  scores[["whole genome"]] <- unlist(scores, use.names = FALSE)
  cat("Summary of best scores:\n")
  t <- t(sapply(scores, FUN = function(x) {
    c(summary(x), count = length(x))
  }))
  
  print(t)
}
