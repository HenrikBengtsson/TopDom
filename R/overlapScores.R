#' Calculates Overlap Scores between Two Sets of Topological Domains
#' 
#' @param a,b Topological domain sets 'A' and 'B', where set 'A' contains
#' D_A domains and set 'B' D_B domains.
#'
#' @param \ldots ...
#'
#' @value
#' Returns a data frame with D_A rows.
#' 
#' @seealso [TopDom].
#'
#' @export
overlapScores <- function(a, b, ...) {
  stopifnot(inherits(a, "TopDom"), inherits(b, "TopDom"))

  domains_A <- a$domain
  domains_B <- b$domain
  
  chromosomes <- unique(domains_A$chr)
  scores <- vector("list", length = length(chromosomes))
  names(scores) <- chromosomes
  for (chr in chromosomes) {
    message(sprintf("Chromosome %s ...", chr))
    doms_A <- domains_A[domains_A$chr == chr, ]
    doms_B <- domains_B[domains_B$chr == chr, ]
    scores[[chr]] <- overlapScoresOneChromosome(doms_A, doms_B, ...)
    message(sprintf("Chromosome %s ... done", chr))
  }

  scores
}


overlapScoresOneChromosome <- function(doms_A, doms_B, ...) {
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
    message(sprintf("TD #%d of %d ...", ii, length(idxs_td)))
    td <- doms_A[idx_td, ]
    print(td)
    
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
    
    message(sprintf("TD #%d of %d ... done", ii, length(idxs_td)))
  } ## for (ii ...)

  list(best_scores = best_scores, best_sets = best_sets)
} ## overlapScores()
