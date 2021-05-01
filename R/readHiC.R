#' Reads Hi-C Contact Data from File
#' 
#' @param file The pathname of a normalize Hi-C contact matrix file
#' stored as a whitespace-delimited file.  See below for details.
#' Also a gzip-compressed file can be used.
#'
#' @param chr,binSize If the file contains a count matrix without bin
#' annotation, the latter is created from these parameters.
#'
#' @param \ldots Arguments passed to [utils::read.table()] as-is.
#'
#' @param debug If `TRUE`, debug output is produced.
#' 
#' @return A list with elements \code{bins} (an N-by-4 data.frame) and
#' \code{counts} (N-by-N matrix).
#'
#' @section Format of HiC contact-matrix file:
#' The contact-matrix file should be a whitespace-delimited text file with
#' neither row names nor column names.  The content should be a N-by-(3+N)
#' table where the first three columns correspond to `chr` (string),
#' `from.coord` (integer position), and `to.coord` (integer position).
#' These column defines the genomic location of the N Hi-C bins (in order).
#' The last N columns should contain normalized contact counts (float) such
#' that element (r,3+c) in this table corresponds to count (r,c) in the
#' normalized contact matrix.
#'
#' If an N-by-(4+N) table, then the first column is assumed to contain an
#' `id` (integer), and everything else as above.
#'
#' Example:
#' \preformatted{
#'   chr10       0   40000  0 0 0 0 ...
#'   chr10   40000   80000  0 0 0 0 ...
#'   chr10   80000  120000  0 0 0 0 ...
#'   chr10  120000  160000  0 0 0 0 ...
#'   ...
#' }
#'
#' @example incl/readHiC.R
#' 
#' @seealso [TopDom].
#'
#' @importFrom utils file_test read.table
#' @export
readHiC <- function(file, chr = NULL, binSize = NULL, ..., debug = getOption("TopDom.debug", FALSE)) {
  stopifnot(file_test("-f", file))
  stopifnot(is.logical(debug), length(debug) == 1L, !is.na(debug))

  if (debug) {
    mcat("#########################################################################")
    mcat("Step 0 : File Read")
    mcat("#########################################################################")
  }
  
  if (!is.null(chr)) {
    chr <- as.character(chr)
    stopifnot(is.character(chr), length(chr) == 1, !is.na(chr))
    binSize <- as.integer(binSize)
    stopifnot(is.integer(binSize), length(binSize) == 1, !is.na(binSize), binSize >= 1)

    args <- list(..., comment.char = "", na.strings = "", quote = "",
                      stringsAsFactors = FALSE)
    bins <- args$bins                        
    args$bins <- NULL
    args <- args[unique(names(args))]
    argsT <- c(list(file, header = FALSE, nrows = 1L), args)
    first <- do.call(read.table, args = argsT)
    if (debug) mcat("  -- reading ", length(first), "-by-", length(first), " count matrix")
    ## Assert that it's a count matrix
    is.numeric <- unlist(lapply(first, FUN = is.numeric), use.names = FALSE)
    stopifnot(all(is.numeric))
    
    if (!is.null(bins)) {
      if (any(bins < 1)) {
        stop("Argument 'bins' specifies non-positive bin indices")
      } else if (any(bins > length(first))) {
        stop(sprintf("Argument 'bins' specifies bin indices out of range [1,%d]", length(first)))
      }
      colClasses <- rep("NULL", times = length(first))
      colClasses[bins] <- "numeric"
    } else {
      colClasses <- rep("numeric", times = length(first))
    }
    argsT <- c(list(file, colClasses = colClasses, header = FALSE), args)
    matrix.data <- do.call(read.table, args = argsT)
    colnames(matrix.data) <- NULL

    if (!is.null(bins)) {
      matrix.data <- matrix.data[bins, , drop = FALSE]
    }

    ## N-by-N count matrix (from file content)
    matrix.data <- as.matrix(matrix.data)
    dimnames(matrix.data) <- NULL
    stopifnot(nrow(matrix.data) == ncol(matrix.data))
    n_bins <- length(first)

    from.coord <- seq(from = 0, by = binSize, length.out = n_bins)
    to.coord   <- seq(from = binSize, by = binSize, length.out = n_bins)
    if (is.null(bins)) {
      stopifnot(n_bins == nrow(matrix.data))
      id <- seq_len(n_bins)
    } else {
      n_bins <- length(bins)
      id <- bins
      from.coord <- from.coord[bins]
      to.coord   <- to.coord[bins]
    }
    
    ## Bin annotation from (chr, binSize)
    bins <- data.frame(
      id               = id,
      chr              = chr,
      from.coord       = from.coord,
      to.coord         = to.coord,
      stringsAsFactors = FALSE
    )
  } else {
    args <- list(..., stringsAsFactors = FALSE)
    args <- args[unique(names(args))]
    argsT <- c(list(file, header = FALSE), args)
    matdf <- do.call(read.table, args = argsT)
    n_bins <- nrow(matdf)
    if (ncol(matdf) - n_bins == 3) {
      colnames(matdf) <- c("chr", "from.coord", "to.coord")
    } else if (ncol(matdf) - n_bins == 4) {
      colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
    } else {
      stop("Unknown format of count-matrix file: ", sQuote(file))
    }

    ## Bin annotation (from file content)
    bins <- data.frame(
      id         = seq_len(n_bins),
      chr        = matdf[["chr"]],
      from.coord = matdf[["from.coord"]],
      to.coord   = matdf[["to.coord"]],
      stringsAsFactors = FALSE
    )

    ## N-by-N count matrix (from file content)
    matdf <- matdf[, (ncol(matdf) - n_bins + 1):ncol(matdf)]
    matrix.data <- as.matrix(matdf)
    rm(list = "matdf")
  }

  stopifnot(is.numeric(matrix.data),
            is.matrix(matrix.data),
            nrow(matrix.data) == ncol(matrix.data),
            nrow(matrix.data) == n_bins)

  if (debug) mcat("-- Done!")
  if (debug) mcat("Step 0 : Done!")

  structure(list(bins = bins, counts = matrix.data), class = "TopDomData")
}
