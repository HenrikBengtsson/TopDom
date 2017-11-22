#' Identify topological domains for given Hi-C contact matrix
#' 
#' @param matrix.file string, matrixFile Address,
#' - Format = {chromosome, bin start, bin end, N numbers normalized value}
#' - N * (N + 3), where N is the number of bins
#' 
#' @param window.size integer, number of bins to extend.
#' 
#' @param outFile string, binSignal file address to write
#' 
#' @param statFilter logical, ...
#'
#' @return A named list with elements `binSignal`, `domain`, and `bed`.
#' 
#' @author Hanjun Shin, Harris Lazaris, and Gangqing Hu.
#' R package, help and code refactoring by Henrik Bengtsson.
#' 
#' @importFrom utils read.table write.table
#' @export
TopDom <- function(matrix.file, window.size, outFile = NULL, statFilter = TRUE) {
  mcat("#########################################################################")
  mcat("Step 0 : File Read")
  mcat("#########################################################################")
  window.size <- as.numeric(window.size)
  matdf <- read.table(matrix.file, header = FALSE)

  if (ncol(matdf) - nrow(matdf) == 3) {
    colnames(matdf) <- c("chr", "from.coord", "to.coord")
  } else if (ncol(matdf) - nrow(matdf) == 4) {
    colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
  } else {
    stop("Unknown type of matrix file: ", sQuote(matrix.file))
  }
  n_bins <- nrow(matdf)
  mean.cf <- rep(0, times = n_bins)
  pvalue <- rep(1, times = n_bins)

  local.ext <- rep(-0.5, times = n_bins)

  bins <- data.frame(
    id = seq_len(n_bins),
    chr = matdf[, "chr"],
    from.coord = matdf[, "from.coord"],
    to.coord = matdf[, "to.coord"]
  )
  matrix.data <- as.matrix(matdf[, (ncol(matdf) - nrow(matdf) + 1):ncol(matdf)])

  mcat("-- Done!")
  mcat("Step 0 : Done!")


  mcat("#########################################################################")
  mcat("Step 1 : Generating binSignals by computing bin-level contact frequencies")
  mcat("#########################################################################")
  ptm <- proc.time()
  for (i in seq_len(n_bins)) {
    diamond <- Get.Diamond.Matrix(mat.data = matrix.data, i = i, size = window.size)
    mean.cf[i] <- mean(diamond)
  }

  eltm <- proc.time() - ptm
  mcat("Step 1 Running Time : ", eltm[3])
  mcat("Step 1 : Done!")

  mcat("#########################################################################")
  mcat("Step 2 : Detect TD boundaries based on binSignals")
  mcat("#########################################################################")

  ptm <- proc.time()
  # gap.idx <- Which.Gap.Region(matrix.data = matrix.data)
  # gap.idx <- Which.Gap.Region2(mean.cf)
  gap.idx <- Which.Gap.Region2(matrix.data = matrix.data, w = window.size)

  proc.regions <- Which.process.region(rmv.idx = gap.idx, n_bins = n_bins, min.size = 3)

  # mcat(proc.regions)

  for (i in seq_len(nrow(proc.regions))) {
    start <- proc.regions[i, "start"]
    end <- proc.regions[i, "end"]

    mcat("Process Regions from ", start, " to ", end)

    local.ext[start:end] <- Detect.Local.Extreme(x = mean.cf[start:end])
  }

  eltm <- proc.time() - ptm
  mcat("Step 2 Running Time : ", eltm[3])
  mcat("Step 2 : Done!")

  if (statFilter) {
    mcat("#########################################################################")
    mcat("Step 3 : Statistical Filtering of false positive TD boundaries")
    mcat("#########################################################################")

    ptm <- proc.time()
    mcat("-- Matrix Scaling....")
    scale.matrix.data <- matrix.data
    for (i in seq_len(2 * window.size)) {
      # diag(scale.matrix.data[, i:n_bins]) <- scale(diag(matrix.data[, i:n_bins]))
      scale.matrix.data[seq(from = 1 + (n_bins * i), to = n_bins * n_bins, by = 1 + n_bins)] <- scale(matrix.data[seq(from = 1 + (n_bins * i), to = n_bins * n_bins, by = 1 + n_bins)])
    }

    mcat("-- Compute p-values by Wilcox Ranksum Test")
    for (i in seq_len(nrow(proc.regions))) {
      start <- proc.regions[i, "start"]
      end <- proc.regions[i, "end"]

      mcat("Process Regions from ", start, " to ", end)

      pvalue[start:end] <- Get.Pvalue(matrix.data = scale.matrix.data[start:end, start:end], size = window.size, scale = 1.0)
    }
    mcat("-- Done!")

    mcat("-- Filtering False Positives")
    local.ext[intersect(union(which(local.ext == -1), which(local.ext == -1)), which(pvalue < 0.05))] <- -2
    local.ext[which(local.ext == -1)] <- 0
    local.ext[which(local.ext == -2)] <- -1
    mcat("-- Done!")

    eltm <- proc.time() - ptm
    mcat("Step 3 Running Time : ", eltm[3])
    mcat("Step 3 : Done!")
  } else {
    pvalue <- 0
  }

  domains <- Convert.Bin.To.Domain.TMP(
    bins = bins,
    signal.idx = which(local.ext == -1),
    gap.idx = which(local.ext == -0.5),
    pvalues = pvalue,
    pvalue.cut = 0.05
  )

  bins <- cbind(
    bins,
    local.ext = local.ext,
    mean.cf = mean.cf,
    pvalue = pvalue
  )

  bedform <- domains[, c("chr", "from.coord", "to.coord", "tag")]
  colnames(bedform) <- c("chrom", "chromStart", "chromEnd", "name")

  if (!is.null(outFile)) {
    mcat("#########################################################################")
    mcat("Writing Files")
    mcat("#########################################################################")

    outBinSignal <- paste0(outFile, ".binSignal")
    mcat("binSignal File : ", outBinSignal)
    write.table(bins, file = outBinSignal, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


    outDomain <- paste0(outFile, ".domain")
    mcat("Domain File : ", outDomain)
    write.table(domains, file = outDomain, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

    outBed <- paste0(outFile, ".bed")
    mcat("Bed File : ", outBed)
    write.table(bedform, file = outBed, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }

  mcat("Done!")

  mcat("Job Complete!")
  list(binSignal = bins, domain = domains, bed = bedform)
}

# @fn Get.Diamond.Matrix
# @param mat.data : N by N matrix, where each element indicate contact frequency
# @param i :integer, bin index
# @param size : integer, window size to expand from bin
# @retrun : matrix.
Get.Diamond.Matrix <- function(mat.data, i, size) {
  n_bins <- nrow(mat.data)
  if (i == n_bins) {
    na <- NA_real_
    storage.mode(na) <- storage.mode(mat.data)
    return(na)
  }

  lowerbound <- max(1, i - size + 1)
  upperbound <- min(i + size, n_bins)

  mat.data[lowerbound:i, (i + 1):upperbound]
}

# @fn Which.process.region
# @param rmv.idx : vector of idx, remove index vector
# @param n_bins : total number of bins
# @param min.size : minimum size of bins
# @retrun : data.frame of proc.regions
Which.process.region <- function(rmv.idx, n_bins, min.size = 3) {
  gap.idx <- rmv.idx

  proc.regions <- data.frame(start = numeric(0), end = numeric(0))
  proc.set <- setdiff(seq_len(n_bins), gap.idx)
  n_proc.set <- length(proc.set)

  i <- 1
  while (i < n_proc.set) {
    start <- proc.set[i]
    j <- i + 1

    while (j <= n_proc.set) {
      if (proc.set[j] - proc.set[j - 1] <= 1) {
        j <- j + 1
      } else {
        proc.regions <- rbind(proc.regions, c(start = start, end = proc.set[j - 1]))
        i <- j
        break
      }
    }

    if (j >= n_proc.set) {
      proc.regions <- rbind(proc.regions, c(start = start, end = proc.set[j - 1]))
      break
    }
  }

  colnames(proc.regions) <- c("start", "end")
  proc.regions <- proc.regions[which(abs(proc.regions[, "end"] - proc.regions[, "start"]) >= min.size), ]

  proc.regions
}

# @fn Which.Gap.Region
# @breif version 0.0.1 used
# @param matrix.data : n by n matrix
# @return gap index
Which.Gap.Region <- function(matrix.data) {
  n_bins <- nrow(matrix.data)
  gap <- rep(0, times = n_bins)

  i <- 1
  while (i < n_bins) {
    j <- i + 1
    while (j <= n_bins) {
      if (sum(matrix.data[i:j, i:j]) == 0) {
        gap[i:j] <- -0.5
        j <- j + 1
        # if (j-i > 1) gap[i:j] <- -0.5
        # j <- j+1
      } else {
        break
      }
    }

    i <- j
  }

  idx <- which(gap == -0.5)
  idx
}

# @fn Which.Gap.Region3
# @param matrix.data : n by n matrix
# @return gap index
Which.Gap.Region3 <- function(mean.cf) {
  n_bins <- length(mean.cf)
  gapidx <- which(mean.cf == 0)

  gapidx
}

# @fn Which.Gap.Region2
# @breif version 0.0.2 used
# @param matrix.data : n by n matrix
# @return gap index
Which.Gap.Region2 <- function(matrix.data, w) {
  n_bins <- nrow(matrix.data)
  gap <- rep(0, times = n_bins)

  for (i in seq_len(n_bins)) {
    if (sum(matrix.data[i, max(1, i - w):min(i + w, n_bins)]) == 0) gap[i] <- -0.5
  }

  idx <- which(gap == -0.5)
  idx
}

# @fn Detect.Local.Extreme
# @param x : original signal to find local minima
# @return vector of local extrme, -1 if the index is local minimum, 1 if the index is local maxima, 0 otherwise.
Detect.Local.Extreme <- function(x) {
  n_bins <- length(x)
  ret <- rep(0, times = n_bins)
  x[is.na(x)] <- 0

  if (n_bins <= 3) {
    ret[which.min(x)] <- -1
    ret[which.max(x)] <- 1

    return(ret)
  }
  # Norm##################################################3
  new.point <- Data.Norm(x = seq_len(n_bins), y = x)
  x <- new.point$y
  ##################################################
  cp <- Change.Point(x = seq_len(n_bins), y = x)

  if (length(cp$cp) <= 2) return(ret)
  if (length(cp$cp) == n_bins) return(ret)
  for (i in 2:(length(cp$cp) - 1)) {
    if (x[cp$cp[i]] >= x[cp$cp[i] - 1] && x[cp$cp[i]] >= x[cp$cp[i] + 1]) {
      ret[cp$cp[i]] <- 1
    } else if (x[cp$cp[i]] < x[cp$cp[i] - 1] && x[cp$cp[i]] < x[cp$cp[i] + 1]) ret[cp$cp[i]] <- -1

    min.val <- min(x[cp$cp[i - 1]], x[cp$cp[i]])
    max.val <- max(x[cp$cp[i - 1]], x[cp$cp[i]])

    if (min(x[cp$cp[i - 1]:cp$cp[i]]) < min.val) ret[cp$cp[i - 1] - 1 + which.min(x[cp$cp[i - 1]:cp$cp[i]])] <- -1
    if (max(x[cp$cp[i - 1]:cp$cp[i]]) > max.val) ret[cp$cp[i - 1] - 1 + which.max(x[cp$cp[i - 1]:cp$cp[i]])] <- 1
  }

  ret
}

# @fn Data.Norm
# @param x : x axis vector
# @param x : y axis vector
# @return list of normalized x and y
Data.Norm <- function(x, y) {
  ret.x <- rep(0, times = length(x))
  ret.y <- rep(0, times = length(y))

  ret.x[1] <- x[1]
  ret.y[1] <- y[1]

  diff.x <- diff(x)
  diff.y <- diff(y)

  scale.x <- 1 / mean(abs(diff(x)))
  scale.y <- 1 / mean(abs(diff(y)))

  # mcat(scale.x)
  # mcat(scale.y)

  for (i in 2:length(x)) {
    ret.x[i] <- ret.x[i - 1] + (diff.x[i - 1] * scale.x)
    ret.y[i] <- ret.y[i - 1] + (diff.y[i - 1] * scale.y)
  }

  list(x = ret.x, y = ret.y)
}

# @fn Change.Point
# @param x : x axis vector
# @param x : y axis vector
# @return change point index in x vector,
# Note that the first and the last point will be always change point
Change.Point <- function(x, y) {
  if (length(x) != length(y)) {
    mcat("ERROR : The length of x and y should be the same")
    return(0)
  }

  n_bins <- length(x)
  Fv <- rep(NA_real_, times = n_bins)
  Ev <- rep(NA_real_, times = n_bins)
  cp <- 1

  i <- 1
  Fv[1] <- 0
  while (i < n_bins) {
    j <- i + 1
    Fv[j] <- sqrt((x[j] - x[i]) ^ 2 + (y[j] - y[i]) ^ 2)

    while (j < n_bins) {
      j <- j + 1
      k <- (i + 1):(j - 1)
      Ev[j] <- (sum(abs((y[j] - y[i]) * x[k] - (x[j] - x[i]) * y[k] - (x[i] * y[j]) + (x[j] * y[i]))) / sqrt((x[j] - x[i]) ^ 2 + (y[j] - y[i]) ^ 2))
      Fv[j] <- sqrt((x[j] - x[i]) ^ 2 + (y[j] - y[i]) ^ 2) - (sum(abs((y[j] - y[i]) * x[k] - (x[j] - x[i]) * y[k] - (x[i] * y[j]) + (x[j] * y[i]))) / sqrt((x[j] - x[i]) ^ 2 + (y[j] - y[i]) ^ 2))

      #################################################
      # Not Original Code
      if (is.na(Fv[j]) || is.na(Fv[j - 1])) {
        j <- j - 1
        cp <- c(cp, j)
        break
      }
      #################################################### 3
      if (Fv[j] < Fv[j - 1]) {
        j <- j - 1
        cp <- c(cp, j)
        break
      }
    }
    i <- j
  }

  cp <- c(cp, n_bins)

  list(cp = cp, objF = Fv, errF = Ev)
}

# @fn Get.Pvalue
# @param matrix.data : matrix
# @param size : size to extend
# @param scale : scale parameter if necessary. deprecated parameter
# @return computed p-value vector
#' @importFrom stats wilcox.test
Get.Pvalue <- function(matrix.data, size, scale = 1.0) {
  n_bins <- nrow(matrix.data)
  pvalue <- rep(1, times = n_bins)

  for (i in seq_len(n_bins - 1)) {
    dia <- as.vector(Get.Diamond.Matrix2(matrix.data, i = i, size = size))
    ups <- as.vector(Get.Upstream.Triangle(matrix.data, i = i, size = size))
    downs <- as.vector(Get.Downstream.Triangle(matrix.data, i = i, size = size))

    wil.test <- wilcox.test(x = dia * scale, y = c(ups, downs), alternative = "less", exact = FALSE)
    pvalue[i] <- wil.test$p.value

    # mcat(i, " = ", wil.test$p.value)
  }

  pvalue[is.na(pvalue)] <- 1
  pvalue
}

# @fn Get.Upstream.Triangle
# @param mat.data : matrix data
# @param i : bin index
# @param size : size of window to extend
# @return upstream triangle matrix
Get.Upstream.Triangle <- function(mat.data, i, size) {
  n_bins <- nrow(mat.data)

  lower <- max(1, i - size)
  tmp.mat <- mat.data[lower:i, lower:i]
  tmp.mat[upper.tri(tmp.mat, diag = FALSE)]
}

# @fn Get.Downstream.Triangle
# @param mat.data : matrix data
# @param i : bin index
# @param size : size of window to extend
# @return downstream triangle matrix
Get.Downstream.Triangle <- function(mat.data, i, size) {
  n_bins <- nrow(mat.data)
  if (i == n_bins) {
    na <- NA_real_
    storage.mode(na) <- storage.mode(mat.data)
    return(na)
  }

  upperbound <- min(i + size, n_bins)
  tmp.mat <- mat.data[(i + 1):upperbound, (i + 1):upperbound]
  tmp.mat[upper.tri(tmp.mat, diag = FALSE)]
}

# @fn Get.Diamond.Matrix2
# @param mat.data : matrix data
# @param i : bin index
# @param size : size of window to extend
# @return diamond matrix
Get.Diamond.Matrix2 <- function(mat.data, i, size) {
  n_bins <- nrow(mat.data)
  na <- NA_real_
  storage.mode(na) <- storage.mode(mat.data)
  new.mat <- matrix(rep(na, times = size * size), nrow = size, ncol = size)

  for (k in seq_len(size)) {
    if (i - (k - 1) >= 1 && i < n_bins) {
      lower <- min(i + 1, n_bins)
      upper <- min(i + size, n_bins)

      new.mat[size - (k - 1), seq_len(upper - lower + 1)] <- mat.data[i - (k - 1), lower:upper]
    }
  }

  new.mat
}

# @fn Convert.Bin.To.Domain
# @param bins : bin information
# @param signal.idx : signal index
# @param signal.idx : gap index
# @param pvalues : pvalue vector
# @param pvalue.cut : pvalue threshold
# @return dataframe storing domain information
Convert.Bin.To.Domain <- function(bins, signal.idx, gap.idx, pvalues = NULL, pvalue.cut = NULL) {
  n_bins <- nrow(bins)
  ret <- data.frame(chr = character(0), from.id = numeric(0), from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0), tag = character(0), size = numeric(0))
  levels(x = ret[, "tag"]) <- c("domain", "gap", "boundary")

  rmv.idx <- setdiff(seq_len(n_bins), gap.idx)
  proc.region <- Which.process.region(rmv.idx, n_bins = n_bins, min.size = 0)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  n_procs <- nrow(proc.region)
  gap <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = rep(0, times = n_procs), from.coord = from.coord, to.id = rep(0, times = n_procs), to.coord = rep(0, times = n_procs), tag = rep("gap", times = n_procs), size = rep(0, times = n_procs), stringsAsFactors = FALSE)

  rmv.idx <- union(signal.idx, gap.idx)
  proc.region <- Which.process.region(rmv.idx, n_bins = n_bins, min.size = 0)
  n_procs <- nrow(proc.region)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  domain <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = rep(0, times = n_procs), from.coord = from.coord, to.id = rep(0, times = n_procs), to.coord = rep(0, times = n_procs), tag = rep("domain", times = n_procs), size = rep(0, times = n_procs), stringsAsFactors = FALSE)

  rmv.idx <- setdiff(seq_len(n_bins), signal.idx)
  proc.region <- as.data.frame(Which.process.region(rmv.idx, n_bins = n_bins, min.size = 1))
  n_procs <- nrow(proc.region)
  if (n_procs > 0) {
    from.coord <- bins[proc.region[, "start"] + 1, "from.coord"]
    boundary <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = rep(0, times = n_procs), from.coord = from.coord, to.id = rep(0, times = n_procs), to.coord = rep(0, times = n_procs), tag = rep("boundary", times = n_procs), size = rep(0, times = n_procs), stringsAsFactors = FALSE)
    ret <- rbind(ret, boundary)
  }

  ret <- rbind(gap, domain)
  ret <- ret[order(ret[, 3]), ]

  ret[, "to.coord"] <- c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
  ret[, "from.id"] <- match(ret[, "from.coord"], table = bins[, "from.coord"])
  ret[, "to.id"] <- match(ret[, "to.coord"], table = bins[, "to.coord"])
  ret[, "size"] <- ret[, "to.coord"] - ret[, "from.coord"]

  if (!is.null(pvalues) && !is.null(pvalue.cut)) {
    for (i in seq_len(nrow(ret))) {
      if (ret[i, "tag"] == "domain") {
        domain.bins.idx <- ret[i, "from.id"]:ret[i, "to.id"]
        p.value.constr <- which(pvalues[domain.bins.idx] < pvalue.cut)

        if (length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] <- "boundary"
      }
    }
  }

  new.bdr.set <- data.frame(chr = character(0), from.id = numeric(0), from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0), tag = character(0), size = numeric(0))
  stack.bdr <- data.frame(chr = character(0), from.id = numeric(0), from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0), tag = character(0), size = numeric(0))

  i <- 1
  while (i <= nrow(ret)) {
    if (ret[i, "tag"] == "boundary") {
      stack.bdr <- rbind(stack.bdr, ret[i, ])
    } else if (nrow(stack.bdr) > 0) {
      new.bdr <- data.frame(
        chr = bins[1, "chr"],
        from.id = min(stack.bdr[, "from.id"]),
        from.coord = min(stack.bdr[, "from.coord"]),
        to.id = max(stack.bdr[, "to.id"]),
        to.coord = max(stack.bdr[, "to.coord"]),
        tag = "boundary",
        size = max(stack.bdr[, "to.coord"]) - min(stack.bdr[, "from.coord"])
      )
      new.bdr.set <- rbind(new.bdr.set, new.bdr)
      stack.bdr <- data.frame(chr = character(0), from.id = numeric(0), from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0), tag = character(0), size = numeric(0))
    }

    i <- i + 1
  }

  ret <- rbind(ret[ret[, "tag"] != "boundary", ], new.bdr.set)
  ret <- ret[order(ret[, "to.coord"]), ]

  ret
}


# @fn Convert.Bin.To.Domain
# @param bins : bin information
# @param signal.idx : signal index
# @param signal.idx : gap index
# @param pvalues : pvalue vector
# @param pvalue.cut : pvalue threshold
# @return dataframe storing domain information
Convert.Bin.To.Domain.TMP <- function(bins, signal.idx, gap.idx, pvalues = NULL, pvalue.cut = NULL) {
  n_bins <- nrow(bins)
  ret <- data.frame(chr = character(0), from.id = numeric(0), from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0), tag = character(0), size = numeric(0))
  levels(x = ret[, "tag"]) <- c("domain", "gap", "boundary")

  rmv.idx <- setdiff(seq_len(n_bins), gap.idx)
  proc.region <- Which.process.region(rmv.idx, n_bins = n_bins, min.size = 0)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  n_procs <- nrow(proc.region)
  gap <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = rep(0, times = n_procs), from.coord = from.coord, to.id = rep(0, times = n_procs), to.coord = rep(0, times = n_procs), tag = rep("gap", times = n_procs), size = rep(0, times = n_procs), stringsAsFactors = FALSE)

  rmv.idx <- union(signal.idx, gap.idx)
  proc.region <- Which.process.region(rmv.idx, n_bins = n_bins, min.size = 0)
  n_procs <- nrow(proc.region)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  domain <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = rep(0, times = n_procs), from.coord = from.coord, to.id = rep(0, times = n_procs), to.coord = rep(0, times = n_procs), tag = rep("domain", times = n_procs), size = rep(0, times = n_procs), stringsAsFactors = FALSE)

  rmv.idx <- setdiff(seq_len(n_bins), signal.idx)
  proc.region <- as.data.frame(Which.process.region(rmv.idx, n_bins = n_bins, min.size = 1))
  n_procs <- nrow(proc.region)
  if (n_procs > 0) {
    from.coord <- bins[proc.region[, "start"] + 1, "from.coord"]
    boundary <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = rep(0, times = n_procs), from.coord = from.coord, to.id = rep(0, times = n_procs), to.coord = rep(0, times = n_procs), tag = rep("boundary", times = n_procs), size = rep(0, times = n_procs), stringsAsFactors = FALSE)
    ret <- rbind(ret, boundary)
  }

  ret <- rbind(gap, domain)
  ret <- ret[order(ret[, 3]), ]

  ret[, "to.coord"] <- c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
  ret[, "from.id"] <- match(ret[, "from.coord"], table = bins[, "from.coord"])
  ret[, "to.id"] <- match(ret[, "to.coord"], table = bins[, "to.coord"])
  ret[, "size"] <- ret[, "to.coord"] - ret[, "from.coord"]

  if (!is.null(pvalues) && !is.null(pvalue.cut)) {
    for (i in seq_len(nrow(ret))) {
      if (ret[i, "tag"] == "domain") {
        domain.bins.idx <- ret[i, "from.id"]:ret[i, "to.id"]
        p.value.constr <- which(pvalues[domain.bins.idx] < pvalue.cut)

        if (length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] <- "boundary"
      }
    }
  }

  ret
}


#' @importFrom utils capture.output
mcat <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "", file = stderr())
}
