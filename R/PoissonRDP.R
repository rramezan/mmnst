#' Poisson Recursive Dyadic Partitioning
#'
#'
#' Calculates the piecewise constant intensity function using penalized likelihood recursive dyadic partitioning method in a Poisson model.
#'
#' @param sig a numeric vector determining the raw data, e.g., the heights (counts) of the histograms to be smoothed.
#' @param gamma a scalar determining the penalty factor in the penalized likelihood method.
#'
#' @return A list containing two numeric vectors:
#' 1. `est`: The estimated intensity function, represented as the RDP piecewise-constant estimate \eqn{c(t)}.
#' 2. `splitvec`: A binary vector of length \eqn{2^J - 1}, where $2^J$ is the length of the input `sig`. Its entries correspond to the internal nodes of the complete dyadic partition tree, ordered breadth-first: starting with the root node, followed by nodes at successively finer scales from left to right within each level.
#' An entry `splitvec[i] = 1` indicates that the corresponding node is retained as a split, so that its two child intervals remain separate. An entry `splitvec[i] = 0` indicates that the node does not split.
#'
#' @source This code is an R translation of the [MATLAB code](http://math.bu.edu/people/kolaczyk/software/msglmcode.zip) from Kolaczyk and Nowak (2005) adjusting for piecewise constant intensity function; i.e., m = 0 in the original MATLAB code.
#'
#' @references Kolaczyk, E.D. and Nowak, R.D. (2004). Multiscale likelihood analysis and complexity penalized estimation, *The Annals of Statistics*, **32**(2), 500-527. doi: 10.1214/009053604000000076.
#'
#' Kolaczyk, E.D. and Nowak, R.D. (2005). Multiscale generalized linear models for nonparametric function estimation. *Biometrika*,  **92**(1), 119–133. doi: 10.1093/biomet/92.1.119.
#'
#' @export

PoissonRDP <- function(sig, gamma) {
  # This is a translation and modification (computationally faster version) of the MATLAB code from Nowak and Kolaczyk (2005)
  # adjusting for piecewise constant, not polynomial, intensity function; i.e., m = 0 in the original MATLAB code
  # note that the length of the vector sig must be a power of 2

  if (sum(sig) < 1) {
    # MODIFICATION: Return a list, including an all-zero hereditary split vector.
    return(
      list(
        est = rep(0, length(sig)),
        splitvec = integer(max(length(sig) - 1, 0))
      )
    )
  }

  n <- length(sig)

  if (log2(n) != round(log2(n))) {
    stop("The length of sig must be a power of 2")
  }

  J <- log2(n)

  zeropadding <- rep(0, 2^ceiling(J) - n)
  sumx <- c(sig, zeropadding)
  n2 <- length(sumx)

  lam <- log(n) * gam

  bestFit <- sig
  bestPL <- sig * log(pmax(bestFit, 1e-50)) - bestFit

  ind_pad <- matrix(1:length(sumx), 1)

  # MODIFICATION: Store split decisions across all scales.
  # The final ordering will be root, then level 1, level 2, ..., level J - 1.
  decor <- integer(0)

  for (j in (J - 1):0) {
    n <- 2^(J - j)
    dim(ind_pad) <- c(n, length(ind_pad) / n)

    sumx <- .Internal(colSums(sumx, 2L, length(sumx) / 2L, FALSE)) # if this line replaces the next one, the code is faster, but it will mess up R Check for CRAN because of the .Internal function.

    bestFit2 <- sumx / n

    pl0 <- sumx * log(pmax(bestFit2, 1e-50)) - bestFit2 * n

    pl1 <- .Internal(colSums(c(bestPL, zeropadding), n, n2 / n, FALSE)) # if this line replaces the next one, the code is faster, but it will mess up R Check for CRAN because of the .Internal function.

    pl1 <- pl1 * 2 / n - lam

    for (k in which(pl1[1:2^j] <= pl0[1:2^j])) {
      bestFit[ind_pad[, k]] <- bestFit2[k]
    }

    comp <- pl1 > pl0

    # Save split decisions at this scale.
    # At scale j, only the first 2^j entries correspond to valid intervals.
    decor_j <- as.integer(comp[seq_len(2^j)])

    # Prepend each coarser level so decor is in tree order:
    # root, level 1, level 2, ..., finest internal level.
    decor <- c(decor_j, decor)

    plmax <- pl1 * comp + pl0 * !comp

    bestPL <- rep(plmax, rep(n, length(pl0)))[1:length(bestPL)]
  }

  # Enforce hereditary splitting.
  # A node can remain split only if its own split decision and all ancestral
  # split decisions are equal to one.
  splitvec <- integer(length(decor))

  if (length(decor) > 0) {
    splitvec[1] <- decor[1]
    active_parent_splits <- splitvec[1]

    if (J >= 2) {
      for (j in 1:(J - 1)) {
        parent_splits <- rep(active_parent_splits, each = 2L)

        current_level_ind <- 2^j:(2^(j + 1) - 1)

        active_parent_splits <- parent_splits * decor[current_level_ind]

        splitvec[current_level_ind] <- active_parent_splits
      }
    }
  }

  # Return both the estimate and hereditary split indicators.
  return(
    list(
      est = bestFit,
      splitvec = splitvec
    )
  )
}
