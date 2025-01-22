#' Runs tigress on an expression matrix
#'
#' This function implements Trustful Inference of Gene REgulation with Stability
#' Selection (TIGRESS) algoritm. It is very similar to the function
#' \code{stabilityselection} in the R package \code{tigress}, except that the
#' \code{stabilityselection} function does not allow for changing some arguments
#' to \code{lars::lars} that significantly speed up processing.
#'
#' @param x A gene expression matrix where rows are samples and columns are genes
#' @param y An Nx1 vector of gene expression values where N = number of genes
#' @param nsteps_tigress Optional. The number of times the data should be split and sampled
#' @param nsteps_lars Optional. The number of steps that lars should use
#' @param alpha Optional. When sampling, weights are randomly uniformly generated
#' in the interval [alpha, 1].
#' @param ... Optional. Other arguments to pass to \code{lars::lars}
#'
#' @return A named vector of co-expression values of gene Y to columns of X
#' @export
tigress <- function(x,
                    y,
                    nsteps_tigress = 100,
                    nsteps_lars = 5,
                    alpha = 0.2,
                    ...) {
  n_samples <- length(y)
  n_genes <- ncol(x)
  halfsize <- floor(n_samples / 2)

  freq <- matrix(0, nsteps_lars, n_genes)

  for (i in 1:nsteps_tigress) {
    xs <- t(t(x) * stats::runif(n_genes, alpha, 1))

    indexVec <- sample(1:n_samples, n_samples)
    i1 = indexVec[1:halfsize]
    i2 = indexVec[(halfsize + 1):n_samples]

    result1 <- lars::lars(x = xs[i1, ],
                          y = y[i1],
                          type = "lar",
                          max.steps = nsteps_lars,
                          normalize = FALSE,
                          use.Gram = FALSE,
                          ...)
    freq <- freq + abs(sign(result1$beta[2:(nsteps_lars + 1), ]))

    result2 <- lars::lars(x = xs[i2, ],
                          y = y[i2],
                          type = "lar",
                          max.steps = nsteps_lars,
                          normalize = FALSE,
                          use.Gram = FALSE,
                          ...)
    freq <- freq + abs(sign(result2$beta[2:(nsteps_lars + 1), ]))
  }

  # Uses 2 * nsteps_tigress because of the two additions in each loop
  freq <- freq / (2 * nsteps_tigress)

  # colMeans() here is equivalent to scoring = "area" in tigress::stabilityselection
  return(colMeans(freq))
}
