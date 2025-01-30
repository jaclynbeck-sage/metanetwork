#' Runs variational Bayes spike regression across a gene expression matrix
#'
#' This function wraps variable Bayes spike regression, using a single gene's
#' expression across all samples as the response variable, and the expression of
#' the other genes as the input.
#'
#' @inheritParams glmnetIC
#' @param fdr Optional. FDR threshold cut off for edge determination.)
#' @param ... Other parameters accepted by \code{vbsr::vbsr()}
#'
#' @return A named list where "Z" is a network matrix constructed from the
#'   Z-values of the fit, and "2Z" is a network matrix where values with a
#'   p-value less than the FDR threshold are re-calculated. TODO
#' @export
vbsrWrapper <- function(x, y, fdr = 0.05, ...) {
  # JB TODO why are we using Z values and not beta values?
  result <- vbsr::vbsr(y = y, X = x, ...)$z
  result_2z <- result

  # JB TODO vbsr returns a list including $pval, but it is slightly different
  # than this calculation. Which is correct?
  pval <- stats::pchisq(result^2, 1, lower.tail = F)
  thres <- fdrThres(pval, fdr = fdr)

  # JB TODO it's possible that everything BUT where pval > thres is supposed to
  # be zeroed out?
  if (sum(pval < thres) > 0) {
    newz <- fastlm_z(y, x[, pval < thres])
    result_2z[pval < thres] <- newz
  }

  networks <- rbind(result, result_2z)
  colnames(networks) <- colnames(x)
  rownames(networks) <- c("Z", "2Z")

  return(networks)
}


#' FDR Threshold
#'
#' This function calculates a corrected FDR threshold for p-values, which is
#' corrected based on the number of p-values.
#'
#' @param pval A vector of uncorrected p-values.
#' @param fdr Optional. Desired FDR cutoff. (Default = 0.05)
#'
#' @return A corrected threshold for p-value cutoff
fdrThres <- function(pval, fdr = 0.05) {
  n <- length(pval)
  comp <- sort(pval) < ((fdr / n) * (1:n))

  if (min(pval) < (fdr / n)) {
    w1 <- which(!comp)[1]
    return((fdr / n) * w1)
  }
  else {
    return(fdr / n)
  }
}
