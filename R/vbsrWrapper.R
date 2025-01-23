#' Runs variational Bayes spike regression across a gene expression matrix
#'
#' This function wraps variable bays spike regression of a genes expression across
#' a matrix of genes expressed in the same samples. Returns both the original
#' output and the output with a 2Z cutoff.
#'
#' @inheritParams glmnetIC
#' @param fdr Optional. FDR threshold cut off for edge determination. (Default = 0.05)
#' @param ... Other parameters accepted by `vbsr::vbsr()`
#'
#' @return A network matrix constructed from the Z-values of the fit
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
#' This function applies a user FDR threshold to input p-values
#'
#' @param pval Required. A vector of uncorrected P-Values.
#' @param fdr Optional. desired FDR cutoff. (Default = 0.05)
#' as y
#' @return Corrected PValues TODO this is wrong...
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
