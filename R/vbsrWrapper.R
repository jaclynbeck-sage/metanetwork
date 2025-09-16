#' Runs variational Bayes spike regression across a gene expression matrix
#'
#' This function wraps variable Bayes spike regression, using a single gene's
#' expression across all samples as the response variable, and the expression of
#' the other genes as the input.
#'
#' @inheritParams glmnetIC
#' @param fdr Optional. FDR threshold cut off for edge determination.)
#' @param ... Other parameters accepted by [vbsr::vbsr]
#'
#' @return A named list where "Z" is a network matrix constructed from the
#'   Z-values of the fit, and "2Z" is a network matrix where values from "Z"
#'   with an adjusted p-value less than the FDR threshold are re-calculated, and
#'   all other values are set to 0. P-values are adjusted using the
#'   Benjamini-Hochberg correction.
#' @export
vbsrWrapper <- function(x, y, fdr = 0.05, ...) {
  # JB TODO why are we using Z values and not beta values?
  result <- vbsr::vbsr(y = y, X = x, ...)$z
  result_2z <- rep(0, length(result))

  # JB TODO vbsr returns a list including $pval, but it is slightly different
  # than this calculation. Which is correct?
  pval <- stats::pchisq(result^2, 1, lower.tail = F)
  pval <- stats::p.adjust(pval, method = "BH", n = length(pval))

  if (sum(pval < fdr) > 0) {
    newz <- fastlm_z(y, x[, pval < fdr])
    result_2z[pval < fdr] <- newz
  }

  networks <- rbind(result, result_2z)
  colnames(networks) <- colnames(x)
  rownames(networks) <- c("Z", "2Z")

  return(networks)
}
