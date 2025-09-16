#' Best AIC or BIC solution for Lasso or Ridge regression
#'
#' This function runs a glmnet() function with lasso or ridge regression and
#' returns the results with the best AIC or BIC. This function is called by both
#' [lassoIC] and [ridgeIC].
#'
#' For gaussian distributions,
#' `AIC = 2k + n*ln(sigma^2)` and
#' `BIC = ln(n)*k + n*ln(sigma^2)`,
#' where `k` is the number of parameters, `n` is the number of
#' observations, and `sigma^2 = sum[(y - residuals)^2] / n`.
#' `ln(sigma^2)` represents the maximum log-likelihood estimate.
#'
#' See https://en.wikipedia.org/wiki/Normal_distribution#Log-likelihood and
#' https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares
#'
#' @param x A design matrix where rows are samples and columns are genes. Can be
#'   in sparse or dense matrix format.
#' @param y The response variable (a single gene), an nx1 matrix where n = the
#'   number of rows in x (the number of samples).
#' @param alpha Optional. If alpha = 1, glmnet will run lasso regression. If
#'   alpha = 0, it will run ridge regression. If alpha > 0 and < 1, the
#'   regression will be a mixture between the two.
#' @param ... Optional. Other arguments accepted by [glmnet]
#'
#' @return a 2 x n_genes matrix, where row 1 is the coefficients for the
#'   solution with the best AIC, and row 2 is for the best BIC. The rows are
#'   named "AIC" and "BIC".
#'
#' @export
glmnetIC <- function(x, y, alpha = 1, ...) {
  res <- glmnet::glmnet(x = x, y = y, alpha = alpha, family = "gaussian", ...)

  # Log likelihood -- sum[(y - resid)^2)] / n can be re-written as mean([(y-resid)^2])
  resid <- y - cbind(1, x) %*% stats::coef(res)
  error <- Matrix::colMeans(resid^2)
  log_like <- log(error)

  n_obs <- nrow(x)
  n_param <- res$df

  # Degrees of freedom for ridge regression require calculation of non-zero
  # eigenvalues of t(X) * X. The eigenvalues should be calculated as
  # svd(t(x) %*% x)$d, however to avoid large matrix multiplication we calculate
  # it as svd(x)$d^2, which is roughly equivalent but much faster to calculate.
  if (alpha == 0) {
    eigen <- svd(x)$d^2
    n_param <- sapply(res$lambda, function(lambda) {
      return(sum(eigen / (eigen + lambda)))
    })
  }

  aic <- n_obs * log_like + n_param * 2
  bic <- n_obs * log_like + n_param * log(n_obs)

  # res$beta is the same as coef(res) without the intercept, which we assume to
  # be extremely close to zero.
  best_res <- rbind(res$beta[, which.min(aic)],
                    res$beta[, which.min(bic)])
  rownames(best_res) <- c("AIC", "BIC")

  return(best_res)
}


#' Best AIC or BIC solution for Lasso regression
#'
#' This function is a wrapper for [glmnetIC] that runs a glmnet() function with
#' lasso regression and returns the results with the best AIC or BIC.
#'
#' @inheritParams glmnetIC
#'
#' @inherit glmnetIC return
#'
#' @seealso [glmnetIC]
#'
#' @export
lassoIC <- function(x, y, ...) {
  return(glmnetIC(x, y, alpha = 1, ...))
}


#' Best AIC or BIC solution for Ridge Regression
#'
#' This function is a wrapper for [glmnetIC] that runs a glmnet() function with
#' ridge regression and returns the results with the best AIC or BIC.
#'
#' @inheritParams glmnetIC
#' @inherit glmnetIC return
#'
#' @seealso [glmnetIC]
#'
#' @export
ridgeIC <- function(x, y, ...) {
  return(glmnetIC(x, y, alpha = 0, ...))
}
