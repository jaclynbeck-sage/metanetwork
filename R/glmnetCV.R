#' Best Cross-validation Solutions for Lasso or Ridge Regression
#'
#' This function runs a glmnet() with lasso or ridge regression and returns the
#' results from the fit for:
#' 1. The value of lambda that gives the minimum mean cross-validation error, and,
#' 2. The largest value of lambda such that error is within 1 standard error of the minimum.
#'
#' This function is called by both `lassoCV()` and `ridgeCV()`
#'
#' @inheritParams glmnetIC
#' @param folds Optional. The number of cross validation folds to partition the
#' data into. (Default = 10). Corresponds to the `nfolds` argument of cv.glmnet.
#'
#' @return a 2 x n_genes matrix, where row 1 is the coefficients for the solution
#' with the minimum CVM (#1 above), and row 2 is for the best SE1 (#2 above). The
#' rows are named "lambda.min" and "lambda.1se".
#'
#' @export
glmnetCV <- function(x, y, alpha = 1, folds = 10, ...) {
  res <- glmnet::cv.glmnet(x = x, y = y, nfolds = folds, alpha = alpha, ...)
  fit <- res$glmnet.fit

  best_res <- rbind(fit$beta[, which(res$lambda == res$lambda.min)],
                    fit$beta[, which(res$lambda == res$lambda.1se)])
  rownames(best_res) <- c("lambda.min", "lambda.1se")

  return(best_res)
}


#' Best Cross-validation Solutions for Lasso regression
#'
#' This function is a wrapper for `glmnetCV` that runs a glmnet() function with
#' lasso regression and returns the results with the best cross-validation
#' solutions.
#'
#' @inheritParams glmnetCV
#'
#' @inherit glmnetCV return
#'
#' @seealso [glmnetCV]
#'
#' @export
lassoCV <- function(x, y, ...) {
  return(glmnetCV(x, y, alpha = 1, ...))
}


#' Best Cross-validation Solutions for Ridge Regression
#'
#' This function is a wrapper for `glmnetCV` that runs a glmnet() function with
#' ridge regression and returns the results with the best cross-validation
#' solutions.
#'
#' @inheritParams glmnetCV
#' @inherit glmnetCV return
#'
#' @seealso [glmnetCV]
#'
#' @export
ridgeCV <- function(x, y, ...) {
  return(glmnetCV(x, y, alpha = 0, ...))
}
