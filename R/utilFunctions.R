#' Fast Linear Modeling
#'
#' This function returns results from a fast linear model TODO
#'
#' @param y Required. A vector of response values
#' @param X Required. A matrix of the same number of observations/rows
#' as y
#' @return A named list containing betahat and ginv
fastlm <- function(y, X) {
  ginv <- solve(t(X) %*% X)
  Xhat <- ginv %*% t(X)
  betahat <- Xhat %*% y

  return(list(betahat = betahat,
              ginv = ginv))
}

#' Z value for fastlm
#'
#' TODO
#'
#' @inheritParams fastlm
#'
#' @return A vector of Z-values
fastlm_z <- function(y, X) {
  # If x is a single vector, this makes it an n x 1 matrix
  X <- as.matrix(X)
  n1 <- nrow(X)
  X <- cbind(rep(1, n1, X))

  results <- fastlm(y, X)

  sig <- (mean((y - X %*% results$betahat)^2)) * ((n1) / (n1 - ncol(X)))
  zval <- results$betahat / (sqrt(sig * (diag(results$ginv))))

  return(zval[-1])
}


#' Fast Linear Modeling BIC
#'
#' This function deploys matrix operations to calculate a model BIC given a vector
#' of model coefficients. TODO
#'
#' @param x Optional. A numeric vector or matrix of model coefficients. If not set x
#' becomes a vector of integers the from 1 to length(y). (?)
#' @param y A numeric vector of model coefficients. If not set x
#' becomes a vector of integers the from 1 to length(y). (?)
#' @param correction A vector of correction factors
#'
#' @return BIC estimate
fastlm_bic <- function(y, x = NULL, correction = 1) {
  if (!is.null(x)) {
    X <- as.matrix(x)
    n1 <- nrow(X)
    X <- cbind(rep(1, n1), X)
  } else {
    n1 <- length(y)
    X <- as.matrix(rep(1, n1))
  }

  results <- fastlm(y, X)

  sig <- mean((y - X %*% results$betahat)^2)

  # Calculate and return the BIC TODO verify this is correct
  return(n1 * (log(sig) + 1 + log(2 * pi)) + (ncol(X) + 1) * log(n1 * correction))
}
