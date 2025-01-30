#' Fast Linear Modeling
#'
#' This function returns results from a fast linear model TODO I think this
#' assumes that ncol(X) < nrow(X) but that's never checked for in the code that
#' calls this function
#'
#' This function solves \code{Ax = b} for non-square matrix A as:
#'
#' \code{x = Ahat * b},
#'
#' where \code{Ahat = (t(A) * A)^-1 * t(A)}.
#'
#' Here, A = X and b = y.
#'
#' @param y A vector of response values
#' @param X A matrix of the same number of observations/rows as y
#' @return A named list containing: "betahat" = coefficients, and "ginv" = the
#'   inverse of t(X) %*% X
fastlm <- function(y, X) {
  ginv <- solve(t(X) %*% X)
  Xhat <- ginv %*% t(X)
  betahat <- Xhat %*% y

  return(list(betahat = betahat,
              ginv = ginv))
}

#' Z value for fastlm
#'
#' This function runs \code{fastlm} and returns the Z-values associated with the
#' solution.
#'
#' @inheritParams fastlm
#'
#' @return A vector of Z-values
fastlm_z <- function(y, X) {
  # If x is a single vector, this makes it an n x 1 matrix
  X <- as.matrix(X)
  n1 <- nrow(X)
  X <- cbind(rep(1, n1), X)

  results <- fastlm(y, X)

  sig <- (mean((y - X %*% results$betahat)^2)) * ((n1) / (n1 - ncol(X)))
  zval <- results$betahat / (sqrt(sig * (diag(results$ginv))))

  return(zval[-1])
}


#' Fast Linear Modeling BIC
#'
#' This function runs \code{fastlm} and calculates the BIC score based on the
#' solution. TODO
#'
#' @param x Optional. A numeric vector or matrix of model coefficients. If not set x
#' becomes a vector of integers the from 1 to length(y). (?)
#' @param y A numeric vector of model coefficients. If not set x
#' becomes a vector of integers the from 1 to length(y). (?)
#' @param correction A vector of correction factors, or a single numeric value.
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


#' Load CSV File using data.table
#'
#' This function loads data from a CSV file using \code{data.table::fread},
#' which is much faster than \code{read.table} or \code{read.csv}. The data can
#' optionally be coerced to a matrix or left as a data.frame.
#'
#' @param filename The path to the file to load
#' @param return_matrix Optional. If \code{TRUE}, the loaded data object will be
#' coerced to a matrix. If \code{FALSE}, the data object will be a \code{data.frame}.
#' @param ... Optional, other arguments to pass to \code{fread()}.
#'
#' @returns either a matrix or a data.frame, depending on the value of \code{return_matrix}
loadCSVFile <- function(filename, return_matrix = TRUE, ...) {
  object <- data.table::fread(file = filename, sep = ",", ...)
  object <- tibble::column_to_rownames(object, var = colnames(object)[1])

  if (return_matrix) {
    object <- data.matrix(object)
  }

  return(object)
}


#' Write CSV File
#'
#' Convenience wrapper for \code{data.table::fwrite}, which is faster than
#' \code{write.csv}.
#'
#' @param object The object to write to disk, which should be of a type supported
#' by \code{fread}. This function assumes that the columns and rows are named.
#' @param filename The path and name of the file where the object should be saved
#' @param ... Optional. Additional arguments to \code{fread}.
#'
#' @returns Nothing
writeCSVFile <- function(object, filename, ...) {
  data.table::fwrite(object,
                     file = filename,
                     sep = ",",
                     row.names = TRUE,
                     col.names = TRUE,
                     ...)
}


#' Write upper triangular matrix
#'
#' Convenience wrapper for \code{writeCSVFile} which converts a matrix to an
#' upper triangular matrix before writing to a CSV file.
#'
#' @inheritParams writeCSVFile
#' @param object A square matrix or other object coercible to a matrix.
#'
#' @returns Nothing
writeUpperTri <- function(object, filename, ...) {
  writeCSVFile(object * upper.tri(object), filename, ...)
}
