#' Winsorize a Gene Expression Matrix
#'
#' Winsorizes a given gene expression matrix.
#'
#' The top and bottom 1% of values for each gene will be replaced with the .99
#' and .01 quantile values for that gene, respectively. All NA values for the
#' gene are replaced with the mean expression for that gene. The data is then
#' scaled and centered on 0.
#'
#' @param x Gene expression matrix where rows are genes and columns are samples
#'
#' @return A winsorized expression matrix, which has been transposed so that
#' rows are samples and columns are genes
#'
#' @export
winsorizeData <- function(x) {
  winsorize <- function(x, per = .99) {
    up <- stats::quantile(x, per, na.rm = T)
    low <- stats::quantile(x, 1 - per, na.rm = T)
    x[x >= up] <- up
    x[x <= low] <- low
    return(x)
  }

  replaceNaMean <- function(x) {
    if (sum(is.na(x)) > 0) {
      y <- x
      y[is.na(x)] <- mean(x, na.rm = T)
      return(y)
    } else {
      return(x)
    }
  }

  x <- t(x)
  x <- apply(x, 2, winsorize)
  x <- apply(x, 2, replaceNaMean)
  x <- scale(x)
  return(x)
}
