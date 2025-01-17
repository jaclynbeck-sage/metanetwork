#' This Function Computes a Network BIC Curve
#'
#' (?) TODO
#'
#' @param network A network object
#' @param exprData Expression matrix where rows are samples and columns are genes
#' @param maxEdges Optional. Maximum number of edges to use.
#'
#' @return  A list object containing values of a sparse network, the best fit
#' minimum BIC, and TODO
#' @export
#'
computeBICcurve <- function(network, exprData, maxEdges = 2e5) {
  maxEdges <- min(maxEdges, round((nrow(exprData) * ncol(exprData)) / 20))
  cat('Maximum network edges:', maxEdges, '\n')

  network <- abs(data.matrix(network))
  upper_net <- network[upper.tri(network)]

  network[lower.tri(network, diag = TRUE)] <- 0

  # fsort is much faster than sort(), but only if decreasing = FALSE, so we
  # sort ascending and set the threshold to <end> - maxEdges to capture the
  # <maxEdges> largest edges
  thresh_ind <- max(length(upper_net) - maxEdges, 0)
  thresVal <- data.table::fsort(upper_net)[thresh_ind]

  # It's possible for the threshold to be 0 -- set threshold to the smallest
  # non-zero value if so.
  if (thresVal == 0) {
    thresVal <- min(upper_net[upper_net > 0])
  }

  cat('Edge threshold:', thresVal, '\n')

  edgeList <- which(network >= thresVal, arr.ind = TRUE)
  edgeList <- as.data.frame(cbind(edgeList, network[network >= thresVal]))

  colnames(edgeList) <- c('node1', 'node2', 'weight')
  rownames(edgeList) <- paste0('e', 1:nrow(edgeList))

  edgeList <- edgeList[order(edgeList$weight, decreasing = TRUE),]

  bicPath <- covarianceSelectionMBPath(exprData,
                                       rankedEdges = edgeList)

  network <- network >= edgeList$weight[which.min(bicPath$bic)]

  return(list(network = Matrix::Matrix(network, sparse = T),
              bicMin = min(bicPath$bic),
              bicPath = bicPath$bic))
}
