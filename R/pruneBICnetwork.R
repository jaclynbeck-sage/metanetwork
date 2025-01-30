#' Prune a network using a BIC Curve
#'
#' This function gets a list of edges from the input network, creates a new
#' empty network, iteratively adds edges to the new network, and computes the BIC
#' after adding each edge. The final network will be the set of edges that
#' produce the smallest BIC.
#'
#' @param network A matrix containing network data. The matrix should be square
#'   and either upper triangular or symmetric. Can be a sparse or dense matrix.
#' @param exprData Expression matrix where rows are samples and columns are
#'   genes
#' @param maxEdges Optional. Maximum number of edges to use.
#'
#' @return  A named list containing: "network" = a sparse matrix where values
#'   are TRUE if that edge is in the list of edges that produce the smallest
#'   BIC, otherwise FALSE. "bicMin" = the smallest BIC value, "bicPath" = a
#'   vector containing the BIC values at each iteration.
#'
#' @export
pruneBICnetwork <- function(network, exprData, maxEdges = 2e5) {
  # Create sparse upper-triangle matrix
  network <- abs(Matrix::triu(network))
  network <- Matrix::Matrix(network, sparse = TRUE)

  # fsort is much faster than sort(), but only if decreasing = FALSE, so we
  # sort ascending and set the threshold to <end> - maxEdges to capture the
  # <maxEdges> largest edges. network@x is the non-zero values of network.
  thresh_ind <- max(Matrix::nnzero(network) - maxEdges, 0)
  thresVal <- data.table::fsort(network@x)[thresh_ind]

  # It's possible for the threshold to be 0 -- set threshold to the smallest
  # non-zero value if so.
  if (thresVal == 0) {
    thresVal <- min(network@x)
  }

  cat('Edge threshold:', thresVal, '\n')

  edgeList <- as.data.frame(Matrix::mat2triplet(network))
  edgeList <- subset(edgeList, edgeList$x > thresVal)

  colnames(edgeList) <- c('node1', 'node2', 'weight')

  edgeList <- edgeList[order(edgeList$weight, decreasing = TRUE),]

  bicPath <- computeBICpath(exprData, rankedEdges = edgeList)

  network <- network >= edgeList$weight[which.min(bicPath$bic)]

  return(list(network = Matrix::Matrix(network, sparse = TRUE),
              bicMin = min(bicPath$bic),
              bicPath = bicPath$bic))
}
