#' This Function Computes a Network BIC Curve
#'
#' (?) TODO
#'
#' @param network A network object
#' @param exprData Expression matrix where rows are samples and columns are genes
#' @param maxEdges Optional. Maximum number of edges to use (Default = NULL)
#'
#' @return  A list object containing values of a sparse network, the best fit
#' minimum BIC, and TODO
#' @export
#'
computeBICcurve <- function(network, exprData, maxEdges = NULL) {
  if (is.null(maxEdges)) {
    maxEdges <- round((nrow(exprData) * ncol(exprData)) / 20)
  }

  maxEdges <- min(maxEdges, round((nrow(exprData) * ncol(exprData)) / 20))
  cat('maxEdges:', maxEdges, '\n')

  network <- data.matrix(network)
  upper_net <- network[which(upper.tri(network))]
  upper_net <- abs(upper_net)

  network[(lower.tri(network))] <- 0
  diag(network) <- 0

  thresVal <- sort(upper_net, decreasing = T)[min(maxEdges, length(upper_net))]
  cat('threshold:', thresVal, '\n')

  # Add in check for zero edges -- set threshold to the smallest non-zero value
  # if so.
  if (thresVal == 0) {
    thresVal <- min(upper_net[which(upper_net > 0)])
  }

  network <- abs(network)
  edgeList <- which(network >= thresVal, arr.ind = TRUE)
  edval <- network[which(network >= thresVal)]

  edgeList <- cbind(edgeList, edval)
  colnames(edgeList) <- c('node1', 'node2', 'weight')
  rownames(edgeList) <- paste0('e', 1:nrow(edgeList))

  edgeList <- data.frame(edgeList, stringsAsFactors = F)
  edgeList <- dplyr::arrange(edgeList, dplyr::desc(.data$weight))

  bicPath <- covarianceSelectionMBPath(data.matrix(exprData),
                                       rankedEdges = edgeList[, 1:2],
                                       startI = 1)
  bicPath2 <- NA

  network <- network >= edgeList$weight[which.min(bicPath$bic)]
  return(list(network = Matrix::Matrix(network, sparse = T),
              bicMin = min(bicPath$bic),
              bicPath = bicPath$bic))
}
