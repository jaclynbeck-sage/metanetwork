#' Ranks Consensus networks
#'
#' This function ranks a list of consensus network objects and returns the best
#' rank network. TODO
#'
#' @param networks Required. A list object containing an individual network as
#' a list entry.
#'
#' @return The best rank consensus network.
#'
#' @export
rankConsensus <- function(networks) {
  aggregateRankFunction <- function(network, upperTriIndices) {
    collapsedEdgeSet <- network[upperTriIndices]
    collapsedRank <- rank(-abs(collapsedEdgeSet), ties.method = 'min')
    collapsedRank <-  bit64::as.integer64(collapsedRank)
    return(collapsedRank)
  }

  upperTriIndices <- which(upper.tri(networks[[1]]))
  aggregateRank <- rep(0, length(upperTriIndices))

  for (i in 1:length(networks)) {
    aggregateRank <- aggregateRank + aggregateRankFunction(networks[[i]],
                                                           upperTriIndices)
  }

  aggregateRank <- -aggregateRank

  finalRank <- bit64::rank.integer64(aggregateRank)
  finalRank <- finalRank / max(finalRank)

  network <- matrix(0, nrow(networks[[1]]), ncol(networks[[1]]))
  colnames(network) <- colnames(networks[[1]])
  rownames(network) <- rownames(networks[[1]])

  network[upperTriIndices] <- finalRank

  return(network)
}
