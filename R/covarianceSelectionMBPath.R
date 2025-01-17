#' Covariance based on Neighborhood Selection
#'
#' This function defines the covariance neighborhood between gene-gene edges
#' in an expression matrix. TODO
#'
#' @param X An expression matrix where rows are samples and columns are genes.
#' @param rankedEdges An edge list in the form of a matrix where column one is
#'   gene one and column 2 is gene two
#' @param startI Optional. Start at the first edge in `rankedEdges` (Default = 1)
#'
#' @return A list object of containing the BIC estimate, bicNeighborhood (?) , neighborhoods (?), flag (?)
#' @export
covarianceSelectionMBPath = function(X, rankedEdges, startI = 1) {
  nedges <- nrow(rankedEdges)
  bic <- rep(0, nedges)

  neighborhoods <- vector('list', ncol(X))
  names(neighborhoods) <- colnames(X)

  bicNeighborhood <- apply(X, 2, fastlm_bic, correction = ncol(X))
  names(bicNeighborhood) <- colnames(X)
  bicCurrent <- sum(bicNeighborhood, na.rm = TRUE)

  flag <- c()

  for (count in 1:nedges) {
    if (count %% 1000 == 0) {
      cat('Count:', count, 'BIC:', bicCurrent, '\n')
    }
    gene1 <- colnames(X)[rankedEdges[count, 1]]
    gene2 <- colnames(X)[rankedEdges[count, 2]]

    neighborhoods[[gene1]] <- c(neighborhoods[[gene1]], gene2)
    neighborhoods[[gene2]] <- c(neighborhoods[[gene2]], gene1)

    # If starting on an edge other than 1, this re-calculates bicNeighborhood
    # with the existing set of neighborhood genes
    if (count == startI) {
      for (i in 1:ncol(X)) {
        bicNeighborhood[i] <- fastlm_bic(X[, i], X[, neighborhoods[[i]]],
                                        correction = ncol(X))
      }
      bicCurrent <- sum(bicNeighborhood, na.rm = TRUE)

    } else if (count > startI) {
      bicGene1 <- NA
      bicGene2 <- NA

      try(bicGene1 <- fastlm_bic(X[, gene1], X[, neighborhoods[[gene1]]],
                                correction = ncol(X)))
      try(bicGene2 <- fastlm_bic(X[, gene2], X[, neighborhoods[[gene2]]],
                                correction = ncol(X)))

      if (!is.na(bicGene1)) {
        bicNeighborhood[gene1] <- bicGene1

      } else {
        flag <- rbind(flag, c(gene1, count))
      }

      if (!is.na(bicGene2)) {
        bicNeighborhood[gene2] <- bicGene2

      } else {
        flag <- rbind(flag, c(gene2, count))
      }

      bicCurrent <- sum(bicNeighborhood, na.rm = TRUE)
    }

    bic[count] <- bicCurrent
  }

  return(list(bic = bic,
              bicNeighborhood = bicNeighborhood,
              neighborhoods = neighborhoods,
              flag = flag))
}
