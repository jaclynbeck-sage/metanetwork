#' Rank Consensus Networks
#'
#' This function ranks edges in list of network objects and combines the ranks
#' to create a consensus network.
#'
#' For each network in the list, its edges are ranked relative to each other.
#' The ranks from all networks are added together, and those values are then
#' ranked to create a final "consensus" of each edge. To conserve memory,
#' this function uses a list of files to work on networks one at a time, rather
#' than using a list of all networks loaded into memory at once.
#'
#' @param exprData A matrix of the gene expression data used to generate the
#'   networks. Rows must be samples and columns must be genes.
#' @param network_files A vector or list of file paths pointing to each
#'   individual network CSV that should be in the consensus. The contents of
#'   each file should be a square matrix with the columns and rows labeled with
#'   the gene names.
#'
#' @return A new NxN matrix, where N = the number of genes in \code{exprData},
#'   that represents the consensus rank of each edge in the network.
#'
#' @export
rankConsensus <- function(exprData, network_files) {
  gene_names <- colnames(exprData)

  aggregateRank <- 0

  # Sum the ranks for each network into one vector
  for (net_file in network_files) {
    cat("Ranking", basename(net_file), "...\n")
    network <- loadCSVFile(net_file)

    # Ensure it's a square matrix with no NAs
    stopifnot(ncol(network) == nrow(network),
              all(!is.na(network)))

    # Some network matrices might not have all genes in them, so we have to add missing genes back
    if (!(all(gene_names %in% colnames(network)))) {
      missing <- setdiff(gene_names, colnames(network))
      network <- cbind(network,
                       matrix(0, nrow = nrow(network), ncol = length(missing),
                              dimnames = list(rownames(network), missing)))
      network <- rbind(network,
                       matrix(0, nrow = length(missing), ncol = ncol(network),
                              dimnames = list(missing, colnames(network))))
    }

    stopifnot(ncol(network) == length(gene_names),
              all(colnames(network) %in% gene_names))

    # Ensure it's in the same order as gene_names
    network <- network[gene_names, gene_names]
    network <- abs(network[upper.tri(network)])

    # Rank network edges in ascending order: The highest weights in the network
    # get the highest rank. network[upper.tri(network)] collapses the upper
    # triangle into a vector so lower triangle & diagonal zeros are ignored.
    # Using "dense" for ties.method and subtracting 1 ensures that 0s in the
    # network get a rank of 0, and non-zero values in sparser networks don't get
    # assigned ranks that cluster at the top 99%. frankv() is significantly
    # faster than rank().
    collapsedRank <- data.table::frankv(network, ties.method = "dense") - 1

    # Avoid integer overflows by normalizing to 1
    collapsedRank <- collapsedRank / max(collapsedRank)

    aggregateRank <- aggregateRank + collapsedRank
    rm(network, collapsedRank)
  }

  # Rank the aggregated ranks to get a final consensus
  finalRank <- data.table::frankv(aggregateRank)
  finalRank <- finalRank / max(finalRank)

  network <- matrix(0,
                    nrow = length(gene_names),
                    ncol = length(gene_names),
                    dimnames = list(gene_names, gene_names))

  network[upper.tri(network)] <- finalRank

  return(network)
}
