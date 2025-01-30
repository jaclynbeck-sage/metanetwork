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
#' @param network_files A vector or list of file paths pointing to each individual
#' network CSV that should be in the consensus. The contents of each file should
#' be a square matrix with the columns and rows labeled with the gene names. All
#' matrices must have the same sets of genes and be the same dimensions.
#'
#' @return A new NxN matrix that represents the consensus rank of each edge in
#' the network.
#'
#' @export
rankConsensus <- function(network_files) {
  # Load the first row of the first network file to get the expected dimensions
  # and column/row names
  gene_names <- data.table::fread(network_files[1], nrows = 1)
  gene_names <- colnames(gene_names)[-1]

  aggregateRank <- 0

  # Sum the ranks for each network into one vector
  for (net_file in network_files) {
    cat("Ranking", basename(net_file), "...\n")
    network <- loadCSVFile(net_file)

    # Ensure it's a square matrix with the same genes as the first network in the list
    stopifnot(ncol(network) == nrow(network))
    stopifnot(ncol(network) == length(gene_names))
    stopifnot(all(colnames(network) %in% gene_names))
    stopifnot(all(!is.na(network)))

    # Ensure it's in the same order as the first network
    network <- network[gene_names, gene_names]
    network <- abs(network[upper.tri(network)])

    # Rank network edges in descending order. network[upper.tri(network)]
    # collapses the upper triangle into a vector so lower triangle & diagonal
    # zeros are ignored. frankv() is significantly faster than rank().
    collapsedRank <- data.table::frankv(network, order = -1)
    collapsedRank <- bit64::as.integer64(collapsedRank)

    aggregateRank <- aggregateRank + collapsedRank
    rm(network, collapsedRank)
  }

  # Rank the aggregated ranks in descending order and normalize to max value --
  # largest aggregate rank should have the smallest value, smallest aggregate
  # rank should have the largest
  finalRank <- data.table::frankv(aggregateRank, order = -1)
  finalRank <- finalRank / max(finalRank)

  network <- matrix(0,
                    nrow = length(gene_names),
                    ncol = length(gene_names),
                    dimnames = list(gene_names, gene_names))

  network[upper.tri(network)] <- finalRank

  return(network)
}
