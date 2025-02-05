#' This function builds the Consensus Network from the component network
#'
#' This function builds a consensus co-expression network.
#'
#' TODO add consensus by weight instead of rank?
#' TODO better save format for bicNetworks
#'
#' @param network_files A vector of file paths pointing to each network matrix.
#' @param exprData A matrix of the gene expression data used to generate the
#'   networks. Rows must be samples and columns must be genes.
#' @param outputpath Path to the folder where the consensus network should be
#'   saved.
#' @param max_edges Optional. The maximum number of network edges to use when
#'   computing the BIC curve.
#'
#' @return Nothing. Saves the rank consensus network to
#'   \code{<outputpath>/rankConsensusNetwork.csv} and saves the BICNetwork
#'   object to \code{<outputpath>/bicNetworks.rds}
#'
#' @export
buildConsensus <- function(network_files,
                           exprData,
                           outputpath,
                           max_edges = 2e5) {
  ranked_network <- rankConsensus(exprData, network_files)

  writeCSVFile(ranked_network,
               filename = file.path(outputpath, "rankConsensusNetwork.csv"))

  bicNetworks <- pruneBICnetwork(ranked_network,
                                 exprData = exprData,
                                 maxEdges = max_edges)

  saveRDS(bicNetworks, file = file.path(outputpath, "bicNetworks.rds"))
}
