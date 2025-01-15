#' This function builds the Consensus Network from the component network
#'
#' This function builds a consensus co-expression network. TODO
#'
#' @param outputpath Local directory to load files from and save files to.
#' @param network_files A vector of file paths pointing to each network matrix.
#' @param fileName Required. The file path to the gene expression file. The gene
#' expression matrix should have rows as samples and columns as genes.
#' TODO input the matrix as an argument instead of reading it in here.
#'
#' @export
#' @return Saves a rank consensus network to `outputpath` and saves the BICNetwork
#'  object to `outputpath`
#'
buildConsensus = function(outputpath, network_files, fileName) {
  loadNetwork <- function(file) {
    network <- data.table::fread(file = file, sep = ",",
                                 header = TRUE, data.table = FALSE)
    network <- tibble::column_to_rownames(network, var = colnames(network)[1])

    return(network)
  }

  # JB TODO this needs a LOT of memory, think of something more efficient
  networks <- lapply(network_files, loadNetwork)
  networks <- lapply(networks, data.matrix)

  ranked_network <- rankConsensus(networks)

  data.table::fwrite(ranked_network,
                     file = file.path(outputpath, "rankConsensusNetwork.csv"),
                     sep = ",",
                     row.names = TRUE,
                     col.names = TRUE)

  dataSet <- data.table::fread(file = fileName, sep = ",",
                               header = TRUE, data.table = FALSE)
  dataSet <- tibble::column_to_rownames(dataSet, var = colnames(dataSet)[1])
  dataSet <- as.matrix(dataSet)

  bicNetworks <- computeBICcurve(ranked_network, dataSet, maxEdges = 2e5)

  saveRDS(bicNetworks, file = file.path(outputpath, "bicNetworks.rds"))
}
