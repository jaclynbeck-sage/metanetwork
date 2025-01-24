#' Find Modules with Megena Clustering
#'
#' This function constructs a network adjacency matrix using MEGENA.
#'
#' @param data A gene expression matrix where rows are samples and columns are
#'   genes. The matrix will be transposed to work with MEGENA.
#' @param n_cores Optional. The number of cores/threads to use in parallel
#'   computation. If n_cores = 1, MEGENA will not run in parallel.
#' @param log_file_path Optional. The folder path to where log files should be
#'   stored. Log files capture any output during parallel execution. If omitted,
#'   log files will be stored in the working directory.
#' @param doPerm Optional. Number of permutations for calculating FDRs for all
#'   correlation pairs. (Default = 10)
#' @param ... Optional. Additional arguments accepted by
#'   \code{MEGENA::calculate_correlation} or \code{MEGENA::calculate.PFN}.
#'
#' @return An NxN matrix where N is the number of genes used in the network,
#'   which may be less than the number of genes in \code{data}.
#'
#' @export
constructNetwork.megena <- function(data,
                                    n_cores = 1,
                                    log_file_path = ".",
                                    doPerm = 10,
                                    ...) {
  # Correlation data needed for calculate.PFN
  ijw <- MEGENA::calculate.correlation(t(data),
                                       doPerm = doPerm,
                                       output.corTable = FALSE,
                                       output.permFDR = FALSE,
                                       ...)

  # calculate.PFN uses foreach to do parallel processing and requires the
  # cluster to be created ahead of time
  clust <- NULL
  if (n_cores > 1) {
    log_file <- file.path(log_file_path, "megena_log.txt")
    clust <- parallel::makeCluster(n_cores, outfile = log_file)
    doParallel::registerDoParallel(clust)
  }

  edge_list <- MEGENA::calculate.PFN(ijw,
                                     doPar = n_cores > 1,
                                     num.cores = n_cores,
                                     keep.track = FALSE,
                                     ...)

  if (!is.null(clust)) {
    parallel::stopCluster(clust)
  }

  # Turn the edge list into an adjacency matrix
  all_genes <- unique(c(edge_list$row, edge_list$col))
  network <- matrix(0, nrow = length(all_genes), ncol = length(all_genes),
                    dimnames = list(all_genes, all_genes))

  for (ind in 1:nrow(edge_list)) {
    edge <- edge_list[ind, ]
    network[edge$row, edge$col] <- edge$weight
    network[edge$col, edge$row] <- edge$weight
  }

  return(network)
}
