#' Wrapper for WGCNA
#'
#' Runs \code{WGCNA::adjacency} and \code{WGCNA::TOMsimilarity} over a gene
#' expression matrix.
#'
#' @param data A gene expression matrix where rows are samples and columns are
#'   genes.
#' @param n_cores Optional. How many cores/threads to allow when running WGCNA.
#'   The default value of 1 means "do not thread".
#' @param ... Optional. Additional arguments accepted by either
#'   \code{WGCNA::pickSoftThreshold}, \code{WGCNA::adjacency}, or
#'   \code{WGCNA::TOMsimilarity}.
#'
#' @return a named list of two network matrices, where "adjacency" is the output
#'   of \code{WGCNA::adjacency} and "TOM" is the output of
#'   \code{WGCNA::TOMsimilarity}.
#'
#' @export
constructNetwork.wgcna <- function(data, n_cores = 1, ...) {
  if (n_cores > 1) {
    WGCNA::enableWGCNAThreads(nThreads = n_cores)
  }

  res <- R.utils::doCall(WGCNA::pickSoftThreshold,
                         data = data,
                         args = list(...),
                         .ignoreUnusedArgs = TRUE)
  print(paste("Power estimate:", res$powerEstimate))

  network_adj <- R.utils::doCall(WGCNA::adjacency,
                                 datExpr = data,
                                 power = res$powerEstimate,
                                 args = list(...),
                                 .ignoreUnusedArgs = TRUE)

  network_tom <- R.utils::doCall(WGCNA::TOMsimilarity,
                                 adjMat = network_adj,
                                 args = list(...),
                                 .ignoreUnusedArgs = TRUE)

  colnames(network_tom) <- colnames(network_adj)
  rownames(network_tom) <- colnames(network_adj)

  if (n_cores > 1) {
    WGCNA::disableWGCNAThreads()
  }

  return(list(adjacency = network_adj,
              TOM = network_tom))
}
