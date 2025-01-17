#' Find Global Modularity (Q)
#'
#' This function calculates global modularity (Q) of a upper triangular
#' adjacency matrix using \code{igraph::modularity}.
#'
#' @inheritParams compute.LocalModularity
#'
#' @return Q = modularity index.
#'
#' @export compute.Modularity
compute.Modularity <- function(adj, mod) {
  # Error functions
  if (!inherits(adj, "matrix")) {
    stop("Adjacency matrix should be of class matrix")
  }

  if (nrow(adj) != ncol(adj)) {
    stop("Adjacency matrix should be symmetric")
  }

  if (!all(adj[lower.tri(adj)] == 0)) {
    stop("Adjacency matrix should be upper triangular")
  }

  if (ncol(mod) != 3) {
    stop("Module label matrix should be a nx3 data frame")
  }

  # Get modules
  modules <- mod$moduleNumber + 1
  names(modules) <- mod$Gene.ID

  # Convert to igraph graph object
  g <- igraph::graph_from_adjacency_matrix(adj,
                                           mode = "upper",
                                           weighted = TRUE,
                                           diag = FALSE)
  Q <- igraph::modularity(g, membership = modules)

  return(Q)
}
