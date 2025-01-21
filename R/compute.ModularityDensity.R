#' Find Global Modularity (Qds)
#'
#' This function calculates modularity density. TODO needs review.
#' TODO this might be combinable with compute.Modularity to avoid calling
#' some of the igraph functions twice.
#'
#' @inheritParams compute.LocalModularity
#'
#' @return Qds = module density.
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export compute.ModularityDensity
compute.ModularityDensity <- function(adj, mod) {
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

  # Convert lsparseNetwork upper adj matrix to graph
  g <- igraph::graph_from_adjacency_matrix(adj,
                                           mode = "upper",
                                           weighted = TRUE,
                                           diag = FALSE)
  rownames(mod) <- mod$Gene.ID
  igraph::V(g)$moduleNumber <- mod[igraph::V(g)$name, "moduleNumber"]

  # Get number of edges between communities

  # JB TODO theoretically modules don't overlap, so this code unnecessarily loops
  # through all module comparisons. I think it's just effectively counting the
  # size of each module?
  edge.comm <- foreach(ci = unique(igraph::V(g)$moduleNumber),
                       .packages = c("foreach"),
                       .combine = cbind) %dopar% {
    foreach(cj = unique(igraph::V(g)$moduleNumber), .combine = c) %dopar% {
      gi <- igraph::induced_subgraph(g,
                                     vids = which(igraph::V(g)$moduleNumber == ci))
      gj <- igraph::induced_subgraph(g,
                                     vids = which(igraph::V(g)$moduleNumber == cj))
      igraph::ecount(igraph::intersection(gi, gj))
    }
  }

  edge.comm <- data.frame(edge.comm)
  rownames(edge.comm) <- unique(igraph::V(g)$moduleNumber)
  colnames(edge.comm) <- unique(igraph::V(g)$moduleNumber)

  # Get size of each modules
  mod.sz <- table(igraph::V(g)$moduleNumber)

  # Calculate local modularity
  E <- sum(edge.comm, na.rm = TRUE)

  # JB TODO Ecc.dcc is always zero if the non-diagonal of edge.comm is all 0
  Qds <- foreach(ci = rownames(edge.comm), .combine = c) %dopar% {
    if (edge.comm[ci, ci] != 0) {
      Ein <- edge.comm[ci, ci]
      Eout <- sum(edge.comm[ci, ], na.rm = TRUE) - Ein

      dc <- 2 * Ein / (mod.sz[ci] * (mod.sz[ci] - 1))

      Ecc.dcc <- 0

      for (cj in rownames(edge.comm)) {
        if (ci != cj) {
          Ecc.dcc <- Ecc.dcc + edge.comm[ci, cj]^2 / (mod.sz[ci] * mod.sz[cj])
        }
      }
      Ecc.dcc <- Ecc.dcc / (2 * E)

      ((Ein / E) * dc) - ((2 * Ein + Eout) * dc / (2 * E))^2 - Ecc.dcc
    }
  }

  Qds <- sum(Qds, na.rm = TRUE)

  return(Qds)
}
