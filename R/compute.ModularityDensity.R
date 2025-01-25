#' Find Modularity Density (Qds)
#'
#' This function calculates modularity density according to Botta and Genio 2016.
#'
#' @param g An \code{igraph} graph of the network
#' @param mod A named vector where names are genes and values are module membership
#'
#' @return Qds = module density.
#'
#' @references Federico Botta, Charo I. del Genio. Finding network communities
#' using modularity density. https://arxiv.org/abs/1612.07297
#'
#' @export
compute.ModularityDensity <- function(g, mod) {
  # Module names need to be characters, not interpretable as numbers
  igraph::V(g)$moduleNumber <- paste0("mod.", mod[igraph::V(g)$name])

  # Get number of edges within and between communities. The diagonal will
  # contain the number of internal edges in each module, and the non-diagonals
  # will contain the number of edges between module i and module j.
  edge.comm <- sapply(unique(igraph::V(g)$moduleNumber), function(ci) {
    vi <- which(igraph::V(g)$moduleNumber == ci)
    edges1 <- igraph::incident_edges(g, v = vi)
    edges1 <- do.call(c, edges1)

    lengths <- sapply(unique(igraph::V(g)$moduleNumber), function(cj) {
      # Internal edges
      if (ci == cj) {
        sg <- igraph::induced_subgraph(g, vids = vi)
        return(igraph::ecount(sg))

      } else {
        # Edges between ci and cj
        vj <- which(igraph::V(g)$moduleNumber == cj)
        edges2 <- igraph::incident_edges(g, v = vj)
        edges2 <- do.call(c, edges2)

        return(length(igraph::intersection(edges1, edges2)))
      }
    })

    return(lengths)
  })

  # Number of genes in each module ("Nc" in the paper)
  N <- table(igraph::V(g)$moduleNumber)

  # Total number of edges -- "m" in paper
  m <- igraph::ecount(g)

  Qds <- sapply(rownames(edge.comm), function(ci) {
    mc <- edge.comm[ci, ci]  # "mc" in paper
    ec <- sum(edge.comm[ci, ], na.rm = TRUE) - mc  # "ec" in paper

    # pc in the paper
    pc <- 2 * mc / (N[ci] * (N[ci] - 1))

    between <- 0

    # sum of [(edges between ci and cj)^2 / (2m * Ni * Nj)] in paper
    for (cj in setdiff(rownames(edge.comm), ci)) {
      between <- between + edge.comm[ci, cj]^2 / (2 * m * N[ci] * N[cj])
    }

    ((mc / m) * pc) - ((2 * m + ec) * pc / (2 * m))^2 - between
  })

  Qds <- sum(Qds, na.rm = TRUE)

  return(Qds)
}
