#' Find Modularity Density (Qds)
#'
#' This function calculates modularity density according to Chen, Nguyen, and
#' Szymanski 2013.
#'
#' @param g An \code{igraph} graph of the network
#' @param mod A named vector where names are genes and values are module membership
#'
#' @return Qds = modularity density.
#'
#' @references Mingming Chen, Tommy Nguyen, Boleslaw K. Szymanski, 2013. A New
#'   Metric for Quality of Network Community Structure.
#'   https://doi.org/10.48550/arXiv.1507.04308
#'
#' @export
compute.ModularityDensity <- function(g, mod) {
  # Get number of edges within and between communities
  edge.comm <- community_edges(g, mod)

  # Number of genes in each module ("|c_i|" in the paper)
  c_i <- table(mod)

  # Total number of edges -- "|E|" in paper
  n_edges <- igraph::ecount(g)

  Ein <- diag(edge.comm)
  Eout <- colSums(edge.comm) - Ein
  d_ci <- 2 * Ein / (c_i * (c_i - 1))

  # Terms one and two of the summation can be calculated as vectors
  term1 <- Ein * d_ci / n_edges
  term2 <- ((2 * Ein + Eout) * d_ci / (2 * n_edges))^2

  # The third term (sum of (|E_cicj| / (2|E|)) * d_cicj) in paper) needs to
  # loop through all clusters
  term3 <- sapply(rownames(edge.comm), function(i) {
    sapply(setdiff(rownames(edge.comm), i), function(j) {
      d_cicj <- edge.comm[i, j] / (c_i[i] * c_i[j])
      return(edge.comm[i, j] * d_cicj / (2 * n_edges))
    }) |>
      sum()
  })

  Qds <- sum(term1 - term2 - term3)

  return(Qds)
}


#' Community Edges
#'
#' Get the number of edges within and between communities.
#'
#' @param g An \code{igraph} graph
#' @param mod A vector where names are genes and values are module membership
#'
#' @returns An NxN matrix, where N is the number of modules/communities. The
#'   diagonal will contain the number of internal edges in each module, and the
#'   non-diagonals will contain the number of edges between module i and module
#'   j.
#' @export
community_edges <- function(g, mod) {
  # Module names need to be characters to show up as row/column names in edge.comm
  igraph::V(g)$moduleNumber <- as.character(mod[igraph::V(g)$name])

  all_modules <- sort(unique(igraph::V(g)$moduleNumber))

  # For each community ci, get the within-community edges and the edges between
  # ci and each other community.
  edge.comm <- sapply(all_modules, function(ci) {
    # All edges involving nodes in this module
    vi <- which(igraph::V(g)$moduleNumber == ci)
    edges1 <- igraph::incident_edges(g, v = vi)
    edges1 <- unique(do.call(c, edges1))

    lengths <- sapply(all_modules, function(cj) {
      # Internal edges
      if (ci == cj) {
        sg <- igraph::induced_subgraph(g, vids = vi)
        return(igraph::ecount(sg))

      } else {
        # Edges between ci and cj
        vj <- which(igraph::V(g)$moduleNumber == cj)
        edges2 <- igraph::incident_edges(g, v = vj)
        edges2 <- unique(do.call(c, edges2))

        return(length(igraph::intersection(edges1, edges2)))
      }
    })

    return(lengths)
  })

  edge.comm <- edge.comm[sort(rownames(edge.comm)), sort(colnames(edge.comm))]
  return(edge.comm)
}
