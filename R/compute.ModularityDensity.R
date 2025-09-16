#' Find Modularity Density (Qds)
#'
#' This function calculates modularity density according to Chen, Nguyen, and
#' Szymanski 2013.
#'
#' @param g An `igraph` graph of the network
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
  results <- community_edges(g, mod)
  edge.comm <- results$edge.comm

  # Number of genes in each module ("|c_i|" in the paper)
  c_i <- table(results$renamed_mod)

  stopifnot(all(colnames(edge.comm) == names(c_i)))

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
#' @param g An `igraph` graph
#' @param mod A vector where names are genes and values are module membership
#'
#' @returns A named list with names "edge.comm" and "renamed_mod". "edge.comm" =
#'   an NxN matrix, where N is the number of modules/communities. The diagonal
#'   will contain the number of internal edges in each module, and the
#'   non-diagonals will contain the number of edges between module i and module
#'   j. "renamed_mod" = `mod`, if `mod` is already a character vector, otherwise
#'   all values of `mod` have an "m" pasted in front to force it to be a
#'   character vector. This avoids ambiguity between module numbers and indexing
#'   into arrays.
#' @export
community_edges <- function(g, mod) {
  # Module names need to be characters to avoid ambiguity between module numbers
  # and indexing
  if (is.numeric(mod)) {
    genes <- names(mod)
    mod <- paste0("m", mod)
    names(mod) <- genes
  }

  igraph::V(g)$moduleNumber <- mod[igraph::V(g)$name]

  all_modules <- sort(unique(mod))

  edge.comm <- matrix(0, nrow = length(all_modules), ncol = length(all_modules),
                      dimnames = list(all_modules, all_modules))

  # For each community ci, get the within-community edges and the edges between
  # ci and each other community.
  for (ind1 in 1:length(all_modules)) {
    ci <- all_modules[ind1]

    # All edges involving nodes in this module
    vi <- which(igraph::V(g)$moduleNumber == ci)
    edges1 <- igraph::incident_edges(g, v = vi)
    edges1 <- unique(do.call(c, edges1))

    for (ind2 in ind1:length(all_modules)) {
      cj <- all_modules[ind2]

      # Internal edges
      if (ci == cj) {
        sg <- igraph::induced_subgraph(g, vids = vi)
        edge.comm[ci, cj] <- igraph::ecount(sg)

      } else {
        # Edges between ci and cj
        vj <- which(igraph::V(g)$moduleNumber == cj)
        edges2 <- igraph::incident_edges(g, v = vj)
        edges2 <- unique(do.call(c, edges2))

        edge_len <- length(igraph::intersection(edges1, edges2))
        edge.comm[ci, cj] <- edge_len
        edge.comm[cj, ci] <- edge_len
      }
    }
  }

  return(list(edge.comm = edge.comm, renamed_mod = mod))
}
