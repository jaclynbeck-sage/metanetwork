#' Find Modules TODO
#' TODO WGCNA and any other algorithms used in construction that also do modules
#'
#' @param adj An n x n upper triangular adjacency matrix where "n" is the number
#'   of genes.
#' @param method Which method to use to find modules. Current options are:
#'   fast_greedy, infomap, label_prop, leading_eigen, link_communities, louvain,
#'   megena, spinglass, walktrap
#' @param nperm Optional. Number of permutations on the gene ordering.
#' @param min.module.size Optional. Integer between 1 and n genes.
#' @param n_cores Optional. For "megena" only, how many cores to use when
#'   computing in parallel. If n_cores = 1, megena will not compute in parallel.
#' @param ... Optional. Additional arguments to be passed to the individual
#'   clustering functions
#'
#' @return An n x 2 data frame with columns for "gene" and "module", where n is
#'   the number of genes
#'
#' @export
findModules <- function(adj, method, nperm = 10, min.module.size = 30, n_cores = 1, ...) {
  if (!inherits(adj, "matrix")) {
    adj <- data.matrix(adj)
  }

  if (nrow(adj) != ncol(adj)) {
    stop("Adjacency matrix should be square")
  }

  if (!all(adj[lower.tri(adj)] == 0)) {
    stop("Adjacency matrix should be upper triangular")
  }

  # Make the upper-triangular adjacency matrix symmetric
  adj <- adj + t(adj)

  set.seed(nperm) # TODO better seed

  # Compute modules by permuting the labels nperm times
  # TODO parallel?
  all.modules <- lapply(1:nperm, function(i, adj, min.module.size) {
    # Permute gene ordering
    ind <- sample(1:nrow(adj), nrow(adj), replace = FALSE)
    adj1 <- adj[ind, ind]

    # Convert to an igraph object
    g <- igraph::graph_from_adjacency_matrix(adj1,
                                             mode = "undirected",
                                             weighted = TRUE,
                                             diag = FALSE)

    # Find modules TODO switch statement
    mod <- findModules.igraphWrapper(g, method, min.module.size, n_cores, ...)

    # Compute modularity and modularity density

    # Mod needs to start at 1, not 0 for this function
    Q <- igraph::modularity(g, membership = mod + 1)

    Qds <- compute.ModularityDensity(g, mod)

    return(list(mod = mod, Q = Q, Qds = Qds))
  }, adj, min.module.size)

  # Find the best module based on highest Q and Qds
  Q <- sapply(all.modules, "[[", "Q")
  Qds <- sapply(all.modules, "[[", "Qds")

  rank <- rank(Q) + rank(Qds)
  ind <- which.max(rank)

  mod <- all.modules[[ind]]$mod

  mod_df <- data.frame(gene = names(mod), module = mod)
  return(mod_df)
}
