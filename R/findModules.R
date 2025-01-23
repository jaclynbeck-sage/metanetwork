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
#'
#' @return GeneModules = n x 3 data frame with column names as Gene.ID,
#'   moduleNumber, and moduleLabel.
#'
#' @importFrom magrittr %>%
#' @export
findModules <- function(adj, method, nperm = 10, min.module.size = 30) {
  if (!inherits(adj, "matrix")) {
    stop("Adjacency matrix should be of class matrix")
  }

  if (nrow(adj) != ncol(adj)) {
    stop("Adjacency matrix should be symmetric")
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

    # Find modules TODO switch statement
    if (method == "megena") {
      # TODO
    } else {
      mod <- findModules.igraphWrapper(adj1, method, min.module.size)
    }

    # Compute local and global modularity
    adj1[lower.tri(adj1)] <- 0
    Q <- compute.Modularity(adj1, mod)
    Qds <- compute.ModularityDensity(adj1, mod)

    return(list(mod = mod, Q = Q, Qds = Qds))
  }, adj, min.module.size)

  # Find the best module based on Q and Qds
  tmp <- plyr::ldply(all.modules, function(x) {
    data.frame(Q = x$Q, Qds = x$Qds)
  }) %>%
    dplyr::mutate(r = base::rank(.data$Q) + base::rank(.data$Qds))
  ind <- which.max(tmp$r)

  mod <- all.modules[[ind]]$mod

  return(mod)
}
