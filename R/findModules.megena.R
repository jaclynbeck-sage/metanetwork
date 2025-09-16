#' Find Modules with Megena Clustering
#'
#' This function finds modules from a network adjacency matrix using MEGENA.
#'
#' @param g An `igraph` graph
#' @param n_cores Optional. The number of cores to use for parallel computing.
#'   If n_cores = 1, MEGENA will not compute in parallel.
#' @param alpha.cut Optional. Resolution cut-off for cutting the cluster
#'   dendrogram.
#' @param ... Optional. Additional arguments for [MEGENA::do.MEGENA] or
#'   [MEGENA::get.union.cut].
#'
#' @return A named vector where the names are genes and the values are cluster
#'   assignments
#'
#' @export
findModules.megena <- function(g,
                               n_cores = 1,
                               alpha.cut = 1,
                               ...) {
  result <- MEGENA::do.MEGENA(g, doPar = n_cores > 1, num.cores = n_cores, ...)

  # Get clusters such that nodes exist in only one cluster
  clusters <- MEGENA::get.union.cut(result$module.output,
                                    alpha.cut = alpha.cut,
                                    output.plot = FALSE,
                                    ...)

  membership <- lapply(names(clusters), function(clust_name) {
    mem <- rep(clust_name, length(clusters[[clust_name]]))
    names(mem) <- clusters[[clust_name]]
    mem
  })

  return(unlist(membership))
}
