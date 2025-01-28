#' Find Modules using igraph
#'
#' This function finds modules with one of igraph's cluster functions, as
#' specified in the arguments.
#'
#' TODO set random seed
#'
#' @param g An \code{igraph} graph
#' @param method Which method to use to find modules. Current options supported
#'   by this wrapper are: fast_greedy, infomap, label_prop, leading_eigen,
#'   linkcommunities, louvain, spinglass, walktrap
#' @param min.module.size Optional. How many genes need to be in a module for it
#'   to be considered valid.
#' @param n_cores Optional. For "megena" only, how many cores to use when
#'   computing in parallel. If n_cores = 1, megena will not compute in parallel.
#' @param ... Optional. Additional arguments to pass through to the cluster
#'   algorithm.
#'
#' @return a named vector where names are genes and values are cluster assignments
#'
#' @references Fast greedy method: http://arxiv.org/abs/cond-mat/0408187
#' @references Infomap algorithm: http://arxiv.org/abs/physics/0512106
#' @references Label propagation algorithm: Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm to detect community structures in large-scale networks. Phys Rev E 76, 036106. (2007)
#' @references Leading eigen algorithm: MEJ Newman: Finding community structure using the eigenvectors of matrices, Physical Review E 74 036104, 2006
#' @references Link communities: TODO
#' @references Louvain method: Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre: Fast unfolding of communities in large networks. J. Stat. Mech. (2008) P10008
#' @references Spinglass algorithm: http://arxiv.org/abs/cond-mat/0603718
#' @references Walktrap algorithm: http://arxiv.org/abs/physics/0512106
#'
#' @export
findModules.igraphWrapper <- function(g,
                                      method,
                                      min.module.size = 30,
                                      n_cores = 1,
                                      ...) {
  mod <- switch(method,
    fast_greedy = igraph::cluster_fast_greedy(g, ...),
    infomap = igraph::cluster_infomap(g, ...),
    label_prop = igraph::cluster_label_prop(g, ...),
    leading_eigen = igraph::cluster_leading_eigen(g, ...),
    linkcommunities = linkcommunities_wrapper(g, ...),
    louvain = igraph::cluster_louvain(g, ...),
    megena = findModules.megena(g, n_cores = n_cores, ...),
    spinglass = spinglass_wrapper(g, min.module.size, ...),
    walktrap = igraph::cluster_walktrap(g, ...),
    # Default: unrecognized method
    NULL
  )

  if (is.null(mod)) {
    message(paste("Unrecognized algorithm name", method))
    return(NULL)
  }

  # The spinglass, linkcommunities, and megena wrappers return a vector, but
  # other algorithms return a "communities" object. The membership vector needs
  # to be extracted.
  if (inherits(mod, "communities")) {
    mod <- igraph::membership(mod)
  }

  # Reassign modules smaller than min.module.size to module 0
  mod_sizes <- table(mod)
  excluded <- names(mod_sizes)[mod_sizes < min.module.size]

  mod[mod %in% excluded] <- 0

  # Add back any genes that are missing from vector
  gene_names <- names(igraph::V(g))
  missing <- setdiff(gene_names, names(mod))

  if (length(missing) > 0) {
    missing_vec <- rep(0, length(missing))
    names(missing_vec) <- missing

    mod <- c(mod, missing_vec)
  }

  mod <- mod[gene_names]

  # Re-set module numbers to be consecutive after filtering. The -1 is to make
  # it start at 0 instead of 1. Using as.numeric erases vector names, so we
  # put them back.
  mod <- as.numeric(factor(mod)) - 1
  names(mod) <- gene_names

  return(mod)
}


#' Spinglass wrapper
#'
#' This function wraps \code{igraph::cluster_spinglass} so it can be run on
#' all connected sub-graphs of the main graph.
#'
#' \code{igraph::cluster_spinglass} only works on connected graphs, and the main
#' graph may have components that are not connected to other components. This
#' function divides the main graph into sub-graphs that are connected, runs
#' \code{cluster_spinglass} on each sub-graph, and concatenates the membership
#' results into one vector.
#'
#' @param g An \code{igraph} graph
#' @param min.module.size Optional. The minimum number of genes per module,
#'   which is used here to avoid running spinglass on connected graphs with
#'   fewer than \code{min.module.size} genes in them.
#' @param ... Optional. Other arguments accepted by
#'   \code{igraph::cluster_spinglass}.
#'
#' @returns a named vector where the names are gene names and the values are
#'   module membership.
#' @export
spinglass_wrapper <- function(g, min.module.size = 30, ...) {
  # Decompose into connected components with at least min.module.size components
  sg_list <- igraph::decompose(g, min.vertices = min.module.size)

  mod_list <- lapply(sg_list, igraph::cluster_spinglass, ...)

  # All results start with module "1", so we can't directly combine them or
  # there will be genes from different sub-graphs assigned to module 1 even
  # though they aren't in the same module. We need to add the number of modules
  # that exist in the list before the current module, so that each module has a
  # unique number across all sub-graphs.
  n_mods <- c(0, cumsum(sapply(mod_list, length)))

  memberships <- lapply(1:length(mod_list), function(mod_ind) {
    igraph::membership(mod_list[[mod_ind]]) + n_mods[mod_ind]
  })

  # Concatenate into one named vector
  return(do.call(c, memberships))
}


#' Linkcommunities wrapper
#'
#' This function wraps \code{linkcomm::linkcommunities}.
#'
#' \code{linkcomm::linkcommunities} requires an edge list instead of a graph,
#' and returns a custom object, so this function creates the edge list and then
#' extracts the cluster membership from the returned object.
#'
#' @param g An \code{igraph} graph
#' @param min.module.size Optional. The minimum number of genes per module,
#'   which is used here to avoid assigning genes to clusters with fewer than
#'   \code{min.module.size} nodes in them.
#' @param ... Optional. Other arguments accepted by
#'   \code{linkcomm::linkcommunities}.
#'
#' @returns a named vector where the names are gene names and the values are
#'   module membership.
#' @export
linkcommunities_wrapper <- function(g, min.module.size = 30, ...) {
  elist <- igraph::as_edgelist(g)

  comm <- linkcomm::getLinkCommunities(elist, ...)

  # Genes can be in multiple clusters, so we assign each gene to the cluster
  # where it has the most edges. The edges data frame is in order by ascending
  # cluster number, so if there are ties in size, the cluster with the smaller
  # ID will always be chosen.

  # We only work with "node1", so we need to make sure that nodes that only
  # exist in "node2" get included. We do this by reversing the edges so that
  # edges are bi-directional.
  edges <- rbind(
    comm$edges,
    data.frame(
      node1 = comm$edges$node2,
      node2 = comm$edges$node1,
      cluster = comm$edges$cluster
    )
  )
  edges$clusterSize <- comm$clustsizes[edges$cluster]
  edges <- subset(edges, edges$clusterSize >= min.module.size)

  # Assign genes to the cluster where they have the most edges.
  node_vec <- sapply(unique(edges$node1), function(gene) {
    clusts <- subset(edges, edges$node1 == gene)
    n_edges <- table(clusts$cluster)
    mem <- names(n_edges)[which.max(n_edges)]
    return(as.numeric(mem))
  })

  return(node_vec)
}
