#' Find Modules using igraph
#'
#' This function finds modules with one of igraph's cluster functions, as
#' specified in the arguments.
#'
#' TODO set random seed
#'
#' @param adj An adjacency matrix which should be n x n, where n is the number
#'   of genes. The matrix should be square and symmetric across the diagonal.
#' @param method Which method to use to find modules. Current options supported
#'   by this wrapper are: fast_greedy, infomap, label_prop, leading_eigen,
#'   link_communities, louvain, spinglass, walktrap
#' @param min.module.size Optional. Integer between 1 and n genes.
#' @param ... Optional. Additional arguments to pass through to the cluster
#'   algorithm.
#'
#' @return GeneModules = n x 3 data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
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
findModules.igraphWrapper <- function(adj, method, min.module.size = 30, ...) {
  # Convert network matrix to igraph graph object
  g <- igraph::graph_from_adjacency_matrix(adj,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)

  mod <- switch(method,
    fast_greedy = igraph::cluster_fast_greedy(g, ...),
    infomap = igraph::cluster_infomap(g, ...),
    label_prop = igraph::cluster_label_prop(g, ...),
    leading_eigen = igraph::cluster_leading_eigen(g, ...),
    link_communities = linkcommunities_wrapper(g, ...),
    louvain = igraph::cluster_louvain(g, ...),
    spinglass = spinglass_wrapper(g, min.module.size, ...),
    walktrap = igraph::cluster_walktrap(g, ...),
    # Default: unrecognized method
    NULL
  )

  if (is.null(mod)) {
    # TODO
  }

  # The spinglass and linkcommunities wrappers return a vector, but other
  # algorithms return a "communities" object. The membership vector needs to be
  # extracted.
  if (inherits(mod, "communities")) {
    mod <- igraph::membership(mod)
  }

  geneModules <- data.frame(Gene.ID = names(mod),
                            moduleNumber = as.numeric(mod))

  # Reassign modules smaller than min.module.size to module 0
  mod_sizes <- table(geneModules$moduleNumber)
  excluded <- names(mod_sizes)[mod_sizes < min.module.size]

  geneModules$moduleNumber[geneModules$moduleNumber %in% excluded] <- 0

  # Add back any genes that are missing from the geneModules data frame
  missing <- setdiff(colnames(adj), geneModules$Gene.ID)

  if (length(missing) > 0) {
    geneModules <- rbind(geneModules,
                         data.frame(Gene.ID = missing, moduleNumber = 0))
  }

  # Re-set module numbers to be consecutive after filtering. The -1 is to make
  # it start at 0 instead of 1.
  geneModules$moduleNumber <- as.numeric(factor(geneModules$moduleNumber)) - 1

  # Change cluster number to color labels
  geneModules$moduleLabel <- WGCNA::labels2colors(geneModules$moduleNumber)

  return(geneModules)
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
#' @param ... Optional. Other arguments accepted by
#'   \code{linkcomm::linkcommunities}.
#'
#' @returns a named vector where the names are gene names and the values are
#'   module membership.
linkcommunities_wrapper <- function(g, ...) {
  elist <- igraph::as_edgelist(g)

  comm <- linkcomm::getLinkCommunities(elist, ...)

  # Extract cluster information for each node / gene
  nodes <- comm$nodeclusters
  nodes$clusterSize <- comm$clustsizes[nodes$cluster]

  # Genes can be in multiple clusters, so we assign each gene to the largest
  # cluster it belongs to. The nodeclusters data frame is in order by ascending
  # cluster number, so if there are ties in size, the cluster with the smaller
  # ID will always be chosen.
  node_vec <- sapply(unique(nodes$node), function(gene) {
    clusts <- subset(nodes, nodes$node == gene)
    largest <- clusts$cluster[which.max(clusts$clusterSize)]
    return(as.numeric(largest))
  })

  return(node_vec)
}
