#' Find Modules using igraph
#'
#' This function finds modules with one of igraph's cluster functions, as
#' specified in the arguments.
#'
#' @param adj An adjacency matrix which should be n x n, where n is the number
#'   of genes
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
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @references Fast greedy method: http://arxiv.org/abs/cond-mat/0408187
#' @references Infomap algorithm: http://arxiv.org/abs/physics/0512106
#' @references Label propagation algorithm: Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm to detect community structures in large-scale networks. Phys Rev E 76, 036106. (2007)
#' @references Louvain method: Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre: Fast unfolding of communities in large networks. J. Stat. Mech. (2008) P10008
#' @references Spinglass algorithm: http://arxiv.org/abs/cond-mat/0603718
#' @references Walktrap algorithm: http://arxiv.org/abs/physics/0512106
#'
#' @export
findModules.igraphWrapper <- function(adj, method, min.module.size, ...) {
  # Convert network matrix to igraph graph object
  g <- igraph::graph_from_adjacency_matrix(adj,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)

  mod <- switch(method,
    fast_greedy = igraph::cluster_fast_greedy(g, ...),
    infomap = igraph::cluster_infomap(g, ...),
    label_prop = igraph::cluster_label_prop(g, ...),
    # Leading eigenvectors TODO
    leading_eigen = c(), # TODO
    # Link communities TODO
    link_communities = c(), # TODO
    louvain = igraph::cluster_louvain(g, ...),
    spinglass = c(), # TODO
    walktrap = igraph::cluster_walktrap(g, ...),
    # Default: unrecognized method
    NULL
  )

  if (is.null(mod)) {
    # TODO
  }

  # Get individual clusters from the igraph community object
  geneModules <- igraph::membership(mod) %>%
    unclass() %>%
    as.data.frame() %>%
    dplyr::rename(moduleNumber = ".")

  geneModules$Gene.ID <- rownames(geneModules)

  # Rename modules with size less than min module size to 0
  filteredModules <- geneModules %>%
    dplyr::group_by(.data$moduleNumber) %>%
    dplyr::summarise(counts = dplyr::n()) %>%
    dplyr::filter(.data$counts >= min.module.size)

  excluded <- !(geneModules$moduleNumber %in% filteredModules$moduleNumber)
  geneModules$moduleNumber[excluded] <- 0

  # Change cluster number to color labels
  geneModules$moduleLabel <- WGCNA::labels2colors(geneModules$moduleNumber)

  return(geneModules[, c("Gene.ID", "moduleNumber", "moduleLabel")])
}
