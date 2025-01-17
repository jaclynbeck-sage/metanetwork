#' Find Modules Fast and Greedy
#'
#' This function finds modules with igraph::cluster_fast_greedy()
#'
#' @inheritParams findModules.CFinder
#'
#' @return GeneModules = n x 3 data frame with column names as Gene.ID,
#' moduleNumber, and moduleLabel.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
findModules.fast_greedy.once <- function(adj, min.module.size) {
  # Convert lsparseNetwork to igraph graph object
  g <- igraph::graph_from_adjacency_matrix(adj,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)

  # Get modules using fast greedy method (http://arxiv.org/abs/cond-mat/0408187)
  mod <- igraph::cluster_fast_greedy(g)

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
