#' Find Module Quality
#'
#' This function finds quality metrics of all modules
#'
#' TODO more description
#'
#' @inheritParams compute.ModularityDensity
#'
#' @return metrics = data frame of module quality metrics, including the number
#'   of internal and external edges in each module, and the internal density,
#'   contraction, expansion, conductance, fitness, and average degree of each
#'   module.
#'
#' @export
compute.ModuleQualityMetric <- function(g, mod) {
  edge.comm <- community_edges(g, mod)
  internal <- as.numeric(diag(edge.comm))
  external <- as.numeric(rowSums(edge.comm, na.rm = TRUE)) - internal

  # Get size of each module
  N <- as.numeric(table(mod))

  metrics <- data.frame(moduleNumber = rownames(edge.comm),
                        # Internal edges in each community
                        internal.edges = internal,
                        # Boundary Edges
                        boundary.edges = external,
                        # Internal density of each community
                        internal.density = 2 * internal / (N * (N - 1)),
                        # Contraction
                        contraction = 2 * internal / N,
                        # Expansion
                        expansion = external / N,
                        # Conductance
                        conductance = external / (2 * internal + external),
                        # Fitness function
                        fitness = internal / (internal + external),
                        # Average modularity degree
                        avg.degree = (2 * internal - external) / N)

  return(metrics)
}
