#' Runs Consensus Clustering Algorithm
#'
#' A modified parallel version of code imported from
#' https://bioconductor.org/packages/ConsensusClusterPlus 1.70.0. Taken out
#' of findModules.consensusCluster.R
#'
#' @inheritParams findModules.consensusCluster
#'
#' @param rep_ind Optional. The index of this random sub-sample, used as an ID
#'   and in setting the random seed.
#' @param main.dist.obj A \code{dist} object containing distances between each
#'   feature/gene in the original data. The object must have gene names as
#'   labels.
#' @param kGrid Optional. A vector of k values to use for clustering.
#' @param seed Optional. The random seed to use. This value will be added to
#'   \code{rep_ind} to create a unique seed value for each random sample, for
#'   reproducibility.
#'
#' @return A named list where names are "k_" plus each \code{k} in \code{kGrid},
#'   and each item is a vector of cluster assignments generated from that
#'   \code{k} value. The vector's names are genes and the values are cluster
#'   numbers.
#'
#' @export
run.consensus.cluster <- function(rep_ind = 0,
                                  main.dist.obj,
                                  kGrid = 2:5,
                                  pGenes = 0.8,
                                  clusterAlg = "hclust",
                                  hclust_method = "average",
                                  distance_metric = "pearson",
                                  verbose = FALSE,
                                  seed = 101) {
  set.seed(seed + rep_ind)

  # Sample genes
  gene_names <- labels(main.dist.obj)
  sampleN <- floor(length(gene_names) * pGenes)
  sample_x <- sort(sample(gene_names, sampleN, replace = FALSE))

  # Sub-sample the distance object
  this_dist <- as.matrix(main.dist.obj)[sample_x, sample_x]
  this_dist <- stats::as.dist(this_dist)

  # Cluster samples using hclust (hier. clustering)
  this_cluster <- NA
  if (clusterAlg == "hclust") {
    this_cluster <- stats::hclust(this_dist, method = hclust_method)
  }

  # Run through all values of k and identify members
  cls_k <- lapply(kGrid, function(k, verbose, clusterAlg, this_dist,
                                  this_cluster, sample_x) {
    if (verbose) {
      message(paste("subsample", rep_ind, "\tk =", k))
    }

    this_assignment <- NA
    if (clusterAlg == "hclust") {
      # Prune to k for hclust
      this_assignment <- stats::cutree(this_cluster, k)

    } else if (clusterAlg == "kmeans") {
      this_assignment <- stats::kmeans(this_dist, k, iter.max = 100)$cluster

    } else if (clusterAlg == "pam") {
      # TODO fpc::pamk(this_dist, krange=k) is equivalent, could get rid of
      # one of fpc or cluster library dependency. cluster::pam is faster.
      this_assignment <- cluster::pam(this_dist, k, metric = distance_metric,
                                      cluster.only = TRUE)
    } else {
      # Optional clusterArg Hook.
      this_assignment <- get(clusterAlg)(this_dist, k)
    }

    names(this_assignment) <- sample_x
    return(this_assignment)
  },
  verbose, clusterAlg, this_dist, this_cluster, sample_x)

  names(cls_k) <- paste0("k_", kGrid)

  return(cls_k)
}
