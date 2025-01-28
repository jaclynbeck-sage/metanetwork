#' Runs Consensus Clustering Algorithm
#'
#' A modified parallel version of code imported from
#' https://bioconductor.org/packages/ConsensusClusterPlus 1.70.0. Taken out
#' of findModules.consensusCluster.R
#'
#' @param d A matrix where rows are samples and columns are genes
#' @param kGrid Optional. A vector of k values to use for clustering.
#' @param repCount Optional. Replicate count for random sampling.
#' @param pGenes Optional. Proportion of genes to sample on each rep
#' @param clusterAlg Optional. One of three cluster algorithms: "hclust",
#' "kmeans", or "pam".
#' @param hclust_method Optional. The "method" argument of \code{hclust}
#' @param distance Optional. Sample distance metrics: "pearson", "spearman", or "euclidean".
#' @param verbose Optional. When TRUE, prints messages to the screen to indicate progress.
#' @param corUse Optional. The "use" argument of the \code{cor} function
#'
#' @return A named list containing "areaUnderCDF" = a numerical vector of areas,
#' and "consensus.matrix" = a list of consensus matrices. Both the vector and
#' the inner list are named with the "k" value that was used to create them.
#'
#' @export
run.consensus.cluster <- function(d,
                                  kGrid = 2:5,
                                  repCount = 10,
                                  pGenes = 0.8,
                                  hclust_method = "average",
                                  distance = "pearson",
                                  clusterAlg = "hclust",
                                  verbose = FALSE,
                                  corUse = "everything") {
  main.dist.obj <- NULL

  if (clusterAlg != "kmeans") {
    if (distance == "pearson" | distance == "spearman") {
      main.dist.obj <- stats::as.dist(1 - stats::cor(d, method = distance,
                                                     use = corUse))

    } else if (inherits(try(get(distance), silent = T), "function")) {
      main.dist.obj <- get(distance)(t(d))

    } else {
      main.dist.obj <- stats::dist(t(d), method = distance)
    }
  }

  cls <- lapply(1:repCount, function(i, verbose, d, kGrid, pGenes, main.dist.obj,
                                     clusterAlg, distance, corUse, hclust_method) {
    if (verbose) {
      message(paste("random subsample", i))
    }

    # Sample genes
    sampleN <- floor(ncol(d) * pGenes)
    sample_x <- colnames(d)[sort(sample(ncol(d), sampleN, replace = FALSE))]

    # Compute distance (if not supplied)
    this_dist <- NA

    if (!is.null(main.dist.obj)) {
      this_dist <- as.matrix(main.dist.obj)[sample_x, sample_x]
      this_dist <- stats::as.dist(this_dist)

    } else {
      this_dist <- d[, sample_x]
    }

    # Cluster samples using hclust (hier. clustering)
    this_cluster <- NA
    if (clusterAlg == "hclust") {
      this_cluster <- stats::hclust(this_dist, method = hclust_method)
    }

    # Run through all values of k and identify members
    cls <- lapply(kGrid, function(k, verbose, clusterAlg, this_dist,
                                  this_cluster, sample_x) {
      if (verbose) {
        message(paste("  k =", k))
      }

      this_assignment <- NA
      if (clusterAlg == "hclust") {
        # Prune to k for hclust
        this_assignment <- stats::cutree(this_cluster, k)

      } else if (clusterAlg == "kmeans") {
        this_assignment <- stats::kmeans(t(this_dist), k, iter.max = 100)$cluster

      } else if (clusterAlg == "pam") {
        this_assignment <- cluster::pam(this_dist, k, metric = distance,
                                        cluster.only = TRUE)
      } else {
        # Optional clusterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }

      names(this_assignment) <- sample_x
      return(this_assignment)
    },
    verbose, clusterAlg, this_dist, this_cluster, sample_x)

    names(cls) <- kGrid

    return(cls)
  },
  verbose, d, kGrid, pGenes, main.dist.obj, clusterAlg, distance, corUse, hclust_method)

  # TODO this should move to findModules.consensusCluster
  # Compute consensus fraction and area under the cdf curve
  areaK <- rep(0, length(cls[[1]]))
  names(areaK) <- names(cls[[1]])
  cns.list <- list()

  for (k in names(areaK)) {
    # Compute consensus matrix
    cns.mtrx <- matrix(0, ncol(d), ncol(d),
                       dimnames = list(colnames(d), colnames(d)))

    counts <- matrix(0, ncol(d), ncol(d),
                     dimnames = list(colnames(d), colnames(d)))

    for (nrep in 1:repCount) {
      res <- cls[[nrep]][[k]]
      sampled <- names(res)
      counts[sampled, sampled] <- counts[sampled, sampled] + 1

      for (clust in unique(res)) {
        clust_genes <- names(res)[res == clust]
        cns.mtrx[clust_genes, clust_genes] <- cns.mtrx[clust_genes, clust_genes] + 1
      }
    }
    cns.mtrx <- cns.mtrx / counts
    cns.mtrx[is.na(cns.mtrx)] <- 0

    cns.list[[k]] <- cns.mtrx

    # Empirical CDF distribution -- create a cumulative distribution function
    # from the matrix, and sum the area under the function curve. Area is
    # estimated from a set of rectangles defined by: height is the CDF at a
    # specific value of cns.mtrx, and width is the difference between that value
    # in cns.mtrx and the next value in sorted order.
    cdf_fn <- stats::ecdf(as.numeric(cns.mtrx))
    vals <- sort(unique(as.numeric(cns.mtrx)))
    height <- cdf_fn(vals)
    width <- diff(vals)

    # The last value in "height" is discarded: the last value in cns.mtrx has
    # no point to the right of it on the graph so there's no width for that point.
    areaK[k] <- sum(height[1:(length(height) - 1)] * width)
  }

  return(list(areaUnderCDF = areaK,
              consensus.matrix = cns.list,
              clusters = cls))
}
