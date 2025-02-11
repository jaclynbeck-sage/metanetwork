#' Finds Consensus Clusters
#'
#' A modified parallel version of code imported from
#' https://bioconductor.org/packages/ConsensusClusterPlus 1.70.0
#'
#' @param d An MxN matrix, where M is the number of clustering methods * the
#'   number of clusters in each method, N is the number of genes, and the values
#'   represent cluster membership: 1 if a gene is in that cluster, 0 if not.
#' @param maxK Optional. Maximum number of clusters to evaluate
#' @param n_subsamples Optional. Number of random sub-samples of the genes to use.
#' @param pGenes Optional. Proportion of genes to sample at each random
#'   sub-sample.
#' @param clusterAlg Optional. One of three cluster algorithms: "hclust",
#'   "kmeans", or "pam".
#' @param hclust_method Optional. The "method" argument of \code{hclust}, if
#'   \code{clusterAlg} is "hclust"
#' @param distance_metric Optional. The metric to use to calculate distances:
#'   "pearson", "spearman", or "euclidean".
#' @param changeCDFArea Optional. Minimum spline distance for seq(2,`maxK`,
#'   length.out = `n_ks`)
#' @param n_ks Optional. Rather than test every \code{k} in the range (2,
#'   \code{maxK}), only \code{n_ks} values in the range will be tested, chosen
#'   with \code{seq(2, maxK, length.out = n_ks)}.
#' @param corUse Optional. The "use" argument of the \code{cor} function. Only
#'   used when \code{clusterAlg} is "hclust" or "pam" and
#'   \code{distance_metric} is "pearson" or "spearman".
#' @param n_cores Optional. Number of parallel cores to use when computing
#'   clusters. If \code{n_cores} = 1, the function will not run in parallel.
#' @param log_file_path Optional. If running in parallel, output will be written
#'   to a log file in the directory specified by \code{log_file_path}. By
#'   default, logs are written to the current working directory.
#' @param verbose Optional. When set to TRUE, prints messages to the screen to
#'   indicate progress. When FALSE, no messages are printed.
#' @param seed Optional. The number to use to set the random seed, for
#'   reproducible results.
#'
#' @return  Final clustered modules. TODO
#'
#' TODO does this really need to be run on partition.adj? Wouldn't it make sense
#' to do the consensus matrix like in run.consensus.cluster using all the different
#' module results instead?
#' @export
findModules.consensusCluster <- function(d,
                                         maxK = 100,
                                         n_subsamples = 100,
                                         pGenes = 0.8,
                                         clusterAlg = "hclust",
                                         hclust_method = "average",
                                         distance_metric = "pearson",
                                         changeCDFArea = 0.001,
                                         n_ks = 20,
                                         corUse = "everything",
                                         n_cores = 1,
                                         log_file_path = ".",
                                         verbose = FALSE,
                                         seed = 101) {
  # Set seed
  set.seed(seed)

  # Run consensus clustering
  kGrid <- seq(2, maxK, length.out = n_ks)
  kGrid <- unique(round(kGrid))

  main.dist.obj <- NULL

  if (verbose) {
    message("Calculating distance...")
  }

  if (clusterAlg != "kmeans") {
    if (distance_metric == "pearson" | distance_metric == "spearman") {
      main.dist.obj <- stats::as.dist(1 - stats::cor(d,
                                                     method = distance_metric,
                                                     use = corUse))

    } else if (inherits(try(get(distance_metric), silent = T), "function")) {
      main.dist.obj <- get(distance_metric)(t(d))

    } else {
      main.dist.obj <- stats::dist(t(d), method = distance_metric)
    }

  } else {
    # For kmeans, the values in "d" are used as distances TODO why?
    main.dist.obj <- stats::as.dist(d)
  }

  if (verbose) {
    message("Clustering with random sub-samples...")
  }

  # Run in parallel
  if (n_cores > 1) {
    log_file <- file.path(log_file_path, "findModules_log.txt")
    clust <- parallel::makeCluster(n_cores, outfile = log_file)

    clusters <- parallel::parLapply(cl = clust,
                                    X = 1:n_subsamples,
                                    fun = run.consensus.cluster,
                                    main.dist.obj = main.dist.obj,
                                    kGrid = kGrid,
                                    pGenes = pGenes,
                                    clusterAlg = clusterAlg,
                                    hclust_method = hclust_method,
                                    distance_metric = distance_metric,
                                    verbose = verbose,
                                    seed = seed)

    parallel::stopCluster(clust)

  } else {
    # Run sequentially
    clusters <- lapply(1:n_subsamples, run.consensus.cluster,
                       main.dist.obj = main.dist.obj,
                       kGrid = kGrid,
                       pGenes = pGenes,
                       clusterAlg = clusterAlg,
                       hclust_method = hclust_method,
                       distance_metric = distance_metric,
                       verbose = verbose,
                       seed = seed)
  }

  # Compute consensus fraction and area under the cdf curve
  areaUnderCDF <- data.frame(k = sapply(clusters[[1]], max),
                             area = 0)

  if (n_cores > 1) {
    log_file <- file.path(log_file_path, "findModules_log.txt")
    clust <- parallel::makeCluster(n_cores, outfile = log_file)

    area <- parallel::parSapply(cl = clust,
                                X = 1:nrow(areaUnderCDF),
                                FUN = function(row_ind, clusters, genes) {
      message(paste("Calculating CDF area for k =", names(clusters[[1]])[row_ind]))
      res <- calculate_areaUnderCDF(lapply(clusters, "[[", row_ind),
                                    genes,
                                    return_consensus = FALSE)$area
      return(c(row_ind, res))
    }, clusters, colnames(d))

    area <- as.data.frame(t(area))
    areaUnderCDF$area[area$V1] <- area$V2

    parallel::stopCluster(clust)

  } else {
    areaUnderCDF$area <- sapply(1:nrow(areaUnderCDF), function(row_ind) {
      message(paste("Calculating CDF area for k =", areaUnderCDF$k[row_ind]))
      calculate_areaUnderCDF(lapply(clusters, "[[", row_ind),
                             colnames(d),
                             return_consensus = FALSE)$area
    })
  }

  # Check if ends are maximum
  fn <- stats::splinefun(areaUnderCDF$k, areaUnderCDF$area)
  ind <- which(diff(fn(2:maxK)) >= changeCDFArea)
  # + 2 because index into diff is 1 less than the index into 2:maxK that
  # produced that item, and add another 1 because (2:maxK)[ind + 1] is ind + 2
  k.final <- ind[length(ind)] + 2

  # k.final might not be one of the k's tested, so we re-run consensus clustering
  # with the final value
  results.final <- run.consensus.cluster(
    rep_ind = 0,
    main.dist.obj = main.dist.obj,
    kGrid = k.final,
    pGenes = pGenes,
    clusterAlg = clusterAlg,
    hclust_method = hclust_method,
    distance_metric = distance_metric,
    verbose = verbose,
    seed = seed + n_subsamples + 1
  )

  final_area <- calculate_areaUnderCDF(lapply(results.final, "[[", 1),
                                       colnames(d),
                                       return_consensus = TRUE)

  # TODO Which one to use? Spot-checking one time, the clustering from hclust is
  # really bad for k = 10 while it looks ok for kmeans
  hc <- stats::hclust(stats::as.dist(1 - final_area$consensus.matrix),
                      method = "average")
  ct <- stats::cutree(hc, k.final)

  # Compute the final clusters
  cluster.final = stats::kmeans(final_area$consensus.matrix,
                                k.final)$cluster

  names(cluster.final) <- colnames(d)

  # TODO module size cutoff
  return(cluster.final)
}


# TODO minimum module size should matter I think? Genes in modules < min module size
# should be put in cluster 0 before calculating area. Otherwise using clusters with
# 1-2 genes in them probably skews the data.
calculate_areaUnderCDF_test <- function(clusters_k, genes, clust) {
  # Compute consensus matrix
  consensus <- parallel::parSapply(cl = clust,
                                   X = genes,
                                   FUN = calculate_consensusMatrix, genes, clusters_k)

  # Empirical CDF distribution -- create a cumulative distribution function
  # from the matrix, and sum the area under the function curve. Area is
  # estimated from a set of rectangles defined by: height is the CDF at a
  # specific value of consensus, and width is the difference between that value
  # in consensus and the next value in sorted order.
  cdf_fn <- stats::ecdf(as.numeric(consensus))
  vals <- sort(unique(as.numeric(consensus)))
  height <- cdf_fn(vals)
  width <- diff(vals)

  # The last value in "height" is discarded: the last value in consensus has
  # no point to the right of it on the graph so there's no width for that point.
  area <- sum(height[1:(length(height) - 1)] * width)

  return(list(area = area, consensus.matrix = consensus))
}


calculate_consensusMatrix <- function(gene, genes, clusters_k) {
  consensus_g <- rep(0, length(genes))
  names(consensus_g) <- genes

  counts_g <- consensus_g

  for (nrep in 1:length(clusters_k)) {
    res <- clusters_k[[nrep]]
    if (gene %in% names(res)) {
      counts_g[names(res)] <- counts_g[names(res)] + 1
      neighbors <- names(res)[res == res[gene]]
      consensus_g[neighbors] <- consensus_g[neighbors] + 1
    }
  }

  consensus_g <- consensus_g / counts_g
  consensus_g[is.na(consensus_g)] <- 0

  return(consensus_g)
}


calculate_areaUnderCDF <- function(clusters_k, genes, return_consensus = TRUE) {
  # Compute consensus matrix
  consensus <- matrix(0, length(genes), length(genes),
                      dimnames = list(genes, genes))

  counts <- matrix(0, length(genes), length(genes),
                   dimnames = list(genes, genes))

  for (nrep in 1:length(clusters_k)) {
    res <- clusters_k[[nrep]]
    clust_sizes <- table(res)
    excluded <- names(clust_sizes)[clust_sizes < 30]
    res[res %in% excluded] <- 0

    sampled <- names(res)
    counts[sampled, sampled] <- counts[sampled, sampled] + 1

    for (clust in unique(res)) {
      clust_genes <- names(res)[res == clust]
      consensus[clust_genes, clust_genes] <- consensus[clust_genes, clust_genes] + 1
    }
  }
  consensus <- consensus / counts
  consensus[is.na(consensus)] <- 0

  rm(counts)

  # Empirical CDF distribution -- create a cumulative distribution function
  # from the matrix, and sum the area under the function curve. Area is
  # estimated from a set of rectangles defined by: height is the CDF at a
  # specific value of consensus, and width is the difference between that value
  # in consensus and the next value in sorted order.
  cdf_fn <- stats::ecdf(as.numeric(consensus))
  vals <- sort(unique(as.numeric(consensus)))
  height <- cdf_fn(vals)
  width <- diff(vals)

  # The last value in "height" is discarded: the last value in consensus has
  # no point to the right of it on the graph so there's no width for that point.
  area <- sum(height[1:(length(height) - 1)] * width)

  if (return_consensus) {
    return(list(area = area, consensus.matrix = consensus))
  } else {
    return(list(area = area, consensus.matrix = NULL))
  }
}
