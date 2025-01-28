#' Finds Consensus Clusters
#'
#' A modified parallel version of code imported from
#' https://bioconductor.org/packages/ConsensusClusterPlus 1.70.0
#'
#' @param d Optional. A matrix where columns=items/samples and rows are features.
#' For example, a gene expression matrix of genes in rows and microarrays in columns.
#' OR ExpressionSet object. (Default = NULL)
#' @param maxK Optional. An  integer value. maximum cluster number to evaluate.
#' (Default = 100)
#' @param reps Optional. An integer value. number of subsamples.  (Default = 100)
#' @param pGenes Optional. A numerical value. proportion of items to sample.
#'  (Default = 0.8)
#' @param clusterAlg Optional. A character value. cluster algorithm. "hc"
#' heirarchical (hclust) or "km" for kmeans. (Default = "kmeans")
#' @param hclust_method Optional. The "method" argument of \code{hclust}
#' @param distance Optional. A character value. sample distance measures:
#' "pearson", "spearman", or "euclidean". (Default = "pearson")
#' @param seed Optional, A numerical value. Sets random seed for reproducible results.
#' @param verbose Optional. A boolean when set to TRUE, prints messages to the
#' screen to indicate progress. This is useful for large datasets.(Default = FALSE)
#' @param changeCDFArea Optional. Minimum spline distance for seq(2,`maxK`,
#' length.out = `nbreaks`) (Default = 0.001)
#' @param corUse Optional. Use all cores avaiable. (Default = "Everything") TODO this is wrong
#' @param nbreaks Optional. Number of breaks to use in
#' seq(2,`maxK`,length.out = `nbreaks`) this becomes the kGrid argument in
#' run.consensus.cluster. (Default = 20)
#'
#' @return  Final clustered modules.
#'
#' TODO does this really need to be run on partition.adj? Wouldn't it make sense
#' to do the consensus matrix like in run.consensus.cluster using all the different
#' module results instead?
#' @export
findModules.consensusCluster <- function(d = NULL,
                                         maxK = 100,
                                         reps = 100,
                                         pGenes = 0.8,
                                         clusterAlg = "kmeans",
                                         hclust_method = "average",
                                         distance = "pearson",
                                         changeCDFArea = 0.001,
                                         nbreaks = 20,
                                         seed = 101,
                                         corUse = "everything",
                                         verbose = F) {
  # Set seed
  set.seed(seed)

  # Run consensus clustering
  kGrid <- seq(2, maxK, length.out = nbreaks)
  kGrid <- unique(round(kGrid))

  results <- run.consensus.cluster(
    d = d,
    kGrid = kGrid,
    repCount = reps,
    pGenes = pGenes,
    hclust_method = hclust_method,
    clusterAlg = clusterAlg,
    distance = distance,
    verbose = verbose,
    corUse = corUse
  )

  # Check if ends are maximum
  fn <- stats::splinefun(names(results$areaUnderCDF), results$areaUnderCDF)
  ind <- which(diff(fn(2:maxK)) >= changeCDFArea)
  # + 2 because index into diff is 1 less than the index into 2:maxK that
  # produced that item, and add another 1 because (2:maxK)[ind + 1] is ind + 2
  k.final <- ind[length(ind)] + 2

  # k.final might not be one of the k's tested, so we re-run consensus clustering
  # with the final value
  results.final <- run.consensus.cluster(
    d = d,
    kGrid = k.final,
    repCount = reps,
    pGenes = pGenes,
    hclust_method = hclust_method,
    clusterAlg = clusterAlg,
    distance = distance,
    verbose = verbose,
    corUse = corUse
  )

  # TODO Which one to use?
  hc <- stats::hclust(stats::as.dist(1 - results.final$consensus.matrix[[1]]),
                      method = "average")
  ct <- stats::cutree(hc, k.final)

  # Compute the final clusters
  cluster.final = stats::kmeans(results.final$consensus.matrix[[1]],
                                k.final)$cluster

  names(cluster.final) <- colnames(d)

  # TODO module size cutoff
  return(cluster.final)
}
