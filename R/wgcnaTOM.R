#' Runs WGCNA TOMsimilarity
#'
#' Runs WGCNA::TOMsimilarity over a gene expression matrix.
#'
#' @inheritParams wgcnaSoftThreshold
#'
#' @return NULL. Writes coexpression network named wgcnaTopologicalOverlapMatrixNetwork.csv
#' to `outputpath`
#'
#' @export
wgcnaTOM <- function(data, outputpath, RsquaredCut = .80, defaultNaPower = 6) {
  res <- WGCNA::pickSoftThreshold(data, RsquaredCut = RsquaredCut, networkType = "signed")
  if (is.na(res$powerEstimate)) {
    res$powerEstimate <- defaultNaPower
  }

  # TODO JB what is this for? Why not use WGCNA's built-in function WGCNA::adjacency?
  network <- abs(stats::cor(t(data), use = "pairwise.complete.obs"))^res$powerEstimate
  # TODO JB why is this multiplied by upper.tri?
  utils::write.csv(network * upper.tri(network),
                   file = file.path(outputpath, "wgcnaSoftThresholdNetwork.csv"))
  gc()
  cat(res$powerEstimate, "\n", sep = "",
      file = file.path(outputpath, "wgcnaPowerEstimate.txt"))
  print(res$powerEstimate)

  cn <- colnames(network)
  network <- WGCNA::TOMsimilarity(network)
  colnames(network) <- cn
  rownames(network) <- cn
  # save(network,file=paste0(outputpath,'result_wgcnaTOM.rda'))

  network <- network * upper.tri(network)
  utils::write.csv(
    network,
    file = file.path(outputpath, "wgcnaTopologicalOverlapMatrixNetwork.csv"),
    quote = F
  )
}
