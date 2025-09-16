#' Finds Modules With Kmeans Clustering
#'
#' Function to get consensus modules from individual partition matrices
#'
#' @param partition.adj An MxN matrix, where M is the number of clustering
#'   methods * the number of clusters in each method, N is the number of genes,
#'   and the values represent cluster membership: 1 if a gene is in that
#'   cluster, 0 if not.
#' @param min.module.size Optional. An integer between 1 and n genes
#' @param usepam Optional. A logical for input into pam based kmeans clustering
#'   to find the number of clusters with the function `fpc::pamk`. If TRUE, pam
#'   is used, otherwise clara (recommended for large datasets with 2,000 or more
#'   observations; dissimilarity matrices can not be used with clara).
#'
#' @return A dataframe of Gene Modules TODO
#'
#' @export
findModules.consensusKmeans <- function(partition.adj,
                                        #maxK = 100, # TODO
                                        min.module.size = 30,
                                        usepam = TRUE) {
  # Error functions
  if (!inherits(partition.adj, "matrix")) {
    partition.adj <- data.matrix(partition.adj)
  }

  # Use pam based kmeans clustering to find the number of clusters
  mod <- fpc::pamk(t(partition.adj), krange = 2:30, usepam = usepam,
                   # both pam and clara
                   keep.data = FALSE,
                   # pam only
                   keep.diss = FALSE,
                   variant = "faster",
                   # clara only
                   samples = 50, pamLike = TRUE)

  mod_final <- mod$pamobject$clustering

  # Reassign modules smaller than min.module.size to module 0
  mod_sizes <- table(mod_final)
  excluded <- names(mod_sizes)[mod_sizes < min.module.size]

  mod_final[mod_final %in% excluded] <- 0

  return(mod_final)
}
