#' Parallel Network Wrapper
#'
#' Wrapper to run network detection in a parallel fashion
#'
#' @param data Required. Expression matrix to be used for network construction,
#' which should have samples as rows and genes as columns.
#' @param regressionFunction Required. The name of the regression function to use,
#' which should exactly match the name of a function in this package. Current
#' options are: lassoIC, lassoCV, ridgeIC, ridgeCV, sparrowZ
#' @param n_cores Optional. The number of parallel cores to use.
#' @param cluster_type Optional. Use "FORK" if running on a Unix system, and
#' "PSOCK" if running on Windows.
#' @param regulatorIndex Optional. A vector of numerical indexes into
#' `colnames(data)` for a subset of genes that should be used in the network,
#' rather than using the full gene set.
#' @param ... TODO
#'
#' @return TODO
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export
parallelNetworkWrapper <- function(data, regressionFunction,
                                   n_cores = 1, cluster_type = "FORK",
                                   regulatorIndex = NULL, ...) {
  data <- as.matrix(data)

  clust <- NULL
  if (n_cores > 1) {
    clust <- parallel::makeCluster(n_cores, type = cluster_type,
                                   outfile = file.path(paste0(regressionFunction, "_log.txt")))
    doParallel::registerDoParallel(clust)
  }

  genes_use <- colnames(data)

  # Subset to just regulator genes if defined
  if (!is.null(regulatorIndex)) {
    genes_use <- colnames(data)[regulatorIndex]
  }

  gene_number <- NULL # Necessary to pass R CMD check

  results <- foreach(gene_number = 1:ncol(data)) %dopar% {
    if (gene_number %% 100 == 0) {
      cat(paste0("Gene ", gene_number, "\n"))
    }

    res <- NA
    set.seed(gene_number)

    gene_query <- colnames(data)[gene_number]

    fxnArgs <- list()
    fxnArgs$y <- as.matrix(data[, gene_query])
    fxnArgs$x <- as.matrix(data[, setdiff(genes_use, gene_query)])

    if (regressionFunction %in% c('sparrowZ', 'sparrow2Z')) {
      fxnArgs$n_orderings <- 12
    }

    try(res <- do.call(regressionFunction,
                       c(fxnArgs, list(...))),
        silent = TRUE)

    if (!any(is.na(res))) {
      # Add 0-entries for the query gene, which are missing in the returned results
      res <- cbind(res, 0)
      colnames(res)[ncol(res)] <- gene_query
    } else {
      res <- matrix(0, ncol = length(genes_use),
                    dimnames = list("", genes_use))
    }

    # Put back in correct order
    res <- res[, genes_use]

    # If res has more than one row (i.e. results from more than one criterion),
    # separate each row into its own list, and ensure it is a 1 x n_genes matrix
    # where the row name is the query gene.
    # If res has only one row, this will be a one-item list.
    res_list <- lapply(rownames(res), function(rname) {
      res_l <- matrix(res[rname, ], nrow = 1, ncol = ncol(res),
                      dimnames = list(gene_query, genes_use))
      res_l
    })
    names(res_list) <- rownames(res)

    return(res_list)
  }

  if (!is.null(clust)) {
    parallel::stopCluster(clust)
  }

  # Results is a list of lists: Each item in the top-level list contains one or
  # more 1-row matrices named after the criterion that generated them (e.g.
  # "AIC", "lambda.min"). All items should have the same set of names, so we use
  # the names from the first item in the list.
  net_types <- names(results[[1]])

  # Extract the matrices for each network type (e.g. "AIC", "lambda.min") from
  # the results list, combine them into one matrix, and write to a file. Then,
  # return the name of the file that was just generated.
  networks <- sapply(net_types, function(net_type) {
    network <- do.call(rbind, lapply(results, "[[", net_type))

    # JB TODO why are we doing this? is it to make the matrix symmetrical?
    network <- network / 2 + t(network) / 2

    return(network)
  })

  names(networks) <- net_types

  return(networks)
}
