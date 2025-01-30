#' Parallel Network Construction Wrapper
#'
#' Wrapper to run network detection in a parallel fashion, using a specified
#' regression function.
#'
#' @param data Expression matrix to be used for network construction, which
#'   should have samples as rows and genes as columns.
#' @param regressionFunction The name of the regression function to use. Current
#'   options are: lassoIC, lassoCV, ridgeIC, ridgeCV, tigress, and vbsr
#' @param log_file_path Optional. The folder path to where log files should be
#' stored. Log files capture any output during parallel execution. If omitted,
#' log files will be stored in the working directory.
#' @param n_cores Optional. The number of parallel cores to use.
#' @param regulatorIndex Optional. A vector of numerical indexes into
#'   \code{colnames(data)} for a subset of genes that should be used in the
#'   network. All genes in \code{data} will be tested against this subset,
#'   rather than against every gene. TODO should this stay?
#' @param ... Optional. Additional arguments that are passed through to the
#'   individual network algorithm.
#'
#' @return a named list, where each item is a network matrix as returned by the
#'   regression function. Regression functions may return multiple matrices or
#'   only one. If there was an error, the list will be empty.
#'
#' @export
constructNetwork.parallelRegression <- function(data, regressionFunction,
                                                log_file_path = ".",
                                                n_cores = 1,
                                                regulatorIndex = NULL,
                                                ...) {
  data <- as.matrix(data)

  clust <- NULL
  if (n_cores > 1) {
    log_file <- file.path(log_file_path, paste0(regressionFunction, "_log.txt"))
    clust <- parallel::makeCluster(n_cores, outfile = log_file)
  }

  genes_use <- colnames(data)

  # Subset to just regulator genes if defined
  if (!is.null(regulatorIndex)) {
    genes_use <- colnames(data)[regulatorIndex]
  }

  results <- parallel::parLapply(cl = clust,
                                 X = 1:ncol(data),
                                 fun = doRegressionFn,
                                 data, genes_use, regressionFunction, ...)

  if (!is.null(clust)) {
    parallel::stopCluster(clust)
  }

  # Results is a list of lists: Each item in the top-level list contains one or
  # more 1-row matrices named after the criterion that generated them (e.g.
  # "AIC", "lambda.min"). All items should have the same set of names, so we use
  # the names from the first item in the list.
  net_types <- names(results[[1]])

  # Extract the matrices for each network type (e.g. "AIC", "lambda.min") from
  # the results list and combine them into one matrix per type.
  networks <- lapply(net_types, function(net_type) {
    network <- do.call(rbind, lapply(results, "[[", net_type))

    # Make the matrix symmetrical
    network <- network / 2 + t(network) / 2
    return(network)
  })

  names(networks) <- net_types

  return(networks)
}


#' Do regression function in parallel
#'
#' This function is called once per gene to calculate the network in parallel.
#'
#' @param gene_number An integer that indicates which gene in \code{data} is
#'   being tested
#' @param data An expression matrix where rows are samples and columns are genes
#' @param genes_use A vector of gene names to include in the network. If this
#'   vector doesn't include all genes in \code{data}, all genes in \code{data}
#'   will be tested against this subset, rather than against every gene.
#' @param regressionFunction The name of the regression function to use. Current
#'   options are: lassoIC, lassoCV, ridgeIC, ridgeCV, vbsr
#' @param ... Optional. Additional arguments to be passed through to the
#'   regression function.
#'
#' @returns a named list, where each list item is a single vector of
#'   coefficients returned by the regression function. Several regression
#'   functions return coefficients for several different solutions, so these are
#'   put in separate list items.
doRegressionFn <- function(gene_number, data, genes_use, regressionFunction, ...) {
  if (gene_number %% 100 == 0) {
    cat(paste0("Gene ", gene_number, "\n"))
  }

  res <- NA
  set.seed(gene_number)

  gene_query <- colnames(data)[gene_number]

  fxnArgs <- list()
  fxnArgs$y <- as.matrix(data[, gene_query])
  fxnArgs$x <- as.matrix(data[, setdiff(genes_use, gene_query)])

  # Special case: the vbsr function is called "vbsrWrapper" to avoid name
  # collisions with the "vbsr" function in the vbsr package.
  if (regressionFunction == "vbsr") {
    try(res <- do.call(vbsrWrapper, c(fxnArgs, list(...))))

  } else {
    try(res <- do.call(regressionFunction, c(fxnArgs, list(...))))
  }

  # Add 0-entries for the query gene, which are missing in the returned results
  if (!any(is.na(res))) {
    # Some functions return a named vector instead of a matrix.
    if (!inherits(res, "matrix")) {
      res <- c(res, 0)
      names(res)[length(res)] <- gene_query
      res <- res[genes_use]

    } else {
      res <- cbind(res, 0)
      colnames(res)[ncol(res)] <- gene_query
      res <- res[, genes_use]
    }
  } else {
    res <- matrix(0, ncol = length(genes_use),
                  dimnames = list(gene_query, genes_use))
  }

  # If res has more than one row (i.e. results from more than one criterion),
  # separate each row into its own list, and ensure it is a 1 x n_genes matrix
  # where the row name is the query gene.
  if (inherits(res, "matrix") && nrow(res) > 1) {
    res_list <- lapply(rownames(res), function(rname) {
      res_l <- matrix(res[rname, ], nrow = 1, ncol = ncol(res),
                      dimnames = list(gene_query, genes_use))
      res_l
    })

    names(res_list) <- rownames(res)

  } else {
    # If res has only one row or is a vector, this will be a one-item list with
    # a one-row matrix
    res <- matrix(res, nrow = 1, dimnames = list(gene_query, genes_use))
    res_list <- list(res)
    names(res_list) <- regressionFunction
  }

  return(res_list)
}
