#' Construct Network
#'
#' Construct a network matrix or matrices using the specified method.
#'
#' Several methods return a list of multiple network matrices, while some return
#' a single-item list with one matrix. The network(s) will be saved to disk
#' using the provided path and base file name. Current methods available are:
#'
#' \itemize{
#'  \item `c3net`: uses [c3net::c3net]
#'  \item `genie3`: uses [GENIE3::GENIE3]
#'  \item `megena`: uses [MEGENA::calculate.PFN] to rank edges
#'  \item `parmigene`: uses [parmigene::knnmi.all] to calculate a mutual
#'        information matrix, and generates 3 separate networks from
#'        [parmigene::aracne.m], [parmigene::mrnet], and [parmigene::clr].
#'  \item `wgcna`: uses [WGCNA::pickSoftThreshold] to pick a power for the
#'        adjacency matrix, calculates the matrix with [WGCNA::adjacency], and
#'        computes the topological overlap matrix with [WGCNA::TOMsimilarity].
#'        Both the adjacency matrix and the TOM matrix are returned.
#'  \item `lassoCV` and `ridgeCV`: uses [glmnet::cv.glmnet] with lasso or
#'  ridge regression, respectively, to construct the network using
#'  cross-validation. The two networks generated using the best `lambda`
#'  (determined by two different criteria) are returned.
#'  \item `lassoIC` and `ridgeIC`: uses [glmnet::glmnet] with lasso or ridge
#'        regression, respectively, to construct the network. The two networks
#'        generated using the best `lambda` (determined by the best AIC or BIC)
#'        are returned.
#'  \item `tigress`: uses [metanetwork::tigress] to construct the network using
#'        the TIGRESS algorithm. NOTE: there is an R package for this algorithm
#'        called `tigress`, however the package does not parallelize in a
#'        memory-efficient way and takes much longer to run than this package's
#'        implementation.
#'  \item `vbsr`: uses [vbsr::vbsr] to compute the network using variational
#'        Bayes spike regression.
#' }
#'
#' @param data A matrix, or an object that can be coerced to a matrix,
#'   containing gene expression values. Rows should be samples and columns
#'   should be genes.
#' @param method_name The name of the method to use. Current accepted values are
#'   "c3net", "genie3", "megena", "parmigene", "wgcna", "lassoCV", "lassoIC",
#'   "ridgeCV", "ridgeIC", "tigress", and "vbsr".
#' @param n_cores Optional. The number of cores to use for algorithms that can
#'   use threads. A value of 1 (default) will result in no threading.
#' @param save_to_disk Optional. If `TRUE`, the network will be saved as a CSV
#'   file with the location and name specified by the `output_filepath` and
#'   `output_filename_base` arguments.
#' @param output_filepath Optional, only used if `save_to_disk` is `TRUE` or if
#'   running one of the parallel regression algorithms. The path to the folder
#'   where results should be stored, if saving to disk. This is also the path
#'   where a log file will be stored for parallel execution. If omitted, results
#'   and logs will be stored in the working directory where the code is
#'   executed.
#' @param output_filename_base Optional, only used if `save_to_disk` is `TRUE`.
#'   The base name of the output file(s) without any extension. In cases where
#'   the network method returns a single matrix, the matrix will be stored at
#'   `<output_filepath>/<output_filename_base>.csv`. When a network method
#'   returns more than one matrix, the name of each matrix will be appended to
#'   `output_filename_base` in the file name, e.g.
#'   `<output_filename_base>_AIC.csv` and `<output_filename_base>_BIC.csv` for a
#'   method that returns a list of two matrices named "AIC" and "BIC".
#' @param ... Optional, additional arguments to pass through to the individual
#'   algorithm function call(s).
#'
#' @returns a single matrix, or a named list where each item is a matrix, for
#'   algorithms that return more than one network. If there was an error, the
#'   list will be empty.
#'
#' @export
#'
#' @examples
#' # Run "parmigene", adding the "k" argument used in knnmi.all()
#' data <- matrix(rnorm(50000), ncol = 1000)
#' network_list <- constructNetwork(data, method_name = "parmigene", k = 7)
constructNetwork <- function(data,
                             method_name,
                             n_cores = 1,
                             save_to_disk = FALSE,
                             output_filepath = ".",
                             output_filename_base = "network",
                             ...) {
  # Ensure data is a matrix
  data <- as.matrix(data)

  networks <- switch(method_name,
    c3net = c3net::c3net(t(data), ...),
    genie3 = GENIE3::GENIE3(t(data), nCores = n_cores, ...),
    megena = constructNetwork.megena(data, n_cores = n_cores,
                                     log_file_path = output_filepath, ...),
    parmigene = constructNetwork.parmigene(data, ...),
    wgcna = constructNetwork.wgcna(data, n_cores = n_cores, ...),
    # These algorithms all run from constructNetwork.parallelRegression()
    lassoCV = ,
    lassoIC = ,
    ridgeCV = ,
    ridgeIC = ,
    tigress = ,
    vbsr = constructNetwork.parallelRegression(data,
                                               regressionFunction = method_name,
                                               n_cores = n_cores,
                                               log_file_path = output_filepath,
                                               ...),
    # Default: unrecognized algorithm returns NULL
    NULL
  )

  if (is.null(networks)) {
    message(paste("Unrecognized algorithm name", method_name))
    return(list())
  }

  # For GENIE: Make sure all rows and columns of the network matrices are in the same order
  # as 'data' and that it's symmetrical across the diagonal
  if (method_name == "genie3") {
    networks <- networks[colnames(data), colnames(data)]
    networks <- networks / 2 + t(networks) / 2
  }

  if (save_to_disk) {
    output_prefix <- file.path(output_filepath, output_filename_base)

    if (inherits(networks, "matrix")) {
      output_filename <- paste0(output_prefix, ".csv")
      writeUpperTri(networks, output_filename)

    } else if (is.list(networks)) {
      if (length(networks) == 1) {
        output_filename <- paste0(output_prefix, ".csv")
        writeUpperTri(networks[[1]], output_filename)

      } else {
        for (net_name in names(networks)) {
          output_filename <- paste0(output_prefix, "_", net_name, ".csv")
          writeUpperTri(networks[[net_name]], output_filename)
        }
      }
    }
  }

  return(networks)
}
