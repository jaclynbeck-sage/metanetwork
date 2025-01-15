#' Construct Network
#'
#' Construct a network matrix or matrices using the specified method.
#'
#' Several methods return a list of multiple network matrices, while some return
#' a single-item list with one matrix. The network(s) will be saved to disk
#' using the provided path and base file name. Current methods available are:
#'
#' \itemize{
#'  \item \code{c3net}: uses \code{c3net::c3net()}
#'  \item \code{genie3}: uses \code{GENIE3::GENIE3()}
#'  \item \code{mrnet}: uses \code{parmigene::knnmi.all} to calculate a mutual
#'        information matrix, and generates 3 separate networks from
#'        \code{parmigene::aracne.m}, \code{parmigene::mrnet}, and
#'        \code{parmigene::clr}.
#'  \item \code{tigress}: uses \code{tigress::tigress}
#'  \item \code{wgcna}: uses \code{WGCNA::pickSoftThreshold}` to pick a power
#'        for the adjacency matrix, calculates the matrix with
#'        \code{WGCNA::adjacency}, and computes the topological overlap matrix
#'        with \code{WGCNA::TOMsimilarity}. Both the adjacency matrix and the
#'        TOM matrix are returned.
#'  \item \code{lassoCV} and \code{ridgeCV}: uses \code{glmnet::cv.glmnet} with
#'        lasso or ridge regression, respectively, to construct the network using
#'        cross-validation. The two networks generated using the best \code{lambda}
#'        (determined by two different criteria) are returned.
#'  \item \code{lassoIC} and \code{ridgeIC}: uses \code{glmnet::glmnet} with
#'        lasso or ridge regression, respectively, to construct the network. The
#'        two networks generated using the best \code{lambda} (determined by the
#'        best AIC or BIC) are returned.
#'  \item \code{vbsr}: uses \code{vbsr::vbsr} to compute the network using
#'        variational Bayes spike regression.
#' }
#'
#' @param data A matrix, or an object that can be coerced to a matrix,
#'   containing gene expression values. Rows should be samples and columns
#'   should be genes.
#' @param method_name The name of the method to use. Current accepted
#'   values are "c3net", "genie3", "mrnet", "tigress", "wgcna", "lassoCV",
#'   "lassoIC", "ridgeCV", "ridgeIC", and "vbsr".
#' @param n_cores Optional. The number of cores to use for algorithms that can
#'   use threads. A value of 1 (default) will result in no threading.
#' @param save_to_disk Optional. If `TRUE`, the network will be saved as a CSV
#'   file with the location and name specified by the `output_filepath` and
#'   `output_filename_base` arguments.
#' @param output_filepath Optional, only used if `save_to_disk` is `TRUE`. The
#'   path to the folder where results should be stored. If omitted, results will
#'   be stored in the working directory where the code is executed.
#' @param output_filename_base Optional, only used if `save_to_disk` is `TRUE`.
#'   The base name of the output file(s) without any extension. In cases where
#'   the network method returns a single matrix, the matrix will be stored at
#'   `<output_filepath>/<output_filename_base>.csv`. When a network method
#'   returns more than one matrix, the name of each matrix will be appended to
#'   `output_filename_base` in the file name, e.g.
#'   `<output_filename_base>_AIC.csv` and `<output_filename_base>_BIC.csv` for a
#'   method that returns a list of two matrices named "AIC" and "BIC".
#'   Default: "network"
#' @param ... Optional, additional arguments to pass through to the individual
#'   algorithm function call(s).
#'
#' @returns a named list, where each item is a network matrix. If there was an
#'   error, the list will be empty.
#'
#' @export
#'
#' @examples
#' # Run "mrnet", adding the "k" argument used in knnmi.all()
#' data <- matrix(rnorm(50000), ncol = 1000)
#' network_list <- constructNetwork(data, method_name = "mrnet", k = 7)
constructNetwork <- function(data, method_name, n_cores = 1,
                             save_to_disk = FALSE, output_filepath = ".",
                             output_filename_base = "network", ...) {
  # Ensure data is a matrix
  data <- as.matrix(data)

  networks <- switch(method_name,
    c3net = c3net::c3net(t(data), ...),
    mrnet = mrnetWrapper(data, ...),
    wgcna = wgcnaWrapper(data, n_cores = n_cores, ...),
    genie3 = GENIE3::GENIE3(t(data), nCores = n_cores, ...),
    tigress = tigress::tigress(data, allsteps = FALSE, usemulticore = n_cores > 1, ...),
    # These algorithms all run from parallelNetworkWrapper()
    lassoCV = ,
    lassoIC = ,
    ridgeCV = ,
    ridgeIC = ,
    vbsr = parallelNetworkWrapper(data,
                                  n_cores = n_cores,
                                  regressionFunction = method_name,
                                  ...),
    # Default: unrecognized algorithm returns NULL
    NULL
  )

  if (is.null(networks)) {
    message(paste("Unrecognized algorithm name", method_name))
    return(list())
  }

  if (save_to_disk) {
    write_upper_tri <- function(mat, filename) {
      data.table::fwrite(mat * upper.tri(mat),
                         file = filename,
                         sep = ",",
                         row.names = TRUE,
                         col.names = TRUE)
    }

    output_prefix <- file.path(output_filepath, output_filename_base)

    if (inherits(networks, "matrix")) {
      output_filename <- paste0(output_prefix, ".csv")
      write_upper_tri(networks, output_filename)

    } else if (is.list(networks)) {
      for (net_name in names(networks)) {
        output_filename <- paste0(output_prefix, "_", net_name, ".csv")
        write_upper_tri(networks[[net_name]], output_filename)
      }
    }
  }

  return(networks)
}
