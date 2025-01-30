#' Wrapper for MRNET, ARACNE, and CLR from \code{parmigene}
#'
#' This function uses \code{parmigene} to calculate network matrices using the
#' MRNET, ARACNE, and CLR algorithms.
#'
#' A mutual information matrix is calculated with \code{parmigene::knnmi.all},
#' and then three network matrices are calculated from this mutual information
#' matrix using MRNET, ARACNE, and CLR. The results from all three algorithms
#' are returned as a named list.
#'
#' @param data A gene expression matrix with rows as sample IDs and columns as
#'   Gene or feature IDs. The matrix will be transposed to work with
#'   \code{parmigene}.
#' @param ... Optional. Additional arguments accepted by either
#'   \code{parmigene::knnmi.all} or \code{parmigene::aracne.m}. The \code{clr}
#'   and \code{mrnet} do not accept additional arguments.
#'
#' @return a named list of three matrices, where "aracne" is the output of
#'   \code{parmigene::aracne.m}, "clr" is the output of \code{parmigene::clr},
#'   and "mrnet" is the output of \code{parmigene::mrnet}.
#'
#' @export
constructNetwork.parmigene <- function(data, ...) {
  mi <- R.utils::doCall(parmigene::knnmi.all,
                        mat = t(data),
                        args = list(...),
                        .ignoreUnusedArgs = TRUE)

  arac <- R.utils::doCall(parmigene::aracne.m,
                          mi = mi,
                          args = list(...),
                          .ignoreUnusedArgs = TRUE)

  networks <- list(aracne = arac,
                   clr = parmigene::clr(mi),
                   mrnet = parmigene::mrnet(mi))

  return(networks)
}
