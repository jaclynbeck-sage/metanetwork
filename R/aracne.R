#' This function applies ARACNE on the data
#'
#' This function takes in a gene expression matrix and applies the ARACNE gene
#' co-expression network analysis framework to find coexpressed gene pairs in
#' the matrix. The ARACNE framework is also installed from within the package
#' located in the `inst/tools/` directory. For more information on the ARACNE
#' framework see:
#' <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7>
#' and <https://califano.c2b2.columbia.edu/aracne>
#'
#' @param data Required. The gene expression matrix with genes as rows and
#' samples as columns.
#' @param pval Optional. Cutoff p-value to determine a coexpressed edge. (Default = 0.05)
#' @param outputpath Required. The path the resulting network should be saved to.
#' @param na_fill Optional. Value to replace `NA` values with. Ideally, a large
#' negative number or use min(data). If NULL, NA values will be removed
#' entirely. (Default = NULL)
#' @param tool_storage_loc Required. Provides the directory to temporarily store
#' the ARACNE files and package.
#'
#' @return the filename containing the network data. The name will be either
#'   `outputpath/aracneThresholdNetwork.csv` (if `pval` < 1) or
#'   `outputpath/aracneNetwork.csv` if `pval` is set to 1.
#' @export
aracne <- function(data, pval = NULL, outputpath, na_fill = NULL, tool_storage_loc) {
  aracne_path <- installAracne(tool_storage_loc = tool_storage_loc)

  if (is.null(pval)) {
    # JB TODO why is this divided? If the user puts in 0.05 instead of NULL it doesn't get divided.
    pval <- 0.05 / choose(nrow(data), 2)
  }

  if (is.null(na_fill)) {
    data <- stats::na.omit(data)
  } else{
    data[is.na(data)] <- na_fill
  }

  old_wd <- getwd()
  setwd(aracne_path)

  dataMatrix <- cbind(rownames(data), rownames(data), data)
  colnames(dataMatrix) <- c('name1', 'name2', colnames(data))

  data.table::fwrite(dataMatrix,
                     file = "dataMatrix.tsv",
                     sep = "\t",
                     row.names = FALSE,
                     col.names = TRUE,
                     quote = FALSE)

  command_string <- paste0(
    "./aracne2 -i dataMatrix.tsv -a adaptive_partitioning -p ",
    pval,
    " -o result.out"
  )

  system(command_string)

  result <- readLines("result.out")
  network <- matrix(0, nrow(data), nrow(data))
  rownames(network) <- rownames(data)
  colnames(network) <- rownames(data)

  fun <- function(x, ref) {
    vec <- strsplit(x, '\t')[[1]]
    model <- list(gene = vec[1],
                  vec = rep(0, length(ref)),
                  keepGene = vec[-1][1:length(vec[-1]) %% 2 == 1],
                  weights = vec[-1][1:length(vec[-1]) %% 2 == 0])

    names(model$vec) <- ref
    model$vec[model$keepGene] <- model$weights
    model$vec <- as.numeric(model$vec)
    return(model)
  }

  resultFormat <- lapply(result[-c(1:17)], fun, ref = rownames(data))
  for (i in 1:length(resultFormat)) {
    network[resultFormat[[i]]$gene, ] <- resultFormat[[i]]$vec
  }

  setwd(old_wd)

  if (pval == 1) {
    fileName <- file.path(outputpath, 'aracneNetwork.csv')
  } else{
    fileName <- file.path(outputpath, 'aracneThresholdNetwork.csv')
  }

  data.table::fwrite(network * upper.tri(network),
                     file = fileName,
                     sep = ",",
                     row.names = TRUE,
                     col.names = TRUE)
  return(fileName)
}
