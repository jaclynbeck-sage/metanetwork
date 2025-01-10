#' This function installs the ARACNE system Function
#'
#' This function automatically installs ARACNE from `inst/aracne`.
#'
#' @param tool_storage_loc Required. Provides the directory to temporarily store
#' the ARACNE files and package.
#'
#' @export
#' @return the file path to the installation folder for ARACNE
#'
installAracne <- function(tool_storage_loc) {
  arc_path <- system.file("aracne", "ARACNE.src.tar.gz", package = "metanetwork")

  # Unzip Aracne into the temp files directory
  system(paste0('tar -xzvf ', arc_path, ' -C ', tool_storage_loc))
  install_path <- file.path(tool_storage_loc, "ARACNE")

  # Change to Aracne directory and install
  old_wd <- getwd()
  setwd(install_path)
  message(paste("Installing ARACNE to", install_path))

  system('make')

  setwd(old_wd)

  return(install_path)
}
