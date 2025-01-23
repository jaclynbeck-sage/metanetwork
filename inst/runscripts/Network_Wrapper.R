## TODOs from Jaclyn:
#   * There needs to be some way to tell which orientation the data is in
#     (gene x sample or sample x gene), because it needs to be gene x sample for
#     winsorize() but sample x gene for everything else. This isn't documented
#     in a clear way anywhere.
#   * Provenance needs to be streamlined
#   * Synapse upload should be optional
#   * Parallel cluster needs to be stopped on error

# Obtaining the data - From User --------------------------------------------

option_list <- list(optparse::make_option(c("-u", "--synapse_authToken"),
                                           type = "character",
                                           action = "store",
                                           help = "Synapse auth token"),
                    optparse::make_option(c("-c", "--config_file"),
                                           type = "character",
                                           action = "store",
                                           help = "Path to the complete config file"))
req_args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
req_args$config_file <- "inst/config/network-construction/construction_template.yml"

# Obtaining the data - From Synapse --------------------------------------------

# Setting up the config file
config <- config::get(file = req_args$config_file)

if (!dir.exists(config$input_profile$temp_storage_loc)) {
  dir.create(config$input_profile$temp_storage_loc, recursive = TRUE)
}
if (!dir.exists(config$output_profile$output_path)) {
  dir.create(config$output_profile$output_path, recursive = TRUE)
}

# Log in to Synapse
synapser::synLogin(authToken = req_args$synapse_authToken)

# Data
synID_input <- config$input_profile$expr_matrix_synid
data <- synapser::synGet(synID_input,
                         downloadLocation = config$input_profile$temp_storage_loc,
                         ifcollision = "overwrite.local")

# Performing the analysis ------------------------------------------------------

# Much faster than using reader or read.csv
data <- data.table::fread(file = data$path, sep = ",",
                          header = TRUE, data.table = FALSE)
data <- tibble::column_to_rownames(data, var = colnames(data)[1])

if (is.null(config$input_profile$na_fill)) {
  print("Data not normalized for missing values. Ignore if using mrnet method.")
} else {
  if (config$input_profile$na_fill == "Winsorize") {
    data <- metanetwork::winsorizeData(data)
  }
}

net_methods <- config$input_profile$network_methods

for (method in net_methods) {
  message(paste0("Running method ", method, "..."))

  algorithm_args <- if (method %in% names(config$algorithm_arguments)) {
    config$algorithm_arguments[[method]]
  } else {
    list()
  }

  tmp <- do.call(metanetwork::constructNetwork, args = c(
    list(data = data,
         method_name = method,
         save_to_disk = TRUE,
         output_filepath = config$output_profile$output_path,
         output_filename_base = method,
         n_cores = config$n_parallel_cores),
    algorithm_args # "..." in the constructNetwork function
  ))
}
