library(synapser, quietly = TRUE)
#library(metanetwork, quietly = TRUE) # TODO temp
library(githubr, quietly = TRUE)
library(optparse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(tibble, quietly = TRUE)

## TODOs from Jaclyn:
#   * There needs to be some way to tell which orientation the data is in
#     (gene x sample or sample x gene), because it needs to be gene x sample for
#     winsorize() but sample x gene for everything else. This isn't documented
#     in a clear way anywhere.
#   * Provenance needs to be streamlined
#   * Synapse upload should be optional

# Obtaining the data - From User --------------------------------------------

option_list <- list(make_option(c("-u", "--synapse_authToken"),
                                type = "character",
                                action = "store",
                                help = "Synapse auth token"),
                    make_option(c("-c", "--config_file"),
                                type = "character",
                                action = "store",
                                help = "Path to the complete config file"))
req_args <- parse_args(OptionParser(option_list = option_list))

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

# Provenance -------------------------------------------------------------------

# Pull Git Provenance
thisRepo <- NULL
thisFile <- NULL

try(thisRepo <- githubr::getRepo(
  repository = config$provenance$code_annotations$repository,
  ref = config$provenance$code_annotations$ref,
  refName = config$provenance$code_annotations$ref_name
),
silent = TRUE)

try(thisFile <- githubr::getPermlink(
  repository = thisRepo,
  repositoryPath = config$provenance$code_annotations$repository_path
),
silent = TRUE)

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

net_methods <- config$input_profile$network_method

medium_algorithms <- c("lassoAIC", "lassoBIC", "lassoCV1se", "lassoCVmin",
                       "ridgeAIC", "ridgeBIC", "ridgeCV1se", "ridgeCVmin",
                       "sparrowZ", "sparrow2Z")

heavy_algorithms <- c("genie3", "tigress")

for (method in net_methods) {
  message(paste0("Running method ", method, "..."))

  output_files <- if (method == "c3net") {
    c3netWrapper(data,
                 outputpath = config$output_profile$output_path,
                 c3net_alpha = config$input_profile$c3net_alpha,
                 debug_save = config$output_profile$debug_save)

  } else if (method == "mrnet") {
    mrnetWrapper(data,
                 outputpath = config$output_profile$output_path,
                 mrnet_k = config$input_profile$mrnet_k,
                 debug_save = config$output_profile$debug_save)

  } else if (method == "wgcna") {
    wgcnaTOM(data,
             outputpath = config$output_profile$output_path,
             RsquaredCut = config$input_profile$wgcna_RsquaredCut,
             defaultPower = config$input_profile$wgcna_defaultPower,
             ncores = config$computing_specs$light_ncores,
             debug_save = config$output_profile$debug_save)

  } else if (method %in% medium_algorithms) {
    parallelNetworkWrapper(data,
                           nodes = config$computing_specs$medium_ncores,
                           regressionFunction = method,
                           outputpath = config$output_profile$output_path)

  } else if (method %in% heavy_algorithms) {
    parallelNetworkWrapper(data,
                           nodes = config$computing_specs$heavy_ncores,
                           regressionFunction = method,
                           outputpath = config$output_profile$output_path)
  } else {
    message(paste("Unrecognized algorithm name", method))
    NULL
  }

  ## Upload to Synapse

  dataFolder <- synapser::Folder(method, parent = config$output_profile$project_id)
  dataFolder <- synapser::synStore(dataFolder)

  config_file <- synapser::File(path = req_args$config_file,
                                parentId = dataFolder$id)

  Config_OBJ <- synapser::synStore(config_file,
                                   used = config$input_profile$expr_matrix_synid,
                                   executed = thisFile,
                                   forceVersion = FALSE)

  syn_config <- Config_OBJ$id

  for (file in output_files) {
    network_file <- synapser::File(path = file, parentId = dataFolder$id)

    Network_OBJ <- synapser::synStore(network_file,
                                      used = c(config$input_profile$expr_matrix_synid, syn_config),
                                      executed = thisFile,
                                      forceVersion = FALSE)
  }
}
