# Obtaining the data - From Synapse --------------------------------------------

# Setting up the cofig file
config_file <- "inst/config/consensus_template.yml"
config <- config::get(file = config_file)

if (!dir.exists(config$input_profile$temp_storage_loc)) {
  dir.create(config$input_profile$temp_storage_loc, recursive = TRUE)
}
if (!dir.exists(config$output_profile$output_path)) {
  dir.create(config$output_profile$output_path, recursive = TRUE)
}

# Log in to Synapse
synapser::synLogin()

input_file <- synapser::synGet(config$input_profile$expr_matrix_synid,
                               downloadLocation = config$input_profile$temp_storage_loc,
                               ifcollision = "overwrite.local")

fileName <- input_file$path

outputpath <- config$output_profile$output_path

# TODO temporary
network_files <- list.files("~/meta_out/construction", pattern = ".csv", full.names = TRUE)

data <- loadCSVFile(fileName)

message("Building Consensus Networks")

metanetwork::buildConsensus(network_files = network_files,
                            exprData = data,
                            outputpath = outputpath)
