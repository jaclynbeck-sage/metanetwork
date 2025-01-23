# Obtaining the data - From User --------------------------------------------

option_list <- list(
  optparse::make_option(c("-u", "--synapse_authToken"),
    type = "character",
    action = "store",
    help = "Synapse auth token"
  ),
  optparse::make_option(c("-c", "--config_file"),
    type = "character",
    action = "store",
    help = "Path to the complete config file"
  )
)
req_args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
req_args$config_file <- "inst/config/network-module/module_template.yml"

# Obtaining the data - From Synapse --------------------------------------------

# Setting up the cofig file
config <- config::get(file = req_args$config_file)

if (!dir.exists(config$temp_storage_loc)) {
  dir.create(config$temp_storage_loc, recursive = TRUE)
}
if (!dir.exists(config$output_path)) {
  dir.create(config$output_path, recursive = TRUE)
}

# Log in to Synapse
synapser::synLogin(authToken = req_args$synapse_authToken)

consensus_file <- synapser::synGet(config$consensus_net_synid,
                                   downloadLocation = config$temp_storage_loc,
                                   ifcollision = "overwrite.local")

bic_file <- synapser::synGet(config$bic_file_synid,
                             downloadLocation = config$temp_storage_loc,
                             ifcollision = "overwrite.local")

#### Get input data from synapse and formulate adjacency matrix ####
# Get bicNetworks.rda
bicNetworks <- readRDS(bic_file$path)

writeLines(paste("Total number of edges", sum(bicNetworks$network@x)))

# Get rank consensus network for weights
rank.cons <- loadCSVFile(consensus_file$path)

# Formulate adjacency matrix
adj <- rank.cons
adj[!as.matrix(bicNetworks$network)] <- 0

rm(bicNetworks, rank.cons)
gc()

#### Compute modules using specified algorithm ####

# TODO list of algorithms from config?

# Get a specific algorithm

# JB TODO remove?
# Running CF Finder - Under comments since it is not working wiht lisencing
# cf_loc = synGet('syn7806853',downloadLocation = '/home/sage')
# temp_command <- paste0("unzip ",cf_loc$path," -d /home/sage/")
# system(temp_command)
# CFinder = metanetwork::findModules.CFinder(adj, '/home/sage/CFinder-2.0.6--1448/', nperm = 3, min.module.size = 30)

# TODO nperm and min.module.size should be configurable
algorithms <- c("fast_greedy", "infomap", "label_prop", "linkcommunities", "louvain", "megena", "spinglass", "walktrap")
results <- lapply(algorithms, function(alg) {
  message(paste0("Running method ", alg, "..."))
  res <- findModules(adj, method = alg, nperm = 3, min.module.size = 30)

  writeCSVFile(res, file.path(config$output_path, paste0(alg, ".csv")))
  return(res)
})

names(results) <- algorithms

# megena
# Data
synID_input <- config$input_synid
data <- synapser::synGet(config$expr_matrix_synid,
  downloadLocation = config$temp_storage_loc,
  ifcollision = "overwrite.local"
)
exprData <- loadCSVFile(data)

# TODO JB megena::calculate.correlation finds 0 things that pass FDR.cutoff on rosmap data... that doesn't seem right
megena <- findModules.megena(
  exprData,
  method = "pearson",
  FDR.cutoff = 0.05,
  module.pval = 0.05,
  hub.pval = 0.05,
  doPar = TRUE
)
megena["algorithms"] <- "megena"
cat("Completed MEGENA algorithm \n")


# Speakeasy
# speakeasy_r = metanetwork::findModules.speakeasy(adj)
# speakeasy_r['algorithms'] = 'speakeasy_r'
