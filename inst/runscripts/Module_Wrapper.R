# Obtaining the data - From User --------------------------------------------

config_file <- "inst/config/module_template.yml"

# Obtaining the data - From Synapse --------------------------------------------

# Setting up the cofig file
config <- config::get(file = config_file)

if (!dir.exists(config$temp_storage_loc)) {
  dir.create(config$temp_storage_loc, recursive = TRUE)
}
if (!dir.exists(config$output_path)) {
  dir.create(config$output_path, recursive = TRUE)
}

# Log in to Synapse
synapser::synLogin()

consensus_file <- synapser::synGet(config$consensus_net_synid,
                                   downloadLocation = config$temp_storage_loc,
                                   ifcollision = "overwrite.local")

bic_file <- synapser::synGet(config$bic_file_synid,
                             downloadLocation = config$temp_storage_loc,
                             ifcollision = "overwrite.local")

#### Get input data from synapse and formulate adjacency matrix ####
# Get bicNetworks.rda
bicNetworks <- readRDS(bic_file$path)

writeLines(paste("Total number of edges:", sum(bicNetworks$network@x)))

# Get rank consensus network for weights
# TODO loadCSVFile isn't exported in the package
rank.cons <- metanetwork::loadCSVFile(consensus_file$path)

# Formulate adjacency matrix
adj <- rank.cons
adj[!as.matrix(bicNetworks$network)] <- 0
adj <- adj * upper.tri(adj)

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
algorithms <- c("fast_greedy", "infomap", "label_prop", "linkcommunities",
                "louvain", "megena", "spinglass", "walktrap")
# TODO megena crashed at "Calculating distance metric and similarity...":
# Error in sample.int(length(x), size, replace, prob) : invalid first argument
# I think it was the second permutation
results <- lapply(algorithms, function(alg) {
  message(paste0("Running method ", alg, "..."))
  # TODO algorithm args from config
  res <- findModules(adj, method = alg, nperm = 3,
                                  min.module.size = 30, n_cores = config$n_cores)

  writeCSVFile(res, file.path(config$output_path, paste0(alg, ".csv")))
  return(res)
})

# Speakeasy
# speakeasy_r = metanetwork::findModules.speakeasy(adj)
# speakeasy_r['algorithms'] = 'speakeasy_r'

partition.adj <- results
names(partition.adj) <- algorithms

partition.adj <- mapply(function(mod, method) {
  mod <- as.data.frame(mod)
  mod$module <- paste0(method, ".", mod$module)
  mod$gene <- rownames(mod)
  return(mod)
}, partition.adj, names(partition.adj), SIMPLIFY = FALSE)

partition.adj <- do.call(rbind, partition.adj)
partition.adj$value <- 1
partition.adj <- tidyr::pivot_wider(partition.adj, names_from = "gene",
                                    values_from = "value", values_fill = 0)

partition.adj <- as.data.frame(partition.adj)
rownames(partition.adj) <- partition.adj$module
partition.adj$module <- NULL

mod <- metanetwork::findModules.consensusCluster(d = partition.adj,
                                                 maxK = 100,
                                                 n_subsamples = 100,
                                                 pGenes = 0.8,
                                                 clusterAlg = "hclust",
                                                 hclust_method = "average",
                                                 distance_metric = "pearson",
                                                 changeCDFArea = 0.001,
                                                 n_ks = 20,
                                                 corUse = "everything",
                                                 n_cores = config$n_cores,
                                                 log_file_path = "~/meta_out/modules",
                                                 verbose = TRUE,
                                                 seed = 101)
