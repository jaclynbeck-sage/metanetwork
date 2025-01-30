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

writeLines(paste("Total number of edges", sum(bicNetworks$network@x)))

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
algorithms <- c("fast_greedy", "infomap", "label_prop", "linkcommunities", "louvain", "megena", "spinglass", "walktrap")
results <- lapply(algorithms, function(alg) {
  message(paste0("Running method ", alg, "..."))
  # TODO algorithm args from config, pass on n_cores too
  res <- metanetwork::findModules(adj, method = alg, nperm = 3, min.module.size = 30, n_cores = 15)

  writeCSVFile(res, file.path(config$output_path, paste0(alg, ".csv")))
  return(res)
})

# Speakeasy
# speakeasy_r = metanetwork::findModules.speakeasy(adj)
# speakeasy_r['algorithms'] = 'speakeasy_r'

partition.adj <- results
names(partition.adj) <- algorithms

partition.adj <- mapply(function(mod, method) {
  mod = mod %>%
    as.data.frame() %>%
    dplyr::select(gene, module) %>%
    dplyr::mutate(value = 1, module = paste0(method, '.', module)) %>%
    tidyr::spread(module, value)
}, partition.adj, names(partition.adj), SIMPLIFY = FALSE) %>%
  plyr::join_all(type = "full")

partition.adj[is.na(partition.adj)] <- 0
rownames(partition.adj) <- partition.adj$gene
partition.adj$gene <- NULL

# Randomise gene order
set.seed(101)
partition.adj <- partition.adj[sample(1:nrow(partition.adj), nrow(partition.adj)), ]

partition.adj <- t(partition.adj)

mod <- metanetwork::findModules.consensusCluster(d = partition.adj,
                                                 maxK = 100,
                                                 reps = 50,
                                                 pGenes = 0.8,
                                                 clusterAlg = "hclust",
                                                 hclust_method = "average",
                                                 distance = "pearson",
                                                 changeCDFArea = 0.001,
                                                 nbreaks = 10,
                                                 seed = 1,
                                                 corUse = "everything",
                                                 verbose = TRUE,
                                                 useParallelFlag = TRUE)
