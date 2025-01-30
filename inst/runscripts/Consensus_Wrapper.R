# Obtaining the data - From Synapse --------------------------------------------

# Setting up the cofig file
config_file <- "inst/config/network-consensus/consensus_template.yml"
config <- config::get(file = config_file)

if (!dir.exists(config$input_profile$temp_storage_loc)) {
  dir.create(config$input_profile$temp_storage_loc, recursive = TRUE)
}
if (!dir.exists(config$output_profile$output_path)) {
  dir.create(config$output_profile$output_path, recursive = TRUE)
}

# Log in to Synapse
synapser::synLogin(authToken = req_args$synapse_authToken)

input_file <- synapser::synGet(config$input_profile$expr_matrix_synid,
                               downloadLocation = config$input_profile$temp_storage_loc,
                               ifcollision = "overwrite.local")

fileName <- input_file$path

outputpath <- config$output_profile$output_path

#child_obj <- synapser::synGetChildren(as.character(networkFolderId), includeTypes = list("folder"))
#child_list <- as.list(child_obj)
#child_names <- c()
#for (k in 1:length(child_list)) {
#  temp_name <- child_list[k]
#  temp_name <- unlist(temp_name)
#  nn <- temp_name["name"]
#  nn <- as.character(nn)
#  child_names <- c(child_names, nn)
#}

#out_list <- list()
#print(length(child_names))
#for (ent in 1:length(child_names)) {
#  print("st1")
#  temp_l <- synapser::synGetChildren(child_list[[ent]]$id)
#  temp_list <- temp_l$asList()
#  print(temp_list)
#  temp_names <- c()
#  temp_ids <- c()
#  for (k in 1:length(temp_list)) {
#    temp_name <- temp_list[k]
#    temp_name <- unlist(temp_name)
#    nn <- temp_name["name"]
#    nn <- as.character(nn)
#    temp_names <- c(temp_names, nn)
#    nn <- as.character(temp_name["id"])
#    temp_ids <- c(temp_ids, nn)
#  }
#
#  if (child_names[ent] == "wgcna") {
#    temp_name_search <- c(
#      "wgcna Topological Overlap Matrix Network",
#      "wgcna Soft Threshold Network"
#    )
#    temp_names_ids <- c(
#      grep(temp_name_search[1], temp_names),
#      grep(temp_name_search[2], temp_names)
#    )
#  } else {
#    temp_name_search <- paste0(child_names[[ent]], pattern_id)
#    temp_names_ids <- grep(temp_name_search, temp_names)
#  }
#
#  cat(paste0(temp_names[temp_names_ids], " \n"))
#  id_t <- temp_ids[temp_names_ids]
#  if (length(id_t) == 0) {
#    cat(paste0("No match for ", temp_name_search, " \n"))
#  } else {
#    for (id in id_t) {
#      temp <- synapser::synGet(id, downloadLocation = outputpath)
#      out_list <- append(out_list, temp)
#    }
#  }
#}

# TODO temporary
network_files <- list.files("~/meta_out/ros_network_synapse", pattern = ".csv", full.names = TRUE)
network_files <- network_files[!grepl("rankConsensus", network_files)]

data <- loadCSVFile(fileName)

message("Building Consensus Networks")

buildConsensus(network_files = network_files,
               exprData = data,
               outputpath = outputpath)


# Obtaining the data - For provenance --------------------------------------------

#activity <- synapser::synGet(config$input_profile$project_id)

#dataFolder <- synapser::Folder("Consensus", parent = config$input_profile$project_id)
#dataFolder <- synapser::synStore(dataFolder)

# file <- File(path = paste0(outputpath,'rankConsensusNetwork.csv'), parent = dataFolder)
# file <- synStore(file)
# file2 <- File(path = paste0(outputpath,'bicNetworks.rda'), parent = dataFolder)
# file2 <- synStore(file2)

#thisRepo <- NULL
#thisFile <- NULL

#try(
#  thisRepo <- githubr::getRepo(
#    repository = config$provenance$code_annotations$repository,
#    ref = config$provenance$code_annotations$ref,
#    refName = config$provenance$code_annotations$ref_name
#  ),
#  silent = TRUE
#)
#try(
#  thisFile <- githubr::getPermlink(
#    repository = thisRepo,
#    repositoryPath = config$provenance$code_annotations$repository_path
#  ),
#  silent = TRUE
#)
# Push Config File to synaopse
#ENRICH_OBJ <- synapser::synStore(
#  synapser::File(
#    path = req_args$config_file,
#    name = "Consensus Config",
#    parentId = dataFolder$properties$id
#  ),
#  used = config$input_profile$input_synid,
#  activityName = config$provenance$activity_name,
#  executed = thisFile,
#  activityDescription = config$provenance$activity_description
#)
#config_syn <- ENRICH_OBJ$properties$id
#
#ENRICH_OBJ <- synapser::synStore(
#  synapser::File(
#    path = paste0(outputpath, "bicNetworks.rda"),
#    name = "bicNetworks",
#    parentId = dataFolder$properties$id
#  ),
#  used = c(config$input_profile$input_synid, config_syn),
#  activityName = config$provenance$activity_name,
#  executed = thisFile,
#  activityDescription = config$provenance$activity_description
#)

#ENRICH_OBJ <- synapser::synStore(
#  synapser::File(
#    path = paste0(outputpath, "rankConsensusNetwork.csv"),
#    name = "RankConsensusNetwork",
#    parentId = dataFolder$properties$id
#  ),
#  used = c(config$input_profile$input_synid, config_syn),
#  activityName = config$provenance$activity_name,
#  executed = thisFile,
#  activityDescription = config$provenance$activity_description
#)
