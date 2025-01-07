# library(dplyr, quietly = TRUE)
# library(glmnet, quietly = TRUE)
# library(randomForest, quietly = TRUE)
# library(Hmisc, quietly = TRUE)
# library(lars, quietly = TRUE)
# library(WGCNA, quietly = TRUE)
library(synapser, quietly = TRUE)
library(metanetwork, quietly = TRUE)
# library(githubr, quietly = TRUE)
# library(c3net, quietly = TRUE)
library(optparse, quietly = TRUE)
# library(data.table, quietly = TRUE)
# library(parmigene, quietly = TRUE)
library(reader, quietly = TRUE)
library(Rmpi)
# library(parallel)
# library(doParallel)

## TODOs from Jaclyn:
#   * file paths need a / at the end with this code, it's not robust
#   * No parallel interface is registered for 'light' algorithms
#   * Commented out libraries above that weren't needed to run WGCNA. Will
#     un-comment as they are needed for more network functions, then delete
#     unneeded ones.
#   * pval_wgcna in config is never used

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

# Setting up the cofig file
# Sys.setenv(R_CONFIG_ACTIVE = "default") # TODO is this necessary?
config <- config::get(file = req_args$config_file)

if (!dir.exists(config$input_profile$temp_storage_loc)) {
  dir.create(config$input_profile$temp_storage_loc, recursive = TRUE)
}
if (!dir.exists(config$output_profile$output_path)) {
  dir.create(config$output_profile$output_path, recursive = TRUE)
}

# Log in to Synapse
synLogin(authToken = req_args$synapse_authToken)

# Data
synID_input <- config$input_profile$input_synid
data <- synGet(synID_input,
               downloadLocation = config$input_profile$temp_storage_loc,
               ifcollision = "overwrite.local")

# Registering the parallel clusters
if (config$computing_specs$medium_ncores > 0) {
  nslaves <- config$computing_specs$medium_ncores
  mpi.spawn.Rslaves(nslaves = nslaves, hosts = NULL)
}
if (config$computing_specs$heavy_ncores > 0) {
  nslaves <- config$computing_specs$heavy_ncores
  mpi.spawn.Rslaves(nslaves = nslaves, hosts = NULL)
}

# Performing the analysis -------------------------------------------------

net_methods <- config$input_profile$network_method
data <- reader::reader(data$path)

if (is.null(config$input_profile$na_fill)) {
  print("Data not normalized for missing values. Ignore if using mrnet method.")
} else {
  if (config$input_profile$na_fill == "Winsorize") {
    data <- metanetwork::winsorizeData(data)
  }
}

for (method in net_methods) {
  # Assuming we have more methods - not developing for now
  switch(
    method,
    "c3net" = c3netWrapper(data,
                           pval = config$input_profile$p_val_c3net,
                           outputpath = config$output_profile$output_path),
    "mrnet" = mrnetWrapper(data,
                           pval = config$input_profile$p_val_mrnet,
                           outputpath = config$output_profile$output_path,
                           tool_storage_loc = config$input_profile$temp_storage_loc),
    "wgcna" = wgcnaTOM(data,
                       outputpath = config$output_profile$output_path,
                       RsquaredCut = config$input_profile$rsquaredCut,
                       defaultNaPower = config$input_profile$defaultnaPower),
    # TODO all of these use the exact same arguments except for medium vs heavy cores,
    # this could be condensed
    "lassoAIC" = mpiWrapper(data,
                            nodes = config$computing_specs$medium_ncores,
                            pathv = NULL,
                            regressionFunction = method,
                            outputpath = config$output_profile$output_path),
    "lassoBIC" = mpiWrapper(data,
                            nodes = config$computing_specs$medium_ncores,
                            pathv = NULL,
                            regressionFunction = method,
                            outputpath = config$output_profile$output_path),
    "lassoCV1se" = mpiWrapper(data,
                              nodes = config$computing_specs$medium_ncores,
                              pathv = NULL,
                              regressionFunction = method,
                              outputpath = config$output_profile$output_path),
    "lassoCVmin" = mpiWrapper(data,
                              nodes = config$computing_specs$medium_ncores,
                              pathv = NULL,
                              regressionFunction = method,
                              outputpath = config$output_profile$output_path),
    "ridgeAIC" = mpiWrapper(data,
                            nodes = config$computing_specs$medium_ncores,
                            pathv = NULL,
                            regressionFunction = method,
                            outputpath = config$output_profile$output_path),
    "ridgeBIC" = mpiWrapper(data,
                            nodes = config$computing_specs$medium_ncores,
                            pathv = NULL,
                            regressionFunction = method,
                            outputpath = config$output_profile$output_path),
    "ridgeCV1se" = mpiWrapper(data,
                              nodes = config$computing_specs$medium_ncores,
                              pathv = NULL,
                              regressionFunction = method,
                              outputpath = config$output_profile$output_path),
    "ridgeCVmin" = mpiWrapper(data,
                              nodes = config$computing_specs$medium_ncores,
                              pathv = NULL,
                              regressionFunction = method,
                              outputpath = config$output_profile$output_path),
    "sparrowZ" = mpiWrapper(data,
                            nodes = config$computing_specs$medium_ncores,
                            pathv = NULL,
                            regressionFunction = method,
                            outputpath = config$output_profile$output_path),
    "sparrow2Z" = mpiWrapper(data,
                             nodes = config$computing_specs$medium_ncores,
                             pathv = NULL,
                             regressionFunction = method,
                             outputpath = config$output_profile$output_path),
    "genie3" = mpiWrapper(data,
                          nodes = config$computing_specs$heavy_ncores,
                          pathv = NULL,
                          regressionFunction = method,
                          outputpath = config$output_profile$output_path),
    "tigress" = mpiWrapper(data,
                           nodes = config$computing_specs$heavy_ncores,
                           pathv = NULL,
                           regressionFunction = method,
                           outputpath = config$output_profile$output_path)
  )
}

if (config$computing_specs$heavy_ncores > 0) {
  mpi.close.Rslaves()
}

if (config$computing_specs$medium_ncores > 0) {
  mpi.close.Rslaves()
}


# Obtaining the data - For provenance --------------------------------------------

all.annotations <- synGetAnnotations(config$input_profile$input_synid)

# JB TODO this doesn't actually work -- it returns NULL, not a list of things
checkAnnotations <- function(annotations, config) {
  annot_default <- list(
    dataType = NULL,
    resourceType = NULL,
    metadataType = NULL,
    isModelSystem = NULL,
    isMultiSpecimen = NULL,
    fileFormat = NULL,
    grant = NULL,
    species = NULL,
    organ = NULL,
    tissue = NULL,
    study = NULL,
    consortium = NULL,
    assay = NULL
  )
  for (item in names(annot_default)) {
    if (!is.null(config$provenance$annotations$item)) {
      annot_default$item <- config$provenance$annotations$item[[1]]
    } else if (!is.null(annotations$item)) {
      annot_default$item <- annotations$item[[1]]
    }
  }
}

all.annotations <- checkAnnotations(all.annotations, config)

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

## Store the file
dataFolder <- Folder(method, parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)

#### - push the config file req_args$config_file

config_file <- synapser::File(path = req_args$config_file,
                              name = basename(req_args$config_file),
                              parentId = dataFolder$id)

Config_OBJ <- synapser::synStore(
  config_file,
  used = config$input_profile$input_synid,
  activityName = config$provenance$activity_name,
  executed = thisFile,
  activityDescription = config$provenance$activity_description,
  forceVersion = FALSE
)

syn_config <- Config_OBJ$id

####

# JB TODO I think this chunk of code needs to be INSIDE the for loop? 'method' is only defined there...
output_files <- list.files(config$output_profile$output_path,
                           pattern = method,
                           full.names = TRUE)

# JB TODO this is a little over-complicated
if (method == "wgcna") {
  filePath <- gsub("//", "/", output_files)
  syn_name <- gsub("\\.txt", "", gsub("\\.csv", "", basename(filePath)))
  syn_name <- gsub("wgcna", "wgcna ", gsub("Network", " Network", syn_name))
  syn_name <- gsub("Power", "Power ", gsub("Soft", "Soft ", syn_name))
  syn_name <- gsub("Overlap", " Overlap ", syn_name)
} else {
  filePath <- gsub("//", "/", output_files)
  syn_name <- config$output_profile$output_name
}

for (file in filePath) {
  network_file <- synapser::File(path = file,
                                 name = syn_name[which(filePath == file)],
                                 parentId = dataFolder$id)
  
  ENRICH_OBJ <- synapser::synStore(
    network_file,
    used = c(config$input_profile$input_synid, syn_config),
    activityName = config$provenance$activity_name,
    executed = thisFile,
    activityDescription = config$provenance$activity_description,
    forceVersion = FALSE
  )

  # TODO synSetAnnotations doesn't exist. Is this necessary?
  # synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
}

if ((!is.na(config$computing_specs$heavy_ncores)) ||
    (!is.na(config$computing_specs$medium_ncores))) {
  Rmpi::mpi.quit(save = "no")
}
