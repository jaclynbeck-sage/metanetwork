library(synapser, quietly = TRUE)
#library(metanetwork, quietly = TRUE) # TODO temp
library(githubr, quietly = TRUE)
library(optparse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(tibble, quietly = TRUE)
#library(Rmpi) # JB TODO Trying to remove RPMI and rely on "parallel" to simplify things for users
#library(parallel)
#library(doParallel)

## TODOs from Jaclyn:
#   * No parallel interface is registered for 'light' algorithms

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

# TODO JB do we need light cluster thread setup?
if (config$computing_specs$light_ncores > 0) {
  #nslaves <- config$computing_specs$light_ncores
  #mpi.spawn.Rslaves(nslaves = nslaves, hosts = NULL)
  #cl <- parallel::makeCluster(config$computing_specs$light_ncores, type = "FORK", outfile = "")
}

# Registering the parallel clusters
if (config$computing_specs$medium_ncores > 0) {
  #nslaves <- config$computing_specs$medium_ncores
  #mpi.spawn.Rslaves(nslaves = nslaves, hosts = NULL)
}
if (config$computing_specs$heavy_ncores > 0) {
  #nslaves <- config$computing_specs$heavy_ncores
  #mpi.spawn.Rslaves(nslaves = nslaves, hosts = NULL)
}

# Performing the analysis -------------------------------------------------

net_methods <- config$input_profile$network_method

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

# JB TODO I actually don't think this needs to be a for loop at all, it looks
# like the intention is that net_methods is only one method at a time? Or that
# the intention was to allow multiple methods but the upload/provenance code
# was never moved into the for loop.
# JB TODO add debug_save arguments once that's put in the config files
for (method in net_methods) {
  # Assuming we have more methods - not developing for now # JB TODO I don't know what this comment means
  switch(
    method,
    "c3net" = c3netWrapper(data,
                           outputpath = config$output_profile$output_path,
                           c3net_alpha = config$input_profile$c3net_alpha,
                           debug_save = FALSE),
    "mrnet" = mrnetWrapper(data,
                           pval = config$input_profile$p_val_mrnet,
                           outputpath = config$output_profile$output_path,
                           tool_storage_loc = config$input_profile$temp_storage_loc),
    "wgcna" = wgcnaTOM(data,
                       outputpath = config$output_profile$output_path,
                       RsquaredCut = config$input_profile$wgcna_RsquaredCut,
                       defaultPower = config$input_profile$wgcna_defaultPower,
                       ncores = config$computing_specs$light_ncores,
                       debug_save = FALSE),
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
  #mpi.close.Rslaves()
}

if (config$computing_specs$medium_ncores > 0) {
  #mpi.close.Rslaves()
}


# Obtaining the data - For provenance --------------------------------------------

all.annotations <- synGetAnnotations(config$input_profile$input_synid)

# TODO this function is broken, calling annotations$item isn't subbing in the
# value of 'item'
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
    if (!is.null(config$provenance$annotations[[item]])) {
      annot_default[[item]] <- config$provenance$annotations[[item]][[1]]
    } else if (!is.null(annotations[[item]])) {
      annot_default[[item]] <- annotations[[item]][[1]]
    }
  }
  annot_default
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
# JB TODO use the vector of file names output by the functions inside the for loop
output_files <- list.files(config$output_profile$output_path,
                           pattern = method,
                           full.names = TRUE)

# JB TODO this is a little over-complicated
#if (method == "wgcna") {
#  filePath <- gsub("//", "/", output_files)
#  syn_name <- gsub("\\.txt", "", gsub("\\.csv", "", basename(filePath)))
#  syn_name <- gsub("wgcna", "wgcna ", gsub("Network", " Network", syn_name))
#  syn_name <- gsub("Power", "Power ", gsub("Soft", "Soft ", syn_name))
#  syn_name <- gsub("Overlap", " Overlap ", syn_name)
#} else {
#  filePath <- gsub("//", "/", output_files)
#  syn_name <- config$output_profile$output_name
#}

for (file in output_files) {
  network_file <- synapser::File(path = file,
                                 #name = syn_name[which(filePath == file)],
                                 parentId = dataFolder$id)

  ENRICH_OBJ <- synapser::synStore(
    network_file,
    used = c(config$input_profile$input_synid, syn_config),
    activityName = config$provenance$activity_name,
    #executed = thisFile, # TODO temporary
    activityDescription = config$provenance$activity_description,
    forceVersion = FALSE
  )

  # TODO synSetAnnotations doesn't exist. Is this necessary?
  # synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
}

if ((!is.na(config$computing_specs$heavy_ncores)) ||
    (!is.na(config$computing_specs$medium_ncores))) {
  #Rmpi::mpi.quit(save = "no")
}
