## Main Script: Validation and Consistency Checks for Data Paper
# Purpose: This script performs various checks to validate the integrity and consistency 
#          of species occurrence data, environmental and thematic variables, and species 
#          distribution models (SDMs).

# Define the directory
 #setwd("path/to/project")

# Source functions
 source("Scripts/functions_checks.R")

## Check species consistency ##

 validate_species("./sp_occurrence/trees/global", 
                  "./sp_occurrence/trees/regional")


## Check variables ##

 check_rasters("./env_variables/regional")
 check_rasters("./env_variables/global")
 check_rasters("./thematic/connectivity")
 check_rasters("./thematic/vegetation_types")
 check_rasters("./thematic/biodiversity")

## Check sdms ##

# SDMs for shrubs
 species_path <- "./sp_occurrence/shrubs/global"
 species_list <- list.files(species_path, full.names = FALSE, recursive = FALSE)
 species_list <- tools::file_path_sans_ext(species_list)

 results_shrubs <- check_sdm_files(species_list, 
                                   "./SDMs/shrubs")

# SDMs for trees (multiply models)
 species_path <- "./sp_occurrence/trees/global"
 species_list <- list.files(species_path, full.names = FALSE, recursive = FALSE)
 species_list <- tools::file_path_sans_ext(species_list)

 results_trees_multiply <- check_sdm_files(species_list, 
                                   "./SDMs/trees/multiply")

# SDMs for trees (covariate models)
 results_trees_covariate <- check_sdm_files(species_list, 
                                   "./SDMs/trees/covariate")

