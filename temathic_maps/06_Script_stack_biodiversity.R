## stack sdms

 setwd("D:/UAM_DATAPAPER")

# libraries
library(terra)

## Binarization trees covariate 1km ##
 base_species_dir <- "./SDMs/trees/covariate"
 threshold_dir <- "./SDMs_org/Arboles_1km_covariate/Covariate/Projections"
 output_dir <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"

 if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
 }

# List spp
 species_list <- list.dirs(base_species_dir, full.names = FALSE, recursive = FALSE)
 species_list <- setdiff(species_list, "Tamarix.canariensis") # rm T.canariensis

# Iterar
# Función para binarizar mapas basados en un threshold para múltiples especies
binarize_species_maps <- function(species_list, base_species_dir, threshold_dir, output_dir, pattern = "Current\\.tif$") {
  for (species in species_list) {
    # species <- species_list[1] #44 Pistacia.terebinthus, 32 Juniperus sabina
    # threshold
    threshold_file <- file.path(base_species_dir, species, "values", paste0(species, "_ensemble.csv"))
  
    if (!file.exists(threshold_file)) {
      cat("Threshold file not found for species:", species, "\n")
      next
    }
  
    threshold_data <- read.csv(threshold_file)
    threshold <- threshold_data$cutoff[threshold_data$metric.eval == "TSS"]
  
    if (length(threshold) == 0) {
      cat("No TSS threshold found for species:", species, "\n")
      next
    }
  
    # raster continuo
    raster_file <- list.files(threshold_dir, 
                              pattern = paste0(species, ".*", pattern), 
                              full.names = TRUE)
  
    if (length(raster_file) == 0) {
      cat("Continuous raster not found for species:", species, "\n")
      next
    }
  
    r <- rast(raster_file[1])
  
    # Binarizar el raster
    values(r) <- ifelse(is.na(values(r)), NA, ifelse(values(r) >= threshold, 1, 0))
  
    # save
    pattern_suffix <- gsub(".*\\*|\\\\|\\.tif\\$", "", pattern)
    output_file <- file.path(output_dir, paste0(species, ".", pattern_suffix, ".bin.TSS.tif"))
    writeRaster(r, output_file, overwrite = TRUE)
  
    cat("Binarized map saved for species:", species, "\n")
  }
}

# Run
# pattern <- ".*Current\\.tif$"
# pattern <- ".*IPSL_CM6A_LR_2070_SSP126\\.tif$"
# pattern <- ".*IPSL_CM6A_LR_2070_SSP585\\.tif$"
# pattern <- ".*MRI_ESM2_0_2070_SSP126\\.tif$"
 pattern <- ".*MRI_ESM2_0_2070_SSP585\\.tif$"

binarize_species_maps(species_list, 
                      base_species_dir, 
                      threshold_dir, 
                      output_dir, 
                      pattern = pattern)


## Reproject Binarizations from latlon to utm
# ! ! ! FUntion at 3_Script_Prepare_sdms_str_reproject_rename.R

# trees 1km resolution
 ref_file <- "./SDMs/shrubs/Acer.granatense/projections/Acer.granatense.Current.tif"
 base_folder <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"
 equalize_projections(base_folder, path_reference_tif = ref_file)



## Stack binary maps shrubs and trees covariate at 1km ##

# Función to stack rasters with patron
sum_rasters_by_pattern <- function(base_path, pattern, exclude_patterns = NULL) {
  raster_files <- list.files(base_path, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  if (!is.null(exclude_patterns)) {
    exclude_regex <- paste(exclude_patterns, collapse = "|")  # Combina los patrones en una sola expresión regular
    raster_files <- raster_files[!grepl(exclude_regex, raster_files)]
  }

  if (length(raster_files) == 0) {
    stop("No raster files found matching the pattern: ", pattern)
  }
  
  raster_stack <- rast(raster_files)
  raster_sum <- sum(raster_stack, na.rm = TRUE)
  
  cat("\nDirectory:", base_path, " | Stacked files: ", length(raster_files), "| Pattern: ", pattern, "\n")

  return(raster_sum)
}


# Run
# Current.bin.TSS.tif
# shrubs
 base_path <- "./SDMs/shrubs"
 current_tss_sum_shrubs <- sum_rasters_by_pattern(base_path, 
                                    pattern = "Current\\.bin\\.TSS\\.tif$")

#trees
 base_path <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"
 exclude <- c("Tamarix\\.canariensis")
 current_tss_sum_trees <- sum_rasters_by_pattern(base_path, 
                                    patter = "Current\\.bin\\.TSS\\.tif$")

 raster_sum_current <- current_tss_sum_shrubs + current_tss_sum_trees
 names(raster_sum_current) <- "spp_richness"
 #writeRaster(raster_sum_current, "./thematic/biodiversity/woody/current/spp_richness.tif", overwrite = TRUE)


# IPSL_CM6A_LR_2070_SSP126.bin.TSS.tif
# shrubs
 base_path <- "./SDMs/shrubs"
 tss_sum_shrubs <- sum_rasters_by_pattern(base_path, 
                                    pattern = "IPSL_CM6A_LR_2070_SSP126\\.bin\\.TSS\\.tif$")

#trees
 base_path <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"
 tss_sum_trees <- sum_rasters_by_pattern(base_path, 
                                    patter = "IPSL_CM6A_LR_2070_SSP126\\.bin\\.TSS\\.tif$")

 raster_sum_ipsl126 <- tss_sum_shrubs + tss_sum_trees
 names(raster_sum_ipsl126) <- "spp_richness"
 #writeRaster(raster_sum_ipsl126, "./thematic/biodiversity/woody/IPSL_CM6A_LR_2070_SSP126/spp_richness.tif", overwrite = TRUE)


# IPSL_CM6A_LR_2070_SSP585.bin.TSS.tif
# shrubs
 base_path <- "./SDMs/shrubs"
 tss_sum_shrubs <- sum_rasters_by_pattern(base_path, 
                                    pattern = "IPSL_CM6A_LR_2070_SSP585\\.bin\\.TSS\\.tif$")

#trees
 base_path <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"
 tss_sum_trees <- sum_rasters_by_pattern(base_path, 
                                    patter = "IPSL_CM6A_LR_2070_SSP585\\.bin\\.TSS\\.tif$")

 raster_sum_ipsl585 <- tss_sum_shrubs + tss_sum_trees
 names(raster_sum_ipsl585) <- "spp_richness"
 #writeRaster(raster_sum_ipsl585, "./thematic/biodiversity/woody/IPSL_CM6A_LR_2070_SSP585/spp_richness.tif", overwrite = TRUE)



# MRI_ESM2_0_2070_SSP126.bin.TSS.tif
# shrubs
 base_path <- "./SDMs/shrubs"
 tss_sum_shrubs <- sum_rasters_by_pattern(base_path, 
                                    pattern = "MRI_ESM2_0_2070_SSP126\\.bin\\.TSS\\.tif$")

#trees
 base_path <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"
 tss_sum_trees <- sum_rasters_by_pattern(base_path, 
                                    patter = "MRI_ESM2_0_2070_SSP126\\.bin\\.TSS\\.tif$")

 raster_sum_mri126 <- tss_sum_shrubs + tss_sum_trees
 names(raster_sum_mri126) <- "spp_richness"
 #writeRaster(raster_sum_mri126, "./thematic/biodiversity/woody/MRI_ESM2_0_2070_SSP126/spp_richness.tif", overwrite = TRUE)


# MRI_ESM2_0_2070_SSP585.TSS.tif
# shrubs
 base_path <- "./SDMs/shrubs"
 tss_sum_shrubs <- sum_rasters_by_pattern(base_path, 
                                    pattern = "MRI_ESM2_0_2070_SSP585\\.bin\\.TSS\\.tif$")

#trees
 base_path <- "./SDMs_org/Arboles_1km_covariate/Covariate/Binarizations"
 tss_sum_trees <- sum_rasters_by_pattern(base_path, 
                                    patter = "MRI_ESM2_0_2070_SSP585\\.bin\\.TSS\\.tif$")

 raster_sum_mri585 <- tss_sum_shrubs + tss_sum_trees
 names(raster_sum_mri585) <- "spp_richness"
 #writeRaster(raster_sum_mri585, "./thematic/biodiversity/woody/MRI_ESM2_0_2070_SSP585/spp_richness.tif", overwrite = TRUE)


### Reproject woody
## Referencia (bio1)
# ref <- rast("./env_variables/global/climatic/current/bio1.tif")
#
# files <- list.files("./thematic/biodiversity/woody", full.names = TRUE, recursive = TRUE)
#
## it
# for (file in files) {
#   # file <- files[1]
#   raster <- rast(file)
#  
#   raster <- project(raster, ref)
#   raster  <- round(raster)
#
#   writeRaster(raster, file, overwrite = TRUE)
#  
#   cat("Processed and saved:", file, "\n")
# }


## Rename biodiversity/protected.
 ref <- rast("./env_variables/global/climatic/current/bio1.tif")

 files <- list.files("./thematic/biodiversity/protected", full.names = TRUE, recursive = TRUE)

# it
 for (file in files) {
   # file <- files[1]
   raster <- rast(file)
  
   raster <- project(raster, ref)
   raster  <- round(raster)
   names(raster) <- "protected_spp_richness"

   writeRaster(raster, file, overwrite = TRUE)
  
   cat("Processed and saved:", file, "\n")
 }

 #raster2 <- rast(files[1]) 


# Rename output file protected
 files <- list.files("./thematic/biodiversity/protected", full.names = TRUE, recursive = TRUE)

# it
 for (file in files) {
   # file <- files[1]
   raster <- rast(file)
  
   new_file <- file.path(dirname(file), "protected_spp_richness.tif")
  
   writeRaster(raster, new_file, overwrite = TRUE)
  
   file.remove(file)
  
   cat("Saved new file and deleted original:", new_file, "\n")
 }

 files <- list.files("./thematic/biodiversity/protected", full.names = TRUE, recursive = TRUE)
 #raster2 <- rast(files[1]) 



