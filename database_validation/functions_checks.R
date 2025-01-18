## Check species consistency ##

# Function to validate species at global and regional level
validate_species <- function(global_folder, regional_folder) {
  # List files
  global_files <- list.files(global_folder, pattern = "\\.csv$", full.names = TRUE)
  regional_files <- list.files(regional_folder, pattern = "\\.csv$", full.names = TRUE)

  # Extract species names
  extract_species <- function(files) {
    tools::file_path_sans_ext(basename(files))
  }

  global_species <- extract_species(global_files)
  regional_species <- extract_species(regional_files)

  # missing species?
  missing_in_global <- setdiff(regional_species, global_species)
  missing_in_regional <- setdiff(global_species, regional_species)

  # 
  if (length(missing_in_global) == 0 && length(missing_in_regional) == 0) {
    cat("All species match between global and regional folders.\n")
  } else {
    if (length(missing_in_global) > 0) {
      cat("Species in regional but missing in global:\n")
      print(missing_in_global)
    }
    if (length(missing_in_regional) > 0) {
      cat("Species in global but missing in regional:\n")
      print(missing_in_regional)
    }
  }
}


## Check variables ##

# Function to check raster resolution, projection, and extent
check_rasters <- function(base_folder, output_file = NULL) {
  tif_files <- list.files(base_folder, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  
  if (length(tif_files) == 0) {
    stop("No .tif files found in the specified folder.")
  }

  # initialize output redirection
  if (!is.null(output_file)) {
    sink(output_file)
  } 
  
  cat("Number of .tif files being checked:", length(tif_files), "\n")

  # Reference file
  ref_raster <- terra::rast(tif_files[1])
  ref_proj <- terra::crs(ref_raster)
  ref_res <- terra::res(ref_raster)
  ref_ext <- as.vector(terra::ext(ref_raster))
  
  #
  proj_mismatches <- c()
  proj_mismatch_details <- list()
  res_mismatches <- list()
  ext_mismatches <- list()
  
  for (file in tif_files) {
    raster <- terra::rast(file)
    current_proj <- terra::crs(raster)
    current_res <- terra::res(raster)
    current_ext <- as.vector(terra::ext(raster))
    
    # check projection
    if (!identical(current_proj, ref_proj)) {
      proj_mismatches <- c(proj_mismatches, file)
      proj_mismatch_details[[file]] <- current_proj
    }
    
    # check resolution
    if (!isTRUE(all.equal(current_res, ref_res))) {
      res_mismatches[[file]] <- current_res
    }
    
    # check extent
    if (!isTRUE(all.equal(current_ext, ref_ext))) {
      ext_mismatches[[file]] <- current_ext
    }
  }
  
  # print
  cat("\nProjection check:\n")
  if (length(proj_mismatches) == 0) {
    cat("  All files have the same projection:\n")
    cat("  ", ref_proj, "\n")
  } else {
    cat("  Reference Projection:\n  ", ref_proj, "\n")
    cat("  Files with different projections:\n")
    for (file in names(proj_mismatch_details)) {
      cat("  ", file, ":\n    ", proj_mismatch_details[[file]], "\n")
    }
  }
  
  cat("\nResolution check:\n")
  if (length(res_mismatches) == 0) {
    cat("  All files have the same resolution:\n")
    cat("  ", ref_res, "\n")
  } else {
    cat("  Reference Resolution:\n  ", ref_res, "\n")
    cat("  Files with different resolutions:\n")
    for (file in names(res_mismatches)) {
      cat("  ", file, ":", res_mismatches[[file]], "\n")
    }
  }
  
  cat("\nExtent check:\n")
  if (length(ext_mismatches) == 0) {
    cat("  All files have the same extent:\n")
    cat("  ", paste(ref_ext, collapse = ", "), "\n")
  } else {
    cat("  Reference Extent:\n  ", paste(ref_ext, collapse = ", "), "\n")
    cat("  Files with different extents:\n")
    for (file in names(ext_mismatches)) {
      cat("  ", file, ":", paste(ext_mismatches[[file]], collapse = ", "), "\n")
    }
  }

  if (!is.null(output_file)) {
    sink()
  }

}


## Check sdms ##

# Function to check sdms
check_sdm_files <- function(species_list, base_path, output_file = NULL) {
  if (!dir.exists(base_path)) {
    stop("The specified base_path does not exist: ", base_path)
  }

  # initialize output redirection
  if (!is.null(output_file)) {
    sink(output_file)
  } 

  results <- list()
  
  for (species in species_list) {
    projections_path <- file.path(base_path, species, "projections")
    values_path <- file.path(base_path, species, "values")
    
    # required files in "projections"
    projections_required <- c(
      paste0(species, ".Current.bin.ROC.tif"),
      paste0(species, ".Current.bin.TSS.tif"),
      paste0(species, ".Current.tif"),
      paste0(species, ".IPSL_CM6A_LR_2070_SSP126.tif"),
      paste0(species, ".IPSL_CM6A_LR_2070_SSP126.bin.ROC.tif"),
      paste0(species, ".IPSL_CM6A_LR_2070_SSP126.bin.TSS.tif"),
      paste0(species, ".IPSL_CM6A_LR_2070_SSP585.tif"),
      paste0(species, ".IPSL_CM6A_LR_2070_SSP585.bin.ROC.tif"),
      paste0(species, ".IPSL_CM6A_LR_2070_SSP585.bin.TSS.tif"),
      paste0(species, ".MRI_ESM2_0_2070_SSP126.tif"),
      paste0(species, ".MRI_ESM2_0_2070_SSP126.bin.ROC.tif"),
      paste0(species, ".MRI_ESM2_0_2070_SSP126.bin.TSS.tif"),
      paste0(species, ".MRI_ESM2_0_2070_SSP585.tif"),
      paste0(species, ".MRI_ESM2_0_2070_SSP585.bin.ROC.tif"),
      paste0(species, ".MRI_ESM2_0_2070_SSP585.bin.TSS.tif")
    )

    if (!grepl("multiply", projections_path, ignore.case = TRUE)) {
      emcv_files <- c(
        paste0(species, ".EMcv.tif"),
        paste0(species, ".IPSL_CM6A_LR_2070_SSP126.EMcv.tif"),
        paste0(species, ".IPSL_CM6A_LR_2070_SSP585.EMcv.tif"),
        paste0(species, ".MRI_ESM2_0_2070_SSP126.EMcv.tif"),
        paste0(species, ".MRI_ESM2_0_2070_SSP585.EMcv.tif")
       )
      projections_required <- c(projections_required, emcv_files)
    }
    
    # required files in "values"
    values_required <- c(
      paste0(species, "_ensemble.csv")
    )

    if (!grepl("multiply", values_path, ignore.case = TRUE)) {
      add_files <- c(
        paste0(species, "_indvar.csv"),
        paste0(species, "_nbestreplicates.csv"),
        paste0(species, "_replica.csv"),
        paste0(species, ".variables.csv")
       )
      values_required <- c(values_required, add_files)
    }
    
    # check projections
    if (dir.exists(projections_path)) {
      projections_present <- list.files(projections_path, full.names = FALSE)
      projections_missing <- setdiff(projections_required, projections_present)
    } else {
      projections_missing <- projections_required
    }
    
    # check values
    if (dir.exists(values_path)) {
      values_present <- list.files(values_path, full.names = FALSE)
      values_missing <- setdiff(values_required, values_present)
    } else {
      values_missing <- values_required
    }
    
    # save
    results[[species]] <- list(
      projections_missing = projections_missing,
      values_missing = values_missing
    )

    # print
    cat("\nSpecies:", species, "\n")
    if (length(projections_missing) > 0) {
      cat("  Directory proyections: missing files:\n")
      cat(paste("    -", projections_missing, collapse = "\n"), "\n")
    } else {
      cat("  Directory proyections: All good!.\n")
    }
    if (length(values_missing) > 0) {
      cat("  Directory values: missing files:\n")
      cat(paste("    -", values_missing, collapse = "\n"), "\n")
    } else {
      cat("  Directory values: All good!.\n")
    }
  }
  
  # stop output redirection 
  if (!is.null(output_file)) {
    sink()
  }
  
  return(results)
}




