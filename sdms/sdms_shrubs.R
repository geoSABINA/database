

# Loading necessary packages
list.of.packages <- c(
  "terra",
  "covsel",
  "biomod2",
  "ecospat",
  "fs",
  "tidyverse",
  "gam",
  "plyr",
  "furrr",
  "sabinaNSDM")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

# Load the packages
for (package in list.of.packages) {
  library(package, character.only = TRUE)
}

VariablesPath <- paste0(getwd(), "/data/covariates1km") 
SpeciesFilePath <- paste0(getwd(), "/data/species/shrubs/global_regional") 

SpeciesDataGlobal <- readRDS(paste0(SpeciesFilePath,"/global_and_regional_species.rds"))
SpeciesName.list <- unique(SpeciesDataGlobal$species) 
head(SpeciesName.list)

sabina.multiple.species <- function(SpeciesName) {
  tryCatch({ 
    library("terra")
    library("biomod2")
    library("covsel")
    library("fs")
    library("ecospat")
    library("tidyverse")
    library("gam")
    library("plyr")
    library("sabinaNSDM")
   
    #Spp data
    SpeciesDataGlobal <- readRDS(paste0(SpeciesFilePath,"/global_and_regional_species.rds"))
    SpeciesDataGlobal <- as.data.frame(SpeciesDataGlobal)
    spp.data.global <- SpeciesDataGlobal %>%
      filter(species == SpeciesName)
    spp.data.global <- spp.data.global[,2:3]
    names(spp.data.global) <- c("x","y")
    head(spp.data.global)

    #Vars
    expl.var.global <- rast(paste0(VariablesPath,"/global/Current.tif"))

    #New.env
    new.env.files <- list.files(paste0(VariablesPath, "/regional"), pattern = "\\.tif$", full.names = TRUE)
    rm_current <- grep("Current.tif", new.env.files)
    new.env.files <- new.env.files[-rm_current]
    
    new.env <- list()
    
    for(file in new.env.files) {
      rast <- terra::rast(file)  
      new.env[[basename(file)]] <- rast
    }

    nsdm_input <-NSDM.InputData(SpeciesName=SpeciesName,
                                   spp.data.global=spp.data.global, 
                                   spp.data.regional=spp.data.global, 
                                   expl.var.global=expl.var.global, 
                                   expl.var.regional=expl.var.global,
                                   new.env=new.env,
                                   new.env.names=NULL, 
                                   Background.Global=NULL, 
                                   Background.Regional=NULL)
  
    nsdm_finput <- NSDM.FormattingData(nsdm_input,
                                        nPoints = 10000,
                                        Min.Dist.Global = "resolution",
                                        Min.Dist.Regional = "resolution",
                                        Background.method="random",
                                        save.output = TRUE)
    
    nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
                                              maxncov.Global = 7,
                                              maxncov.Regional = 7,
                                              corcut = 0.7,
                                              algorithms = c("glm","gam","rf"),
                                              ClimaticVariablesBands = NULL,
                                              save.output = TRUE)
    
    nsdm_global <- NSDM.Global(
      nsdm_selvars,
      algorithms = c("GBM", "RF", "GLM"),
      CV.nb.rep = 10,
      CV.perc = 0.8,
      save.output = TRUE,
      metric.select.thresh = 0.8,
      rm.biomod.folder = TRUE
    )
    
    
  }, error = function(err) {
    message("Error processing species ", SpeciesName, ": ", conditionMessage(err))
    return(list(error = err))
  })
}

plan(multisession, workers = 1)
sabina.Results <- future_map(SpeciesName.list, safely(sabina.multiple.species), .progress = TRUE)
 

print(sabina.Results[[1]]) # to see possible errors

sabina.Results <- transpose(sabina.Results) # list-of-lists "inside-out" 

