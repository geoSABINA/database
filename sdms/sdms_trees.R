
# Setting the working directory
setwd("/Users/RubenGMateo1/MacHD2/sabinaNSDM/Web/R")
setwd("/home/r_mateo/SABINA/R")

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

VariablesPath <- paste0(getwd(), "/data/covariates") 
SpeciesFilePath <- paste0(getwd(), "/data/species/trees") 

SpeciesName.list <- list.files(paste0(SpeciesFilePath,"/global"))
SpeciesName.list <-  gsub(".csv","", SpeciesName.list, fixed=TRUE) 
# Uncomment the next line for testing parallelization with five species
# SpeciesName.list <- SpeciesName.list[40:60]
SpeciesName.list <- SpeciesName.list[1:40]
head(SpeciesName.list)

# Uncomment the next line for testing parallelization with a single file
# SpeciesName <- SpeciesName.list[51]
# SpeciesName <- SpeciesName.list[1]
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
    spp.data.global <- read.csv(paste0(SpeciesFilePath,"/global/", SpeciesName, ".csv"))
    spp.data.regional <- read.csv(paste0(SpeciesFilePath,"/regional/", SpeciesName, ".csv"))
    names(spp.data.global) <- c("x","y")
    names(spp.data.regional) <- c("x","y")
    
    #Vars
    expl.var.global <- rast(paste0(VariablesPath,"/global/Current.tif"))
    expl.var.regional <- rast(paste0(VariablesPath,"/regional/Current.tif"))
    
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
                                   spp.data.regional=spp.data.regional, 
                                   expl.var.global=expl.var.global, 
                                   expl.var.regional=expl.var.regional,
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
                                              maxncov.Global = 5,
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
      #CustomModelOptions = opt.b,  
      save.output = TRUE,
      metric.select.thresh = 0.8,
      rm.biomod.folder = TRUE
    )
    
    nsdm_regional <- NSDM.Regional(
      nsdm_selvars,
      algorithms = c("GBM", "RF", "GLM"),
      CV.nb.rep = 10,
      CV.perc = 0.8,
      #CustomModelOptions = NULL,  
      metric.select.thresh = 0.7,
      save.output = TRUE,
      rm.biomod.folder = TRUE
    )
    
    nsdm_covariate <- NSDM.Covariate(
      nsdm_global,
      algorithms = c("GBM", "RF", "GLM"),
      CV.nb.rep = 10,
      CV.perc = 0.8,
      rm.corr = TRUE,
      #CustomModelOptions = NULL,
      metric.select.thresh = 0.7,
      rm.biomod.folder = TRUE,
      save.output = TRUE
    )
    
    nsdm_multiply <- NSDM.Multiply(
      nsdm_global,
      nsdm_regional,
      method = c("Geometric"),
      rescale = FALSE,
      save.output = TRUE
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

