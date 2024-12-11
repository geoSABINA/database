####################################################
## Calculate corridors effective distance (Least cost path)
####################################################

# Nodes: forest patches (from CORINE) for a species with a 5km median dispersal distance (minimum area 51.02 ha, min distance 370m)
# Resistance from CORINE, 200m resolution
# Extent: Peninsular Spain divided by vegetation types (TCE)

######################
#Starting parameters
######################

library(raster) 
library(gdistance) 
library(sf)

crs_utm<-CRS("+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

scenarios<-c("Present","IPSL_CM6A_LR_2070_SSP126","IPSL_CM6A_LR_2070_SSP585","MRI_ESM2_0_2070_SSP126","MRI_ESM2_0_2070_SSP585")
TCEs<-data.frame(1:6,c("A","C","E","H","S","T"),c("High-mountain vegetation",  "Deciduous","Sclerophylluous", 
                                                  "Hiperxerophilous", "Subsclerophylluous","Mountain conifers")) #vegetation types
colnames(TCEs)<-c("number","code","name")
disp_dist<-c(1,3,5,10,30) #median dispersal distance in km
resolution<-200


#########
#RUN LCP - EFFECTIVE DISTANCE with gdistance
#########

res<-readRDS(paste0("R/data/resistance/resist_",resolution,"_utm.rds"))
res[which(values(res)<1)]<-NA
con<-1/res
time_table <- data.frame(matrix(ncol=4,nrow=0)) 
colnames(time_table)<-c("Process","Scenario","TCE","Time (min)")

if (!file.exists(paste0("R/data/resistance/transition_",resolution,".rds"))) {
  tr <- transition(con, transitionFunction=mean, directions=8) 
  tr <- geoCorrection(tr, type="c",multpl=FALSE, scl=FALSE) 
  saveRDS(tr,"R/data/resistance/transition_200.rds")
} else {tr<-readRDS(paste0("R/data/resistance/transition_",resolution,".rds"))}
crs(tr)<-crs_utm
gc()

for (i in 1:length(scenarios)) {
  scenario_i<-scenarios[i]
  scenario_pepa_i<-scenarios_pepa[i]
  print(scenario_pepa_i) 

  for (t in TCEs$number) {
    t1<-Sys.time()
    print(TCEs$name[t])
    print(paste("Starting time:",t1))
    patches<-readRDS(paste0("R/data/patches/patches_",t,"_",scenario_i,".rds"))   
    print(paste("Ner of patches", nrow(patches)))
    xy_SP<-SpatialPoints(coords = cbind(patches$x,patches$y),proj4string = crs_utm)
    xy_SP2<-SpatialPoints(coords = cbind(patches$x,patches$y),proj4string = crs_laea)
    
    saveRDS(xy_SP,paste0("R/data/nodes/nodes_",t,"_",scenario_i,".rds"))
    
    dist.mat<-gdistance::costDistance(tr,xy_SP,xy_SP)
    print(dim(dist.mat))
    
    dist.table<-as.data.frame(dist.mat)
    colnames(dist.table)<-rownames(dist.table)
    saveRDS(dist.table,paste0("R/results/LCP/eff_distance/eff_distance_",t,"_",scenario_i,".rds"))
    
    t2<-Sys.time()
    time<-difftime(t2,t1,units='mins')
    print(time)
    time_table<-rbind(time_table,c("LC effective distance",scenario_i,TCEs$name[t],time))
    saveRDS(time_table,"R/results/LCP/eff_distance/processing_time.rds")
    gc()
  }}  


