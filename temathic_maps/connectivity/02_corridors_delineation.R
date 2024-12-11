ve###############################
## Delineate corridors with Least Cost Path
#############################

# Nodes: forest patches (from CORINE) for a species with a 5km median dispersal distance (minimum area 51.02 ha, min distance 370m)
# Resistance from CORINE, 200m resolution
# Extent: Peninsular Spain divided by vegetation types (TCE)


########################
## Strating paramenters
########################

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

max_distance<-300000 #300km


#########
#RUN LCP - CALCULATE PATHS
#########

res<-readRDS(paste0("R/data/resistance/resist_",resolution,"_utm.rds"))
res[which(values(res)<1)]<-NA
con<-1/res

if (!file.exists(paste0("R/data/resistance/transition_",resolution,".rds"))) {
  tr <- transition(con, transitionFunction=mean, directions=8) 
  tr <- geoCorrection(tr, multpl=FALSE, scl=FALSE) #type=c for LCP, for random and RSP needs first c and then r. The SCL scalin should be on FALSE to obtain the real resistance values
  saveRDS(tr,"R/data/resistance/transition_200.rds")
} else {tr<-readRDS(paste0("R/data/resistance/transition_",resolution,".rds"))}
crs(tr)<-crs_utm

gc()

time_table <- data.frame(matrix(ncol=4,nrow=0)) 
colnames(time_table)<-c("Process","Scenario","TCE","Time (min)")
for (i in 1:length(scenarios)) {
  scenario_i<-scenarios[i]
  scenario_pepa_i<-scenarios_pepa[i]
  print(scenario_pepa_i) 
  
  for (t in TCEs$number) {
    t1<-Sys.time()
    print(TCEs$name[t])
    
    patches<-readRDS(paste0("R/data/patches/patches_",t,"_",scenario_i,".rds"))   
    print(paste("Ner of patches", nrow(patches)))
    xy_SP<-readRDS(paste0("R/data/nodes/nodes_",t,"_",scenario_i,".rds"))
    rm(patches)
    
    LCP_merged<-list()
    saveRDS(LCP_merged,paste0("R/data/LCP_merged_",t,"_",scenario_i,"_",0))
    
    kk<-list.files("R/data",pattern = "LCP_merged_")
    kk<-kk[grep(paste0(t,"_",scenario_i,"_"),kk)]
    if (length(kk)!=0) {
      library(stringr)
      kk <- as.numeric(str_extract(kk, "\\d+$"))
      first<-max(kk)+1
    } else {first<-1}
    
    #LCP only for patches closer than max dist of 300km
    for (x in first:(length(xy_SP)-1)) {
      print(paste(x,"of",(length(xy_SP)-1)))
      origin_i<-xy_SP[x,]
      goalsID <- (x + 1):length(xy_SP)
      goal_i<-xy_SP[goalsID,]
      
      #remove the goal points that are too far away, more than maximum distance
      distances <- spDistsN1(as.matrix(coordinates(goal_i)), as.matrix(coordinates(origin_i)), longlat = FALSE)
      close_goals <- which(distances <= max_distance)
      
      if (length(close_goals) > 0) {
        goal_i<- goal_i[close_goals,]
        close_goals_ID <- goalsID[close_goals]
        
        table <- data.frame(Source = rep(x, length(goal_i)), Destination = close_goals_ID)
        LCP_i<-gdistance::shortestPath(x=tr,origin=origin_i,goal=goal_i,output= "SpatialLines")
        LCPdf_i <- SpatialLinesDataFrame(sl = LCP_i, data = table, match.ID = FALSE)
        rm(LCP_i)
        
        if (length(xy_SP) > 1700) { #for big networks chunk the work
          if (x!=1) {
            LCP_merged<-readRDS(paste0("R/data/LCP_merged_",t,"_",scenario_i,"_",x-1))
          }
          LCP_merged[[x]]<-LCPdf_i
          rm(LCPdf_i)
          gc()
          saveRDS(LCP_merged,paste0("R/data/LCP_merged_",t,"_",scenario_i,"_",x))
          file.remove(paste0("R/data/LCP_merged_",t,"_",scenario_i,"_",x-1))
          rm(LCP_merged)
        } else {
          LCP_merged[[x]]<-LCPdf_i
          rm(LCPdf_i)
        }
      }}
    
    #Join results of all source nodes
    if (length(xy_SP) > 1700) {LCP_merged<-readRDS(paste0("R/data/LCP_merged_",t,"_",scenario_i,"_",(length(xy_SP)-1)))}
    LCP <- do.call(rbind, LCP_merged) #convert the list of sp objects to an sp object
    LCP@proj4string<-crs_utm
    
    saveRDS(LCP,paste0("R/results/LCP/paths/LCP_",t,"_",scenario_i,".rds"))
    LCP_sf<-st_as_sf(LCP,crs=crs_utm)
    LCP_sf <- st_simplify(LCP_sf, dTolerance = 0.001)
    saveRDS(LCP_sf,paste0("R/results/LCP/paths/LCP_sf_",t,"_",scenario_i,".rds"))
    print("saved")
    
    t2<-Sys.time()
    time<-difftime(t2,t1,units='mins') 
    print(time)
    time_table<-rbind(time_table,c("LCP",scenario_i,TCEs$name[t],time))
    saveRDS(time_table,"R/results/LCP/paths/processing_time.rds")
    gc()
  }}  




