####################################################
## Calculate patches importance (dPC) and overall PC and ECA
####################################################


# Nodes: forest patches (from CORINE) for a species with a 5km median dispersal distance (minimum area 51.02 ha, min distance 370m)
# Resistance from CORINE, 200m resolution
# Extent: Peninsular Spain divided by vegetation types (TCE)
# Connectivity metrics: conefor habitat availability metrics (PC,ECA,dPC) 


######################
#Starting parameters
######################


library(Makurhini)
library(sf)
library(parallel)

crs_utm<-CRS("+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

scenarios<-c("Present","IPSL_CM6A_LR_2070_SSP126","IPSL_CM6A_LR_2070_SSP585","MRI_ESM2_0_2070_SSP126","MRI_ESM2_0_2070_SSP585")
TCEs<-data.frame(1:6,c("A","C","E","H","S","T"),c("High-mountain vegetation",  "Deciduous","Sclerophylluous", 
                                                  "Hiperxerophilous", "Subsclerophylluous","Mountain conifers")) #vegetation types
colnames(TCEs)<-c("number","code","name")
resolution<-200

# convert dispersal distance to effective dispersal distances
res<-readRDS(paste0("R/data/resistance/resist_",resolution,"_utm.rds"))
res[which(raster::values(res)<1)]<-NA
mean_resistance<-mean(na.omit(raster::values(res)))
disp_dist<-c(1,3,5,10,30) #median dispersal distance in km
distance_thres<-disp_dist*1000*mean_resistance

##########################
# overall PC #with Makurhini functions
#########################

time_table <- data.frame(matrix(ncol=4,nrow=0)) 
colnames(time_table)<-c("Process","Scenario","TCE","Time (min)")
n_cores<- detectCores()-1

for (i in 1:length(scenarios)) {
  scenario_i<-scenarios[i]
  scenario_pepa_i<-scenarios_pepa[i]
  print(scenario_pepa_i) 
  
  for (t in TCEs$number) {
    t1<-Sys.time()
    print(TCEs$name[t])
    print(paste("Starting time:",t1))

#Charge input nodes and distances for TCE t and scenario i:
    if (file.exists(paste0("conefor/inputs/nodes_",t,"_",scenario_i,".txt")) ) {
      nodes<-read.table(paste0("conefor/inputs/nodes_",t,"_",scenario_i,".txt"))
      colnames(nodes)<-c("id","attribute")
      
      distances<-readRDS(paste0("R/results/LCP/eff_distance/eff_distance_",t,"_",scenario_i,".rds"))
      distance_matrix<-as.matrix(distances)
      colnames(distance_matrix)<-nodes$id
      row.names(distance_matrix)<-nodes$id
      
      for (d in 1:length(disp_dist)) {
        
        distance_thres_d<-distance_thres[d]
        disp_dist_j<-disp_dist[d]
        print(paste(disp_dist_j, "km"))
        overall<-cbind(0,0)
        
        if (!file.exists(paste0("conefor/outputs/only_overall_",TCE_i,"_",scenario_i,"_",disp_dist_j,"km.rds"))){
          
          PCnodes <- MK_dPCIIC(nodes = nodes, area_unit = "ha",
                               restoration=NULL,
                               distance = distance_matrix,
                               metric = "PC", probability = 0.5,
                               overall = TRUE, onlyoverall = TRUE,
                               distance_thresholds = distance_thres_d,
                               parallel=n_cores) #parallel NULL if no many nodes, or number of cores to paralelize
          
          saveRDS(PCnodes,paste0("conefor/outputs/only_overall_",t,"_",scenario_i,"_",disp_dist_j,"km.rds"))
    }    
}
}
  }
  }
