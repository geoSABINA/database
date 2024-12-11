####################################################
## Calculate corridors importance (dPC) and overall PC and ECA
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
library(raster)
library(future)
library(furrr)

crs_utm<-CRS("+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

scenarios<-c("Present","IPSL_CM6A_LR_2070_SSP126","IPSL_CM6A_LR_2070_SSP585","MRI_ESM2_0_2070_SSP126","MRI_ESM2_0_2070_SSP585")
TCEs<-data.frame(1:6,c("A","C","E","H","S","T"),c("High-mountain vegetation",  "Deciduous","Sclerophylluous", 
                                                  "Hiperxerophilous", "Subsclerophylluous","Mountain conifers")) #vegetation types
colnames(TCEs)<-c("number","code","name")

# convert dispersal distance to effective dispersal distances
resolution<-200
res<-readRDS(paste0("R/data/resistance/resist_",resolution,"_utm.rds"))
res[which(values(res)<1)]<-NA
mean_resistance<-mean(na.omit(values(res)))
disp_dist<-c(1,3,5,10,30) #median dispersal distance in km
distance_thres<-disp_dist*1000*mean_resistance

##########################
# LINKS dPC, and overal PC #with modified Makurhini functions
#########################

time_table <- data.frame(matrix(ncol=4,nrow=0)) 
colnames(time_table)<-c("Process","Scenario","TCE","Time (min)")
n_cores<- detectCores()-1

array_table<-data.frame(array=1:(length(scenarios)*length(TCEs$number)*length(distance_thres)),scenario=rep(1:length(scenarios),each=(length(TCEs$number)*length(distance_thres))),TCE= rep(rep(TCEs$number,each=length(distance_thres)),length(scenarios)),distance=rep(distance_thres,(length(scenarios)*length(TCEs$number))))

for (arrayID in 1:nrow(array_table)) {
  
  i<-array_table$scenario[arrayID]
  t<-array_table$TCE[arrayID]  
  d<-array_table$distance[arrayID] 

  scenario_i<-scenarios[i]
  scenario_pepa_i<-scenarios_pepa[i]
  print(scenario_pepa_i) 
  
    t1<-Sys.time()
    print(TCEs$name[t])
    print(paste("Starting time:",t1))
    
    #Charge input nodes and distances for TCE t and scenario i:
    nodes<-read.table(paste0("conefor/inputs/nodes_",t,"_",scenario_i,".txt"))
    colnames(nodes)<-c("id","attribute")
    
    print(paste("Effectifve distance:",d))
    dist<-d/1000/mean_resistance
    print(paste("Distance:",dist, "km"))
    
    #Charge input effective distances for TCE t and scenario i:
    distances<-readRDS(paste0("R/results/LCP/eff_distance/eff_distance_",t,"_",scenario_i,".rds"))
    distance_matrix<-as.matrix(distances)
    colnames(distance_matrix)<-nodes$id
    row.names(distance_matrix)<-nodes$id
    
    #remove too big distances
    kk<--log(0.5)/d
    prob_min<-0.001 
    dist_max<--log(prob_min)/kk
    dist_max

      PCconnections <- dPC_connections(nodes = nodes, area_unit = "ha",
                                       restoration=NULL, 
                                       dist_max=dist_max, #give them 0 dPC
                                       distance = distance_matrix,
                                       metric = "PC", probability = 0.5,
                                       overall = TRUE, onlyoverall = FALSE,
                                       distance_thresholds = d,  parallel=n_cores,
                                       delta_connections=TRUE)
      saveRDS(PCconnections,paste0("conefor/outputs/dPClinks_",t,"_",scenario_i,"_",dist,"km.rds"))
    
    
    t2<-Sys.time()
    time<-difftime(t2,t1,units='mins')
    print(paste("Time for calculating overal PC, and links dPC:",time))
    time_table<-rbind(time_table,c("dPC links",scenario_i,dist,TCEs$name[t],time,"UNIL cluster","Makurhini"))
    write.csv(time_table,paste0("conefor/outputs/processing_time_dPClinks_",t,"_",i,".csv"))
    
}
