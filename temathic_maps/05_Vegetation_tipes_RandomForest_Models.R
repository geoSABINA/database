

######################################################################################
## randomForest models for vegetation types (TCE) as a function of climate variables
######################################################################################


# Vegetation types (TCE) data source:  Forest Map of Spain 1:200,000 (Ruiz de la Torre J., 1990)
# Bioclimatic variables: BIO1 - BIO19 Earth's Surface Areas dataset (Karger et al, 2017)
# Extent: Peninsular Spain 


######################
#Starting parameters
######################

library(randomForest)
library(covTestR)
library(terra)
library(raster)


setwd("C:/WisdomDir/")
Clim_data<-read.csv2("Clim.csv")
Clim_data$TCE <-as.factor(Clim_data$TCE)


### Variables

#BIO1 = Annual Mean Temperature                 T
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) Trange
#BIO3 = Isothermality (BIO2/BIO7) (×100)            termic
#BIO4 = Temperature Seasonality (standard deviation ×100) stdT
#BIO5 = Max Temperature of Warmest Month              TmaxC
#BIO6 = Min Temperature of Coldest Month              TminF
#BIO7 = Temperature Annual Range (BIO5-BIO6)          TmaxRange
#BIO8 = Mean Temperature of Wettest Quarter           TmMesHum
#BIO9 = Mean Temperature of Driest Quarter            TmMesSec
#BIO10 = Mean Temperature of Warmest Quarter          TmC
#BIO11 = Mean Temperature of Coldest Quarter          TmF
#BIO12 = Annual Precipitation                         P
#BIO13 = Precipitation of Wettest Month               PMesHum
#BIO14 = Precipitation of Driest Month                PMesSec
#BIO15 = Precipitation Seasonality (Coefficient of Variation) CvarP
#BIO16 = Precipitation of Wettest Quarter             P4Sec
#BIO17 = Precipitation of Driest Quarter              P4Hum
#BIO18 = Precipitation of Warmest Quarter             P4C
#BIO19 = Precipitation of Coldest Quarter             P4F

# TCE = Vegetation type according to Forest Map of Spain 1:200.000
  # "A": high-mountain vegetation; "C": deciduous; "T": mountain conifers;
  # "S": Subsclerophylluous; "E": sclerophylluous; "H": hyperxerophilous.
  # Sample size : 7136 


#############################
# Train-Test Validation Split
#############################


set.seed(222)
ind <- sample(2, nrow(CLim_data), replace = TRUE, prob = c(0.7, 0.3))
trainCl <- CLim_data[ind==1,]
testCl <- CLim_data[ind==2,]


#####################
# randomForest model 
#####################

#### 1.Fitting de basic model

rfCl<-randomForest(TCE~., data=trainCl, ntree=500)

#### 2.Model parameters tuning

mtry <- tuneRF(trainCl[-1],trainCl$TCE, ntreeTry=500,
               stepFactor=1.5,improve=0.1, trace=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]

#### 3.Optimized model

rf2Cl<-randomForest(TCE~., data=trainCl, mtry=best.m, importance=TRUE, ntree=500)

### 4.Importance of variables 

rf2Cl$importance
varImpPlot(rf2Cl)

#### 5.Model predictions 

pred1=predict(rf2Cl,type = "prob")

predNewCl<- predict(rf2Cl, newdata = testCl, type = "prob")
mCl<-colnames(predNewCl)[apply(predNewCl, 1, which.max)]


#############################################
# Mapping TCE predictions for peninsular Spain 
#############################################

setwd("C:/WisdomDir/ClimaticRaster")
rlist <- list.files("C:/WisdomDir/ClimaticRaster", "rds$")
x<-lapply(rlist, readRDS)
vars <- stack(x)
names(vars)<-substring(names(vars), first = 40) 

rb<-c("#FF00FF","#CAE1FF","#FFFF00","#FF7F24","#B3EE3A","#5F9EA0")
Rrf2Cl<-predict(vars,rf2Cl, type='response', progress='window')
plot(Rrf2Cl, legend = FALSE, col =rb, main="Present" )

legend("bottomright",legend= c("High-mountain", "Deciduous", "Sclerophylluous", "Hyperxerophilous", 
       "Subsclerophylluous", "Mountain conifers"),
       fill=c("#FF00FF","#CAE1FF","#FFFF00","#FF7F24","#B3EE3A","#5F9EA0"),cex =0.3,bty = "n" )
       
saveRDS(Rrf2Cl,file="Present.rds")



