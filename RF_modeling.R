


####################################################################################################
####################################################################################################
## 1. Train the random forest model: modeling soil NPK as a repsonse of soil properties and control/farmers' current yield using random forest model
####################################################################################################
## variables for RF

RFvars <- c("country", "aluminium_extractable_0_20","aluminium_extractable_20_50", "bedrock_depth_0_200",
  "bulk_density_0_20", "bulk_density_20_50", "calcium_extractable_0_20", "calcium_extractable_20_50", 
  "carbon_organic_0_20", "carbon_organic_20_50", "carbon_total_0_20", "carbon_total_20_50",
  "cation_exchange_capacity_0_20", "cation_exchange_capacity_20_50", "clay_content_0_20", "clay_content_20_50",
  "iron_extractable_0_20", "iron_extractable_20_50", "magnesium_extractable_0_20", "magnesium_extractable_20_50",
  "nitrogen_total_0_20", "nitrogen_total_20_50", "ph_0_20", "ph_20_50", "phosphorous_extractable_0_20", 
  "phosphorous_extractable_20_50" ,"potassium_extractable_0_20", "potassium_extractable_20_50","sand_content_0_20",
  "sand_content_20_50", "silt_content_0_20", "silt_content_20_50", "stone_content_0_20", "stone_content_20_50", 
  "sulphur_extractable_0_20", "sulphur_extractable_20_50", "texture_class_0_20", "texture_class_20_50","slope_angle", 
  "zinc_extractable_0_20", "zinc_extractable_20_50", "soilOM_0_20", "soilOM_20_50", "percentSOM_0_20", 
  "percentSOM_20_50", "wp_0_20", "wp_20_50", "FC_0_20", "FC_20_50", "sws_0_20", "sws_20_50", "soilN","soilP","soilK", "CON")

### read data
randomF_Data <- readRDS("randomF_Data.RDS")


### split data to training and validation set 
splitdataACI <- function(probs=c(0.7,0.3), soil_BLUPSData){
  set.seed(444)
  trianData <- NULL
  testDatA <- NULL
  for(country in unique(soil_BLUPSData$country)){
    Country_data <- droplevels(soil_BLUPSData[soil_BLUPSData$country == country,])
    indexCountry <- sample(2, nrow(Country_data), replace=TRUE, prob=probs)## where control yield is used as a covariate
    trainDataCountry <- Country_data[indexCountry == 1, ]
    testDataCountry <- Country_data[indexCountry == 2, ]
    trianData <- rbind(trianData, trainDataCountry) 
    testDatA <- rbind(testDatA, testDataCountry)
  }
  return(list(trianData, testDatA))
}

traintestData_BLUP <- splitdataACI(soil_BLUPSData = randomF_Data)
trainData_BLUP <- traintestData_BLUP[[1]]
testData_BLUP <- traintestData_BLUP[[2]]

Ndata_Train_BLUP <- subset(trainData_BLUP, select=-c(soilP, soilK))
Pdata_Train_BLUP <- subset(trainData_BLUP, select=-c(soilN, soilK))
Kdata_Train_BLUP <- subset(trainData_BLUP, select=-c(soilN, soilP))

Ndata_Valid_BLUP <- subset(testData_BLUP, select=-c(soilP, soilK))
Pdata_Valid_BLUP <- subset(testData_BLUP, select=-c(soilN, soilK))
Kdata_Valid_BLUP <- subset(testData_BLUP, select=-c(soilN, soilP))

### run RF
require(randomForest)
require(gtools)
## Coustome control parameter 
#custom <- trainControl(method="repeatedcv", number=10, repeats=5, verboseIter=TRUE)
custom <- trainControl(method="oob", number=10)

##soil N
set.seed(773)
Ndata_Train_BLUP <- subset(Ndata_Train_BLUP, select=-c(CON, bedrock_depth_0_200))
RF_N_B <- randomForest(soilN ~ ., Ndata_Train_BLUP, importance=TRUE, ntree=1000)
Ndata_Valid_BLUP$predN <- predict(RF_N_B, Ndata_Valid_BLUP)
sst <- sum((Ndata_Valid_BLUP$soilN - mean(Ndata_Valid_BLUP$soilN))^2)
sse <- sum((Ndata_Valid_BLUP$soilN - Ndata_Valid_BLUP$predN)^2)
rsq <- 1 - sse / sst
rsq

varImpPlot(RF_N_B)
varimportance <- data.frame(RF_N_B$importance)
varimportance$vars <- rownames(varimportance)
varimportance <- varimportance[order(varimportance$X.IncMSE, decreasing = TRUE), ]
rownames(varimportance) <- NULL
varimportance$Importance <- c(1:nrow(varimportance))

ggn <-ggplot(Ndata_Valid_BLUP, aes(soilN, predN, col=country)) +
        geom_point()+
        geom_abline(slope = 1, intercept = 0) +
        xlim(0,150) + ylim(0,150)+
        xlab("soil N from QUEFTS optimization")+
        ylab("perdicted soil N from GIS data")+
        theme_bw()

###soil P 
set.seed(773) 
Pdata_Train_BLUP <- subset(Pdata_Train_BLUP, select=-c(CON, bedrock_depth_0_200))
RF_P1_B <- randomForest(soilP ~ ., Pdata_Train_BLUP , importance=TRUE, ntree=1000)
Pdata_Valid_BLUP$predP <- predict(RF_P1_B, Pdata_Valid_BLUP)
sst <- sum((Pdata_Valid_BLUP$soilP - mean(Pdata_Valid_BLUP$soilP))^2)
sse <- sum((Pdata_Valid_BLUP$soilP - Pdata_Valid_BLUP$predP)^2)
rsq <- 1 - sse / sst
rsq

ggp <- ggplot(Pdata_Valid_BLUP, aes(soilP, predP, col=country)) +
          geom_point()+
          geom_abline(slope = 1, intercept = 0) +
          xlim(0,100) + ylim(0,100)+
          xlab("soil P from QUEFTS optimization")+
          ylab("perdicted soil P from GIS data")+
          theme_bw()




## soil K 
set.seed(773)
Kdata_Train_BLUP <- subset(Kdata_Train_BLUP, select=-c(CONclass, bedrock_depth_0_200))
RF_K_B <- randomForest(soilK ~ ., Kdata_Train_BLUP , importance=TRUE, ntree=1000)
Kdata_Valid_BLUP$predK <- predict(RF_K_B, Kdata_Valid_BLUP)
sst <- sum((Kdata_Valid_BLUP$soilK - mean(Kdata_Valid_BLUP$soilK))^2)
sse <- sum((Kdata_Valid_BLUP$soilK - Kdata_Valid_BLUP$predK)^2)
rsq <- 1 - sse / sst
rsq

ggk <- ggplot(Kdata_Valid_BLUP, aes(soilK, predK, col=country)) +
          geom_point()+
          geom_abline(slope = 1, intercept = 0) +
          xlim(0,200) + ylim(0,200)+
          xlab("soil K from QUEFTS optimization")+
          ylab("perdicted soil K from GIS data")+
          theme_bw()





####################################################################################################
####################################################################################################
## 2. get the soil data for the whole of Burundi and get soilNPK
####################################################################################################               
################################################################################################################################
################################################################################################################################
require(lpSolve)
require(tidyr)
require(randomForest)
require(gtools)
require(caret)

## read Burundi soils data
Bu_soils <- readRDS("Burundi_soil.RDS")

### source teh QUEFTS function in reverse QUEFTS folder
source("QUEFTS_Functions.R")
str(Bu_soils)

Bu_soils$country <- "Bu"
Bu_soils$soilN <- NA
Bu_soils$soilP <- NA
Bu_soils$soilK <- NA


## soil NPK will be estimated for the whole area assuming the five FCY levels
Bu_soils$index <- c(1:nrow(Bu_soils))
Bu_soils_FCY1 <- Bu_soils
Bu_soils_FCY2 <- Bu_soils
Bu_soils_FCY3 <- Bu_soils
Bu_soils_FCY4 <- Bu_soils
Bu_soils_FCY5 <- Bu_soils

Bu_soils_FCY1$CONclass <- "class1"
Bu_soils_FCY2$CONclass <- "class2"
Bu_soils_FCY3$CONclass <- "class3"
Bu_soils_FCY4$CONclass <- "class4"
Bu_soils_FCY5$CONclass <- "class5"


### ACAI, RW and Bu data as a training set
trainData <- randomF_Data
trainData$CONclass <- as.factor(ifelse(trainData$CON < 7.5, "class1",
                                      ifelse(trainData$CON >= 7.5 & trainData$CON < 15, "class2",
                                             ifelse(trainData$CON >= 15 & trainData$CON < 22.5, "class3",
                                                    ifelse(trainData$CON >= 22.5 & trainData$CON < 33, "class4", "class5")))))


trainData <- trainData[!is.na(trainData$slope_angle), ]
trainData <- subset(trainData, select = -c(CON, bedrock_depth_0_200,stone_content_0_20, stone_content_20_50 ))

######### getting soil NPK RW ############

Bu_FCY1_soilNPK <- get_SoilNPK_AOI(trainData=trainData, AIOdata=Bu_soils_FCY1, FCY = "FCY1")

Bu_FCY2_soilNPK <- get_SoilNPK_AOI(trainData=trainData, AIOdata=Bu_soils_FCY2, FCY = "FCY2")

Bu_FCY3_soilNPK <- get_SoilNPK_AOI(trainData=trainData, AIOdata=Bu_soils_FCY3, FCY = "FCY3")

Bu_FCY4_soilNPK <- get_SoilNPK_AOI(trainData=trainData, AIOdata=Bu_soils_FCY4, FCY = "FCY4")

Bu_FCY5_soilNPK <- get_SoilNPK_AOI(trainData=trainData, AIOdata=Bu_soils_FCY5, FCY = "FCY5")

BY_soilNPK <- rbind(Bu_FCY5_soilNPK, Bu_FCY4_soilNPK, Bu_FCY3_soilNPK, Bu_FCY2_soilNPK, Bu_FCY1_soilNPK)
BY_soilNPK$location <- paste(BY_soilNPK$lat, BY_soilNPK$lon, sep="_")


