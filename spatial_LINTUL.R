
###########################################################################################
### Data source and run LINTUL
###########################################################################################

source('LINTUL_functions.R')

## define the country and the year LINTUL is going to run 
Country_year <- rbind(data.frame(country="Ghana", Q3Year = 2016, Q2Year = 2012, Q1Year = 2019))

## define planting and harvest dates as the count of number of days in the year.
## for cassava we assume it can have weekly planting throughout the year and 
## can stay on the field up to 15 months = 455 growing days. Weekly harvest between 8 
## and 15 months is also assumed.    

Planting_Harvest_GH <- data.frame(startDate=seq(1, 365, 7), 
                                    endDate = (seq(1, 365, 7) + 454), 
                                    weekNr = seq(1:53),
                                    months = c(rep("Jan", 5), rep("Feb", 4), rep("Mar", 5), 
                                               rep("Apr", 4), rep("May", 4), rep("Jun", 5), 
                                               rep("Jul", 4), rep("Aug", 5), rep("Sep", 4), 
                                               rep("Oct",4), rep("Nov",5), rep("Dec", 4)))


##GH  2016 read soil, rainfall and solar data 
Solar_GH_16 <- readRDS("Solar_GH_Data.RDS")
RF_GH_16 <- readRDS("Rainfall_GH_Data.RDS")
Soil_GH <- readRDS("Soil_GH_Data.RDS")


Solar_GH_16$location <- paste(Solar_GH_16$lat, Solar_GH_16$lon, sep="_")
RF_GH_16$location <- paste(RF_GH_16$lat, RF_GH_16$long, sep="_")
Soil_GH$location <- paste(Soil_GH$lat, Soil_GH$lon, sep="_")


GH_WLY <- runLINTUL_server(year=2016,  country = "Ghana", 
                 SoilData = Soil_GH, RainData = RF_GH_16, Solar = Solar_GH_16,
                 Planting_Harvest_days = Planting_Harvest_GH)


str(GH_WLY)

