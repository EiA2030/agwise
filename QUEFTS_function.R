
nutrient_uptake <- function(S1=NA, S2=NA, d1=NA, a1=NA, d2=NA, a2=NA, r1=NA, r2=NA) {
  # N, P and K uptakes based on QUEFTS
  uptakeX_givenY = S1 - 0.25 * (S1 - r1 - (S2 - r2) * (a2 / d1))^2 / ((S2 - r2) * (d2 / a1 - a2 / d1))
  uptakeX_givenY[S1 > r1 + ((S2 - r2) * (2 * d2 / a1 - a2 / d1))] <- (r1 + (S2 - r2) * (d2 / a1))[S1 > r1 + ((S2 - r2) * (2 * d2 / a1 - a2 / d1))]
  uptakeX_givenY[S1 < r1 + ((S2 - r2) * a2 / d1)] <- S1[S1 < r1 + ((S2 - r2) * a2 / d1)]
  # Nutrient uptake given availability of other nutrient
  
  return(uptakeX_givenY)
  
}


######################################################################
######################################################################


yield_nutrients_combined <- function(U1=NA, d1=NA, a1=NA, Y2A=NA, Y2D=NA, Y3D=NA, r1=NA){
  # Determine which nutrient limited to force yield being the lowest.
  YxD = pmin(Y2D,Y3D)
  # If the uptake of one of the nutrients, and therefore the yield associated with that nutrient, is zero the overall yield is also zero.
  Y12 <- Y2A + (2 * (YxD - Y2A) * (U1 - r1 - Y2A / d1)) / (YxD / a1 - Y2A / d1) -
    (YxD - Y2A) * (U1 - r1 - Y2A / d1)^2 / (YxD / a1 - Y2A / d1)^2
  Y12[U1 == 0 || YxD == 0] <- 0
  
  # Return the calculated yield based on the uptake of nutrients 1 and 2
  return(Y12)
}



######################################################################
######################################################################


water_dependent_nutrient_uptake <- function(S1=NA, WLY=NA, d1=NA, a1=NA, r1=NA) {
  uptakeX_givenWater <- S1 - 0.25 * (S1 - r1 - WLY/d1)^2 / (WLY / a1 - WLY / d1)
  uptakeX_givenWater[S1 > r1 + 2*WLY/a1 - WLY/d1] <- (WLY / a1)[S1 > r1 + 2*WLY/a1 - WLY/d1]
  uptakeX_givenWater[S1 < r1 + WLY / d1] <- S1[S1 < r1 + WLY / d1]
  
  return(uptakeX_givenWater)
}


######################################################################
## if the crop parameters are known this function can be ised
NUE <- function(HI, CmaxNroots, CminNroots, CmaxNtops, CminNtops, CmaxProots, CminProots, CmaxPtops, CminPtops, 
                CmaxKroots, CminKroots, CmaxKtops, CminKtops){

  aN = round(1000 * HI/(HI * CmaxNroots + (1 - HI) * CmaxNtops), digits=0)
  dN = round(1000 * HI/(HI * CminNroots + (1 - HI) * CminNtops), digits=0)
  
  aP = round(1000 * HI/(HI * CmaxProots + (1 - HI) * CmaxPtops), digits=0)
  dP = round(1000 * HI/(HI * CminProots + (1 - HI) * CminPtops), digits=0)
  
  aK = round(1000 * HI/(HI * CmaxKroots + (1 - HI) * CmaxKtops), digits=0)
  dK = round(1000 * HI/(HI * CminKroots + (1 - HI) * CminKtops), digits=0)
  
  return(data.frame(aN=aN, dN=dN,aP=aP,dP=dP,aK=aK,dK=dK))
  
}



######################################################################
######################################################################


## Maize: aN=29	; dN=74	; aP=95	; dP=476; aK=38	; dK=143	; rN=5	; rP=0.4	, rK=2
## Cassava: aN=41, dN=96, rN=0, aP=233, dP=588, rP=0, aK=34, dK=161, rK=0)
##Yield_S <- function(SN, SP, SK, WLY, aN=41, dN=96, rN=0, aP=233, dP=588, rP=0, aK=34, dK=161, rK=0){ ## cassava
Yield_S <- function(SN, SP, SK, WLY, crop){
  if(crop == "Cassava"){
    aN=41; dN=96; rN=0; aP=233; dP=588; rP=0; aK=34; dK=161; rK=0
  }else if (crop == "Maize"){
    aN=29	; dN=74	; aP=95	; dP=476; aK=38	; dK=143	; rN=5	; rP=0.4	; rK=2
  }else if (crop == "Potato"){
    aN=175	; dN=435	; aP=625	; dP=4348; aK=109	; dK=513	; rN=0	; rP=0	; rK=0
  }
  
  ## uptake of one nutrient conditioned by availablility of one other nutrient or water, rule of minimum
  UNP <- nutrient_uptake(S1 = SN, S2 = SP, d1 = dN, a1 = aN, d2 = dP, a2 = aP, r1 = rN, r2 = rP)
  UNK <- nutrient_uptake(S1 = SN, S2 = SK, d1 = dN, a1 = aN, d2 = dK, a2 = aK, r1 = rN, r2 = rK)
  UNW <- water_dependent_nutrient_uptake(S1 = SN, WLY = WLY, d1 = dN, a1 = aN, r1 = rN)
  UN <- min(UNP, UNK, UNW)
  
  UPN <- nutrient_uptake(S1 = SP, S2 = SN, d1 = dP, a1 = aP, d2 = dN, a2 = aN, r1 = rP, r2 = rN)
  UPK <- nutrient_uptake(S1 = SP, S2 = SK, d1 = dP, a1 = aP, d2 = dK, a2 = aK, r1 = rP, r2 = rK)
  UPW <- water_dependent_nutrient_uptake(S1 = SP, WLY = WLY, d1 = dP, a1 = aP, r1 = rP)
  UP <- min(UPN, UPK, UPW)				
  
  UKN <- nutrient_uptake(S1 = SK, S2 = SN, d1 = dK, a1 = aK, d2 = dN, a2 = aN, r1 = rK, r2 = rN)
  UKP <- nutrient_uptake(S1 = SK, S2 = SP, d1 = dK, a1 = aK, d2 = dP, a2 = aP, r1 = rK, r2 = rP)
  UKW <- water_dependent_nutrient_uptake(S1 = SK, WLY = WLY, d1 = dK, a1 = aK, r1 = rK)
  UK <- min(UKN, UKP, UKW)		
  
  ## yield base don uptake of Nitrogen at dilution or accumulation
  YNA <- max((UN - rN), 0) * aN
  YND <- max((UN - rN), 0) * dN
  YPA <- max((UP - rP), 0) * aP
  YPD <- max((UP - rP), 0) * dP
  YKA <- max((UK - rK), 0) * aK
  YKD <- max((UK - rK), 0) * dK	
  
  ## yield as defined by uptake of one nutrient conditioned to max and min yield based on secnd nutrient not to be constrined by a third nutrient
  YNP <- yield_nutrients_combined(U1 = UN, d1 = dN, a1 = aN, Y2A = YPA, Y2D = YPD, Y3D = YKD, r1 = rN)
  YNK <- yield_nutrients_combined(U1 = UN, d1 = dN, a1 = aN, Y2A = YKA, Y2D = YKD, Y3D = YPD, r1 = rN)
  YPN <- yield_nutrients_combined(U1 = UP, d1 = dP, a1 = aP, Y2A = YNA, Y2D = YND, Y3D = YKD, r1 = rP)
  YPK <- yield_nutrients_combined(U1 = UP, d1 = dP, a1 = aP, Y2A = YKA, Y2D = YKD, Y3D = YND, r1 = rP)
  YKN <- yield_nutrients_combined(U1 = UK, d1 = dK, a1 = aK, Y2A = YNA, Y2D = YND, Y3D = YPD, r1 = rK)
  YKP <- yield_nutrients_combined(U1 = UK, d1 = dK, a1 = aK, Y2A = YPA, Y2D = YPD, Y3D = YND, r1 = rK)
  
  # Make sure the nutrient limited yields do not exceed the maximum possible yield = WLY
  YNPc <- min(c(YNP, YND, YPD, YKD, WLY))
  YNKc <- min(c(YNK, YND, YPD, YKD, WLY))
  YPNc <- min(c(YPN, YND, YPD, YKD, WLY))
  YPKc <- min(c(YPK, YND, YPD, YKD, WLY))
  YKNc <- min(c(YKN, YND, YPD, YKD, WLY))
  YKPc <- min(c(YKP, YND, YPD, YKD, WLY))
  
  #Final estimate
  YEc <- mean(c(YNPc, YNKc, YPNc, YPKc, YKNc, YKPc))
  return(YEc)				
}

######################################################################
######################################################################



## get yield estimates using QUEFTS and giving in SN, SP, SK and WLY.  get TSS for the estimates of yield and NOT blup. 
## SN, SP and SK are supposed to be the sum of what the soil originally has and added fertilizer. 
## added fertilizer is provided as it is defined in NOT trials and INS is the one to be defined by iterating on while controling the output
## by the squared yield difference between simulated and measured from NOT  (Blups of NOT)
optim_INS <- function(INS, addedFertilizer, YieldMatrix, WLY, crop){	
  
  if(crop == "Cassava"){
    RecoveryFraction <- data.frame(rec_N = 0.5,	rec_P = 0.21,	rec_K = 0.49)
  }else if (crop == "Maize"){
    RecoveryFraction <- data.frame(rec_N = 0.5,	rec_P = 0.1,	rec_K = 0.5)
  }else if (crop == "Potato"){
    RecoveryFraction <- data.frame(rec_N = 0.41, rec_P = 0.25,	rec_K = 0.65)
  }
  
  SoilN <- INS[1]
  SoilP <- INS[2]
  SoilK <- INS[3]
  
  #Yield_Control <- Yield_S(SN = addedFertilizer[1,1] + SoilN , SP = addedFertilizer[1,2] + SoilP, SK = addedFertilizer[1,3] + SoilK, WLY = WLY)
  Yield_pk <- Yield_S(SN = (addedFertilizer[2,1]*RecoveryFraction$rec_N) + SoilN , 
                      SP = (addedFertilizer[2,2]*RecoveryFraction$rec_P) + SoilP, 
                      SK = (addedFertilizer[2,3]*RecoveryFraction$rec_K) + SoilK, WLY = WLY, crop=crop)
  
  Yield_nk <- Yield_S(SN = (addedFertilizer[3,1]*RecoveryFraction$rec_N) + SoilN ,
                      SP = (addedFertilizer[3,2]*RecoveryFraction$rec_P) + SoilP, 
                      SK = (addedFertilizer[3,3]*RecoveryFraction$rec_K) + SoilK, WLY = WLY, crop=crop)
  
  Yield_np <- Yield_S(SN = (addedFertilizer[4,1]*RecoveryFraction$rec_N) + SoilN , 
                      SP = (addedFertilizer[4,2]*RecoveryFraction$rec_P) + SoilP, 
                      SK = (addedFertilizer[4,3]*RecoveryFraction$rec_K) + SoilK, WLY = WLY, crop=crop)

  
  #control_SS <- (Yield_Control - YieldMatrix$NOT_yield[1])^2
  pk_SS <- (Yield_pk - YieldMatrix$NOT_yield[2])^2
  nk_SS <- (Yield_nk - YieldMatrix$NOT_yield[3])^2
  np_SS <- (Yield_np - YieldMatrix$NOT_yield[4])^2
  #halfnpk_SS <- (Yield_halfnpk - YieldMatrix$NOT_yield[5])^2
  
  TSS <- sum( pk_SS, nk_SS, np_SS)	
  return(TSS)	
}

######################################################################
######################################################################


# PD = "2018-03-01"; HD = "2019-05-31"; lat = 10.024; lon = 4.025; country = "NG"; cassUW = 1000; cassUP = 12000; maxInv = 72000;
# SoilData = SoilGridData_NG , userName = "acai cassava", userPhoneCC = 254, userPhoneNr = 702974480,userEmail = "acai.akilimo@gmail.com",
# cassPD = "roots"
#' 
#' @param dss 
#' @returnType 
#' @return 
#' 
#' @author Meklit
#' @export
getsupply <- function(dss){
  supply <- data.frame(lat=dss$lat, long=dss$long, rel_N=dss$rel_N, rel_P=dss$rel_P, rel_K=dss$rel_K, SN=dss$soilN, SP=dss$soilP, SK=dss$soilK, water_limited_yield = dss$water_limited_yield,
                       aN=dss$aN, dN=dss$dN, aP=dss$aP, dP=dss$dP, aK=dss$aK, dK=dss$dK, rN=dss$rN, rP=dss$rP, rK=dss$rK, max_yield=dss$max_yield,  tolerance=dss$tolerance,
                       WLY = dss$water_limited_yield)
}

######################################################################
######################################################################


actual_uptake_tool <- function(df){ 
  UNP <- nutrient_uptake(S1 = df$SN, S2 = df$SP, d1 = df$dN, a1 = df$aN, d2 = df$dP, a2 = df$aP, r1 = df$rN, r2 = df$rP)
  UNK <- nutrient_uptake(S1 = df$SN, S2 = df$SK, d1 = df$dN, a1 = df$aN, d2 = df$dK, a2 = df$aK, r1 = df$rN, r2 = df$rK)
  UNW <- water_dependent_nutrient_uptake(S1 = df$SN, WLY = df$WLY, d1 = df$dN, a1 = df$aN, r1 = df$rN)
  df$UN <- pmin(UNP, UNK, UNW)
  
  UPN <- nutrient_uptake(S1 = df$SP, S2 = df$SN, d1 = df$dP, a1 = df$aP, d2 = df$dN, a2 = df$aN, r1 = df$rP, r2 = df$rN)
  UPK <- nutrient_uptake(S1 = df$SP, S2 = df$SK, d1 = df$dP, a1 = df$aP, d2 = df$dK, a2 = df$aK, r1 = df$rP, r2 = df$rK)
  UPW <- water_dependent_nutrient_uptake(S1 = df$SP, WLY = df$WLY, d1 = df$dP, a1 = df$aP, r1 = df$rP)
  df$UP <- pmin(UPN, UPK, UPW)
  
  UKN <- nutrient_uptake(S1 = df$SK, S2 = df$SN, d1 = df$dK, a1 = df$aK, d2 = df$dN, a2 = df$aN, r1 = df$rK, r2 = df$rN)
  UKP <- nutrient_uptake(S1 = df$SK, S2 = df$SP, d1 = df$dK, a1 = df$aK, d2 = df$dP, a2 = df$aP, r1 = df$rK, r2 = df$rP)
  UKW <- water_dependent_nutrient_uptake(S1 = df$SK, WLY = df$WLY, d1 = df$dK, a1 = df$aK, r1 = df$rK)
  df$UK <- pmin(UKN, UKP, UKW)  
  return(df)
}

######################################################################
######################################################################

max_min_yields_tools <- function(dss){ 
  
  dss$YNA <- pmax((dss$UN - dss$rN), 0) * dss$aN
  dss$YND <- pmax((dss$UN - dss$rN), 0) * dss$dN
  dss$YPA <- pmax((dss$UP - dss$rP), 0) * dss$aP
  dss$YPD <- pmax((dss$UP - dss$rP), 0) * dss$dP
  dss$YKA <- pmax((dss$UK - dss$rK), 0) * dss$aK
  dss$YKD <- pmax((dss$UK - dss$rK), 0) * dss$dK
  
  return(dss)
}

######################################################################
######################################################################

quefts_tools <- function(supply_wly){
  # Actual uptake of nutrients.
  supply_wly <- actual_uptake_tool(supply_wly)
  
  # Maximum and minimum yields, depending on maximum accumulation and dilution.
  supply_wly <- max_min_yields_tools(supply_wly)
  
  # Final yield based on the combinations of nutrient uptake and minimum + maximum yields.
  supply_wly$TargetYield <- final_yield_tools(supply_wly)
  
  return(supply_wly)
}

######################################################################
######################################################################

final_yield_tools <- function(df){
  YNP <- yield_nutrients_combined(U1 = df$UN, d1 = df$dN, a1 = df$aN, Y2A = df$YPA, Y2D = df$YPD, Y3D = df$YKD, r1 = df$rN)
  YNK <- yield_nutrients_combined(U1 = df$UN, d1 = df$dN, a1 = df$aN, Y2A = df$YKA, Y2D = df$YKD, Y3D = df$YPD, r1 = df$rN)
  YPN <- yield_nutrients_combined(U1 = df$UP, d1 = df$dP, a1 = df$aP, Y2A = df$YNA, Y2D = df$YND, Y3D = df$YKD, r1 = df$rP)
  YPK <- yield_nutrients_combined(U1 = df$UP, d1 = df$dP, a1 = df$aP, Y2A = df$YKA, Y2D = df$YKD, Y3D = df$YND, r1 = df$rP)
  YKN <- yield_nutrients_combined(U1 = df$UK, d1 = df$dK, a1 = df$aK, Y2A = df$YNA, Y2D = df$YND, Y3D = df$YPD, r1 = df$rK)
  YKP <- yield_nutrients_combined(U1 = df$UK, d1 = df$dK, a1 = df$aK, Y2A = df$YPA, Y2D = df$YPD, Y3D = df$YND, r1 = df$rK)
  
  
  # Make sure the nutrient limited yields do not exceed the maximum possible yield = WLY
  YNPc <- pmin(YNP, df$YND, df$YPD, df$YKD, df$WLY)
  YNKc <- pmin(YNK, df$YND, df$YPD, df$YKD, df$WLY)
  YPNc <- pmin(YPN, df$YND, df$YPD, df$YKD, df$WLY)
  YPKc <- pmin(YPK, df$YND, df$YPD, df$YKD, df$WLY)
  YKNc <- pmin(YKN, df$YND, df$YPD, df$YKD, df$WLY)
  YKPc <- pmin(YKP, df$YND, df$YPD, df$YKD, df$WLY)
  
  
  #Final estimate
  YEc <- (YNPc + YNKc + YPNc + YPKc + YKNc + YKPc) / 6
  
  YEc[YEc > df$WLY] <- df$WLY[YEc > df$WLY]
  
  
  return(YEc)
}


######################################################################
######################################################################
#' @description fucntion to run teh reverse QUEFTS to obtain apparnet soil indegenous NPK supply, given observed yield and rate of fertilizer used
#' @param ds_GPS is a data frame with dry matter/grain yield for the NPK, NP, NK, PK and contreol yield
#' @param rateN N in the fertilizer used as full NPK, NP and NK treatment
#' @param rateP P in the fertilizer used as full NPK, NP and PK treatment
#' @param rateK K in the fertilizer used as full NPK, NK and PK treatment
#' @param halfrate trueif other treatment tested in NOT setup with differnt N, P and K rate
#' @param hrN N rate for the halfrate treatment
#' @param hrP P rate for the halfrate treatment
#' @param hrK K rate for the halfrate treatment
#' @return ds_GPS with columns for soil N, P and K and estimated yield for the control treatment generated using estiamted soil NPK. 
Reverse_QUEFTS <- function(ds_GPS, rateN, rateP, rateK, crop, halfrate=FALSE, hrN=0, hrP=0, hrK=0){
  soilNPK <- NULL
  for(i in 1:nrow(ds_GPS)){
    print(i)
    oneTrial <- ds_GPS[i,]	
    # WLY <- oneTrial$NPK/3*1000 # convert to kg.ha
    WLY <- oneTrial$DM_NPK # for cassava this is NPK yield in kg/ha dry root
    
    if(halfrate == TRUE){
      addedFertilizer <- as.data.frame(matrix(nrow=5,ncol=3))
      colnames(addedFertilizer) <- c("FN", "FP", "FK")
      addedFertilizer[,1] <- c(0, 0, rateN, rateN, hrN)
      addedFertilizer[,2] <- c(0, rateP, 0, rateP, hrP)
      addedFertilizer[,3] <- c(0, rateK, rateK, 0, hrK)	
      
      YieldMatrix <- as.data.frame(matrix(nrow=5,ncol=2))
      colnames(YieldMatrix) <- c("treatment", "NOT_yield")
      YieldMatrix[,1] <- c("Control", "PK", "NK", "NP", "halfNPK")
      YieldMatrix[,2] <- c(oneTrial$DM_CON, oneTrial$DM_PK, oneTrial$DM_NK, oneTrial$DM_NP, oneTrial$DM_half_NPK)
      
    }else{
      
      addedFertilizer <- as.data.frame(matrix(nrow=4,ncol=3))
      colnames(addedFertilizer) <- c("FN", "FP", "FK")
      addedFertilizer[,1] <- c(0, 0, rateN, rateN)
      addedFertilizer[,2] <- c(0, rateP, 0, rateP)
      addedFertilizer[,3] <- c(0, rateK, rateK, 0)	
      
      YieldMatrix <- as.data.frame(matrix(nrow=4,ncol=2))
      colnames(YieldMatrix) <- c("treatment", "NOT_yield")
      YieldMatrix[,1] <- c("Control", "PK", "NK", "NP")
      YieldMatrix[,2] <- c(oneTrial$DM_CON, oneTrial$DM_PK, oneTrial$DM_NK, oneTrial$DM_NP)
    }
    
    ## optimizing by min TSS. 
    soil_NPK <- optim(par=c(0,0,0), optim_INS, lower=c(0.1, 0.1, 0.1), method = "L-BFGS-B", addedFertilizer=addedFertilizer, YieldMatrix=YieldMatrix, WLY=WLY, crop=crop)$par	
    
   
    if(!is.na(oneTrial$DM_CON)){
      Yield_control <- Yield_S(SN = addedFertilizer[1,1] + soil_NPK[1] , SP = addedFertilizer[1,2] + soil_NPK[2], 
                               SK = addedFertilizer[1,3] + soil_NPK[3], WLY = WLY, crop=crop)
      oneTrial$Control_observed <- oneTrial$DM_CON
      oneTrial$Control_estimated <-  Yield_control
    }else{
      oneTrial$Control_observed <- NA
      oneTrial$Control_estimated <-  NA
    }
    
    if(halfrate == TRUE){
      Yield_halfnpk <- Yield_S(SN = addedFertilizer[5,1] + soil_NPK[1] , SP = addedFertilizer[5,2] + soil_NPK[2], 
                               SK = addedFertilizer[5,3] + soil_NPK[3], WLY = WLY, crop=crop)
      oneTrial$halfNPK_estimated <-  Yield_halfnpk
      
    }
    
    
    snpk <- oneTrial
    snpk$soilN <- soil_NPK[1]
    snpk$soilP <- soil_NPK[2]
    snpk$soilK <- soil_NPK[3]
    soilNPK <- rbind(soilNPK, snpk)
  }
  
  # 
  
  return(soilNPK)
}


######################################################################
######################################################################

#' The soil NPK as obtained from randdom forest model
#' @param zone, is used to define the varieties and HI to get NUE. Lake zone (the default) is proxy for Mkobozi and E & S zone is Kiroba  
#' @param Queft_Input_Data: per lat and long, crop param and soil param, water limited yield, fertlizer recommendation
#' @return 
#' 
#' @author Meklit
#' @export
QUEFTS_WLY_CY <- function(SoilData=SoilData, country=country, wlyd=wlyd, crop){	
  #wly_plDate <- wly_data[wly_data$plantingDate == pl_Date, c("lat", "long", "wly_KgHa")]
  wly_plDate <- wlyd[,  c("lat", "long", "water_limited_yield")]
  
  # colnames(wly_plDate) <- c("lat", "long", "water_limited_yield")	
  Quefts_Input_Data_wly <- merge(SoilData, wly_plDate, by=c("lat", "long"))
  
  if(crop == "cassava"){
    crop_param <- data.frame(aN = 41, dN = 96, aP = 233, dP = 588, aK = 34, dK = 161, rec_N = 0.5, rec_P = 0.21, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=Quefts_Input_Data_wly$water_limited_yield, tolerance=0.01,
                             lat=Quefts_Input_Data_wly$lat, long=Quefts_Input_Data_wly$long )
  }else if (crop == "potato"){
    crop_param <- data.frame(aN = 175, dN = 435, aP = 625, dP = 4348, aK = 109, dK = 513, rec_N = 0.41, rec_P = 0.25, rec_K = 0.65,
                             rN=0, rP=0, rK=0, max_yield=Quefts_Input_Data_wly$water_limited_yield, tolerance=0.01,
                             lat=Quefts_Input_Data_wly$lat, long=Quefts_Input_Data_wly$long)
  }else if (crop == "maize"){
    crop_param <- data.frame(aN = 29, dN = 74, aP = 95, dP = 476, aK = 38, dK = 143, rec_N = 0.5, rec_P = 0.1, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=Quefts_Input_Data_wly$water_limited_yield, tolerance=0.01,
                             lat=Quefts_Input_Data_wly$lat, long=Quefts_Input_Data_wly$long)
  }else if (crop == "soybean"){
    crop_param <- data.frame(aN = 13, dN = 22, aP = 100, dP = 286, aK = 43, dK = 71, rec_N = 0.5, rec_P = 0.21, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=Quefts_Input_Data_wly$water_limited_yield, tolerance=0.01,
                             lat=Quefts_Input_Data_wly$lat, long=Quefts_Input_Data_wly$long)
  }else if (crop == "commonbean"){
    crop_param <- data.frame(aN = 13, dN = 22, aP = 21, dP = 177, aK = 48, dK = 125, rec_N = 0.5, rec_P = 0.21, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=Quefts_Input_Data_wly$water_limited_yield, tolerance=0.01,
                             lat=Quefts_Input_Data_wly$lat, long=Quefts_Input_Data_wly$long)
    
  }
 
  Quefts_Input_Data_wly$rel_N <- 1
  Quefts_Input_Data_wly$rel_P <- Quefts_Input_Data_wly$soilP / Quefts_Input_Data_wly$soilN
  Quefts_Input_Data_wly$rel_K <- Quefts_Input_Data_wly$soilK / Quefts_Input_Data_wly$soilN
  
  
  ## 1. get soil nutrient supply
  Queft_Input_Data_Var <- cbind(Quefts_Input_Data_wly, crop_param)	
  supply <- getsupply(Queft_Input_Data_Var) ## to get yield at zero input level
  names(supply)
  
  ## 2. Current yield: 
  actualUptake <- actual_uptake_tool(supply)
  names(actualUptake)
  minmax_Yield <-  max_min_yields_tools(actualUptake)
  names(minmax_Yield)
  Current_Yield <- ddply(minmax_Yield,.(lat, long), final_yield_tools)## yield at zero input
  colnames(Current_Yield) <- c("lat", "long", "CurrentYield")
  Yield_Fertilizer <- merge(wly_plDate, Current_Yield, by=c("lat", "long"))
  return(Yield_Fertilizer)
}

######################################################################
######################################################################
######################################################################
######################################################################

#' @description Function to convert root DM yield into root fresh matter yield (RFY).Function to predict root FM yield based on date of harvest and country, using data from gravimetric starch measurements conducted across ACAI trials.
#' @param HD: harvest date (Date format)
#' @param RDY: root dry matter yield (user's units)
#' @param country = c("NG", "TZ")
#' @return RFY: root fresh yield in the same units as root DM yield input
getRFY <- function(HD, 
                   RDY, 
                   country = c("NG", "TZ")){
  

  # d  <- as.numeric(strftime(HD, format = "%j"))
  d <- HD
  fd <- read.csv("fd2.csv") #data.frame with day of the year (dayNr = [1..366]) and %DM (DMCont = [0..100], by country)
  DC <- merge(data.frame(dayNr=d), fd[fd$country=="NG",], sort=FALSE)$DMCont
  RFY <- RDY / DC * 100
  
  return(RFY)
  
}
######################################################################
######################################################################



#' @description Function to convert root FM yield into root dry matter yield (RDY): user define CY in FM in ton/ha, QUEFTS require Cy in DM kg/ha
#' @param HD: harvest date (Date format)
#' @param RFY: root fresh matter yield (user's units)
#' @param country = c("NG", "TZ")
getRDY <- function(HD, RFY, country){
  if(HD > 366){
    HD <- HD - 366
  }
  d <- HD
  fd <- read.csv("fd2.csv") 
  DC <- merge(data.frame(dayNr=d), fd[fd$country==country,], sort=FALSE)$DMCont
  RDY <- (RFY * DC)/100
  return(RDY)
}

######################################################################
######################################################################


#' computes target yield in tonnes/ha from a given NPK rate
#' @param QID a data frame containing soil NPK, WLY (kg/ha dry wt.), 
#' @param rec recomended NPK rate
#' @returnType 
#' @return target yield in ton/ha dry matter
#' 
#' @author Meklit
#' @export
QUEFTS1_Pedotransfer <- function(QID, rec, crop){		
  QID$WLY <- QID$water_limited_yield
  
  
  if(crop == "cassava"){
    crop_param <- data.frame(aN = 41, dN = 96, aP = 233, dP = 588, aK = 34, dK = 161, rec_N = 0.5, rec_P = 0.21, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=QID$WLY, tolerance=0.01, lat=QID$lat, long=QID$lon )
  }else if (crop == "potato"){
    crop_param <- data.frame(aN = 175, dN = 435, aP = 625, dP = 4348, aK = 109, dK = 513, rec_N = 0.41, rec_P = 0.25, rec_K = 0.65,
                             rN=0, rP=0, rK=0, max_yield=QID$WLY, tolerance=0.01, lat=QID$lat, long=QID$lon)
  }else if (crop == "maize"){
    crop_param <- data.frame(aN = 29, dN = 74, aP = 95, dP = 476, aK = 38, dK = 143, rec_N = 0.5, rec_P = 0.1, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=QID$WLY, tolerance=0.01, lat=QID$lat, long=QID$lon)
  }else if (crop == "soybean"){
    crop_param <- data.frame(aN = 13, dN = 22, aP = 100, dP = 286, aK = 43, dK = 71, rec_N = 0.5, rec_P = 0.21, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=QID$WLY, tolerance=0.01, lat=QID$lat, long=QID$lon)
  }else if (crop == "commonbean"){
    crop_param <- data.frame(aN = 13, dN = 22, aP = 21, dP = 177, aK = 48, dK = 125, rec_N = 0.5, rec_P = 0.21, rec_K = 0.5,
                             rN=0, rP=0, rK=0, max_yield=QID$WLY, tolerance=0.01, lat=QID$lat, long=QID$lon)
    
  }
  
  Queft_Input_Data_Var1 <- cbind(QID, crop_param)
  indata <- Queft_Input_Data_Var1[,c("lat","long" ,"WLY","aN", "dN", "aP", "dP","aK","dK", "rN", "rP", "rK", "soilN", "soilP", "soilK","max_yield", "tolerance")]
  
  N_rate <- rec[1]
  P_rate <- rec[2]
  K_rate <- rec[3]
  
  TargetYield_from_NPK <- NPK_TargetYield_forOutput(NutrUse_soilNPK=indata, N_rate, P_rate, K_rate)
  
  return(TargetYield_from_NPK$TargetYield)
}


#####################################################################
######################################################################
######################################################################




#' using the output of function "NPK_TargetYield_forinput" and a dat frame per lon and lat for intended NPK input
#' this function calculates the yield that can be obtained for intended NPK rate.
#' @param NutrUse_soilNPK 
#' @param NPKdata: needs to be provided
#' @return 
#' 
#' @author Meklit
#' @export
NPK_TargetYield_forOutput <- function(NutrUse_soilNPK, N_rate, P_rate, K_rate){	
  NutrUse_soilNPK$N_rate <- N_rate
  NutrUse_soilNPK$P_rate <- P_rate
  NutrUse_soilNPK$K_rate <- K_rate
  
  ## Supply of nutrients to the crop
  NutrUse_soilNPK$SN <- NutrUse_soilNPK$N_rate + NutrUse_soilNPK$soilN
  NutrUse_soilNPK$SP <- NutrUse_soilNPK$P_rate + NutrUse_soilNPK$soilP
  NutrUse_soilNPK$SK <- NutrUse_soilNPK$K_rate + NutrUse_soilNPK$soilK
  
  ## Actual Uptake of nutrients: crop param + nutrient supply
  tmp <- ddply(NutrUse_soilNPK,.(lat, long), actual_uptake_tool)

  ## max and min yield: actual uptake and crop param. min of N uptake constrianed by availability of P, K and water
  maxminY <- ddply(tmp,.(lat, long), max_min_yields_tools)

  ## final yield: min yield for combined uptake of 2 nutrients assuming the 3rd is not limiting, should be < WLY, and take meanof the six combinations 
  Target_Yield <- ddply(maxminY,.(lat, long), quefts_tools)	
  TY <- data.frame(lat=Target_Yield$lat, lon=Target_Yield$long, TargetYield=Target_Yield$TargetYield)
  
  return(TY)
}

######################################################################

#' Title get CY
#'
#' @param SoilData 
#' @param LINTUL_WLY 
#'
#' @return
#' @export
#'
#' @examples
get_CY <- function(SoilData, LINTUL_WLY, crop, country){
  
  LINTUL_WLY$location <- paste(LINTUL_WLY$lat, LINTUL_WLY$lon, sep="_")
  SoilData$location <- paste(SoilData$lat, SoilData$lon, sep="_")
  SoilData <- SoilData[SoilData$location %in% LINTUL_WLY$location, ]
  LINTUL_WLY <- LINTUL_WLY[LINTUL_WLY$location %in% SoilData$location, ]
  
  ## one planting week
  pl_Date <- unique(LINTUL_WLY$plantingDate)
  pl_cy <- NULL
  for(pld in pl_Date){
    
    wlypd_ll <- LINTUL_WLY[LINTUL_WLY$plantingDate == pld, ]
    wlyHDs <- c("WLY_242", "WLY_249", "WLY_256", "WLY_263", "WLY_270", "WLY_277", "WLY_284",
                "WLY_291", "WLY_298", "WLY_305", "WLY_312", "WLY_319", "WLY_326", "WLY_333",
                "WLY_340", "WLY_347", "WLY_354", "WLY_361", "WLY_368", "WLY_375", "WLY_382",
                "WLY_389", "WLY_396", "WLY_403", "WLY_410", "WLY_417", "WLY_424", "WLY_431",
                "WLY_438", "WLY_445", "WLY_452")
    
    
    ## one harvest week
    pl_hv_cy <- NULL
    for(hd in wlyHDs){
      wlydata2 <- wlypd_ll[, c("lat", "long", "location", "plantingDate", hd)]
      hdays <- as.numeric(as.character(gsub("WLY_", "", hd)))
      colnames(wlydata2) <- c("lat", "long", "location", "pl_Date", "water_limited_yield")
      wlydata2$zone <- country
      wlydata2$daysOnField <- hdays
      wlydata2 <- wlydata2[, c("lat", "long", "water_limited_yield", "location", "pl_Date", "zone", "daysOnField")]
      
      SoilData$long <- SoilData$lon
      
      ## get CY
      Current_Yield <- QUEFTS_WLY_CY(SoilData = SoilData, country = country, wlyd = wlydata2, crop = crop)
      
      wlydata2 <- unique(merge(wlydata2, Current_Yield[, c("lat", "long", "CurrentYield")], by=c("lat", "long")))
      wlydata2$weekNr <- unique(wlypd_ll$weekNr)
      
      pl_hv_cy <- rbind(pl_hv_cy, wlydata2 )
      
    }
    
    pl_cy <- rbind(pl_cy, pl_hv_cy )
  }
  return(pl_cy)
  
}


### get fertilizer types and prices 
fertilizerFunc <- function(NPK201216available = FALSE, NPK201216CostperBag = NA, NPK201216BagWt = 50,
                           ureaavailable = FALSE, ureaCostperBag = NA, ureaBagWt = 50,
                           MOPavailable = FALSE, MOPCostperBag = NA, MOPBagWt = 50,
                           DAPavailable = FALSE, DAPCostperBag = NA, DAPBagWt = 50,
                           NPK201010available = FALSE, NPK201010CostperBag = NA, NPK201010BagWt = 50,
                           NPK151515available = FALSE, NPK151515CostperBag = NA, NPK151515BagWt = 50,
                           TSPavailable = FALSE, TSPCostperBag = NA, TSPBagWt = 50,
                           NPK171717available = FALSE, NPK171717CostperBag = NA, NPK171717BagWt = 50,
                           CANavailable = FALSE, CANCostperBag = NA, CANBagWt = 50,
                           SSPavailable = FALSE, SSPCostperBag = NA, SSPBagWt = 50, 
                           NPK112221available = FALSE, NPK112221CostperBag = NA, NPK112221BagWt = 50, 
                           NPK251010available = FALSE, NPK251010CostperBag = NA, NPK251010BagWt = 50, 
                           NPK152020available = FALSE, NPK152020CostperBag = NA, NPK152020BagWt = 50, 
                           NPK23105available = FALSE, NPK23105CostperBag = NA, NPK23105BagWt = 50, 
                           NPK123017available = FALSE, NPK123017CostperBag = NA, NPK123017BagWt = 50, 
                           country) {
  
  if (country == "NG") { 
    if (ureaCostperBag == 0) { ureaCostperBag <- 7500 } else { ureaCostperBag <- as.numeric(ureaCostperBag) }
    if (MOPCostperBag == 0) { MOPCostperBag <- 13500 }else { MOPCostperBag <- as.numeric(MOPCostperBag) }
    if (DAPCostperBag == 0) { DAPCostperBag <- 13250 }else { DAPCostperBag <- as.numeric(DAPCostperBag) }
    if (NPK201010CostperBag == 0) { NPK201010CostperBag <- 7200 }else { NPK201010CostperBag <- as.numeric(NPK201010CostperBag) }
    if (NPK151515CostperBag == 0) { NPK151515CostperBag <- 8500 }else { NPK151515CostperBag <- as.numeric(NPK151515CostperBag) }
    if (TSPCostperBag == 0) { TSPCostperBag <- 13250 }else { TSPCostperBag <- as.numeric(TSPCostperBag) }
    if (NPK171717CostperBag == 0) { NPK171717CostperBag <- 7800 }else { NPK171717CostperBag <- as.numeric(NPK171717CostperBag) }
    if (CANCostperBag == 0) { CANCostperBag <- 7500 }else { CANCostperBag <- as.numeric(CANCostperBag) }
    if (SSPCostperBag == 0) { SSPCostperBag <- 22364 }else { SSPCostperBag <- as.numeric(SSPCostperBag) }
    if (NPK201216CostperBag == 0) { NPK201216CostperBag <- 8000 }else { NPK201216CostperBag <- as.numeric(NPK201216CostperBag) }
  }else if(country == "TZ"){
    if (ureaCostperBag == 0) { ureaCostperBag <- 58000 }else { ureaCostperBag <- as.numeric(ureaCostperBag) } ##65000
    if (MOPCostperBag == 0) { MOPCostperBag <- 67000 }else { MOPCostperBag <- as.numeric(MOPCostperBag) } ##120000
    if (DAPCostperBag == 0) { DAPCostperBag <- 60000 }else { DAPCostperBag <- as.numeric(DAPCostperBag) } ##85000
    if (NPK201010CostperBag == 0) { NPK201010CostperBag <- 68000 }else { NPK201010CostperBag <- as.numeric(NPK201010CostperBag) }
    if (TSPCostperBag == 0) { TSPCostperBag <- 64000 }else { TSPCostperBag <- as.numeric(TSPCostperBag) } #80000
    if (NPK171717CostperBag == 0) { NPK171717CostperBag <- 63000 }else { NPK171717CostperBag <- as.numeric(NPK171717CostperBag) } ##61000
    if (CANCostperBag == 0) { CANCostperBag <- 53000 }else { CANCostperBag <- as.numeric(CANCostperBag) } #65000
   
  }else if(country == "Rw"){
    if (ureaCostperBag == 0) { ureaCostperBag <- 28200 }else { ureaCostperBag <- as.numeric(ureaCostperBag) }#564*50
    if (MOPCostperBag == 0) { MOPCostperBag <- 27900}else { MOPCostperBag <- as.numeric(MOPCostperBag) }#558*50 
    if (DAPCostperBag == 0) { DAPCostperBag <- 31650 }else { DAPCostperBag <- as.numeric(DAPCostperBag) } #633*50
    if (NPK171717CostperBag == 0) { NPK171717CostperBag <- 35650}else { NPK171717CostperBag <- as.numeric(NPK171717CostperBag) }#713 *50
  }else if(country == "GH"){
    if (ureaCostperBag == 0) { ureaCostperBag <- 142 }else { ureaCostperBag <- as.numeric(ureaCostperBag) }
    if (NPK112221CostperBag == 0) { NPK112221CostperBag <- 160}else { NPK112221CostperBag <- as.numeric(NPK112221CostperBag) }
    if (NPK251010CostperBag == 0) { NPK251010CostperBag <- 150 }else { NPK251010CostperBag <- as.numeric(NPK251010CostperBag) }
    if (NPK152020CostperBag == 0) { NPK152020CostperBag <- 156}else { NPK152020CostperBag <- as.numeric(NPK152020CostperBag) }
    if (NPK201010CostperBag == 0) { NPK201010CostperBag <- 146 }else { NPK201010CostperBag <- as.numeric(NPK201010CostperBag) }
    if (NPK23105CostperBag == 0)  { NPK23105CostperBag <- 146}else { NPK23105CostperBag <- as.numeric(NPK23105CostperBag) } 
    if (NPK123017CostperBag == 0) { NPK123017CostperBag <- 160 }else { NPK123017CostperBag <- as.numeric(NPK123017CostperBag) } 
  }
  
  ureaFert <- data.frame(type = 'Urea', available = ureaavailable, N_cont = 0.46, P_cont = 0, K_cont = 0, costPerBag = ureaCostperBag, bagWeight = ureaBagWt)
  MOPFert <- data.frame(type = 'MOP', available = MOPavailable, N_cont = 0.00, P_cont = 0.00, K_cont = 0.60, costPerBag = MOPCostperBag, bagWeight = MOPBagWt)
  DAPFert <- data.frame(type = 'DAP', available = DAPavailable, N_cont = 0.18, P_cont = 0.20, K_cont = 0.0, costPerBag = DAPCostperBag, bagWeight = DAPBagWt)
  NPK201010Fert <- data.frame(type = 'NPK20_10_10', available = NPK201010available, N_cont = 0.20, P_cont = 0.044, K_cont = 0.083, costPerBag = NPK201010CostperBag, bagWeight = NPK201010BagWt)
  NPK151515Fert <- data.frame(type = 'NPK15_15_15', available = NPK151515available, N_cont = 0.15, P_cont = 0.07, K_cont = 0.125, costPerBag = NPK151515CostperBag, bagWeight = NPK151515BagWt)
  TSPFert <- data.frame(type = 'TSP', available = TSPavailable, N_cont = 0.0, P_cont = 0.2, K_cont = 0.0, costPerBag = TSPCostperBag, bagWeight = TSPBagWt)
  NPK171717Fert <- data.frame(type = 'NPK17_17_17', available = NPK171717available, N_cont = 0.17, P_cont = 0.074, K_cont = 0.15, costPerBag = NPK171717CostperBag, bagWeight = NPK171717BagWt) # TODO get price
  NPK201216Fert <- data.frame(type = 'NPK20_12_26', available = NPK201216available, N_cont = 0.20, P_cont = 0.052, K_cont = 0.216, costPerBag = NPK201216CostperBag, bagWeight = NPK201216BagWt)
  CANFert <- data.frame(type = 'CAN', available = CANavailable, N_cont = 0.27, P_cont = 0.00, K_cont = 0.00, costPerBag = CANCostperBag, bagWeight = CANBagWt) ## not correct value TODO check
  SSPFert <- data.frame(type = 'SSP', available = SSPavailable, N_cont = 0.00, P_cont = 0.15, K_cont = 0.00, costPerBag = SSPCostperBag, bagWeight = SSPBagWt) ## not correct value TODO check
 
  NPK112221Fert <- data.frame(type = 'NPK112221', available = NPK112221available, N_cont = 0.11, P_cont = 0.1, K_cont = 0.17, costPerBag = 160, bagWeight = 50)
  NPK251010Fert <- data.frame(type = 'NPK251010', available = NPK251010available, N_cont = 0.25, P_cont = 0.044, K_cont = 0.083, costPerBag = 150, bagWeight = 50)
  NPK152020Fert <- data.frame(type = 'NPK152020', available = NPK152020available, N_cont = 0.15, P_cont = 0.088, K_cont = 0.166, costPerBag = 156, bagWeight = 50)
  NPK201010Fert <- data.frame(type = 'NPK201010', available = NPK201010available, N_cont = 0.20, P_cont = 0.044, K_cont = 0.083, costPerBag = 146, bagWeight = 50)
  NPK23105Fert <- data.frame(type = 'NPK23105', available = NPK23105available, N_cont = 0.23, P_cont = 0.044, K_cont = 0.0415, costPerBag = 146, bagWeight = 50)
  NPK123017Fert <- data.frame(type = 'NPK123017', available = NPK123017available, N_cont = 0.12, P_cont = 0.132, K_cont = 0.14, costPerBag = 160, bagWeight = 50)
  
  if (country == "NG") {
    fd_cont <- rbind(ureaFert, MOPFert, DAPFert, CANFert, NPK171717Fert, NPK151515Fert, NPK201010Fert, TSPFert, SSPFert, NPK201216Fert)
  }else if (country == "NG") {
    fd_cont <- rbind(ureaFert, MOPFert, DAPFert, CANFert, NPK171717Fert, NPK151515Fert, NPK201010Fert, TSPFert, SSPFert)
  }else if (country == "RW") {
    fd_cont <- rbind(ureaFert, MOPFert, DAPFert, NPK171717Fert)
  }else if (country == "GH") {
    fd_cont <- rbind(ureaFert, NPK112221Fert, NPK251010Fert, NPK152020Fert, NPK201010Fert, NPK23105Fert, NPK123017Fert)
  }
  
  
  fd_cont <- droplevels(fd_cont[fd_cont$available == TRUE,])
  fd_cont$costPerBag <- as.numeric(fd_cont$costPerBag)
  fd_cont$price <- fd_cont$costPerBag / fd_cont$bagWeight
  fd_cont <- subset(fd_cont, select = -c(available))

  return(fd_cont)
}






#' get optimized fertilizer rate, target yield for the recommended rate and net revenue given cost and investment
run_Optim <- function(rootUP, QID, fertilizer, invest, plDatell, lat, lon, areaHa, HD, WLY, DCY, country, crop) {
  
  
  ## input of CY and WLY are in dry wt in KG/ha
  QID$water_limited_yield <- WLY
  initial <- rep(0, nrow(fertilizer))
  lowerST <- rep(0, nrow(fertilizer))
  
  ## both CY and TY should be changed to user land size and should be in ton/ha and fresh wt
  if(country %in% c("NG", "TZ")){
    CY_user <- ((getRFY(HD = HD, RDY = DCY, country = country)) / 1000) * areaHa ## TZ model is extrememly high
    WLY_user <- ((getRFY(HD = HD, RDY = WLY, country = country)) / 1000) * areaHa
  } else if (country == "Rw"){
    CY_user <- DCY * 0.003 * areaHa 
    WLY_user <- WLY * 0.003 * areaHa
  }else if (country == "GH"){
    CY_user <- ((getRFY(HD = HD, RDY = DCY, country = "NG")) / 1000) * areaHa ## TZ model is extrememly high
    WLY_user <- ((getRFY(HD = HD, RDY = WLY, country = "NG")) / 1000) * areaHa
  }
  

  FR <- optim(par = initial, fn = optim_NR, lower = lowerST, method = "L-BFGS-B", control = list(fnscale = -1), rootUP = rootUP,
              QID = QID, CY = DCY, fertilizer = fertilizer, invest = invest, HD = HD, country = country, crop = crop)$par
  
  
  
  if (all(FR == 0)) {
    return(data.frame(lat = lat, lon = lon, plDate= plDatell, N = 0, P = 0, K = 0, WLY = WLY_user, CurrentY = CY_user, TargetY = CY_user, TC = 0, NR = 0))
  }else {
    
    fertilizer$FR <- FR
    
    ## NPK rate for ha of land
    N <- as.vector(FR %*% fertilizer$N_cont)
    P <- as.vector(FR %*% fertilizer$P_cont)
    K <- as.vector(FR %*% fertilizer$K_cont)
    rec <- c(N, P, K)
    
    ## NPK rate for user land size
    NPK_user <- rec * areaHa
    
    ## TY for ha of land
    TY <- QUEFTS1_Pedotransfer(QID, rec, crop=crop)    # Yield possible at recommended NPK in kg/ha dry wt.
    
    if(country %in% c("NG", "TZ")){
      TY_user <- ((getRFY(HD = HD, RDY = TY, country = country)) / 1000) * areaHa
    } else if (country == "Rw"){
      TY_user <- TY * 0.003 * areaHa
    }else if (country == "GH"){
      TY_user <- ((getRFY(HD = HD, RDY = TY, country = "NG")) / 1000) * areaHa
    }
  
    
    ## reporting the recommended fertilizers
    Recomfr <- fertilizer[fertilizer$FR > 0,]
    Recomfr$FR <- round(Recomfr$FR * areaHa, digits=0)
    
    
    ## total cost per ha
    # TC <- (sum(round(FR,digits=0) * fertilizer$price)) * areaHa
    TC <- as.numeric(Recomfr$FR %*% Recomfr$price)
    
    ## net revenue on the users land size
    GR <- (TY_user - CY_user) * rootUP                      # Gross revenue given root up is for fresh wt ton/ha
    TC <- round_any(TC, 100)
    GR <- round_any(GR, 100)
    
    NR <- round(GR - TC, digits = 0)                                                # Net Revenue
    
    Recomfr_wide <- spread(Recomfr[, c('type', 'FR')], type, FR)
    
    d1 <- data.frame(lat = lat, lon = lon, plDate = plDatell, N = NPK_user[1], P = NPK_user[2], K = NPK_user[3],
                     WLY = WLY_user, CurrentY = CY_user, TargetY = TY_user, TC = TC, NR = NR)
    d2 <- cbind(d1, Recomfr_wide)
    row.names(d2) <- NULL
    if (d2$NR <= 0 | d2$TargetY <= d2$CurrentY) {
      fertinfo <- subset(d2, select = c(lat, lon, plDate= plDatell, N, P, K, WLY, CurrentY, TargetY, TC, NR))
      fertinfo$N <- fertinfo$P <- fertinfo$K <- fertinfo$TC <- fertinfo$NR <- 0
      fertinfo$TargetY <- fertinfo$CurrentY
      d2 <- fertinfo
    }
    return(d2)
  }
  
}



#'  Optimize the UREA, TSP and MOP needed to maximize the NR. x1, x2, x3 = Urea, MOP and Nafaka kg/ha.
optim_NR <- function(fertRate, rootUP, QID, CY, fertilizer, invest, HD, country, crop) {
  f_price <- fertilizer$price
  TC <- sum(fertRate * f_price)
  
  ## Kg of Urea, Kg of NPK151515, Kg of NPK201010, Kg of MOP
  
  N <- as.vector(fertRate %*% fertilizer$N_cont)
  P <- as.vector(fertRate %*% fertilizer$P_cont)
  K <- as.vector(fertRate %*% fertilizer$K_cont)
  
  rec <- c(N, P, K)
  
  TotalYield <- QUEFTS1_Pedotransfer(QID, rec, crop=crop)
  
  if(country %in% c("NG", "TZ")){
    AdditionalYield <- (getRFY(HD = HD, RDY = (TotalYield - CY), country = country)) / 1000 ## DM is converted to FW and then from KG/ha to ton/ha
  }else{
    AdditionalYield <- (TotalYield - CY)*3/1000
  }
 
  
  #AdditionalYield <- (TotalYield - CY)*0.003
  PriceYield <- AdditionalYield * rootUP
  NetRev <- PriceYield - TC
  if (!is.na(invest) & TC > invest) { NetRev <- NetRev - (invest - TC)^2 } #penalize NR if costs exceed investment cap
  return(NetRev)
}





### see if profit is > (0.18 * total cost) + total cost
NRabove18Cost <- function(ds, riskAtt) {
  
  #minimal required net revenue increase from fertilizer needed (taking into account risk attitude of user)
  # dNRmin <- ifelse(riskAtt == 0, 2.8, ifelse(riskAtt == 1, 2, 1.2))
  dNRmin <- ifelse(riskAtt == 0, 1.8, ifelse(riskAtt == 1, 1, 0.2))
  
  
  #if(ds$NR <= (ds$TC + (ds$TC * 0.18))){
  if (ds$NR < ds$TC * dNRmin) {
    fertRecom <- subset(ds, select = c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
    fertRecom$N <- fertRecom$P <- fertRecom$K <- fertRecom$TC <- fertRecom$NR <- 0
    fertRecom$TargetY <- fertRecom$CurrentY
    
    onlyFert <- subset(ds, select = -c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
    row.names(onlyFert) <- NULL
    for (j in 1:ncol(onlyFert)) {
      onlyFert[, j] <- 0
    }
    fertRecom <- cbind(fertRecom, onlyFert)
    ds <- fertRecom
  }
  row.names(ds) <- NULL
  return(ds)
}




#' after setting fertilizer recommendation <25 kg/ha Urea, MOP or Nafaka, target yield with the remaining recommended fertilizer is  re-estimated  and
#'  total cost, gross and net revenue are re calcuated.
#' @param rootUP cassava root price
#' @param zone
#' @param wdd has dry wt
#' @param rdd has fresh wt
#' @param fertilizer
#' @author Meklit
#' @export

Rerun_25kgKa_try <- function(rootUP, rdd, fertilizer, QID, onlyFert, country, WLY = WLY, DCY = DCY, HD = HD, areaHa = areaHa) {
  
  
  QID$water_limited_yield <- WLY
  fertilizer <- merge(fertilizer, onlyFert, by = 'type')
  TC <- (sum(fertilizer$price %*% fertilizer$rate)) * areaHa
  TC <- round_any(TC, 100)
  N <- as.vector(fertilizer$rate %*% fertilizer$N_cont)
  P <- as.vector(fertilizer$rate %*% fertilizer$P_cont)
  K <- as.vector(fertilizer$rate %*% fertilizer$K_cont)
  rec <- c(N, P, K)
  
  ## NPK rate for user land size
  NPK_user <- rec * areaHa
  
  TY <- QUEFTS1_Pedotransfer(QID, rec, crop)                    #dry wt yield in kg/ha
  #TY_user  <- ((getRFY(HD = as.Date(HD), RDY = TY, country = country))/1000) * areaHa
  if(country %in% c("NG", "TZ")){
    TY_user <- ((getRFY(HD = HD, RDY = TY, country = country)) / 1000) * areaHa
    CY_user <- ((getRFY(HD = HD, RDY = DCY, country = country)) / 1000) * areaHa
  }else if (country == "Rw"){
    TY_user <- (TY * 3 / 1000) * areaHa
    CY_user <- (DCY * 3 /1000) * areaHa
  }else if (country == "GH"){
    TY_user <- ((getRFY(HD = HD, RDY = TY, country = "NG")) / 1000) * areaHa
    CY_user <- ((getRFY(HD = HD, RDY = DCY, country = "NG")) / 1000) * areaHa
  }
  
  
  
  rdd$CurrentY <- CY_user
  rdd$TargetY <- TY_user
  rdd$TC <- TC
  nr <- (rdd$TargetY - rdd$CurrentY) * rootUP
  nr <- round_any(nr, 100)
  rdd$NR <- nr - rdd$TC
  rdd$N <- NPK_user[1]
  rdd$P <- NPK_user[2]
  rdd$K <- NPK_user[3]
  
  if (rdd$TargetY <= rdd$CurrentY) {
    rdd$N <- rdd$P <- rdd$K <- rdd$TC <- rdd$NR <- 0
    rdd$TargetY <- CY_user
  }
  
  if (rdd$NR <= 0 | rdd$TargetY <= rdd$CurrentY) {
    fertinfo <- subset(rdd, select = c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
    fertinfo$N <- fertinfo$P <- fertinfo$K <- fertinfo$TC <- fertinfo$NR <- 0
    fertinfo$TargetY <- fertinfo$CurrentY
    rdd <- fertinfo
  }
  
  return(rdd)
}




## 
getFRrecommendations <- function(WLY_CYll, SoilData_ll, maxInv, fertilizers, rootUP, areaHa, country, FCY, riskAtt, crop) {
  
  lat <- WLY_CYll$lat
  lon <- WLY_CYll$long
  ## 1. get WLY, CY, fert recom and soil data
  WLY <- WLY_CYll$water_limited_yield ## DM in kg/ha
  DCY <- WLY_CYll$CurrentYield ## DM in kg/ha
  plDatell <- WLY_CYll$pl_Date
  
  
  ## 2. change investment from given areaHa to 1ha
  invest <- (maxInv / areaHa)
  
  ## for NG and TZ to convert dry to fresh matter harvest daty as the day of teh year is required
  WLY_CYll$HD <- WLY_CYll$daysOnField - WLY_CYll$pl_Date
  
  if(WLY_CYll$HD > 365){
    WLY_CYll$HD <- WLY_CYll$HD - 365
  }
  
  HD <- WLY_CYll$HD
  
  ## 3. optimize the fertilizer recommendation for maxInv in local currency and provide expected target yield in kg
  fert_optim <- run_Optim(rootUP = rootUP, QID = SoilData_ll, fertilizer = fertilizers, invest = invest, plDatell,
                          lat = lat, lon = lon, areaHa, HD = HD, DCY = DCY, WLY = WLY, country = country, crop=crop)
  
  if (fert_optim$NR == 0) { ## no fertilizer recommendation
    fertilizer_rates <- NULL
    return(list(rec = fert_optim, fertilizer_rates = fertilizer_rates))
  }else {
    fertinfo <- subset(fert_optim, select = c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
    onlyFert <- subset(fert_optim, select = -c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
    
    ## 4. remove ferilizer application < 25 kg/ha and re run the TY and NR calculation
    RecomperHa <- onlyFert / areaHa
    RecomperHa2 <- gather(RecomperHa, type, rate)
    onlyFert2 <- droplevels(RecomperHa2[RecomperHa2$rate > 25,])
    
    if (nrow(onlyFert2) == 0) { ## if all fertilizer recom < 25 kg/ha all will be set to 0
      fertinfo$N <- fertinfo$P <- fertinfo$K <- fertinfo$NR <- fertinfo$TC <- 0
      fertinfo$TargetY <- fertinfo$CurrentY
      fertilizer_rates <- NULL
      return(list(rec = fertinfo, fertilizer_rates = fertilizer_rates))
    }else if (ncol(onlyFert) == nrow(onlyFert2)) { ## if all fertilizer recom are >= 25 kg/ha they will be kept and only checked for NR >= 18% of invest
      Reset_fert_Cont <- fert_optim
      GPS_fertRecom <- NRabove18Cost(ds = Reset_fert_Cont, riskAtt=riskAtt)
      rec <- subset(GPS_fertRecom, select = c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
      frates <- subset(GPS_fertRecom, select = -c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
      frates2 <- gather(frates, type, rate)
      return(list(rec = rec, fertilizer_rates = frates2))
      
    }else {
      fert25 <- spread(onlyFert2, type, rate) ## when some fertilizer recom are dropped b/c < 25 kg/ha, ty and NR should be recalculated
      fert_optim2 <- cbind(fertinfo, fert25)
      fertilizer <- fertilizers[fertilizers$type %in% onlyFert2$type,]
      Reset_fert_Cont <- Rerun_25kgKa_try(rootUP = rootUP, rdd = fert_optim2, fertilizer = fertilizer, QID = SoilData_ll,
                                          onlyFert = onlyFert2,
                                          country = country, WLY = WLY, DCY = DCY, HD = HD, areaHa = areaHa)
      if (Reset_fert_Cont$NR <= 0) { ## after rerunning after avoiding <25KG/ha fertilizers, if NR <=0
        fertilizer_rates <- NULL
        return(list(rec = Reset_fert_Cont, fertilizer_rates = fertilizer_rates))
      }else {
        GPS_fertRecom <- NRabove18Cost(ds = Reset_fert_Cont,riskAtt=riskAtt)
        rec <- subset(GPS_fertRecom, select = c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
        frates <- subset(GPS_fertRecom, select = -c(lat, lon, plDate, N, P, K, WLY, CurrentY, TargetY, TC, NR))
        frates2 <- gather(frates, type, rate)
        return(list(rec = rec, fertilizer_rates = frates2))
      }
    }
  }
}




#' @description  generting FR recom for every loc, and pl date while running the harvest dates in parlles with several cores.
FR_parallel_byLoc <- function( WLY_CY_FCY_ll, SoilData_ll, maxInv = maxInv, fertilizers=fertilizers, rootUP = rootUP, areaHa=1, country=country,
                               FCY = FCY, riskAtt = riskAtt, crop = crop){
  
  WLY_CY_FCY_ll$water_limited_yield <- round( WLY_CY_FCY_ll$water_limited_yield, digits = 0)
  WLY_CY_FCY_ll$CurrentYield <- round( WLY_CY_FCY_ll$CurrentYield, digits = 0)
  
  
  ## for every planting date
  plDate <- unique( WLY_CY_FCY_ll$pl_Date)
  
  ## For every harvest day: harvest is possible weekly between 8 and 15 months after planting 
  daysOnField <- unique( WLY_CY_FCY_ll$daysOnField )
  
  PlDatr_FR_Recom <- NULL 
  for(pld in plDate){
    print(pld)
    WLY_CY_FCY_ll_pd <-  WLY_CY_FCY_ll[ WLY_CY_FCY_ll$pl_Date == pld, ]
    
    FRrecFCY5 <- parallel::parLapply(parallelCluster, daysOnField,
                                     try(QUEFTS_Parallel(WLY_CY_FCY_ll_pd, SoilData_ll=SoilData_ll, maxInv, fertilizers, rootUP, areaHa, country, FCY, riskAtt, crop)))
    
    FRrecFCY5_df <- ldply(FRrecFCY5, data.frame)
    FRrecFCY5_df$NAME_1 <- unique(SoilData_ll$NAME_1)
    FRrecFCY5_df$NAME_2 <- unique(SoilData_ll$NAME_2)
    PlDatr_FR_Recom <- rbind(PlDatr_FR_Recom, FRrecFCY5_df)
  }
  return(PlDatr_FR_Recom)

}



## in order to run QUEFTS in parlled with defined number of cores. 
QUEFTS_Parallel <- function(WLY_CY_FCY_ll_pd, SoilData_ll, maxInv, fertilizers, rootUP, areaHa, country, FCY, riskAtt, crop){
  force(WLY_CY_FCY_ll_pd)
  force(SoilData_ll)
  force(maxInv)
  force(fertilizers)
  force(rootUP)
  force(areaHa)
  force(country)
  force(FCY)
  force(riskAtt)
  force(crop)
  
  require(plyr)
  
  # source("QUEFTS_function.R", local=TRUE)
  
  worker <- function(hd) {
    FR_forHD(hd, WLY_CY_FCY_ll_pd, SoilData_ll, maxInv, fertilizers, rootUP, areaHa, country, FCY, riskAtt, crop)
  }
  return(worker)
  
}


## is to be used to parallelize across the harvest dates after filtering for lat, lon and pl_Date
FR_forHD <- function(hd, WLY_CY_FCY_ll_pd, SoilData_ll, maxInv, fertilizers, rootUP, areaHa, country, FCY, riskAtt, crop){
  require(plyr)
  require(tidyr)
  WLY_CY_FCY_ll_pd_hd <- WLY_CY_FCY_ll_pd[WLY_CY_FCY_ll_pd$daysOnField == hd, ]
  
  FRadv <- getFRrecommendations(WLY_CY=WLY_CY_FCY_ll_pd_hd, SoilData_ll=SoilData_ll, maxInv=maxInv, 
                                fertilizers=fertilizers, rootUP=rootUP, areaHa=areaHa, country=country, 
                                riskAtt = riskAtt, crop = crop)
  
  FRadv_A <- FRadv$rec
  FRadv_B <- FRadv$fertilizer_rates
  
  if(is.null(FRadv_B)){
    if(country=="Rw"){
      FRadv_A$Urea <- FRadv_A$MOP <- FRadv_A$DAP <- FRadv_A$NPK17_17_17 <- 0
      FRadv_A$FCY <- FCY
    }else if (country == "GH"){
      FRadv_A$Urea <- FRadv_A$NPK112221 <- FRadv_A$NPK251010 <- FRadv_A$NPK152020 <- FRadv_A$NPK201010  <- FRadv_A$NPK23105  <- FRadv_A$NPK123017  <- 0
      FRadv_A$FCY <- FCY
    }else if (country == "Bu"){
      FRadv_A$FOMI_BAGARA <- FRadv_A$FOMI_IMBURA <- FRadv_A$ FOMI_TOTAHAZA <- 0
      FRadv_A$FCY <- FCY
    }
    
  }else{
    FRadv_B <- spread(FRadv_B, type ,rate)
    
    if(country == "Rw"){
      FRadv_A$Urea <- ifelse("Urea" %in% colnames(FRadv_B), FRadv_B$Urea, 0)
      FRadv_A$MOP <- ifelse("MOP" %in% colnames(FRadv_B), FRadv_B$MOP, 0)
      FRadv_A$DAP <- ifelse("DAP" %in% colnames(FRadv_B), FRadv_B$DAP, 0)   
      FRadv_A$NPK17_17_17 <- ifelse("NPK17_17_17" %in% colnames(FRadv_B), FRadv_B$NPK17_17_17, 0)
    }else if(country == "GH"){
      FRadv_A$Urea <- ifelse("Urea" %in% colnames(FRadv_B), FRadv_B$Urea, 0)
      FRadv_A$NPK112221 <- ifelse("NPK112221" %in% colnames(FRadv_B), FRadv_B$NPK112221, 0)
      FRadv_A$NPK251010 <- ifelse("NPK251010" %in% colnames(FRadv_B), FRadv_B$NPK251010, 0)
      FRadv_A$NPK152020 <- ifelse("NPK152020" %in% colnames(FRadv_B), FRadv_B$NPK152020, 0)
      FRadv_A$NPK201010 <- ifelse("NPK201010" %in% colnames(FRadv_B), FRadv_B$NPK201010, 0)
      FRadv_A$NPK23105 <- ifelse("NPK23105" %in% colnames(FRadv_B), FRadv_B$NPK23105, 0)
      FRadv_A$NPK123017 <- ifelse("NPK123017" %in% colnames(FRadv_B), FRadv_B$NPK123017, 0)
      
    }else  if(country == "Bu"){
      FRadv_A$FOMI_BAGARA <- ifelse("FOMI_BAGARA" %in% colnames(FRadv_B), FRadv_B$FOMI_BAGARA, 0)
      FRadv_A$FOMI_IMBURA <- ifelse("FOMI_IMBURA" %in% colnames(FRadv_B), FRadv_B$FOMI_IMBURA, 0)
      FRadv_A$FOMI_TOTAHAZA <- ifelse("FOMI_TOTAHAZA" %in% colnames(FRadv_B), FRadv_B$FOMI_TOTAHAZA, 0) 
    }
    FRadv_A$FCY <- FCY
  }
  
  FRadv_A$weekNr <- WLY_CY_FCY_ll_pd_hd$weekNr
  FRadv_A$daysOnField <- WLY_CY_FCY_ll_pd_hd$daysOnField
  
  return(FRadv_A)
}








#' @description get soil NPK based on the RF model developed using the NOT location QUEFTS generated soil NPK
#'
#' @param trainData is data with soil NPK of teh NOT locations
#' @param AIOdata is area of intersts iSDA soil data
#' @param FCY 
#'
#' @return
#' @export
#'
#' @examples
get_SoilNPK_AOI <- function(trainData, AIOdata, FCY, country){
  AIOdata2 <- AIOdata[, colnames(trainData)]
  trainData$purpose <- "Traindata"
  AIOdata2$purpose <- "Testdata"
  
  trainData$index <- c(1:nrow(trainData))
  AIOdata2$index <- c(1:nrow(AIOdata2))
  AIOdata$index <- c(1:nrow(AIOdata))
  
  TT <- rbind(trainData, AIOdata2)
  
  TT <-TT[!TT$slope_angle == "No data", ]
  TT$slope_angle <- as.numeric(as.character(TT$slope_angle))
 
  
  trainData_FCY <- TT[TT$purpose == "Traindata", ]
  getData_FCY <- TT[TT$purpose == "Testdata", ]
  
  
  
  
  ## link the GPS for the test data 
  getData_FCY <- unique(merge(getData_FCY, AIOdata[, c("lat", "lon", "NAME_1", "NAME_2", "index")], by="index"))
  

  Ndata_Train_FCY <- subset(trainData_FCY, select=-c(soilP, soilK, purpose, index))
  Pdata_Train_FCY <- subset(trainData_FCY, select=-c(soilN, soilK, purpose, index))
  Kdata_Train_FCY <- subset(trainData_FCY, select=-c(soilN, soilP, purpose, index))
  
  Ndata_Valid_FCY <- subset(getData_FCY, select=-c(soilP, soilK, purpose, index))
  Pdata_Valid_FCY <- subset(getData_FCY, select=-c(soilN, soilK, purpose, index))
  Kdata_Valid_FCY <- subset(getData_FCY, select=-c(soilN, soilP, purpose, index))
  
 
  ## Coustome control parameter 
  #custom <- trainControl(method="repeatedcv", number=10, repeats=5, verboseIter=TRUE)
  custom <- trainControl(method="oob", number=10)
  
  ##soil N
  set.seed(773)
  RF_N1 <- randomForest(soilN ~ ., Ndata_Train_FCY , importance=TRUE, ntree=1000)
  Ndata_Valid_FCY$soilN <- predict(RF_N1, Ndata_Valid_FCY)
  summary(Ndata_Valid_FCY$soilN)
  
  ###soil P
  set.seed(773)
  RF_P1 <- randomForest(soilP ~ ., Pdata_Train_FCY , importance=TRUE, ntree=1000)
  Pdata_Valid_FCY$soilP <- predict(RF_P1, Pdata_Valid_FCY)
  summary(Pdata_Valid_FCY$soilP)
  
  ## soil K
  set.seed(773)
  RF_K1 <- randomForest(soilK ~ ., Kdata_Train_FCY , importance=TRUE, ntree=1000)
  Kdata_Valid_FCY$soilK <- predict(RF_K1, Kdata_Valid_FCY)
  summary(Kdata_Valid_FCY$soilK)
  
  
  FCY_soilNP <- merge(Ndata_Valid_FCY, Pdata_Valid_FCY[c("lat", "lon", "NAME_1", "NAME_2", "soilP")], by=c("lat", "lon", "NAME_1", "NAME_2"))
  FCY_soilNPK <- merge(FCY_soilNP, Kdata_Valid_FCY[c("lat", "lon", "NAME_1", "NAME_2", "soilK")], by=c("lat", "lon", "NAME_1", "NAME_2"))
  FCY_soilNPK$FCY <- FCY
  return(FCY_soilNPK)
  
}



