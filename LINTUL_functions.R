
list.of.packages <- c("readr", "dplyr", "deSolve", "plyr", "lubridate", "parallel", "tidyr", "doParallel", "foreach", "gtools") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(gtools)
require(deSolve)   
require(plyr)
require(lubridate) 
require(parallel)
require(tidyr)
require(doParallel)
require(foreach)


#This file contains functions: 
#
# LINTUL2,                    containing the LINTUL code 
# LINTUL2_DEFAUL_PARAMETERS,  containing all default settings
# LINTUL2_iniSTATES,          providing intial state values
# PENMAN,                     computation of the PENMAN EQUATION 
# EVAPTR,                     to compute actual rates of evaporation and transpiration
# GLA,                        to computes daily increase of leaf area index  
# DRUNIR,                     to compute rates of drainage, runoff and irrigation
# plot_years,                 function for plotting multiple years (max 15) in one graph
# get_weather,                function to read CABO weather file
#
# Part of this LINTUL-R code was developed by Mink Zijlstra and adapted to work with the deSolve package.
# This code is based on the original FST based LINTUL-Cassava developed by Ezui KS et Leffelar PA, 2016
# Wageningen University (Tom Schut & Rob van Beuken), 2017
#State <- LINTUL2_CASSAVA_iniSTATES_MTC(param_noIrrigation)
#Time <- seq(STTIME, FINTIM, by = DELT)
#Pars <-  param_noIrrigation

#--------------------------------------Lintul---------------------------------------#
#                                                                                   # 
#-----------------------------------------------------------------------------------#
LINTUL2_CASSAVA <-function(Time, State, Pars, WDATA){
  #intergration with delay for RROOTD
  with(as.list(c(State, Pars)), {
    
    RTRAIN <- WDATA$RAIN[Time]                       # rain rate, mm d-1
    DTEFF  <- max(0, WDATA$DAVTMP[Time] - TBASE)     # effective daily temperature (for crop development a treshold temperature (TBASE) needs to be exceeded)
    RPAR   <- FPAR * WDATA$DTR[Time]                 # PAR MJ m-2 d-1
    
    #determine rates when crop is still growing
    if(TSUM < FINTSUM){
      # Determine water content of rooted soil
      WC  <- 0.001 * WA/ROOTD
      # RTSUM <- DTEFF*ifelse((Time-DOYPL) >= 0, 1, 0) # TSUM after planting
      RTSUM <- DTEFF*ifelse((Time-STTIME) >= 0, 1, 0) # TSUM after planting
      
      # Once the emergence date is reached and enough water is available the crop emerges (1), once the crop is established is does not disappear again (2)
      #     if((WC-WCWP) >= 0 && (TSUM-OPTEMERGTSUM) >= 0 && (TSUM-SPROUTFAIL <=0)) { 	# (1) ## for sprouting correction
      if((WC-WCWP) >= 0 && (TSUM-OPTEMERGTSUM) >= 0) { 	# (1)
        emerg1 <- 1} else { emerg1 <- 0 }
      if(TSUMCROP > 0) {									# (2)
        emerg2 <- 1 } else { emerg2 <- 0}
      # Emergence of the crop is used to calculate the accumulated temperature sum.
      EMERG  <- max(emerg1,emerg2)
      RTSUMCROP <- DTEFF*EMERG
      
      # Root depth
      # If soil water content drops to, or below, wilting point fibrous root growth stops.
      # As long as the crop has not reached its maximum rooting depth and has not started flowering yet, fibrous root growth continues.
      if((ROOTD-ROOTDM) < 0 && (WC-WCWP) >= 0) {
        # The rooting depth (m) is calculated from a maximum rate of change in rooting depth, the emergence of the crop and the constraints mentioned above.
        RROOTD <- RRDMAX * EMERG 
      }else{ 
        RROOTD = 0
      }
      
      EXPLOR <- 1000 * RROOTD * WCFC                   # exploration rate of new soil water layers by root depth growth (mm d-1)
      
      
      RNINTC <- min(RTRAIN, (FRACRNINTC * LAI))                # interception of rain by the canopy (mm d-1)
      
      # Potential evaporation (mm d-1) and transpiration (mm d-1) are calculated according to Penman-Monteith
      PENM   <- penman(WDATA$DAVTMP[Time],WDATA$VP[Time],WDATA$DTR[Time],LAI,WDATA$WN[Time],RNINTC)
      RPTRAN <- PENM$PTRAN
      RPEVAP <- PENM$PEVAP
      
      WCSD <- WCWP * TWCSD
      WCCR <- WCWP + pmax(WCSD-WCWP, (PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP))/1000)
      
      # Actual evaporation (mm d-1) and transpiration (mm d-1)
      EVA  <- evaptr(RPEVAP,RPTRAN,ROOTD,WA,WCAD,WCWP,TWCSD,WCFC,WCWET,WCST,TRANCO,DELT)
      RTRAN <- EVA$TRAN
      REVAP <- EVA$EVAP
      
      # The transpiration reduction factor is defined as the ratio between actual and potential transpiration
      TRANRF = ifelse(RPTRAN <= 0, 1, RTRAN/RPTRAN)
      
      # Drainage (below the root zone; mm d-1), surface water runoff (mm d-1) and irrigation rate (mm d-1)
      DRUNIR    <- drunir(RTRAIN,RNINTC,REVAP,RTRAN,IRRIGF,DRATE,DELT,WA,ROOTD,WCFC,WCST)
      
      # Rate of change of soil water amount (mm d-1)
      RWA <- (RTRAIN + EXPLOR + DRUNIR$IRRIG) - (RNINTC + DRUNIR$RUNOFF + RTRAN + REVAP + DRUNIR$DRAIN)
      WC <- 0.001 * WA/ROOTD
      
      # Light interception (MJ m-2 d-1) and total crop growth rate (g m-2 d-1)
      PARINT <- RPAR * (1 - exp(-K_EXT * LAI))
      LUE    <- LUE_OPT * approx(TTB[,1], TTB[,2], WDATA$DAVTMP[Time])$y
      
      # Dormancy and recovery from dormancy
      if ((WC-WCSD) <= 0 && (LAI - LAI_MIN) <= 0){
        dormancy = 1
      } else {
        dormancy = 0
      }
      
      if ((WC - RECOV * WCCR) >= 0 && (WC - WCWP) >= 0){
        pushdor = 1
      } else {
        pushdor = 0
      }
      
      if (WSO == 0) {
        WSOREDISTFRAC = 1
      } else {
        WSOREDISTFRAC = REDISTSO/WSO
      }
      
      PUSHREDISTEND = max(ifelse((WSOREDISTFRAC-WSOREDISTFRACMAX) >= 0, 1, 0), ifelse((REDISTLVG - WLVGNEWN)>= 0, 1, 0), ifelse((PUSHREDISTSUM - TSUMREDISTMAX) >= 0, 1, 0)) *ifelse(-PUSHREDISTSUM >= 0, 0, 1)
      PUSHREDIST = ifelse((PUSHDORMRECTSUM - DELREDIST) >= 0, 1, 0)* (1 - PUSHREDISTEND)
      PUSHDORMREC = pushdor*ifelse(-DORMTSUM >= 0, 0, 1) * (1 - PUSHREDIST) * ifelse((TSUMCROP - TSUMSBR) >= 0, 1, 0)
      DORMANCY = max(dormancy, PUSHDORMREC) * (1 - PUSHREDIST) * ifelse((TSUMCROP - TSUMSBR) >= 0, 1, 0)
      RDORMTSUM = DTEFF *DORMANCY - (DORMTSUM/DELT) * PUSHREDIST
      RPUSHDORMRECTSUM = DTEFF * PUSHDORMREC - (PUSHDORMRECTSUM/DELT) * (1 - PUSHDORMREC) * (1 - PUSHREDIST)
      RPUSHREDISTSUM = DTEFF * PUSHREDIST - (PUSHREDISTSUM/DELT) * PUSHREDISTEND
      RPUSHREDISTENDTSUM = DTEFF * PUSHREDIST - (PUSHREDISTENDTSUM/DELT) * (1 - PUSHREDISTEND)
      RDORMTIME = DORMANCY
      
      # Dry matter redistribution after dormancy
      RREDISTSO = RRREDISTSO * WSO * PUSHREDIST - (REDISTSO/DELT) *ifelse(-DORMTSUM >= 0, 0, 1)
      RREDISTLVG = SO2LV * RREDISTSO * (1- DORMANCY)
      RREDISTMAINTLOSS = (1 - SO2LV) * RREDISTSO
      
      GTOTAL <- LUE * PARINT * TRANRF * (1 - DORMANCY)
      
      # Relative death rate (d-1) due to aging
      RDRDV = ifelse(TSUMCROPLEAFAGE - TSUMLLIFE >= 0, approx(RDRT[,1], RDRT[,2], WDATA$DAVTMP[Time])$y, 0)
      
      # Relative death rate (d-1) due to self shading
      RDRSH <- RDRSHM * (LAI-LAICR) / LAICR
      if(RDRSH < 0) {
        RDRSH <- 0
      } else if(RDRSH >=RDRSHM) {
        RDRSH <- RDRSHM
      }
      
      RTSUMCROPLEAFAGE <- DTEFF * EMERG - (TSUMCROPLEAFAGE/DELT) * PUSHREDIST
      ENHSHED <- max(ifelse((WC-WCSD) >= 0, 0, 1), ifelse((WC-WCWET) >= 0, 1, 0))*ifelse((TSUMCROPLEAFAGE-FRACTLLFENHSH*TSUMLLIFE) >= 0, 1, 0)
      
      # Relative death rate (d-1) due to severe drought
      RDRSD <- RDRB *ENHSHED
      
      # Effective relative death rate (1; d-1) and the resulting decrease in LAI (2; m2 m-2 d-1) and leaf weight (3; g m-2 d-1)
      RDR   <- max(RDRDV, RDRSH, RDRSD) * ifelse((TSUMCROPLEAFAGE - TSUMLLIFE) >= 0, 1, 0) 	# (1)
      DLAI  <- LAI * RDR * (1 - FASTRANSLSO) * (1 - DORMANCY)  			  # (2)
      
      # Allocation to roots (2), leaves (4), stems (5) and storage organs (6)
      # fractions allocated are modified for water availability (1 and 3)
      FRTMOD <- max(1, 1/(TRANRF+0.5))						    # (1)
      FRT    <- approx(FRTTB[,1],FRTTB[,2],TSUMCROP)$y * FRTMOD		# (2)
      FSHMOD <- (1 - FRT) / (1 - FRT / FRTMOD)				    # (3)
      FLV    <- approx(FLVTB[,1],FLVTB[,2],TSUMCROP)$y * FSHMOD		# (4)
      FST    <- approx(FSTTB[,1],FSTTB[,2],TSUMCROP)$y * FSHMOD		# (5)
      FSO    <- approx(FSOTB[,1],FSOTB[,2],TSUMCROP)$y * FSHMOD		# (6)
      
      # stem cutting growth
      #      WCUTTINGMIN <- WCUTTINGMINPRO * WCUTTINGIP
      WCUTTINGIP <- NCUTTINGS * WCUTTINGUNIT # Estimated here but was a parameter in the original version
      WCUTTINGMIN <- WCUTTINGMINPRO * WCUTTINGIP
      
      
      # Leaf growth and senescence
      FRACSLACROPAGE <- approx(FRACSLATB[,1], FRACSLATB[,2], TSUMCROP)$y
      SLA <- SLA_MAX *FRACSLACROPAGE
      RWSOFASTRANSLSO = WLVG * RDR * FASTRANSLSO * (1 - DORMANCY)
      DLV <- (WLVG * RDR - RWSOFASTRANSLSO) * (1 - DORMANCY)
      RWLVD <- DLV
      
      # Change in biomass (g m-2 d-1) for green leaves (1), stems (2), storage organs (3) and roots (4)
      if (TSUM > OPTEMERGTSUM && WST == 0){
        RWCUTTING <- WCUTTING *(-FST_CUTT - FRT_CUTT - FLV_CUTT - FSO_CUTT)
        RWST <- WCUTTINGIP * FST_CUTT
        RWRT <- WCUTTINGIP * FRT_CUTT
        RWLVG <- WCUTTINGIP * FLV_CUTT
        RWSO <- WCUTTINGIP * FSO_CUTT
      } else if (TSUMCROP > 16.8) {   # A little bit of cheating
        RWCUTTING <- -RDRWCUTTING * WCUTTING * ifelse((WCUTTING-WCUTTINGMIN) >= 0, 1, 0) * TRANRF * EMERG * (1 - DORMANCY)
        RWLVG  <- (abs(GTOTAL)+abs(RWCUTTING)) * FLV - DLV + RREDISTLVG * PUSHREDIST 	# (1)
        RWST   <- (abs(GTOTAL)+abs(RWCUTTING)) * FST			  # (2)
        RWSO   <- (abs(GTOTAL)+abs(RWCUTTING)) * FSO + RWSOFASTRANSLSO - RREDISTSO			  # (3)
        RWRT   <- (abs(GTOTAL)+abs(RWCUTTING)) * FRT			  # (4)
      } else{
        RWCUTTING <- 0
        RWLVG <- 0
        RWST <- 0
        RWSO <- 0
        RWRT <- 0
      }
      RWLV = RWLVG+RWLVD
      WGTOTAL = WLV+WST+WCUTTING+WSO+WRT
      GLV <- FLV * (GTOTAL + abs(RWCUTTING)) + RREDISTLVG * PUSHREDIST
      
      GLAI <- gla(DTEFF,TSUMCROP,LAII,RGRL,DELT,SLA,LAI,GLV,TSUMLA_MIN,TRANRF,WC,WCWP,RWCUTTING,FLV,LAIEXPOEND,DORMANCY)
      
      # Change in LAI (m2 m-2 d-1) due to new growth of leaves
      RLAI <- GLAI - DLAI
      
      
    }else{
      #all plant related rates are set to 0
      RROOTD <- 0
      RWLVG <- 0
      RWA <- 0
      RTSUMCROP <- 0
      RLAI <- 0
      RWLVG <- 0
      RWLVD <- 0
      RWLV <- 0
      RWST <- 0
      RWSO <- 0
      RWRT <- 0
      RTRAN <- 0 
      REVAP <- 0 
      RPTRAN <- 0 
      RPEVAP <- 0
      RGTOTAL <- 0
    }
    return(list(c(RROOTD,RWA,RTSUM, RTSUMCROP, RTSUMCROPLEAFAGE,RDORMTSUM,RPUSHDORMRECTSUM,RPUSHREDISTENDTSUM, RDORMTIME, RWCUTTING, RTRAIN,RPAR,RLAI,RWLVD, RWLV,RWST,RWSO,RWRT, RWLVG,RTRAN,REVAP,RPTRAN,RPEVAP,RREDISTLVG,RREDISTSO,RPUSHREDISTSUM,RWSOFASTRANSLSO)))
    
  })
}

#-----------------------------Initial conditions------------------------------------#
#                                                                                   #
# A number of variables need to be given an initial value. These variables are used #
# by the model before they are calculated during the first time step. During all    #
# subsequent time steps the updated values will be used.                            # 
#                                                                                   #
#-----------------------------------------------------------------------------------#
LINTUL2_CASSAVA_iniSTATES <-function(){
  parvalues <- LINTUL2_CASSAVA_parameters()
  return( c(ROOTD   = parvalues[['ROOTDI']],                                               # m     :    initial rooting depth (at crop emergence)
            WA      = 1000 * parvalues[['ROOTDI']] * parvalues[['WCI']],                                # mm    :    initial soil water amount
            TSUM    = 0,                                                 # deg. C:    temperature sum 
            TSUMCROP =0,
            TSUMCROPLEAFAGE =0,
            DORMTSUM = 0,
            PUSHDORMRECTSUM = 0,
            PUSHREDISTENDTSUM = 0,
            DORMTIME = 0,
            WCUTTING = parvalues[['WCUTTINGUNIT']] * parvalues[['NCUTTINGS']], 
            TRAIN   = 0,                                                 # mm   :    rain sum
            PAR     = 0,                                                 # MJ m-2:    PAR sum
            LAI     = 0,                                                 # m2 m-2:    leaf area index 
            WLVD    = 0,                                                 # g m-2 :    dry weight of dead leaves
            WLV     = 0,
            WST     = 0,                                                 # g m-2 :    dry weight of stems
            WSO     = 0,                                                 # g m-2 :    dry weight of storage organs
            WRT     = 0,                                                 # g m-2 :    dry weight of roots
            WLVG    = 0,
            TRAN    = 0,                                                 # mm    :    actual transpiration
            EVAP    = 0,                                                # mm    :    actual evaporation
            PTRAN   = 0,                                                 # mm    :    potential transpiration
            PEVAP   = 0,                                                 # mm    :    potential evaporation
            REDISTLVG = 0,
            REDISTSO = 0,
            PUSHREDISTSUM = 0,
            WSOFASTRANSLSO =0
  ))
}

LINTUL2_CASSAVA_iniSTATES_MTC <-function(param_noIrrigation){
  parvalues <- param_noIrrigation
  return( c(ROOTD   = parvalues[['ROOTDI']],                                               # m     :    initial rooting depth (at crop emergence)
            WA      = 1000 * parvalues[['ROOTDI']] * parvalues[['WCI']],                                # mm    :    initial soil water amount
            TSUM    = 0,                                                 # deg. C:    temperature sum 
            TSUMCROP =0,
            TSUMCROPLEAFAGE =0,
            DORMTSUM = 0,
            PUSHDORMRECTSUM = 0,
            PUSHREDISTENDTSUM = 0,
            DORMTIME = 0,
            WCUTTING = parvalues[['WCUTTINGUNIT']] * parvalues[['NCUTTINGS']], 
            TRAIN   = 0,                                                 # mm   :    rain sum
            PAR     = 0,                                                 # MJ m-2:    PAR sum
            LAI     = 0,                                                 # m2 m-2:    leaf area index 
            WLVD    = 0,                                                 # g m-2 :    dry weight of dead leaves
            WLV     = 0,
            WST     = 0,                                                 # g m-2 :    dry weight of stems
            WSO     = 0,                                                 # g m-2 :    dry weight of storage organs
            WRT     = 0,                                                 # g m-2 :    dry weight of roots
            WLVG    = 0,
            TRAN    = 0,                                                 # mm    :    actual transpiration
            EVAP    = 0,                                                # mm    :    actual evaporation
            PTRAN   = 0,                                                 # mm    :    potential transpiration
            PEVAP   = 0,                                                 # mm    :    potential evaporation
            REDISTLVG = 0,
            REDISTSO = 0,
            PUSHREDISTSUM = 0,
            WSOFASTRANSLSO =0
  ))
}

#---------------------------------------------------------------------#
# FUNCTION LINTUL2_DEFAUL_PARAMETERS                                               #
# Purpose: Listing the default input parameters for Lintul2                   #
#---------------------------------------------------------------------#
LINTUL2_CASSAVA_DEFAUL_PARAMETERS <-function() {
  #All defaults
  return(c( 
    # LINTUL 2
    ROOTDI       = 0.1,    # M           :     initial rooting depth (at crop emergence)
    WCI          = 0.41,   # m3 m-3      :     initial soil water content 
    SLAI         = 0.017,  # m2 g-1      :     specific leaf area
    LAII         = 0.02603125,# m2 m-2
    ROOTDM       = 0.35,    # m           :     maximum rooting depth
    RRDMAX       = 0.012,  # m  d-1      :     max rate increase of rooting depth
    SLAI         = 0.017,
    WCAD         = 0.01,   # m3 m-3      :     soil water content at air dryness 
    WCWP         = 0.12,   # m3 m-3      :     soil water content at wilting point
    WCFC         = 0.41,   # m3 m-3      :     soil water content at field capacity 
    WCWET        = 0.46,   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    WCST         = 0.52,   # m3 m-3      :     soil water content at full saturation
    FRACRNINTC   = 0.25,   # (-)         :     fraction of rain intercepted
    TRANCO       = 8,      # mm d-1      :     transpiration constant (indicating level of drought tolerance)
    DRATE        = 50,     # mm d-1      :     max drainage rate
    IRRIGF       = 0,      # (-)         :     irrigation rate relative to the rate required to keep soil water status at field capacity
    TBASE        = 15,     # deg. C      :     base temperature 
    #TBASE        = 12,     # deg. C      :     base temperature 
    LUE          = 3.0,    # g MJ-1      :     light use efficiency
    K_EXT        = 0.67,   # (-)         :     extinction coefficient for PAR
    RGRL         = 0.003,  # 1/(deg. C d):     relative growth rate of LAI during exponential growth
    TSUMAN       = 1110,   # deg. C      :     temperature sum at anthesis
    FINTSUM      = 5460, #4320,   # deg. C      :     temperature sum at which simulation stops, 5400 is used to simulate until 15 MAP
    LAICR        = 3.5,     # m2 m-2     :     critical LAI beyond which leaf shedding is stimulated
    RDRSHM       = 0.09,   # d-1         :     max relative death rate of leaves due to shading
    FPAR         = 0.5,    # (-)         :     Fraction PAR, MJ PAR/MJ DTR
    
    
    # LINTUL CASSAVA
    DOYPL        = 113,     # d                        :    daynumber at crop planting. 
    TWCSD        = 1.05,   # (-)                      :    Proportion of WCWP at which water content at severe drought (WCSD) is assumed to be reached (5% above WCWP)
    SLA_MAX      = 0.03,   # m2 g-1                   :    maximum specific leaf area
    LAIEXPOEND   = 0.75,   # m2 m-2                   :    maximum leaf area index for exponential growth phase
    TSUMLA_MIN   = 180,    # deg. C d                 :    Temperature sum accumulation for minimum leaf area
    FRACTLLFENHSH= 0.85,   # (-)                      :    Fraction of leaf life at which enhanced shedding can be induced
    DELREDIST    = 12,     # deg. C d                 :    Delay for redistribution of dry matter
    TSUMRROOTEND = 720,    # deg. C                   :    development time from leaf appearance to leaf senescence
    RECOV        = 0.7,    # (-)                      :    Proportion of critical soil water content above which the crop recovers from drought
    # a. emergence
    OPTEMERGTSUM = 180,    # deg. c                   :    Optimal amount of Tsum accumulated from planting to emergence
    # b. simulation of first branching
    TSUMSBR      = 780,    # deg. C                   :    Temperature sum accumulation to first branching 
    # c. simulation of maturity and harvest
    # d. stem cutting growth
    FST_CUTT     = 0.03,   # gstem gcuttings-1        :    Fraction of initial cutting allocated to stem
    FRT_CUTT     = 0.03,   # groots gcuttings-1       :    Fraction of initial cutting allocated to adventitious roots
    FLV_CUTT     = 0.07,   # gleaves gcuttings-1      :    Fraction of initial cutting allocated to leaves
    FSO_CUTT     = 0,      # gstorage gcuttings-1     :    Fraction of initial cutting allocated to storage roots
    RDRWCUTTING  = 0.017,  # d-1                      :    Relative decrease rate of cutting weight
    WCUTTINGUNIT = 14,     # g                        :    Average weight per cutting 
    WCUTTINGMINPRO= 0.15,  # g m-2                    :    Minimum proportion of the stem cutting weight that remains unchanged
    #  NCUTTINGS    = 1.5625, # m-2                      :    Number of cuttings planted per m2
    # NCUTTINGS    = 1, # m-2                          :    Number of cuttings planted per m2; for TZ 
    NCUTTINGS    = 1.25, # m-2                        :    Number of cuttings planted per m2; for NG MC
    WCUTTINGIP   = 17.5, #21.875,
    # e. leaf area index growth
    # f. total biomass growth
    #LUE_OPT      = 1.5,    # g MJ PAR-1               :    Light use effeciency at optimum growing conditions
    LUE_OPT      = 1.7,    # g MJ PAR-1               :    Light use effeciency at optimum growing conditions
    # g. leaf, stem, fibrous roots and storage roots growth
    RRREDISTSO   = 0.01,    # d-1                      :    Relative rate of redistribution of dry matter from storage roots to leaves
    # h. senescence
    FASTRANSLSO  = 0.45,   # (-)                      :    Proportion of senesced leaf weight translocated to storage roots before shedding of the leaf 
    RDRB         = 0.09,   # d-1                      :    basic relative death rate of the leaves
    TSUMLLIFE    = 1200,   # deg. C d                 :    Temperature sum accumulation to leaf life
    # i. dry matter partitioning
    TSUMSOBULKINIT = 540,  # deg. C                   :    Temperature sum accumulation for start of storage roots bulking
    # j. dormancy and recovery from dormancy
    LAI_MIN      = 0.09,   # m2 m-2                   :    minimum leaf area index
    WSOREDISTFRACMAX = 0.05, # (-)                    :    Maximum proportion of dry matter redistribution from storage roots for the formation of new leaves
    WLVGNEWN     = 10,     # g DM m-2                 :    Minimum amount of new leaves weight produced in the redistribution phase
    TSUMREDISTMAX= 144,    # deg. C                   :    Maximum temperature sum accumulation to indicate the duration of dry matter redistribution
    # k. biomass production upon the recovery from dormancy
    SO2LV        = 0.8,     # g leaves DM g-1 storage  :    Converstion rate of storage organs dry matter to leaf dry matter
    DELT         = 1
    
  ))
}

#Look-up tables oC, FRACSLAB 
FRACSLATB <- matrix(c(0,  0.57, 
                      1440,  0.57, 
                      2880,  0.65, 
                      3864,  1, 
                      4320,  1,
                      5460,  1),ncol=2,byrow=TRUE) # (-)  : table of the fraction of SLA_MAX at different physiological times

#Look-up tables   oC,  RDR
RDRT <- matrix(c(-10,  0.02, 
                 10,  0.02, 
                 15,  0.03, 
                 30,  0.06, 
                 50,  0.06),ncol=2,byrow=TRUE) # (-)  :     table of RDR (Relative Death Rate) as function of temperature

#Look-up tables   oC,  TTB
TTB  <- matrix(c(-10,  0, 
                 15,  0, 
                 25,  1, 
                 29,  1,
                 40,  0, 
                 50,  0),ncol=2,byrow=TRUE) # (-)

#Partitioning    TSUM,  fraction of daily growth to roots ## MC: modified to allow 15 MAP
FRTTB <-matrix(c(   0, 0.11,   
                    540, 0.10,   
                    720, 0.094,   
                    900, 0.01,
                    1488, 0.01,  
                    1980, 0.01,  
                    2676, 0.01,  
                    3854, 0.01,  
                    5460, 0.01),ncol=2,byrow=TRUE) # (-)         :     table of FRT as a function of TSUM
#4320, 0.01),ncol=2,byrow=TRUE) # (-)         :     table of FRT as a function of TSUM

#Partitioning    TSUM,  fraction of daily growth to leaves ## MC: modified to allow 15 MAP
FLVTB <-matrix(c(   0, 0.71,   
                    540, 0.515,   
                    720, 0.393,   
                    900, 0.24,
                    1488, 0.21,  
                    1980, 0.18,
                    2676, 0.13,
                    3864, 0.21,
                    5460, 0.06),ncol=2,byrow=TRUE) # (-)         :     table of FLV as a function of TSUM
# 4320, 0.21),ncol=2,byrow=TRUE) # (-)         :     table of FLV as a function of TSUM

#Partitioning    TSUM,  fraction of daily growth to stem ## MC: modified to allow 15 MAP
FSTTB <-matrix(c(   0, 0.18,   
                    540, 0.385,   
                    720, 0.393,   
                    900, 0.26, 
                    1488, 0.26,  
                    1980, 0.19,  
                    2676, 0.29,
                    3864, 0.29,
                    5460, 0.31),ncol=2,byrow=TRUE) # (-)         :     table of FST as a function of TSUM
# 4320, 0.29),ncol=2,byrow=TRUE) # (-)         :     table of FST as a function of TSUM

#Partitioning    TSUM,  fraction of daily growth to storage organs ## MC: modified to allow 15 MAP
FSOTB <-matrix(c(   0, 0,
                    540, 0,  
                    720, 0.12,  
                    900, 0.49,  
                    1488, 0.52,
                    1980, 0.62,
                    2676, 0.57,
                    3864, 0.49, 
                    5460, 0.62),ncol=2,byrow=TRUE) # (-)         :     table of FSO as a function of TSUM
#      4320, 0.49),ncol=2,byrow=TRUE) # (-)         :     table of FSO as a function of TSUM

#---------------------------------------------------------------------#
# SUBROUTINE PENMAN                                                   #
# Purpose: Computation of the PENMAN EQUATION                         #
#---------------------------------------------------------------------#

penman <-function(DAVTMP,VP,DTR,LAI,WN,RNINTC) {
  
  DTRJM2 <-DTR * 1E6        # J m-2 d-1    :    Daily radiation in Joules 
  BOLTZM <-5.668E-8 	      # J m-1 s-1 K-4:    Stefan-Boltzmann constant 
  LHVAP  <-2.4E6            # J kg-1       :    Latent heat of vaporization 
  PSYCH  <-0.067            # kPa deg. C-1 :    Psychrometric constant
  
  BBRAD  <-BOLTZM * (DAVTMP+273)^4 * 86400                 # J m-2 d-1   :     Black body radiation 
  SVP    <-0.611 * exp(17.4 * DAVTMP / (DAVTMP + 239))     # kPa         :     Saturation vapour pressure
  SLOPE  <-4158.6 * SVP / (DAVTMP + 239)^2                 # kPa dec. C-1:     Change of SVP per degree C
  RLWN   <-BBRAD * pmax(0, 0.55 * (1 - VP / SVP))     # J m-2 d-1   :     Net outgoing long-wave radiation
  WDF    <-2.63 * (1.0 + 0.54 * WN)                  # kg m-2 d-1  :     Wind function in the Penman equation
  
  # Net radiation (J m-2 d-1) for soil (1) and crop (2)
  NRADS  <-DTRJM2 * (1 - 0.15) - RLWN     # (1)
  NRADC  <-DTRJM2 * (1 - 0.25) - RLWN     # (2)
  
  # Radiation terms (J m-2 d-1) of the Penman equation for soil (1) and crop (2)
  PENMRS <-NRADS * SLOPE / (SLOPE + PSYCH)    # (1)
  PENMRC <-NRADC * SLOPE / (SLOPE + PSYCH)    # (2)
  
  # Drying power term (J m-2 d-1) of the Penman equation
  PENMD  <-LHVAP * WDF * (SVP - VP) * PSYCH / (SLOPE + PSYCH)
  
  # Potential evaporation and transpiration are weighed by a factor representing the plant canopy (exp(-0.5 * LAI)).
  PEVAP  <-exp(-0.5 * LAI)  * (PENMRS + PENMD) / LHVAP
  PTRAN  <-(1 - exp(-0.5 * LAI)) * (PENMRC + PENMD) / LHVAP
  PTRAN  <-pmax(0, PTRAN - 0.5 * RNINTC) 
  
  
  PENM = data.frame(cbind(PEVAP,PTRAN))
  
  return(PENM)
}

#---------------------------------------------------------------------#
# SUBROUTINE EVAPTR                                                   #
# Purpose: To compute actual rates of evaporation and transpiration   #
#---------------------------------------------------------------------#

evaptr <-function(PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,TWCSD,WCFC,WCWET,WCST,TRANCO,DELT=1) {
  
  # Soil water content (m3 m-3) and the amount of soil water (mm) at air dryness (AD) and field capacity (FC).     
  WC   <-0.001 * WA / ROOTD
  WAAD <-1000 * WCAD * ROOTD
  WAFC <-1000 * WCFC * ROOTD
  
  # Evaporation is decreased when water content is below field capacity, but continues until WC = WCAD.
  limit.evap <-(WC-WCAD)/(WCFC-WCAD)
  #Ensure to stay within 0-1 range
  limit.evap <- pmin(1,pmax(0,limit.evap))
  EVAP <-PEVAP * limit.evap
  
  WCSD <- WCWP * TWCSD
  WCCR <- WCWP + pmax(WCSD-WCWP, PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP))
  
  # If water content is below the critical soil water content a correction factor is calculated that reduces 
  # transpiration until it stops at WC = WCWP.
  FR <-(WC-WCWP) / (WCCR - WCWP)
  
  # If water content is above the critical soil water content a correction factor is calculated that reduces 
  # transpiration when the crop is hampered by waterlogging (WC > WCWET).
  FRW <- (WCST-WC) / (WCST - WCWET)
  
  #Replace values for wet days with high WC values, above WCCR 
  FR[WC > WCCR] <- FRW[WC > WCCR]
  
  #Ensure to stay within the 0-1 range
  FR=pmin(1,pmax(0,FR))
  
  TRAN <-PTRAN * FR
  
  aux <- EVAP+TRAN
  aux[aux <= 0] <- 1
  
  # A final correction term is calculated to reduce evaporation and transpiration when evapotranspiration exceeds 
  # the amount of water in soil present in excess of air dryness.
  AVAILF <- pmin(1, (WA-WAAD)/(DELT*aux))
  EVA <- data.frame(EVAP = EVAP * AVAILF,
                    TRAN = TRAN * AVAILF)
  return(EVA)
}     

# ---------------------------------------------------------------------#
#  SUBROUTINE GLA                                                      #
#  Purpose: This subroutine computes daily increase of leaf area index #
#           (ha leaf/ ha ground/ d)                                    #
# ---------------------------------------------------------------------#

gla <-function(DTEFF,TSUMCROP,LAII,RGRL,DELT,SLA,LAI,GLV,TSUMLA_MIN,TRANRF,WC,WCWP,RWCUTTING,FLV,LAIEXPOEND,DORMANCY) {
  
  # Growth during maturation stage:
  GLAI <-SLA * GLV * (1-DORMANCY)
  
  # Growth during juvenile stage:
  if(TSUMCROP < TSUMLA_MIN && LAI < LAIEXPOEND) {
    GLAI <-((LAI * (exp(RGRL * DTEFF * DELT) - 1) / DELT) +abs(RWCUTTING)*FLV*SLA) * TRANRF
  }
  
  # Growth at day of seedling emergence:
  if(TSUMCROP > 0 && LAI == 0 && WC > WCWP) {
    GLAI <- LAII / DELT
  }
  
  # Growth before seedling emergence:
  if(TSUMCROP == 0) {
    GLAI <-0
  }
  
  return(GLAI)
}

# ---------------------------------------------------------------------#
#  SUBROUTINE DRUNIR                                                   #
#  Purpose: To compute rates of drainage, runoff and irrigation        #
# ---------------------------------------------------------------------#

drunir <-function(RAIN,RNINTC,EVAP,TRAN,IRRIGF,DRATE,DELT,WA,ROOTD,WCFC,WCST) {
  
  # Soil water content (m3 m-3) and the amount of soil water (mm) at field capacity (FC) and full saturation (ST).     
  WC   <-0.001 * WA / ROOTD
  WAFC <-1000 * WCFC * ROOTD
  WAST <-1000 * WCST * ROOTD
  
  # Drainage below the root zone occurs when the amount of water in the soil exceeds field capacity or when the amount of rainfall 
  # in excess of interception and evapotranspiration fills up soil water above field capacity.
  DRAIN <-(WA-WAFC)/DELT + (RAIN - (RNINTC + EVAP + TRAN))
  
  DRAIN <- ifelse(DRAIN < 0, 0, DRAIN)
  DRAIN <- ifelse(DRAIN > DRATE, DRATE, DRAIN)
  
  # Surface runoff occurs when the amount of soil water exceeds total saturation or when the amount of rainfall 
  # in excess of interception, evapotranspiration and drainage fills up soil water above total saturation.
  RUNOFF = max(0, (WA - WAST) / DELT + (RAIN - (RNINTC + EVAP + TRAN + DRAIN)))
  
  
  # The irrigation rate is the extra amount of water that is needed to keep soil water at a fraction of field capacity that 
  # is defined by setting the parameter IRRIGF. If IRRIGF is set to 1, the soil will be irrigated every timestep to keep the 
  # amount of water in the soil at field capacity. IRRIGF = 0 implies rainfed conditions.
  IRRIG  = IRRIGF * max(0, (WAFC - WA) / DELT - (RAIN - (RNINTC + EVAP + TRAN + DRAIN + RUNOFF)))
  
  
  DRUNIR <- data.frame(DRAIN = DRAIN,RUNOFF=RUNOFF, IRRIG=IRRIG)
  
  return(DRUNIR)
}

#------------------------------------------------------------------------------------------------------#
# help FUNCTION                                                                                        #
# Purpose: creating full output of states rates and relevant auxillary variables                       #
# Adds all rates and 
#------------------------------------------------------------------------------------------------------#
get_results <- function(state_out, year, wdata, STTIME,FINTIM){
  #get only relevant weather data
  WDATA <- wdata[STTIME:FINTIM,]
  parv <- LINTUL2_CASSAVA_parameters()
  
  state_out = data.frame(state_out)
  #Determine rate variables from change in states. 
  dt= state_out[2:nrow(state_out),1] - state_out[1:nrow(state_out)-1,1]
  dstate = (state_out[2: nrow(state_out),   2:ncol(state_out)] 
            - state_out[1:(nrow(state_out)-1),2:ncol(state_out)])
  #Add 0s to first day and determine rates of change
  rate_out = rbind(0,dstate/dt)
  colnames(rate_out)<-c("RROOTD","RWA","RTSUM", "RTSUMCROP", "RTSUMCROPLEAFAGE", "RDORMSUM", "RPUSHDORMRECTSUM", "RPUSHREDISTENDTSUM", "RDORMTIME", "RWCUTTING", "RTRAIN","RPAR","RLAI", 
                        "RWLVD", "RWLV", "RWST","RWSO","RWRT",
                        "RLVG","RTRAN","REVAP","RPTRAN","RPEVAP" ,"RREDISTLVG","RREDISTSO","RPUSHREDISTSUM","RWSOFASTRANSLSO")
  
  # Determine water content of rooted soil
  WC  <- 0.001 * state_out$WA/state_out$ROOTD
  
  RNINTC <- pmin(WDATA$RAIN,0.25 * state_out$LAI)         # interception of rain by the canopy (mm d-1)
  
  # Potential evaporation (mm d-1) and transpiration (mm d-1) are calculated according to Penman-Monteith
  PENM   <- penman(WDATA$DAVTMP,WDATA$VP,WDATA$DTR,state_out$LAI,WDATA$WN,RNINTC)
  
  # Actual evaporation (mm d-1) and transpiration (mm d-1)
  EVA  <- evaptr(PENM$PEVAP,PENM$PTRAN,state_out$ROOTD,state_out$WA,
                 parv[["WCAD"]],parv[["WCWP"]],parv[["WCFC"]],parv[["WCWET"]],parv[["WCST"]],parv[["TRANCO"]],1)
  
  # The transpiration reduction factor is defined as the ratio between actual and potential transpiration
  PENM$PTRAN[PENM$PTRAN== 0] <- 1 
  TRANRF <- rate_out[["RTRAN"]]/rate_out[["RPTRAN"]]
  
  #Determine harvest index
  HI <- state_out[["WSO"]] / (state_out[["WSO"]]+state_out[["WLVG"]]+state_out[["WRT"]]+state_out[["WST"]])
  #convert from g m-2 to t ha-1 with factor 10E-2: t/g=10E-6 * factor m2/ha=10E4 = 0.01
  WSOTHA <- state_out[["WSO"]] * 0.01
  
  
  result <- data.frame(cbind(year, state_out, rate_out,TRANRF,HI,WSOTHA))
  
  return(result) 
}



# provides dry matter in kg/ha
runLINTUL2 <- function(inputData_soil=inputData_soil,inputData_Weather=inputData_Weather, STTIME=STTIME, FINTIM=FINTIM, DELT=1){	
  long <- inputData_Weather$long
  lat <- inputData_Weather$lat	
  YR <- unique(inputData_Weather$Year)[1]
  YR2 <- YR+1
  param_noIrrigation <- LINTUL2_CASSAVA_parameters(irri=FALSE, soiltype="test")
  
  param_noIrrigation[["ROOTDM"]]  <- unique(inputData_soil$ROOTDM)
  param_noIrrigation[["WCFC"]]  <- inputData_soil$WCFC
  param_noIrrigation[["WCWP"]]  <- inputData_soil$WCWP
  param_noIrrigation[["WCST"]]  <- inputData_soil$WCST
  param_noIrrigation[["WCWET"]]  <- inputData_soil$WCWET 
  param_noIrrigation[["WCAD"]]  <- inputData_soil$WCAD
  param_noIrrigation[["DRATE"]]  <- 50
  param_noIrrigation[["STTIME"]] <- STTIME
  require(deSolve)
  #state_wlim1 <- ode(LINTUL2_CASSAVA_iniSTATES_MTC(param_noIrrigation), seq(STTIME, FINTIM, by = DELT), LINTUL2_CASSAVA, param_noIrrigation,  WDATA = inputData_Weather, method = "euler")
  state_wlim1 <- ode(LINTUL2_CASSAVA_iniSTATES_MTC(param_noIrrigation), seq(STTIME, FINTIM, by = DELT), LINTUL2_CASSAVA_MC, param_noIrrigation,  WDATA = inputData_Weather, method = "euler")
  
  var_wlim1 <- get_results(state_wlim1, YR, inputData_Weather, STTIME, FINTIM)
  locData <- unique(data.frame(lat = lat, long = long, plantingDate=STTIME))
  ## get weekly WLY for harvest time from 8 to 12 months
  WLY_estimates <- var_wlim1$WSO
  #WLY_time <- as.data.frame(matrix(nrow=1, ncol=22, data=WLY_estimates[seq(214, 365, 7)]))
  WLY_time <- as.data.frame(matrix(nrow=1, ncol=32, data=WLY_estimates[seq(235, 455, 7)]))
  WLY_time <- WLY_time * 10
  colnames(WLY_time) <- paste("WLY", seq(235, 455, 7), sep="_")
  tsumDD <- data.frame(lat=unique(locData$lat), long=unique(locData$long), Tsumcrop = var_wlim1$TSUMCROP, 
                       TSUM = var_wlim1$TSUM, SproutDate = c(1:455))
  tsumDD <- tsumDD[order(tsumDD$Tsumcrop), ]
  tsumDD <- tsumDD[tsumDD$Tsumcrop > 0, ][1,]
  QWW <- cbind(locData, WLY_time)
  QWW$TSUM <- tsumDD$TSUM
  QWW$TSUMCROP <- tsumDD$Tsumcrop
  QWW$SproutDate <- tsumDD$SproutDate
  return(QWW)
  
}




#' Title
#'
#' @param ll one lat and lon reading, concatnated 
#' @param SoilData 
#' @param Planting_Harvest_days 
#' @param RainData 
#' @param Solar 
#' @param country 
#' @param year 
#'
#' @return a data frame with GPS readings, plantingDate, and weekly WLY estimates between 8 - 15 
#' months after harest. The estimated WLY shuld be read as dry matter in kg/ha. It writes the same 
#' output in /home/akilimo/projects/AKILIMO_Modeling/LINTUL/country/LINTUL_WLY_year
#'
#' @examples  Lintul_oneCoord_WLY(ll="10.075_-1.375", SoilData = SoilData, Planting_Harvest_days= Planting_Harvest_days,
#' RainData = RainData, Solar = Solar, country = country, year=year)
#' @author Meklit
Lintul_oneCoord_WLY <- function(ll, SoilData, Planting_Harvest_days, RainData, Solar, country, year){
  
  lat <- strsplit(ll, "_")[[1]][1]
  long <- strsplit(ll, "_")[[1]][2]
  
  dd <- RainData[RainData$location==ll, ]
  OnePixel_RF <- dd[order(dd$Month, dd$Date), ]
  
  ## if all the arainfall data for a pixel is -9999, leave that pixel out, otherwise 
  ## if there is missing data for rainfall, take the mean of the days before and after
  if(nrow(OnePixel_RF[OnePixel_RF$RAIN >= 0,]) > 300){
    if(any(OnePixel_RF$RAIN <0)){
      OnePixel_RF <- OnePixel_RF[order(OnePixel_RF$Date),	]		
      impute <- NULL
      for(n in 1:nrow(OnePixel_RF)){
        if(OnePixel_RF[n, "RAIN"] < 0){
          psrad <- OnePixel_RF[n,]
        }else{
          psradP <- OnePixel_RF[n-1, ]
          psradL <- OnePixel_RF[n+1, ]				
          psrad <- OnePixel_RF[n, ]
          psrad$RAIN <- mean(OnePixel_RF$RAIN, OnePixel_RF$RAIN)
        }			
        impute <- rbind(impute, psrad)				
      }
      OnePixel_RF <- impute
    }				
    
    
    ### Solar Data
    Solar_coor <- Solar[Solar$location == ll, ]
    Solar_coor$DateYear <- yday(gsub("_", "-", Solar_coor$YYYYMMDD))
    
    for(k in 1:nrow(Solar_coor)){
      Solar_coor$YEAR[k] <- strsplit(Solar_coor$YYYYMMDD[k], "_")[[1]][1]
      Solar_coor$MONTH[k] <- strsplit(Solar_coor$YYYYMMDD[k], "_")[[1]][2]
      Solar_coor$DATE[k] <- strsplit(Solar_coor$YYYYMMDD[k], "_")[[1]][3]
    }
    
    Solar_coor$DATE <- as.numeric(as.character(Solar_coor$DATE))
    Solar_coor$MONTH <- as.numeric(as.character(Solar_coor$MONTH))
    
    Solar_coor <- Solar_coor[, c("lat", "lon", "DATE", "DTR", "TMIN", "TMAX", "VP", "WIND", "DateYear")]
    colnames(Solar_coor) <- c("lat", "long", "Day", "DTR", "TMIN", "TMAX", "VP", "WIND", "Date")
    
    
    ## merge rain and solar data
    Rain_Solar <- merge(OnePixel_RF, Solar_coor, by=c("Date", "Day","long", "lat"))
    #RFSolar$DTR <- RFSolar$DTR / 1000 ## conversion to comply the LINTUL units in the calculation of soalr radiation
    Rain_Solar$DAVTMP <- (Rain_Solar$TMIN + Rain_Solar$TMAX)/2	
    
    RFSolar <- Rain_Solar[, c("Year","VP", "WIND","RAIN","DTR","DAVTMP","Month","Date","lat","long")]
    
    rss1 <- RFSolar
    rss1$Year <- rss1$Year +1
    stdate2year <- max(RFSolar$Date) + 1
    endate2year <- max(RFSolar$Date) + nrow(RFSolar)
    rss1$Date <- c(stdate2year:endate2year)
    
    rss2 <- rss1
    rss2$Year <- rss2$Year +1
    stdate3year <- max(rss1$Date) + 1
    endate3year <- max(rss1$Date) + nrow(rss1)
    rss2$Date <- c(stdate3year:endate3year)
    
    rss3 <- rss2
    rss3$Year <- rss3$Year +1
    stdate4year <- max(rss2$Date) + 1
    endate4year <- max(rss2$Date) + nrow(rss2)
    rss3$Date <- c(stdate4year:endate4year)
    
    RRSS <- rbind(RFSolar, rss1, rss2, rss3)
    colnames(RRSS) <- c("Year", "VP", "WN", "RAIN", "DTR", "DAVTMP", "Month", "Date", "lat", "long")
    
    sdd <- SoilData[SoilData$location==ll,  ]	
    sdd$ROOTDM <- 1.8
    
    
    WLY <- NULL
    for(g in 1:nrow(Planting_Harvest_days)){		
      pldate <- Planting_Harvest_days$startDate[g]
      harvestingdate <- Planting_Harvest_days$endDate[g]	
      STTIME <- pldate
      FINTIM <- harvestingdate
      wly_LL <- NULL
      wly_LL <- runLINTUL2(inputData_Weather = RRSS, inputData_soil = sdd, STTIME = STTIME, FINTIM = FINTIM)      
      
      WLY <- rbind(WLY, wly_LL)			
    }	
   
     return(WLY)
  }
}





LINTUL2_CASSAVA_MC <-function(Time, State, Pars, WDATA){
  Pars <- as.data.frame(as.list(Pars))
  State <- as.data.frame(as.list(State))
  parsSatate <- cbind(Pars, State)
  
  #intergration with delay for RROOTD
  #with(as.list(c(State, Pars)), {
  
  RTRAIN <- WDATA$RAIN[Time]                       # rain rate, mm d-1
  DTEFF  <- max(0, WDATA$DAVTMP[Time] - parsSatate$TBASE)     # effective daily temperature (for crop development a treshold temperature (TBASE) needs to be exceeded)
  RPAR   <- parsSatate$FPAR * WDATA$DTR[Time]                 # PAR MJ m-2 d-1
  
  #determine rates when crop is still growing
  if(parsSatate$TSUM < parsSatate$FINTSUM){
    # Determine water content of rooted soil
    WC  <- 0.001 * parsSatate$WA/parsSatate$ROOTD
    
    #	TAP <- ifelse((Time-parsSatate$DOYPL) >= 0, 1, 0)
    RTSUM <- DTEFF * ifelse((Time-parsSatate$DOYPL) >= 0, 1, 0)
    #RTSUM <- DTEFF*ifelse((Time-parsSatate$DOYPL) >= 0, 1, 0) # TSUM after planting
    #RTSUM <- DTEFF*ifelse((Time-STTIME) >= 0, 1, 0) # TSUM after planting
    
    # Once the emergence date is reached and enough water is available the crop emerges (1), once the crop is established is does not disappear again (2)
    if((WC-parsSatate$WCWP) >= 0 && (parsSatate$TSUM-parsSatate$OPTEMERGTSUM) >= 0) { 	# (1)
      emerg1 <- 1} else { emerg1 <- 0 }
    if(parsSatate$TSUMCROP > 0) {									# (2)
      emerg2 <- 1 } else { emerg2 <- 0}
    # Emergence of the crop is used to calculate the accumulated temperature sum.
    EMERG  <- max(emerg1,emerg2)
    RTSUMCROP <- DTEFF*EMERG
    
    # Root depth
    # If soil water content drops to, or below, wilting point fibrous root growth stops.
    # As long as the crop has not reached its maximum rooting depth and has not started flowering yet, fibrous root growth continues.
    if((parsSatate$ROOTD-parsSatate$ROOTDM) < 0 && (WC-parsSatate$WCWP) >= 0) {
      # The rooting depth (m) is calculated from a maximum rate of change in rooting depth, the emergence of the crop and the constraints mentioned above.
      RROOTD <- parsSatate$RRDMAX * EMERG 
    }else{ 
      RROOTD = 0
    }
    
    EXPLOR <- 1000 * RROOTD * parsSatate$WCFC                   # exploration rate of new soil water layers by root depth growth (mm d-1)
    
    
    RNINTC <- min(RTRAIN, (parsSatate$FRACRNINTC * parsSatate$LAI))                # interception of rain by the canopy (mm d-1)
    
    # Potential evaporation (mm d-1) and transpiration (mm d-1) are calculated according to Penman-Monteith
    PENM   <- penman(WDATA$DAVTMP[Time],WDATA$VP[Time],WDATA$DTR[Time],parsSatate$LAI,WDATA$WN[Time],RNINTC)
    RPTRAN <- PENM$PTRAN
    RPEVAP <- PENM$PEVAP
    
    WCSD <- parsSatate$WCWP * parsSatate$TWCSD
    WCCR <- parsSatate$WCWP + pmax(WCSD-parsSatate$WCWP, (parsSatate$PTRAN/(parsSatate$PTRAN+parsSatate$TRANCO) * (parsSatate$WCFC-parsSatate$WCWP))/1000)
    
    # Actual evaporation (mm d-1) and transpiration (mm d-1)
    EVA  <- evaptr(RPEVAP,RPTRAN,parsSatate$ROOTD,parsSatate$WA,parsSatate$WCAD,parsSatate$WCWP,parsSatate$TWCSD,parsSatate$WCFC,parsSatate$WCWET,parsSatate$WCST,parsSatate$TRANCO, parsSatate$DELT)
    RTRAN <- EVA$TRAN
    REVAP <- EVA$EVAP
    
    # The transpiration reduction factor is defined as the ratio between actual and potential transpiration
    TRANRF = ifelse(RPTRAN <= 0, 1, RTRAN/RPTRAN)
    
    # Drainage (below the root zone; mm d-1), surface water runoff (mm d-1) and irrigation rate (mm d-1)
    DRUNIR    <- drunir(RTRAIN,RNINTC,REVAP,RTRAN,parsSatate$IRRIGF,parsSatate$DRATE,parsSatate$DELT,parsSatate$WA,parsSatate$ROOTD,parsSatate$WCFC,parsSatate$WCST)
    
    # Rate of change of soil water amount (mm d-1)
    RWA <- (RTRAIN + EXPLOR + DRUNIR$IRRIG) - (RNINTC + DRUNIR$RUNOFF + RTRAN + REVAP + DRUNIR$DRAIN)
    WC <- 0.001 * parsSatate$WA/parsSatate$ROOTD
    
    # Light interception (MJ m-2 d-1) and total crop growth rate (g m-2 d-1)
    PARINT <- RPAR * (1 - exp(-parsSatate$K_EXT * parsSatate$LAI))
    LUE    <- parsSatate$LUE_OPT * approx(TTB[,1], TTB[,2], WDATA$DAVTMP[Time])$y
    
    # Dormancy and recovery from dormancy
    if ((WC-WCSD) <= 0 && (parsSatate$LAI - parsSatate$LAI_MIN) <= 0){
      dormancy = 1
    } else {
      dormancy = 0
    }
    
    if ((WC - parsSatate$RECOV * WCCR) >= 0 && (WC - parsSatate$WCWP) >= 0){
      pushdor = 1
    } else {
      pushdor = 0
    }
    
    if (parsSatate$WSO == 0) {
      WSOREDISTFRAC = 1
    } else {
      WSOREDISTFRAC = parsSatate$REDISTSO/parsSatate$WSO
    }
    
    PUSHREDISTEND = max(ifelse((WSOREDISTFRAC-parsSatate$WSOREDISTFRACMAX) >= 0, 1, 0), ifelse((parsSatate$REDISTLVG - parsSatate$WLVGNEWN)>= 0, 1, 0), ifelse((parsSatate$PUSHREDISTSUM - parsSatate$TSUMREDISTMAX) >= 0, 1, 0)) *ifelse(-parsSatate$PUSHREDISTSUM >= 0, 0, 1)
    PUSHREDIST = ifelse((parsSatate$PUSHDORMRECTSUM - parsSatate$DELREDIST) >= 0, 1, 0)* (1 - PUSHREDISTEND)
    PUSHDORMREC = pushdor*ifelse(-parsSatate$DORMTSUM >= 0, 0, 1) * (1 - PUSHREDIST) * ifelse((parsSatate$TSUMCROP - parsSatate$TSUMSBR) >= 0, 1, 0)
    DORMANCY = max(dormancy, PUSHDORMREC) * (1 - PUSHREDIST) * ifelse((parsSatate$TSUMCROP - parsSatate$TSUMSBR) >= 0, 1, 0)
    RDORMTSUM = DTEFF *DORMANCY - (parsSatate$DORMTSUM/parsSatate$DELT) * PUSHREDIST
    RPUSHDORMRECTSUM = DTEFF * PUSHDORMREC - (parsSatate$PUSHDORMRECTSUM/parsSatate$DELT) * (1 - PUSHDORMREC) * (1 - PUSHREDIST)
    RPUSHREDISTSUM = DTEFF * PUSHREDIST - (parsSatate$PUSHREDISTSUM/parsSatate$DELT) * PUSHREDISTEND
    RPUSHREDISTENDTSUM = DTEFF * PUSHREDIST - (parsSatate$PUSHREDISTENDTSUM/parsSatate$DELT) * (1 - PUSHREDISTEND)
    RDORMTIME = DORMANCY
    
    # Dry matter redistribution after dormancy
    RREDISTSO = parsSatate$RRREDISTSO * parsSatate$WSO * PUSHREDIST - (parsSatate$REDISTSO/parsSatate$DELT) *ifelse(-parsSatate$DORMTSUM >= 0, 0, 1)
    RREDISTLVG = parsSatate$SO2LV * RREDISTSO * (1- DORMANCY)
    RREDISTMAINTLOSS = (1 - parsSatate$SO2LV) * RREDISTSO
    
    GTOTAL <- LUE * PARINT * TRANRF * (1 - DORMANCY)
    
    # Relative death rate (d-1) due to aging
    RDRDV = ifelse(parsSatate$TSUMCROPLEAFAGE - parsSatate$TSUMLLIFE >= 0, approx(RDRT[,1], RDRT[,2], WDATA$DAVTMP[Time])$y, 0)
    
    # Relative death rate (d-1) due to self shading
    RDRSH <- parsSatate$RDRSHM * (parsSatate$LAI - parsSatate$LAICR) / parsSatate$LAICR
    if(RDRSH < 0) {
      RDRSH <- 0
    } else if(RDRSH >= parsSatate$RDRSHM) {
      RDRSH <- parsSatate$RDRSHM
    }
    
    RTSUMCROPLEAFAGE <- DTEFF * EMERG - (parsSatate$TSUMCROPLEAFAGE/parsSatate$DELT) * PUSHREDIST
    ENHSHED <- max(ifelse((WC-WCSD) >= 0, 0, 1), ifelse((WC-parsSatate$WCWET) >= 0, 1, 0))*ifelse((parsSatate$TSUMCROPLEAFAGE-parsSatate$FRACTLLFENHSH*parsSatate$TSUMLLIFE) >= 0, 1, 0)
    
    # Relative death rate (d-1) due to severe drought
    RDRSD <- parsSatate$RDRB *ENHSHED
    
    # Effective relative death rate (1; d-1) and the resulting decrease in LAI (2; m2 m-2 d-1) and leaf weight (3; g m-2 d-1)
    RDR   <- max(RDRDV, RDRSH, RDRSD) * ifelse((parsSatate$TSUMCROPLEAFAGE - parsSatate$TSUMLLIFE) >= 0, 1, 0) 	# (1)
    DLAI  <- parsSatate$LAI * RDR * (1 - parsSatate$FASTRANSLSO) * (1 - DORMANCY)  			  # (2)
    
    # Allocation to roots (2), leaves (4), stems (5) and storage organs (6)
    # fractions allocated are modified for water availability (1 and 3)
    FRTMOD <- max(1, 1/(TRANRF+0.5))						    # (1)
    FRT    <- approx(FRTTB[,1],FRTTB[,2],parsSatate$TSUMCROP)$y * FRTMOD		# (2)
    FSHMOD <- (1 - FRT) / (1 - FRT / FRTMOD)				    # (3)
    FLV    <- approx(FLVTB[,1],FLVTB[,2],parsSatate$TSUMCROP)$y * FSHMOD		# (4)
    FST    <- approx(FSTTB[,1],FSTTB[,2],parsSatate$TSUMCROP)$y * FSHMOD		# (5)
    # Error in approx(FSTTB[, 1], FSTTB[, 2], TSUMCROP) : 
    #   object 'TSUMCROP' not found
    
    FSO    <- approx(FSOTB[,1],FSOTB[,2],parsSatate$TSUMCROP)$y * FSHMOD		# (6)
    
    # stem cutting growth
    WCUTTINGMIN <- parsSatate$WCUTTINGMINPRO * parsSatate$WCUTTINGIP
    
    # Leaf growth and senescence
    FRACSLACROPAGE <- approx(FRACSLATB[,1], FRACSLATB[,2], parsSatate$TSUMCROP)$y
    SLA <- parsSatate$SLA_MAX *FRACSLACROPAGE
    RWSOFASTRANSLSO = parsSatate$WLVG * RDR * parsSatate$FASTRANSLSO * (1 - DORMANCY)
    DLV <- (parsSatate$WLVG * RDR - RWSOFASTRANSLSO) * (1 - DORMANCY)
    RWLVD <- DLV
    
    # Change in biomass (g m-2 d-1) for green leaves (1), stems (2), storage organs (3) and roots (4)
    if (parsSatate$TSUM > parsSatate$OPTEMERGTSUM && parsSatate$WST == 0){
      RWCUTTING <- parsSatate$WCUTTING *(-parsSatate$FST_CUTT - parsSatate$FRT_CUTT - parsSatate$FLV_CUTT - parsSatate$FSO_CUTT)
      RWST <- parsSatate$WCUTTINGIP * parsSatate$FST_CUTT
      RWRT <- parsSatate$WCUTTINGIP * parsSatate$FRT_CUTT
      RWLVG <- parsSatate$WCUTTINGIP * parsSatate$FLV_CUTT
      RWSO <- parsSatate$WCUTTINGIP * parsSatate$FSO_CUTT
    } else if (parsSatate$TSUMCROP > 16.8) {   # A little bit of cheating
      RWCUTTING <- -parsSatate$RDRWCUTTING * parsSatate$WCUTTING * ifelse((parsSatate$WCUTTING-WCUTTINGMIN) >= 0, 1, 0) * TRANRF * EMERG * (1 - DORMANCY)
      RWLVG  <- (abs(GTOTAL)+abs(RWCUTTING)) * FLV - DLV + RREDISTLVG * PUSHREDIST 	# (1)
      RWST   <- (abs(GTOTAL)+abs(RWCUTTING)) * FST			  # (2)
      RWSO   <- (abs(GTOTAL)+abs(RWCUTTING)) * FSO + RWSOFASTRANSLSO - RREDISTSO			  # (3)
      RWRT   <- (abs(GTOTAL)+abs(RWCUTTING)) * FRT			  # (4)
    } else{
      RWCUTTING <- 0
      RWLVG <- 0
      RWST <- 0
      RWSO <- 0
      RWRT <- 0
    }
    RWLV = RWLVG+RWLVD
    WGTOTAL = parsSatate$WLV+parsSatate$WST+parsSatate$WCUTTING+parsSatate$WSO+parsSatate$WRT
    GLV <- FLV * (GTOTAL + abs(RWCUTTING)) + RREDISTLVG * PUSHREDIST
    
    GLAI <- gla(DTEFF, parsSatate$TSUMCROP, parsSatate$LAII, parsSatate$RGRL, parsSatate$DELT, SLA, parsSatate$LAI, GLV, parsSatate$TSUMLA_MIN,
                TRANRF, WC, parsSatate$WCWP, RWCUTTING, FLV, parsSatate$LAIEXPOEND,DORMANCY)
    # Error in gla(DTEFF, parsSatate$TSUMCROP, LAII, RGRL, DELT, SLA, parsSatate$LAI,  : 
    #   object 'RGRL' not found
    
    
    # Change in LAI (m2 m-2 d-1) due to new growth of leaves
    RLAI <- GLAI - DLAI
    
    
  }else{
    #all plant related rates are set to 0
    RROOTD <- 0
    RWLVG <- 0
    RWA <- 0
    RTSUMCROP <- 0
    RLAI <- 0
    RWLVG <- 0
    RWLVD <- 0
    RWLV <- 0
    RWST <- 0
    RWSO <- 0
    RWRT <- 0
    RTRAN <- 0 
    REVAP <- 0 
    RPTRAN <- 0 
    RPEVAP <- 0
    RGTOTAL <- 0
    RDORMTSUM <- 0
    RTSUM <- 0
    RTSUMCROPLEAFAGE <- 0
    RPUSHREDISTENDTSUM <- 0
    RPUSHDORMRECTSUM <- 0
    RDORMTIME <- 0
    RWCUTTING <- 0
    RREDISTLVG <- 0
    RREDISTSO <- 0
    RPUSHREDISTSUM <- 0
    RWSOFASTRANSLSO <- 0
  }
  return(list(c(RROOTD,RWA,RTSUM, RTSUMCROP, RTSUMCROPLEAFAGE,RDORMTSUM,RPUSHDORMRECTSUM,RPUSHREDISTENDTSUM, RDORMTIME, RWCUTTING, RTRAIN,RPAR,RLAI,RWLVD, RWLV,RWST,RWSO,RWRT, RWLVG,RTRAN,REVAP,RPTRAN,RPEVAP,RREDISTLVG,RREDISTSO,RPUSHREDISTSUM,RWSOFASTRANSLSO)))
  
  #})
}




#---------------------------------------------------------------------#
# FUNCTION adapted parameters                                               #
# Purpose: Listing the input parameters for Lintul2                   #
#---------------------------------------------------------------------#
LINTUL2_CASSAVA_parameters <-function(irri=FALSE,soiltype="clay") {
  #get the true defaults
  PARAM <- LINTUL2_CASSAVA_DEFAUL_PARAMETERS() 
  
  #change what is desired
  PARAM[["DOYEM"]]  <- 32    # day nr of emergence
  
  if(irri==TRUE){ PARAM[["IRRIGF"]] <- 1 }
  
  if(soiltype=="clay"){
    PARAM[["ROOTDM"]]  <- 1.2    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.08   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.20   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.46   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.49   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.52   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }else  if( soiltype=="sand"){
    PARAM[["ROOTDM"]]  <- 0.6    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.05   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.08   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.16   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.40   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.42   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }else  if( soiltype=="sandPRAC"){
    PARAM[["ROOTDM"]]  <- 0.3    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.04   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.09   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.28   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.32   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.38   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }else  if( soiltype=="Nitisol"){
    PARAM[["ROOTDM"]]  <- 0.9    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.01   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.20   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.36   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.37   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.48   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }else  if( soiltype=="Gleysol"){
    PARAM[["ROOTDM"]]  <- 0.9    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.01   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.15   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.28   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.30   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.39   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }else  if( soiltype=="Ferralsol"){
    PARAM[["ROOTDM"]]  <- 0.9    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.01   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.18   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.34   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.39   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.53   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }else  if( soiltype=="Arenosol"){
    PARAM[["ROOTDM"]]  <- 0.9    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.01   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.04   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.07   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.08   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.11   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  } else if (soiltype=="test"){
    PARAM[["ROOTDM"]]  <- 0.9    # m           :     maximum rooting depth
    PARAM[["WCAD"]]    <- 0.01   # m3 m-3      :     soil water content at air dryness 
    PARAM[["WCWP"]]    <- 0.12   # m3 m-3      :     soil water content at wilting point
    PARAM[["WCFC"]]    <- 0.41   # m3 m-3      :     soil water content at field capacity 
    PARAM[["WCWET"]]   <- 0.46   # m3 m-3      :     critical soil water content for transpiration reduction due to waterlogging
    PARAM[["WCST"]]    <- 0.52   # m3 m-3      :     soil water content at full saturation 
    PARAM[["DRATE"]]   <- 50     # mm d-1      :     max drainage rate
  }
  
  return(PARAM)  
}





#' Title
#'
#' @param year the year for whcih LINTUL is running
#' @param country 
#' @param SoilData the soil data 
#' @param RainData 
#' @param Solar 
#' @param Planting_Harvest_days a data framw with start date, end date (measuring the growing period) the week 
#' number and the month of the start date
#' @param nrCores if the script is running on parllel it refers to the number of cores to be used otherwise it is NA
#'
#' @return it does not return any data but writes out WLY estimate for every pixel in 
#' /home/akilimo/projects/AKILIMO_Modeling/LINTUL/country/LINTUL_WLY_year
#' @author Meklit
runLINTUL_server <- function(year=year, country = country,SoilData, RainData, Solar,
                             Planting_Harvest_days=Planting_Harvest_days){
  
  require(deSolve)   
  require(plyr)
  require(lubridate) 
  require(parallel)
  require(tidyr)
  library(doParallel)
  library(foreach)
  
  
  SoilData <- droplevels(SoilData[complete.cases(SoilData), ])
  SoilData$WCFC <- (SoilData$FC_wAverage)/100
  SoilData$WCWP <- (SoilData$wp_wAverage)/100
  SoilData$WCST <- (SoilData$sws_wAverage)/100
  SoilData$WCWET <- ((SoilData$sws_wAverage) - 6)/100
  SoilData$WCAD <- ((SoilData$wp_wAverage)/2)/100
  SoilData$ROOTDM <- 2
  
  
  RainData <- RainData[RainData$location %in% SoilData$location, ]
  Solar <- Solar[Solar$location %in% SoilData$location, ]
  
  ## removing locations with no rain data
  RainData <- RainData[RainData$RAIN >= 0, ]
  SoilData <- SoilData[SoilData$location %in% RainData$location, ]
  Solar <- Solar[Solar$location %in% RainData$location, ]
  latlong <-  unique(SoilData$location)
  
  SoilData <- SoilData[,c("wp_wAverage", "sws_wAverage", "FC_wAverage","WCFC","WCWP", "WCST", 
                          "WCWET", "WCAD","location")]
  
  names(SoilData) <- c("WiltingPoint","WaterAtSaturation", "FieldCapacity", "WCFC","WCWP","WCST","WCWET",
                       "WCAD","location")
  
  
  print("data sourcing is finished.")
  
 
  WLYres <- NULL
    for(i in 1:length(latlong)){
      ll <- latlong[i]
      WLYres_ll <- Lintul_oneCoord_WLY(ll, SoilData = SoilData, Planting_Harvest_days= Planting_Harvest_days,
                          RainData = RainData, Solar = Solar, country = country, year=year)
      WLYres <- rbind(WLYres, WLYres_ll)
    }
  
   return(WLYres)
  
}

