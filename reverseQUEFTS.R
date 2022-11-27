

######################################################################################
######################################################################################
## Reverse QUEFTS
######################################################################################
## read the example data
RQ_Data <- readRDS("RQ_Data.RDS")

## run reverse QUEFTS 
source("QUEFTS_Functions.R")
QUEFTS_soilNPK <- Reverse_QUEFTS(ds_GPS = RQ_Data, rateN = 150, rateP = 40, rateK = 180, 
                                         crop = "Cassava", halfrate=FALSE, hrN = 75, hrP = 20, hrK = 90)

## check the result
ACAI_CON_BLUP_reversQ <- ggplot(QUEFTS_soilNPK, aes(Control_observed, Control_estimated, col=country)) +
  geom_point() + 
  geom_abline(intercept=0, slope=1) +
  facet_wrap(~country, scale="free") +
  # xlim(0,maxrange) + ylim(0,maxrange)+
  xlab("Control yeild Observed (BLUPS)") +
  ylab(" Control yield estimated through QUEFTS") +
  theme_bw() +
  theme(axis.text  = element_text(size=10), 
        axis.title = element_text(size=15), strip.text = element_text(size=12), 
        legend.text =element_text(size=15), legend.position = "none")


