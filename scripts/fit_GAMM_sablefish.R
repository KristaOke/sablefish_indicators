#=========================================================================================
# Fitting GAMM to sablefish indicator data from ESP
#
#=========================================================================================
#Notes:
#=========================================================================================
#packages
library(tidyverse)
library(mgcv)
library(gamm4)
library(corrplot)
library(cowplot)

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#=============================================================
#### get data ####

#first load whole dataset and look to see if any indicators look like they have nonlinear relationships
#with recruitment

scaled_dat <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)

#ONLY NONCORRELATED-----

noncor_covars3 <- c("Year", "ln_rec",
                    "Spr_ST_SEBS_scaled"  ,
                    "Smr_temp_250m_GOA_scaled"  , 
                    "Spr_chlA_biom_GOA_scaled"  , 
                    "Spr_chlA_biom_SEBS_scaled"  ,  
                    "Spr_chlA_peak_GOA_scaled" ,
                    "Spr_chlA_peak_SEBS_scaled" ,
                    "ann_Copepod_size_EGOA_scaled" ,
                    "YOY_grwth_Middleton_scaled",
                    "Smr_CPUE_juv_ADFG_ln_scaled"  ,      
                    "smr_adult_cond_scaled")


subset_dat <- scaled_dat[,names(scaled_dat) %in% noncor_covars3]

plot_covars <- subset_dat %>% gather(key=type, value=value, -c(Year, ln_rec))

rec.plot <- ggplot(plot_covars, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + geom_smooth() +
  #scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
rec.plot


#plot what's left after removing NAs except 250m temp, too short

nona_covars <- na.omit(subset_dat[,-c(4)])
nona_covars <- nona_covars %>% gather(key=type, value=value, -c(Year, ln_rec))

rec.plotn <- ggplot(nona_covars, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + #geom_smooth(method="lm") +
  #scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type) +
  theme(legend.position = "NA") + ylab("ln(Recruitment)")
rec.plotn

#update with new covars if using
# labels_type <- c("Annual copepod size EGOA",
#                  "Annual copepod size WGOA",
#                  "Summer adult condition",
#                  "Summer condition age4 females",
#                  "Summer CPUE juveniles ADFG 
# large mesh survey",
# "Spring chlorophyll A biomass SEBS",
# "Spring chlorophyll A peak GOA",
# "Spring chlorophyll A peak SEBS",
# "Spring SST SEBS")

rec.plotn <- ggplot(nona_covars, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + #geom_smooth(method="lm") +
  #scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, labeller = labeller(type=labels_type)) +
  theme(legend.position = "NA") + ylab("ln(Recruitment)")
rec.plotn

#Get data=======

# load data that is already z-scored, checked for correlations

train1 <- read.csv(file=paste(wd,"/data/dataset_training1.csv", sep=""), row.names = 1)
train2 <- read.csv(file=paste(wd,"/data/dataset_training2.csv", sep=""), row.names = 1)
train3 <- read.csv(file=paste(wd,"/data/dataset_training3.csv", sep=""), row.names = 1)
train4 <- read.csv(file=paste(wd,"/data/dataset_training4.csv", sep=""), row.names = 1)
train5 <- read.csv(file=paste(wd,"/data/dataset_training5.csv", sep=""), row.names = 1)

testing1 <- read.csv(file=paste(wd,"/data/dataset_testing1.csv", sep=""), row.names = 1)
testing2 <- read.csv(file=paste(wd,"/data/dataset_testing2.csv", sep=""), row.names = 1)
testing3 <- read.csv(file=paste(wd,"/data/dataset_testing3.csv", sep=""), row.names = 1)
testing4 <- read.csv(file=paste(wd,"/data/dataset_testing4.csv", sep=""), row.names = 1)
testing5 <- read.csv(file=paste(wd,"/data/dataset_testing5.csv", sep=""), row.names = 1)

#DATA CONTROL SECTION----

#select covariates

#using ONLY the indicator subset that's been checked for collinearity
#using noncor_covars from process_data for now

#GAM cannot handle missing data, euphasiid and the every other yr GOA indicator were already removed


train1_gam_dat <- train1[,names(train1) %in% noncor_covars3]
train2_gam_dat <- train2[,names(train2) %in% noncor_covars3]
train3_gam_dat <- train3[,names(train3) %in% noncor_covars3]
train4_gam_dat <- train4[,names(train4) %in% noncor_covars3]
train5_gam_dat <- train5[,names(train5) %in% noncor_covars3]

scaled_gam_dat <- scaled_dat[,names(scaled_dat) %in% noncor_covars3]



#new gams-----------

#updated to run on whole dataset Jun 17, 2023

#first try
mod1_1 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
               # s(Smr_temp_250m_GOA_scaled, k=4)  +  #too short    
                s(Spr_chlA_biom_GOA_scaled, k=4)  + 
  s(Spr_chlA_biom_SEBS_scaled, k=4)  +     
  s(Spr_chlA_peak_GOA_scaled, k=4)  +     
  s(Spr_chlA_peak_SEBS_scaled, k=4)  +  
s(ann_Copepod_size_EGOA_scaled, k=4)  +     
s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
s(YOY_grwth_Middleton_scaled, k=4) +
s(smr_adult_cond_scaled, k=4), data=scaled_gam_dat) #more coefficients than data
gam.check(mod1_1) #not good
summary(mod1_1)

#next shortest are chl covars try dropping
mod1_2 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
                # s(Smr_temp_250m_GOA_scaled, k=4)  +  #too short    
               # s(Spr_chlA_biom_GOA_scaled, k=4)  + 
               # s(Spr_chlA_biom_SEBS_scaled, k=4)  +     
               # s(Spr_chlA_peak_GOA_scaled, k=4)  +     
              #  s(Spr_chlA_peak_SEBS_scaled, k=4)  +  
                s(ann_Copepod_size_EGOA_scaled, k=4)  +     
                s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
                s(YOY_grwth_Middleton_scaled, k=4) +
                s(smr_adult_cond_scaled, k=4), data=scaled_gam_dat) #runs once all chl covars removed
gam.check(mod1_2) #not good
summary(mod1_2) #nothing significantly nonlinear, only YOY_grwth has edf over 1

#change all edf < 2 to linear and try again
mod1_3 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                # s(Smr_temp_250m_GOA_scaled, k=4)  +  #too short    
                # s(Spr_chlA_biom_GOA_scaled, k=4)  + 
                # s(Spr_chlA_biom_SEBS_scaled, k=4)  +     
                # s(Spr_chlA_peak_GOA_scaled, k=4)  +     
                #  s(Spr_chlA_peak_SEBS_scaled, k=4)  +  
                ann_Copepod_size_EGOA_scaled  +     
                Smr_CPUE_juv_ADFG_ln_scaled +    
                s(YOY_grwth_Middleton_scaled, k=4) +
                smr_adult_cond_scaled, data=scaled_gam_dat) #runs once all chl covars removed
gam.check(mod1_3) #not good
summary(mod1_3) #none are significantly nonlinear


mod1_4 <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +      
                ann_Copepod_size_EGOA_scaled  +     
                Smr_CPUE_juv_ADFG_ln_scaled +    
                YOY_grwth_Middleton_scaled +
                smr_adult_cond_scaled, data=scaled_gam_dat) #
summary(mod1_4)
anova(mod1_4)

#only ADFG survey is significant!


#run on ONLY 4 longest time series

mod1_5 <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +        
               Smr_CPUE_juv_ADFG_ln_scaled +    
               YOY_grwth_Middleton_scaled +
               smr_adult_cond_scaled, data=scaled_gam_dat) #
summary(mod1_5)
anova(mod1_5) #still only ADFG





