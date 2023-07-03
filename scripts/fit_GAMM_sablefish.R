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


#DATA CONTROL SECTION----

#select covariates

#using ONLY the indicator subset that's been checked for collinearity
#using noncor_covars from process_data for now

#GAM cannot handle missing data, euphasiid and the every other yr GOA indicator were already removed


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
anova(mod1_5) #still only ADFG significant

plot(mod1_5)
ggplot(scaled_gam_dat, aes(Smr_CPUE_juv_ADFG_ln_scaled, ln_rec)) + geom_point() +
  geom_smooth(method="lm")


#LOOCV - all time series=======

#STEP 1 - Loop through training sets and fit models-------

#ONLY 4 longest time series
scaled_loop_dat <- scaled_gam_dat

yrs <- unique(scaled_loop_dat$Year)
output_df <- data.frame(matrix(ncol=3, nrow = length(yrs)))
colnames(output_df) <- c("Year", "observed_ln_recruit", "predicted_ln_recruit")

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  dropped_yr <- scaled_loop_dat[i,]
  output_df$observed_ln_recruit[i] <- dropped_yr$ln_rec
  dropped_yr <- dropped_yr[,!names(dropped_yr) %in% "ln_rec"]
  #fit model
  modtemp <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +      
                  ann_Copepod_size_EGOA_scaled  +     
                  Smr_CPUE_juv_ADFG_ln_scaled +    
                  YOY_grwth_Middleton_scaled +
                  smr_adult_cond_scaled, data=temp_dat)
  #have model predict to missing year
   temp_predict <- predict(modtemp, newdata=dropped_yr)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict[1]
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  geom_point() + geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year))

#STEP 2 - get MSE, MAE, and R2------

#get MSE & MAE------

#these need to be double checked!
GAM_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)

GAM_MAE <- ((sum((abs(output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2), na.rm = TRUE)))/length(output_df$observed_ln_recruit)

obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)




#LOOCV - long time series=======

#STEP 1 - Loop through training sets and fit models-------

#
scaled_loop_dat <- scaled_gam_dat

yrs <- unique(scaled_loop_dat$Year)
output_df <- data.frame(matrix(ncol=3, nrow = length(yrs)))
colnames(output_df) <- c("Year", "observed_ln_recruit", "predicted_ln_recruit")

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  dropped_yr <- scaled_loop_dat[i,]
  output_df$observed_ln_recruit[i] <- dropped_yr$ln_rec
  dropped_yr <- dropped_yr[,!names(dropped_yr) %in% "ln_rec"]
  #fit model
  modtemp <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +        
                  Smr_CPUE_juv_ADFG_ln_scaled +    
                  YOY_grwth_Middleton_scaled +
                  smr_adult_cond_scaled, data=temp_dat)
  #have model predict to missing year
  temp_predict <- predict(modtemp, newdata=dropped_yr)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict[1]
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  geom_point() + geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year))

#STEP 2 - get MSE, MAE, and R2------

#get MSE & MAE------

#these need to be double checked!
GAM_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)

GAM_MAE <- ((sum((abs(output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2), na.rm = TRUE)))/length(output_df$observed_ln_recruit)

obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)





