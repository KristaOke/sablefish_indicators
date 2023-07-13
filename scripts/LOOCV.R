#============================================================================================================================================
# Bespoke LOOCV code that will work across models

#Created by Krista, 2023
#============================================================================================================================================
#Notes: Caret doesn't support all the model types we use it seems, so let's
#code it up ourselves
#============================================================================================================================================

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#Get data =======


scaled <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)

#THIS WILL VARY BY MODEL TYPE

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

scaled_loop_dat <- scaled[,names(scaled) %in% noncor_covars3]



#STEP 1 - Loop through training sets and fit models-------

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
  #have model predict to missing year
 # temp_predict <- predict()
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  geom_point() + geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1)+
  ylim(c(0,5)) + xlim(c(0,5))

#STEP 2 - get MSE, MAE, and R2------

#get MSE & MAE------
library(yardstick)

obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
  geom_point() + geom_smooth(method="lm")

rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

mae <- mae(output_df, truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)






