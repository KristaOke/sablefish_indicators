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

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  #fit model
}

#STEP 2 - get MSE, MAE, and R2------






