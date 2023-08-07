#============================================================================================================================================
# Boosted regression tree

#Created by Krista, 2023
#borrowing heavily from code written by Curry
#============================================================================================================================================
#Notes:
#============================================================================================================================================

require(tidyverse)
require(ggthemes)
require(viridis)
#require(R2jags)
require(corrplot)
require(reshape2)
require(BEST)
# require(ggmcmc)
require(mcmcplots)
require(gbm)
require(dismo)
require(visreg)


#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")


#Get data=======

# load data that is already z-scored, checked for correlations

scaled <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)


#DATA CONTROL SECTION----

#select covariates

#using ONLY the indicator subset that's been checked for collinearity
#using noncor_covars from process_data for now

#BRT CAN handle missing data, but euphasiid and the every other yr GOA indicator were already removed
#during collinearity check


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


scaled_brt_dat <- scaled[,names(scaled) %in% noncor_covars3]

#BRT cannot include NAs in response variable
#these are often at the END of the time series

scaled_brt_dat <- scaled_brt_dat[which(is.na(scaled_brt_dat$ln_rec)==FALSE),]

#leave out 2020 b/c it's a longterm mean
scaled_brt_dat <- scaled_brt_dat[which(scaled_brt_dat$Year<2020),]



# Fit BRT Model ======================================================
# if(fit==TRUE) {

res1 <- "ln_rec"

fit.covars <- names(scaled_brt_dat[,!names(scaled_brt_dat) %in% c("Year", "ln_rec")]) 

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(res1, "~",form.covars)


#fit BRT----

#see goood tutorial on https://afit-r.github.io/tree_based_methods

real.fit.1 <- gbm.step(data=scaled_brt_dat, gbm.y=res1, gbm.x=fit.covars, family='gaussian', 
                       #real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                       tree.complexity = 3, #b/c very small sample
                       learning.rate = 0.001, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8, 
                       n.minobsinnode=3) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
real.fit.1$cv.statistics
gbm.plot(real.fit.1)
gbm.plot.fits(real.fit.1)
gbm.perspec(real.fit.1, x=1, y=9, theta = 45, phi=10) #
summary(real.fit.1)
gbm.simplify(real.fit.1, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(real.fit.1) #no interaction
perf_n <- gbm.perf(real.fit.1)[1] #optimal # trees 710 but that's pretty low

par(mar = c(5, 15, 2, 2))
summary(real.fit.1, las=1)

#painful tuning section ==============

#from real.fit.1
total_dev <- 0.73

#big loop to try all reasonable combos
counter <- 1
learningrates <- c(0.001, 0.0001, 0.00001)
families <- c("gaussian", "laplace")
output_list <- data.frame(matrix(ncol=5, nrow=10*length(learningrates)*length(families)))
colnames(output_list) <- c("deviance_explained", "n_trees", "learning_rate", "tree_complexity", "family") #
for(i in 1:10) {
  #tree complexity level
  for(j in 1:length(learningrates)){
    #learning rate level
    for(k in 1:length(families)){
  # fit model
  mod_boost <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                        data = scaled_brt_dat,
                        family = families[k],                       
                        tree.complexity = i, #
                        learning.rate = learningrates[j], #
                        bag.fraction = 0.8,
                        n.minobsinnode=3,
                        n.folds = 5,
                        max.trees = 50000)
  
  dev_ex_temp <- (total_dev - mod_boost$cv.statistics$deviance.mean)/total_dev
  output_list[counter,1] <- dev_ex_temp #should likely ask n trees and store too
  perf_n_temp <- gbm.perf(mod_boost)[1]
  output_list[counter,2] <- perf_n_temp
  output_list[counter,3] <- learningrates[j]
  output_list[counter,4] <- i
  output_list[counter,5] <- families[k]
  
  counter <- counter + 1
    }  
  }
}

output_list

bestsofar <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                      data = scaled_brt_dat,
                      family = "gaussian",                       
                      tree.complexity = 5, 
                      learning.rate = 0.0001, 
                      bag.fraction = 0.8,
                      n.minobsinnode=3,
                      n.folds = 5,
                      max.trees = 50000)
perf_n <- gbm.perf(bestsofar)[1] 
(total_dev-bestsofar$cv.statistics$deviance.mean)/total_dev

finalbrt <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                                  data = scaled_brt_dat,
                                  family = "gaussian",                       
                                  tree.complexity = 5, 
                                  learning.rate = 0.0001, 
                                  bag.fraction = 1,
                                  n.minobsinnode=3,
                                  n.folds = 5,
                                  max.trees = 50000)

plot(finalbrt$fitted, finalbrt$residuals)
finalbrt$cv.statistics
gbm.plot(finalbrt)
gbm.plot.fits(finalbrt)
gbm.perspec(finalbrt, x=1, y=9, theta = 45, phi=10) #
summary(finalbrt)
gbm.simplify(finalbrt, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(finalbrt) #


dev_ex <- (total_dev-finalbrt$cv.statistics$deviance.mean)/total_dev


#plot rel inf----

summary_out <- summary(finalbrt)

summary_out$var[which(summary_out$var=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
summary_out$var[which(summary_out$var=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
summary_out$var[which(summary_out$var=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
summary_out$var[which(summary_out$var=="smr_adult_cond_scaled")] <- "Summer adult condition"


ggplot(summary_out, aes(rel.inf, var)) + geom_col(fill="blue") + theme_bw() + ylab("Indicator") +
  xlab("Relative influence")


#plot within sample predictions----

withinpred <- predict.gbm(finalbrt, scaled_brt_dat[which(scaled_brt_dat$Year>1980),c(3:12)], n.trees = 8121)
plotwithin <- scaled_brt_dat[which(scaled_brt_dat$Year>1980),]
plotwithin$predicted <- withinpred

ggplot(plotwithin, aes(Year, ln_rec)) + geom_point(aes(col="red")) + 
  geom_line() + geom_point(aes(Year, predicted)) + geom_line(aes(Year, predicted)) + theme_bw()


#from BAS
# Plot Model Predictions vs. Observed ==============================
#pdf(file.path(dir.figs,"Model Fit.pdf"), height=6, width=9)
par(oma=c(1,1,1,1), mar=c(4,4,1,1), mfrow=c(1,2))

# Omit NAs
dat.temp <- plotwithin

plot(x=dat.temp$ln_rec, y=plotwithin$predicted,
     xlab="Observed ln(Recruitment)", ylab="Predicted ln(Recruitment)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=plotwithin$fit, y=plotwithin$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)

# Timeseries
plot(x=dat.temp$Year, y=dat.temp$ln_rec,
     xlab="Year", ylab="ln(Recruitment)", type='l', col=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
grid(lty=3, col='dark gray')
points(x=dat.temp$Year, y=dat.temp$ln_rec,
       pch=21, bg=rgb(1,0,0,alpha=0.5))
lines(x=dat.temp$Year, y=plotwithin$predicted, lwd=3, col=rgb(0,0,1, alpha=0.5))
#points(x=dat.temp$Year, y=plotwithin$predicted,
#      pch=21, bg=rgb(0,1,0,alpha=0.5))

legend('topleft', legend=c("Observed","Predicted"), lty=1, col=c(rgb(1,0,0,alpha=0.5),
                                                                 rgb(0,0,1, alpha=0.5)), bg="white")





#REPEAT ON JUST 4 LONG TIME SERIES SET=====

# Fit BRT Model ======================================================


res1 <- "ln_rec"

fit.covars <- names(scaled_brt_dat[,names(scaled_brt_dat) %in% c("Spr_ST_SEBS_scaled", "Smr_CPUE_juv_ADFG_ln_scaled",
                                                                 "YOY_grwth_Middleton_scaled", "smr_adult_cond_scaled")]) 

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(res1, "~",form.covars)


#see goood tutorial on https://afit-r.github.io/tree_based_methods

real.fit.1 <- gbm.step(data=scaled_brt_dat, gbm.y=res1, gbm.x=fit.covars, family='gaussian', 
                       #real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                       tree.complexity = 3, #b/c very small sample
                       learning.rate = 0.001, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8, 
                       n.minobsinnode=3) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
total_dev <- 0.73

real.fit.1$cv.statistics
gbm.plot(real.fit.1)
gbm.plot.fits(real.fit.1)
gbm.perspec(real.fit.1, x=1, y=3, theta = 45, phi=10) #
summary(real.fit.1)
gbm.simplify(real.fit.1, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(real.fit.1) #THERE IS A SST x ADFG interaction
perf_n <- gbm.perf(real.fit.1)[1] #693 too low

par(mar = c(5, 15, 2, 2))
summary(real.fit.1, las=1)

#painful tuning section ==============

#from real.fit.1
total_dev <- 0.73

#big loop to try all reasonable combos
counter <- 1
learningrates <- c(0.001, 0.0001, 0.00001)
families <- c("gaussian", "laplace")
output_list <- data.frame(matrix(ncol=5, nrow=10*length(learningrates)*length(families)))
colnames(output_list) <- c("deviance_explained", "n_trees", "learning_rate", "tree_complexity", "family") #
for(i in 1:10) {
  #tree complexity level
  for(j in 1:length(learningrates)){
    #learning rate level
    for(k in 1:length(families)){
      # fit model
      mod_boost <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                            data = scaled_brt_dat,
                            family = families[k],                       
                            tree.complexity = i, #
                            learning.rate = learningrates[j], #
                            bag.fraction = 0.8,
                            n.minobsinnode=3,
                            n.folds = 5,
                            max.trees = 50000)
      
      dev_ex_temp <- (total_dev - mod_boost$cv.statistics$deviance.mean)/total_dev
      output_list[counter,1] <- dev_ex_temp #should likely ask n trees and store too
      perf_n_temp <- gbm.perf(mod_boost)[1]
      output_list[counter,2] <- perf_n_temp
      output_list[counter,3] <- learningrates[j]
      output_list[counter,4] <- i
      output_list[counter,5] <- families[k]
      
      counter <- counter + 1
    }  
  }
}

output_list

bestsofar <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                      data = scaled_brt_dat,
                      family = "gaussian",                       
                      tree.complexity = 7, 
                      learning.rate = 0.0001, 
                      bag.fraction = 0.8,
                      n.minobsinnode=3,
                      n.folds = 5,
                      max.trees = 50000)
perf_n <- gbm.perf(bestsofar)[1] #8063
(total_dev-bestsofar$cv.statistics$deviance.mean)/total_dev

finalbrt <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                     data = scaled_brt_dat,
                     family = "gaussian",                       
                     tree.complexity = 7, 
                     learning.rate = 0.0001, 
                     bag.fraction = 1,
                     n.minobsinnode=3,
                     n.folds = 5,
                     max.trees = 50000)


plot(finalbrt$fitted, finalbrt$residuals)
finalbrt$cv.statistics
gbm.plot(finalbrt)
gbm.plot.fits(finalbrt)
gbm.perspec(finalbrt, x=1, y=3, theta = 45, phi=10, z.range=c(0,4)) #
summary(finalbrt)
gbm.simplify(finalbrt, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(finalbrt) #interaction between SST and CPUE ADFG


dev_ex <- (total_dev-finalbrt$cv.statistics$deviance.mean)/total_dev

#plot rel inf----

summary_out <- summary(finalbrt)

summary_out$var[which(summary_out$var=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
summary_out$var[which(summary_out$var=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
summary_out$var[which(summary_out$var=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
summary_out$var[which(summary_out$var=="smr_adult_cond_scaled")] <- "Summer adult condition"


ggplot(summary_out, aes(rel.inf, var)) + geom_col(fill="blue") + theme_bw() + ylab("Indicator") +
  xlab("Relative influence")

#plot within sample predictions----

withinpred <- predict.gbm(finalbrt, scaled_brt_dat[which(scaled_brt_dat$Year>1980),c(3:12)], n.trees = 8121)
plotwithin <- scaled_brt_dat[which(scaled_brt_dat$Year>1980),]
plotwithin$predicted <- withinpred

ggplot(plotwithin, aes(Year, ln_rec)) + geom_point(aes(col="red")) + 
  geom_line() + geom_point(aes(Year, predicted)) + geom_line(aes(Year, predicted)) + theme_bw()


#from BAS
# Plot Model Predictions vs. Observed ==============================
#pdf(file.path(dir.figs,"Model Fit.pdf"), height=6, width=9)
par(oma=c(1,1,1,1), mar=c(4,4,1,1), mfrow=c(1,2))

# Omit NAs
dat.temp <- plotwithin

plot(x=dat.temp$ln_rec, y=plotwithin$predicted,
     xlab="Observed ln(Recruitment)", ylab="Predicted ln(Recruitment)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=plotwithin$fit, y=plotwithin$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)

# Timeseries
plot(x=dat.temp$Year, y=dat.temp$ln_rec,
     xlab="Year", ylab="ln(Recruitment)", type='l', col=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
grid(lty=3, col='dark gray')
points(x=dat.temp$Year, y=dat.temp$ln_rec,
       pch=21, bg=rgb(1,0,0,alpha=0.5))
lines(x=dat.temp$Year, y=plotwithin$predicted, lwd=3, col=rgb(0,0,1, alpha=0.5))
#points(x=dat.temp$Year, y=plotwithin$predicted,
#      pch=21, bg=rgb(0,1,0,alpha=0.5))

legend('topleft', legend=c("Observed","Predicted"), lty=1, col=c(rgb(1,0,0,alpha=0.5),
                                                                 rgb(0,0,1, alpha=0.5)), bg="white")




#LOOCV====

scaled_loop_dat <- scaled_brt_dat
fit.covars <- names(scaled_brt_dat[,!names(scaled_brt_dat) %in% c("Year", "ln_rec")]) 


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
  temp_mod <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
           data = scaled_brt_dat,
           family = "gaussian",                       
           tree.complexity = 5, 
           learning.rate = 0.0001, 
           bag.fraction = 1,
           n.minobsinnode=3,
           n.folds = 5,
           max.trees = 50000)
  temp_perf <- gbm.perf(bestsofar)[1] 
  
  #have model predict to missing year
 temp_predict <-  predict.gbm(temp_mod, dropped_yr, n.trees=temp_perf)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  geom_point() + geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1)+ 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year)) +
  ylim(c(0,5)) + xlim(c(0,5))

#STEP 2 - get MSE, MAE, and R2------
library(yardstick)
#get MSE & MAE------

#these need to be double checked!
# BRT_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)
# 

obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

BRT_rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

BRT_mae <- mae(output_df, truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)

obs_pred_modpost2002 <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df[which(output_df$Year>2001),])
summary(obs_pred_modpost2002)

BRT_rmsepost2002 <- rmse(output_df[which(output_df$Year>2001),], truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

BRT_maepost2002 <- mae(output_df[which(output_df$Year>2001),], truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)


output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
  geom_point() + geom_smooth(method="lm")

write.csv(output_df, file=paste(wd,"/data/BRT_obsvpreds_reduced.csv", sep=""))

output_df <- read.csv(file=paste(wd,"/data/BRT_obsvpreds_reduced.csv", sep=""))

#LOOCV long time series====

scaled_loop_dat <- scaled_brt_dat
fit.covars <- names(scaled_brt_dat[,names(scaled_brt_dat) %in% c("Spr_ST_SEBS_scaled", "Smr_CPUE_juv_ADFG_ln_scaled",
                                                                 "YOY_grwth_Middleton_scaled", "smr_adult_cond_scaled")]) 



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
  # temp_mod <- gbm(formula(form2), data=temp_dat, distribution='gaussian', 
  #      n.trees=10000, interaction.depth=1, n.minobsinnode=1,
  #      bag.fraction = 0.7,
  #      shrinkage=0.001)
  #whole reduced data mod
  temp_mod <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                       data = scaled_brt_dat,
                       family = "gaussian",                       
                       tree.complexity = 5, 
                       learning.rate = 0.0001, 
                       bag.fraction = 1,
                       n.minobsinnode=3,
                       n.folds = 5,
                       max.trees = 50000)
  temp_perf <- gbm.perf(bestsofar)[1] 
  
  #have model predict to missing year
  temp_predict <-  predict.gbm(temp_mod, dropped_yr, n.trees=temp_perf)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1)+ 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year)) +
  ylim(c(0,5)) + xlim(c(0,5)) + theme_bw()

#STEP 2 - get MSE, MAE, and R2------
library(yardstick)
#get MSE & MAE------

#these need to be double checked!
# BRT_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)
# 

obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

BRT_rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

BRT_mae <- mae(output_df, truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)

#compare some years

obs_pred_mod96 <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df[which(output_df$Year>1995),])
summary(obs_pred_mod96)

BRT_rmse96 <- rmse(output_df[which(output_df$Year>1995),], truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

BRT_mae96 <- mae(output_df[which(output_df$Year>1995),], truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)





output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
  geom_point() + geom_smooth(method="lm")

write.csv(output_df, file=paste(wd,"/data/BRT_obsvpreds_long.csv", sep=""))

output_df <- read.csv(file=paste(wd,"/data/BRT_obsvpreds_long.csv", sep=""))

