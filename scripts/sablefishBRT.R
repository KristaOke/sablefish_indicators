

require(tidyverse)
require(ggthemes)
require(viridis)
require(R2jags)
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

#Get data plot data =======

dat <- read.csv(file.path(dir.data, "sablefish_BAS_indicators_2022.csv"))

scaled_dat <- dat %>% #group_by(Year) %>%
  mutate(recruit_scaled=scale(Recruitment),
         Spr_ST_SEBS_scaled=scale(Spring_Temperature_Surface_SEBS_Satellite),
         YOY_grwth_Middleton_scaled=scale(Annual_Sablefish_Growth_YOY_Middleton_Survey),
         Smr_CPUE_juv_ADFG_scaled=scale(Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey),
         spawner_mean_age_scaled=scale(Annual_Sablefish_Mean_Age_Female_Adult_Model),
         spawner_age_evenness_scaled=scale(Annual_Sablefish_Age_Evenness_Female_Adult_Model),
         arrowtooth_biomass_scaled=scale(Annual_Arrowtooth_Biomass_GOA_Model),
         sablefish_bycatch_arrowtooth_fishery_scaled=scale(Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery),
         smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Female_Adult_GOA_Survey))

brtdat <- scaled_dat
brtdat <- brtdat[is.na(brtdat$recruit_scaled)==FALSE,] #remove NAs from residuals
rownames(brtdat) <- brtdat[,1]
brtdat <- brtdat[,-1]


#Combine matrices

# Create Input Data Frame ==========================================
#dat <- data.frame(cbind(t(input.rec),t(input.cov)))


# Fit BRT Model ======================================================
# if(fit==TRUE) {

fit.covars <- names(brtdat[,11:18]) 
#fit.covars <- fit.covars[1:4] #to get things running lets subset out just a few

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(names(brtdat[10]), "~",form.covars)

brtdat <- brtdat[,10:18]
res <- names(brtdat[1])

#NEED TO SUBSET out a training dataset
#train <- brtdat[]

#see paper section on simplifying the predictor set esp imp if some predictors are repetitive and n is small

temp.fit <- gbm(formula(form), data=brtdat, distribution='gaussian', n.trees=1e4, interaction.depth=1, n.minobsinnode=1)

temp.fit.2 <- gbm.step(data=brtdat, gbm.y=res, gbm.x=fit.covars, family='gaussian', n.minobsinnode=1)
gbm.plot(temp.fit.2)
gbm.plot.fits(temp.fit.2)

temp.fit.3 <- gbm(formula(form), data=brtdat, distribution='gaussian', n.minobsinnode=1, n.trees=1e5, shrinkage=0.001, interaction.depth=1)#, interaction.depth=1, cv.folds=1, n.minobsinnode=1)
gbm.plot(temp.fit.3) #not working?

temp.fit.4 <- gbm.step(data=brtdat, gbm.y=res, gbm.x=fit.covars, family='gaussian', 
                                       tree.complexity = 1, #b/c very small sample
                                       learning.rate = 0.005, #slower b/c tc is low and want enough trees
                                       bag.fraction = 0.7) #paper recommends not fewer than 1000 trees
gbm.plot(temp.fit.4)
gbm.plot.fits(temp.fit.4)
summary(temp.fit.4)

temp.fit.5 <- gbm(formula(form), data=brtdat, distribution='gaussian', bag.fraction=0.7, n.trees=1e5, shrinkage=0.001, interaction.depth=1)#, interaction.depth=1, cv.folds=1, n.minobsinnode=1)
gbm.plot(temp.fit.5)
summary(temp.fit.5)
