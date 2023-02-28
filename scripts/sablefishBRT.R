

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


#NEED TO SUBSET out a training dataset
datlen <- length(brtdat$recruit_scaled)
train <- brtdat[sample(nrow(brtdat),(round(datlen*0.8))),]
train <- train[order(row.names(train)),]

testing <- anti_join(brtdat, train)

#save this training set?

rownames(train) <- train[,1]
train <- train[,-1]

rownames(testing) <- testing[,1]
testing <- testing[,-1]


#Combine matrices

# Create Input Data Frame ==========================================
#dat <- data.frame(cbind(t(input.rec),t(input.cov)))


# Fit BRT Model ======================================================
# if(fit==TRUE) {


train <- train[,10:18]
res <- names(train[1])

testing <- testing[,10:18]

fit.covars <- names(train[,c(2:4, 7:9)]) #NO mean age or evenness see below
#fit.covars <- fit.covars[1:4] #to get things running lets subset out just a few

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(names(train[1]), "~",form.covars)


#getting worried about mean age, is it important because it is - by definition - autocorrelated (w 2 yr lag)
#w recruitment?? PLus it is very odd smooth pattern

#try without mean age or evenness




#playing with parameters------

#see paper section on simplifying the predictor set esp imp if some predictors are repetitive and n is small

temp.fit <- gbm(formula(form), data=train, distribution='gaussian', n.trees=1e4, interaction.depth=1, n.minobsinnode=1)

temp.fit.2 <- gbm.step(data=train, gbm.y=res, gbm.x=fit.covars, family='gaussian', n.minobsinnode=1)
gbm.plot(temp.fit.2)
gbm.plot.fits(temp.fit.2)

temp.fit.3 <- gbm(formula(form), data=train, distribution='gaussian', n.minobsinnode=1, n.trees=1e5, shrinkage=0.001, interaction.depth=1)#, interaction.depth=1, cv.folds=1, n.minobsinnode=1)
gbm.plot(temp.fit.3) #not working?

temp.fit.4 <- gbm.step(data=train, gbm.y=res, gbm.x=fit.covars, family='gaussian', 
                                       tree.complexity = 1, #b/c very small sample
                                       learning.rate = 0.005, #slower b/c tc is low and want enough trees
                                       bag.fraction = 0.7) #paper recommends not fewer than 1000 trees
gbm.plot(temp.fit.4)
gbm.plot.fits(temp.fit.4)
summary(temp.fit.4)

temp.fit.5 <- gbm(formula(form), data=train, distribution='gaussian', bag.fraction=0.7, n.trees=1e5, shrinkage=0.001, interaction.depth=1)#, interaction.depth=1, cv.folds=1, n.minobsinnode=1)
gbm.plot(temp.fit.5)
summary(temp.fit.5)


res <- names(brtdat[1])

temp.fit.6 <- gbm.step(data=train, gbm.y=res, gbm.x=fit.covars, family='gaussian', 
                       tree.complexity = 1, #b/c very small sample
                       learning.rate = 0.005, #slower b/c tc is low and want enough trees
                       bag.fraction = 0.7) #paper recommends not fewer than 1000 trees
gbm.plot(temp.fit.6)
gbm.plot.fits(temp.fit.6)
summary(temp.fit.6)


#for real on training dataset-----

#see goood tutorial on https://afit-r.github.io/tree_based_methods

res <- names(train[1])

real.fit.1 <- gbm.step(data=train, gbm.y=res, gbm.x=fit.covars, family='gaussian', 
#real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                       tree.complexity = 1, #b/c very small sample
                       learning.rate = 0.005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8,
                       n.minobsinnode=1) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
gbm.plot(real.fit.1)
gbm.plot.fits(real.fit.1)
gbm.perspec(real.fit.1, x=1, y=2) #not working
summary(real.fit.1)
gbm.simplify(real.fit.1, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(real.fit.1)
perf_n <- gbm.perf(real.fit.1)[1] #optimal # trees 264 but that's pretty low

preds <- predict.gbm(real.fit.1, testing, n.trees = perf_n, type="response")

dev <- calc.deviance(obs=testing$recruit_scaled, pred=testing$predictions, family="gaussian", 
              calc.mean=TRUE)


testing$predictions <- preds
ggplot(testing, aes(predictions, recruit_scaled)) + geom_point() + geom_abline()

sum_squared_errors <- sum((preds-testing$recruit_scaled)^2, na.rm=TRUE)

#can I find a covar that's messing everything up? No
real.fit.2 <- gbm.step(data=train, gbm.y=10, gbm.x=c(11:13,16:18), family='gaussian', 
                       tree.complexity = 1, #b/c very small sample
                       learning.rate = 0.005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8,
                       n.minobsinnode=1,
                       interaction.depth=3) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
gbm.plot(real.fit.2)
gbm.plot.fits(real.fit.2)
gbm.perspec(real.fit.2, x=1, y=2) #not working
summary(real.fit.2)
gbm.simplify(real.fit.2, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(real.fit.2) #interactions continue to not work
