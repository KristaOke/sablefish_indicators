

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

#dat <- read.csv(file.path(dir.data, "sablefish_BAS_indicators_2022.csv"))

# train1 <- read.csv(file=paste(wd,"/data/dataset_training1.csv", sep=""))
# train2 <- read.csv(file=paste(wd,"/data/dataset_training2.csv", sep=""))
# train3 <- read.csv(file=paste(wd,"/data/dataset_training3.csv", sep=""))
# train4 <- read.csv(file=paste(wd,"/data/dataset_training4.csv", sep=""))
# train5 <- read.csv(file=paste(wd,"/data/dataset_training5.csv", sep=""))
# 
# testing1 <- read.csv(file=paste(wd,"/data/dataset_testing1.csv", sep=""))
# testing2 <- read.csv(file=paste(wd,"/data/dataset_testing2.csv", sep=""))
# testing3 <- read.csv(file=paste(wd,"/data/dataset_testing3.csv", sep=""))
# testing4 <- read.csv(file=paste(wd,"/data/dataset_testing4.csv", sep=""))
# testing5 <- read.csv(file=paste(wd,"/data/dataset_testing5.csv", sep=""))
# 
# 
# scaled_dat <- dat %>% #group_by(Year) %>%
#   mutate(recruit_scaled=scale(Recruitment),
#          Spr_ST_SEBS_scaled=scale(Spring_Temperature_Surface_SEBS_Satellite),
#          YOY_grwth_Middleton_scaled=scale(Annual_Sablefish_Growth_YOY_Middleton_Survey),
#          Smr_CPUE_juv_ADFG_scaled=scale(Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey),
#          spawner_mean_age_scaled=scale(Annual_Sablefish_Mean_Age_Female_Adult_Model),
#          spawner_age_evenness_scaled=scale(Annual_Sablefish_Age_Evenness_Female_Adult_Model),
#          arrowtooth_biomass_scaled=scale(Annual_Arrowtooth_Biomass_GOA_Model),
#          sablefish_bycatch_arrowtooth_fishery_scaled=scale(Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery),
#          smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Female_Adult_GOA_Survey))
# 
# brtdat <- scaled_dat
# brtdat <- brtdat[is.na(brtdat$recruit_scaled)==FALSE,] #remove NAs from residuals
# 
# 
# #NEED TO SUBSET out a training dataset
# datlen <- length(brtdat$recruit_scaled)
# train <- brtdat[sample(nrow(brtdat),(round(datlen*0.8))),]
# train <- train[order(row.names(train)),]
# 
# testing <- anti_join(brtdat, train)
# 
# #save this training set?
# 
# rownames(train) <- train[,1]
# train <- train[,-1]
# 
# rownames(testing) <- testing[,1]
# testing <- testing[,-1]


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


#try different interaction depths, compare SSE
#loop isn't working yet
(IntDepthG <- sapply(1:8, function(x) {
  # fit model
  mod_boost <- gbm.step(gbm.y=res, gbm.x=fit.covars,
                      data = train,
                      family = "gaussian",                       
                      tree.complexity = 1, #b/c very small sample
                      learning.rate = 0.005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                      bag.fraction = 0.8,
                      n.minobsinnode=1,
                      n.folds = 5,
                      interaction.depth=x)
  # calculate best value for number of trees
  bestHit <- gbm.perf(mod_boost)[1]
  # calculate the sum of squared errors
  boostSSE <- sum((predict(mod_boost, testing, n.trees = bestHit) - testing$recruit_scaled)^2,
                  na.rm = TRUE)
  return(boostSSE)}))







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



#NEW========================================================================================================

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


train1_brt_dat <- train1[,names(train1) %in% noncor_covars3]
train2_brt_dat <- train2[,names(train2) %in% noncor_covars3]
train3_brt_dat <- train3[,names(train3) %in% noncor_covars3]
train4_brt_dat <- train4[,names(train4) %in% noncor_covars3]
train5_brt_dat <- train5[,names(train5) %in% noncor_covars3]

test1_brt_dat <- testing1[,names(testing1) %in% noncor_covars3]
test2_brt_dat <- testing2[,names(testing2) %in% noncor_covars3]
test3_brt_dat <- testing3[,names(testing3) %in% noncor_covars3]
test4_brt_dat <- testing4[,names(testing4) %in% noncor_covars3]
test5_brt_dat <- testing5[,names(testing5) %in% noncor_covars3]

scaled_brt_dat <- scaled[,names(scaled) %in% noncor_covars3]

#BRT cannot include NAs in response variable
#these are often at the END of the time series

train1_brt_dat <- train1_brt_dat[which(is.na(train1_brt_dat$ln_rec)==FALSE),]
train2_brt_dat <- train2_brt_dat[which(is.na(train2_brt_dat$ln_rec)==FALSE),]
train3_brt_dat <- train3_brt_dat[which(is.na(train3_brt_dat$ln_rec)==FALSE),]
train4_brt_dat <- train4_brt_dat[which(is.na(train4_brt_dat$ln_rec)==FALSE),]
train5_brt_dat <- train5_brt_dat[which(is.na(train5_brt_dat$ln_rec)==FALSE),]

test1_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
test2_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
test3_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
test4_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
test5_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]

scaled_brt_dat <- scaled_brt_dat[which(is.na(scaled_brt_dat$ln_rec)==FALSE),]


# Fit BRT Model ======================================================
# if(fit==TRUE) {

res1 <- "ln_rec"

fit.covars <- names(train1_brt_dat[,!names(train1_brt_dat) %in% c("Year", "ln_rec")]) 

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(res1, "~",form.covars)


#fit BRT on training set 1-----

#see goood tutorial on https://afit-r.github.io/tree_based_methods

real.fit.1 <- gbm.step(data=train1_brt_dat, gbm.y=res1, gbm.x=fit.covars, family='gaussian', 
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
perf_n <- gbm.perf(real.fit.1)[1] #optimal # trees 139 but that's pretty low

par(mar = c(5, 15, 2, 2))
summary(real.fit.1, las=1)

preds <- predict.gbm(real.fit.1, test1_brt_dat, n.trees = perf_n, type="response")

dev <- calc.deviance(obs=test1_brt_dat$ln_rec, pred=test1_brt_dat$predictions, family="gaussian", 
                     calc.mean=TRUE)


test1_brt_dat$predictions <- preds
ggplot(test1_brt_dat, aes(ln_rec, predictions)) + geom_point() + geom_abline() + theme_bw() + ylab("Predicted ln(Recruitment)") +
  xlab("Observed ln(Recruitment)")

sum_squared_errors <- sum((preds-test1_brt_dat$ln_rec)^2, na.rm=TRUE)

#plot within sample predictions----

withinpred <- predict.gbm(real.fit.1, train1_brt_dat[which(train1_brt_dat$Year>1980),c(3:12)])
plotwithin <- train1_brt_dat[which(train1_brt_dat$Year>1980),]
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



#try different interaction depths, compare SSE
(IntDepthG <- sapply(1:8, function(x) {
  # fit model
  mod_boost <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                        data = train1_brt_dat,
                        family = "gaussian",                       
                        tree.complexity = 1, #b/c very small sample
                        learning.rate = 0.00005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                        bag.fraction = 0.8,
                        n.minobsinnode=1,
                        n.folds = 5,
                        interaction.depth=x)
  # calculate best value for number of trees
  bestHit <- gbm.perf(mod_boost)[1]
  # calculate the sum of squared errors
  boostSSE <- sum((predict(mod_boost, test1_brt_dat, n.trees = bestHit) - test1_brt_dat$ln_rec)^2,
                  na.rm = TRUE)
  return(boostSSE)}))

#changing interaction depth doesn't make big difference, depth of 2 is best at 5.264870 but
#not much different than the original sse at 5.266896



#try different folds
(foldsG <- sapply(5:10, function(x) {
  # fit model
  mod_boost <- gbm.step(gbm.y=res1, gbm.x=fit.covars,
                        data = train1_brt_dat,
                        family = "gaussian",                       
                        tree.complexity = 1, #b/c very small sample
                        learning.rate = 0.000005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                        bag.fraction = 0.8,
                        n.minobsinnode=1,
                        n.folds = x,
                        interaction.depth=2)
  # calculate best value for number of trees
  bestHit <- gbm.perf(mod_boost)[1]
  # calculate the sum of squared errors
  boostSSE <- sum((predict(mod_boost, test1_brt_dat, n.trees = bestHit) - test1_brt_dat$ln_rec)^2,
                  na.rm = TRUE)
  return(boostSSE)}))
#under 5 folds won't run
#5 to 10 folds makes almost no difference





#REPEAT ON REDUCED SET=====

# Fit BRT Model ======================================================
# if(fit==TRUE) {

res1 <- "ln_rec"

fit.covars <- names(train1_brt_dat[,names(train1_brt_dat) %in% c("Spr_ST_SEBS_scaled", "Smr_CPUE_juv_ADFG_ln_scaled",
                                                                 "YOY_grwth_Middleton_scaled", "smr_adult_cond_scaled")]) 

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(res1, "~",form.covars)


#fit BRT on training set 1-----

#see goood tutorial on https://afit-r.github.io/tree_based_methods

real.fit.1 <- gbm.step(data=train1_brt_dat, gbm.y=2, gbm.x=c(3,11,10,12), family='gaussian', 
                       #real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                       tree.complexity = 1, #b/c very small sample
                       learning.rate = 0.005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8,
                       n.minobsinnode=1) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
gbm.plot(real.fit.1)
gbm.plot.fits(real.fit.1)
gbm.perspec(real.fit.1, x=1, y=8) #not working
summary(real.fit.1)
gbm.simplify(real.fit.1, n.drops = 5) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(real.fit.1)
perf_n <- gbm.perf(real.fit.1)[1] #limit n trees because overfitting is a concern for BRT
#optimal # trees 143 but that's pretty low

preds <- predict.gbm(real.fit.1, test1_brt_dat, n.trees = perf_n, type="response")

dev <- calc.deviance(obs=test1_brt_dat$ln_rec, pred=test1_brt_dat$predictions, family="gaussian", 
                     calc.mean=TRUE)


test1_brt_dat$predictions <- preds
ggplot(test1_brt_dat, aes(predictions, ln_rec)) + geom_point() + geom_abline()

sum_squared_errors <- sum((preds-test1_brt_dat$ln_rec)^2, na.rm=TRUE)

#refit with smallest n trees from gbm.perf
refit.1 <- gbm(formula(form), data=train1_brt_dat, distribution='gaussian', 
               n.trees=143, interaction.depth=1, n.minobsinnode=1)

gbm.plot(refit.1)
gbm.plot.fits(refit.1)
gbm.perspec(refit.1, x=1, y=8) #not working
summary(refit.1)

plot(refit.1)

refit.1.1 <- gbm.step(data=train1_brt_dat, gbm.y=res1, gbm.x=fit.covars, family='gaussian', 
                       #real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                      # tree.complexity = 1, #b/c very small sample
                      max.trees = 30000,
                       learning.rate = 0.0005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8,
                       n.minobsinnode=1) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
gbm.plot(refit.1.1)
gbm.plot.fits(refit.1.1)
gbm.perspec(refit.1.1, x=1, y=8) #not working
summary(refit.1.1)
gbm.simplify(refit.1.1, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(refit.1.1)


perf_n1.1 <- gbm.perf(refit.1.1)[1] #limit n trees because overfitting is a concern for BRT
#optimal # trees 1018

preds <- predict.gbm(refit.1.1, test1_brt_dat, n.trees = perf_n, type="response")

dev <- calc.deviance(obs=test1_brt_dat$ln_rec, pred=test1_brt_dat$predictions, family="gaussian", 
                     calc.mean=TRUE)


test1_brt_dat$predictions <- preds
ggplot(test1_brt_dat, aes(predictions, ln_rec)) + geom_point() + geom_abline()

sum_squared_errors <- sum((preds-test1_brt_dat$ln_rec)^2, na.rm=TRUE)

refit.1.1.1 <- gbm(formula(form), data=train1_brt_dat, distribution='gaussian', 
               n.trees=1018, interaction.depth=1, n.minobsinnode=1)

gbm.plot(refit.1.1.1)
gbm.plot.fits(refit.1.1.1)
gbm.perspec(refit.1.1.1, x=1, y=8) #not working
summary(refit.1.1.1)

#plot within sample predictions----

withinpred <- predict.gbm(real.fit.1, train1_brt_dat[,c(3,11,10,12)])
plotwithin <- train1_brt_dat
plotwithin$predicted <- withinpred

ggplot(plotwithin, aes(Year, ln_rec)) + geom_point(aes(col="red")) + 
  geom_line() + geom_point(aes(Year, predicted)) + geom_line(aes(Year, predicted)) + theme_bw()
                          

#need to also check if different distributions give better SSE


ggplot(train1_brt_dat, aes(Spr_ST_SEBS_scaled, ln_rec, col=as.factor(Year))) + geom_point()
ggplot(train1_brt_dat, aes(Spr_ST_SEBS_scaled, ln_rec, col=as.factor(Year))) + geom_point() +geom_text(aes(Spr_ST_SEBS_scaled, ln_rec,label=Year))

ggplot(train1_brt_dat, aes(Spr_ST_SEBS_scaled, Smr_CPUE_juv_ADFG_ln_scaled, col=ln_rec)) + geom_point() +geom_text(aes(Spr_ST_SEBS_scaled, Smr_CPUE_juv_ADFG_ln_scaled,label=Year))

ggplot(train1_brt_dat, aes(Smr_CPUE_juv_ADFG_ln_scaled, ln_rec)) + geom_point()

ggplot(train1_brt_dat, aes(Smr_temp_250m_GOA_scaled, ln_rec)) + geom_point()

ggplot(train1_brt_dat, aes(Year, Smr_temp_250m_GOA_scaled)) + geom_point()

ggplot(train1, aes(ann_heatwave_GOA_scaled, Smr_CPUE_juv_ADFG_ln_scaled, col=ln_rec)) + geom_point() +geom_text(aes(ann_heatwave_GOA_scaled, Smr_CPUE_juv_ADFG_ln_scaled,label=Year))

ggplot(train1, aes(ann_heatwave_GOA_scaled, ln_rec)) + geom_point() +geom_text(aes(ann_heatwave_GOA_scaled, ln_rec,label=Year))

ggplot(train1, aes(ann_heatwave_GOA_scaled, Spr_ST_SEBS_scaled, col=ln_rec)) + geom_point()



#REPEAT without shorter time series that aren't likely to be included===========================

#both "Smr_condition_fem_age4_GOA_scaled" and "Smr_temp_250m_GOA_scaled"
#are shorter than others and have low influence

# Fit BRT Model ======================================================
# if(fit==TRUE) {

res1 <- "ln_rec"

fit.covars <- names(train1_brt_dat[,!names(train1_brt_dat) %in% c("Year", "ln_rec", 
                                  "Smr_temp_250m_GOA_scaled", "Smr_condition_fem_age4_GOA_scaled")]) 

form.covars <- paste(fit.covars, collapse=" + ")
form <- paste(res1, "~",form.covars)


#fit BRT on training set 1-----

#see goood tutorial on https://afit-r.github.io/tree_based_methods

real.fit.1 <- gbm.step(data=train1_brt_dat, gbm.y=2, gbm.x=3:12, family='gaussian', 
                       #real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                       tree.complexity = 1, #b/c very small sample
                       learning.rate = 0.005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                       bag.fraction = 0.8,
                       n.minobsinnode=1) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
gbm.plot(real.fit.1)
gbm.plot.fits(real.fit.1)
gbm.perspec(real.fit.1, x=1, y=8) #not working
summary(real.fit.1)
gbm.simplify(real.fit.1, n.drops = 5) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(real.fit.1)
perf_n <- gbm.perf(real.fit.1)[1] #limit n trees because overfitting is a concern for BRT
#optimal # trees 143 but that's pretty low

preds <- predict.gbm(real.fit.1, test1_brt_dat, n.trees = perf_n, type="response")

dev <- calc.deviance(obs=test1_brt_dat$ln_rec, pred=test1_brt_dat$predictions, family="gaussian", 
                     calc.mean=TRUE)


test1_brt_dat$predictions <- preds
ggplot(test1_brt_dat, aes(predictions, ln_rec)) + geom_point() + geom_abline()

sum_squared_errors <- sum((preds-test1_brt_dat$ln_rec)^2, na.rm=TRUE)

#refit with smallest n trees from gbm.perf
refit.1 <- gbm(formula(form), data=train1_brt_dat, distribution='gaussian', 
               n.trees=76, interaction.depth=1, n.minobsinnode=1)

gbm.plot(refit.1)
gbm.plot.fits(refit.1)
gbm.perspec(refit.1, x=1, y=8) #not working
summary(refit.1)

plot(refit.1)

refit.1.1 <- gbm.step(data=train1_brt_dat, gbm.y=res1, gbm.x=fit.covars, family='gaussian', 
                      #real.fit.1 <- gbm.step(data=train, gbm.y=10, gbm.x=11:18, family='gaussian', 
                      # tree.complexity = 1, #b/c very small sample
                      max.trees = 30000,
                      learning.rate = 0.0005, #slower b/c tc is low and want enough trees, paper recommends not fewer than 1000 trees
                      bag.fraction = 0.8,
                      n.minobsinnode=1) #won't run with bag fraction lower than 0.8 for training dataset of 34 years
gbm.plot(refit.1.1)
gbm.plot.fits(refit.1.1)
gbm.perspec(refit.1.1, x=1, y=8) #not working
summary(refit.1.1)
gbm.simplify(refit.1.1, n.drops = 3) #error that nTrain * bag.fraction <= n.minobsinnode` but it is not!!
gbm.interactions(refit.1.1)


perf_n1.1 <- gbm.perf(refit.1.1)[1] #limit n trees because overfitting is a concern for BRT
#optimal # trees 1007

preds <- predict.gbm(refit.1.1, test1_brt_dat, n.trees = perf_n, type="response")

dev <- calc.deviance(obs=test1_brt_dat$ln_rec, pred=test1_brt_dat$predictions, family="gaussian", 
                     calc.mean=TRUE)


test1_brt_dat$predictions <- preds
ggplot(test1_brt_dat, aes(predictions, ln_rec)) + geom_point() + geom_abline()

sum_squared_errors <- sum((preds-test1_brt_dat$ln_rec)^2, na.rm=TRUE)

refit.1.1.1 <- gbm(formula(form), data=train1_brt_dat, distribution='gaussian', 
                   n.trees=1018, interaction.depth=1, n.minobsinnode=1)

gbm.plot(refit.1.1.1)
gbm.plot.fits(refit.1.1.1)
gbm.perspec(refit.1.1.1, x=1, y=8) #not working
summary(refit.1.1.1)











