#============================================================================================================================================
# DFA using BAYES DFA

#Created by Krista, May 2023
#============================================================================================================================================
#Notes: following vignettes at: https://fate-ewi.github.io/bayesdfa/articles/a1_bayesdfa.html
#============================================================================================================================================
library(tidyverse)
library(cowplot)
library(AKesp)
library(corrplot)
library(bayesdfa)
library(rstan)

chains = 3
iter = 1000

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")




#Get data=======

# load data

scaled <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)

#DATA CONTROL SECTION----

#select the z-scored columns, no recruitment

scaled_dfa_dat <- scaled[,c(1,24:43)]

#remove problem covariates=====

#DFA can handle missing data but there is still too much missing in euphasiid data
#remove the mean age and evenness indicators
#also remove arrowtooth biomass b/c it's not converging even in best model
#same reason dropping sebs ST


scaled_dfa_dat <- scaled_dfa_dat[,!names(scaled_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled",
                                                                "Smr_condition_fem_age4_GOA_scaled")]



#manupulate data for DFA========

#check to make sure years are in order!

scaled_dfa_dat$Year==scaled_dfa_dat$Year[order(scaled_dfa_dat$Year)] #SHOULD BE ALL TRUE

#put data into wide format
b.mat1 <- t(as.matrix(scaled_dfa_dat))
colnames(b.mat1) <- b.mat1[1,]
b.mat1 <- b.mat1[-1,]

#fit first model
f1 <- fit_dfa(
  y = b.mat1, num_trends = 1, scale="none", #already z-scored
  iter = iter, chains = chains, thin = 1
)

#check convergence
is_converged(f1, threshold = 1.05)

#rotate trends
r <- rotate_trends(f1)

#plot trend
plot_trends(r) + theme_bw() #if invalid graphics state error try dev.off()

#plot loadings
plot_fitted(f1) + theme_bw()
plot_loadings(r) + theme_bw()

#LOO
loo1 <- loo(f1)
loo1$estimates


#now multiple models at once=======

#let's try 1:5 trends, normal and Student-t distribution, and equal and unequal varianes
m <- find_dfa_trends(
  y = b.mat1, iter = iter,
  kmin = 1, kmax = 5, chains = chains, compare_normal = TRUE,
  variance = c("equal", "unequal"), estimate_trend_ar=TRUE
)

#seems like best is model 13, 3 trends, equal, normal, let's fit that to look at?
names(m)

m$best_model


#check convergence
is_converged(m$best_model, threshold = 1.05)

#rotate trends
r <- rotate_trends(m$best_model)

#plot trend
plot_trends(r) + theme_bw() #if invalid graphics state error try dev.off()

#plot loadings
plot_fitted(m$best_model) + theme_bw()
plot_loadings(r) + theme_bw()



#let bayesdfa do scaling======



unscaled_dfa_dat <- scaled[,c(1:17, 19:22)]

#remove problem covariates=====

#DFA can handle missing data but there is still too much missing in euphasiid data
#remove the mean age and evenness indicators
#also remove arrowtooth biomass b/c it's not converging even in best model
#same reason dropping sebs ST


unscaled_dfa_dat <- unscaled_dfa_dat[,!names(unscaled_dfa_dat) %in% c("Summer_Euphausiid_Abundance_Kodiak_Survey",
                                                                "Annual_Sablefish_Mean_Age_Female_Adult_Model",
                                                                "Annual_Sablefish_Age_Evenness_Female_Adult_Model",
                                                                "Annual_Arrowtooth_Biomass_GOA_Model",
                                                                "Summer_Sablefish_Condition_Female_Age4_GOA_Survey")]



#manupulate data for DFA========

#check to make sure years are in order!

unscaled_dfa_dat$Year==unscaled_dfa_dat$Year[order(unscaled_dfa_dat$Year)] #SHOULD BE ALL TRUE

#put data into wide format
b.mat2 <- t(as.matrix(unscaled_dfa_dat))
colnames(b.mat2) <- b.mat1[1,]
b.mat2 <- b.mat2[-1,]


#now multiple models at once=======

#let's try 1:5 trends, normal and Student-t distribution, and equal and unequal varianes
m_4000_3_noAR <- find_dfa_trends(
  y = b.mat2, iter = 4000, scale="z-score",
  kmin = 1, kmax = 5, chains = chains, compare_normal = TRUE,
  variance = c("equal", "unequal")#, estimate_trend_ar=TRUE
)

#seems like best is model 11, 1 trends, equal, normal, let's fit that to look at?
names(m)

m$best_model


#check convergence
is_converged(m$best_model, threshold = 1.05)

#rotate trends
r <- rotate_trends(m$best_model)

#plot trend
plot_trends(r) + theme_bw() #if invalid graphics state error try dev.off()

#plot loadings
plot_fitted(m$best_model) + theme_bw()
plot_loadings(r) + theme_bw()


