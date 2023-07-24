#=============================================================================== 
#Project Name: SABLEFISH ESP - State-Space Regression - v1
#Creator: Curry James Cunningham, College of Fisheries and Ocean Sciences, UAF
#Date: 7.23.2023
#
#Purpose: To fit a regression 
#
#  1) Read in data
#  2) Craft Stan Input Objects
#  3) Function to Generate Initial Values
#  4) Fit Stan model
#  5) Plot Results
#
#
#===============================================================================
#NOTES:
#  1) Set fit variable to TRUE/FALSE to fit model or read saved .rds model objects
#ideally we would run separate state space models regressing the trends from
#each model against ln_recruitment
#e.g.
# ln_rec ~ model1_state1 
# AND
# ln_rec ~ model2_state1 + model2_state2
#===============================================================================
require(BEST)
require(rstan)
require(tidyverse)
require(ggthemes)
require(viridis)
require(shinystan)
require(beepr)
require(reshape2)
require(dplyr)
require(tidybayes)
require(loo)
require(here)


# CURRENTLY NEEDED WITH NEW R FOR PARALLELIZATION
# if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
#     Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
#   parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
# }

options(mc.cores = parallel::detectCores())

rstan_options(javascript=FALSE)

# Define Workflow Paths ========================================================
wd <- here()
setwd(wd)

dir.output <- file.path(wd, "output")
dir.figs <- file.path(wd, "figs")
dir.stan <- file.path(wd, "stan")
dir.data <- file.path(wd,"data")
# dir.R <- file.path(wd, "R")

# CONTROL SECTION ==============================================================
version <- 1 

# Model Name
mod <- c("model1", "model2")[1]
mod.name <- paste0(mod, "_", "v", version)

# Do we fit the model, or just load saved .rds outputs
fit <- TRUE 

# Use Brood table data as runsize target (TRUE), or use total CE as runsize target (FALSE)
run.target.bt <- TRUE 

# MCMC Parameters
n.chains <- 3
n.iter <- 2e3
n.thin <- 2

# Update figure and output directories
dir.figs <- file.path(dir.figs, mod.name)
dir.create(dir.figs, recursive=TRUE)

dir.output <- file.path(dir.output, mod.name)
dir.create(dir.output, recursive=TRUE)

#  1) Read in data =============================================================

# DFA Trends (predictor)
dat.dfa <- read_csv(file.path(dir.data, "pred_rec_std_aligned.csv"))
head(dat.dfa)

# Sablefish Recruitment Data (response)
dat.rec <- read_csv(file.path(dir.data, "DFA_trends_recruit_data.csv"))
head(dat.rec)

# Join Data together
dat.comb <- dat.rec %>% left_join(dat.dfa) %>% dplyr::select(-ln_rec)
head(dat.comb)

#  2) Craft Stan Input Objects =================================================

# Calculate CV (normal space)
dat.comb <- dat.comb %>% mutate("cv_norm"=std.dev/pred_rec,
                                "sd_ln"=sqrt(log(cv_norm^2 + 1)),
                                "rec_ln"=log(pred_rec))

head(dat.comb)

# Checkup
# cv_norm_test <- sqrt(exp(dat.comb$sd_ln^2)-1) - Passed

#  3) Function to Generate Initial Values ======================================


#  4) Fit Stan model ===========================================================



#  5) Plot Results =============================================================
