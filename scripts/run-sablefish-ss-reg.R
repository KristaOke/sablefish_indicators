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

# 

#  3) Function to Generate Initial Values ======================================
# Initialization Function
init_fn <- function(chain_id=1) {
  list( 
    "scale" = runif(1, 1e3, 3e4),
    "scale_hist" = runif(nYearPM, 1e4, 3e4),
    "muDay" = rnorm(1,prior.muDay.mean,1),
    "sigmaDay" = rnorm(1,prior.sigmaDay.mean,1),
    "muDay_hist" = rnorm(nYearPM,prior.muDay.mean,1),
    "sigmaDay_hist" = rnorm(nYearPM,prior.sigmaDay.mean,1),
    "sigmaOE"=runif(1,0.5,1),
    "sigmaOE_hist"=runif(nYearPM,0.5,1),
    "prop_logis_pm"=runif(1,quantile(ce.props$prop.bay$prop, probs=0.25), 
                          quantile(ce.props$prop.bay$prop, probs=0.25)),
    "int_saa"=runif(n=1, min=min(Robs_saa), max=max(Robs_saa))
  )
}
# init_fn()
# Initial List of Lists for Multiple Chains
init_ll <- lapply(1:n.chains, function(id) init_fn(chain_id = id))

#  4) Fit Stan model ===========================================================
#Fit the model
stan.fit <- NULL
if(fit==TRUE) {
  stan.fit <- stan(file=file.path(dir.stan, paste0("sablefish-ss-reg-v", version, ".stan")),
                   model_name=paste0("sablefish-ss-reg-v", version),
                   data=list(
                     "n_years"=n_years,
                     "rec_ln"=,
                     ""=
                   ),
                   chains=n.chains, iter=n.iter, thin=n.thin,
                   # chains=3, iter=5e3, thin=5,
                   cores=n.chains, verbose=FALSE,
                   seed=101,
                   control = list(adapt_delta = 0.99),
                   init=init_ll)
  # Save Model Fit
  saveRDS(stan.fit, file.path(dir.output, paste0(,".rds")))

}else {
  stan.fit <- readRDS(file.path(dir.output, paste0(, ".rds")))
}


pars <- rstan::extract(stan.fit)


#  5) Plot Results =============================================================