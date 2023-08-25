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

# Determine number of DFA trends to fit
if(mod=="model1") {
  n_trends <- 1  
}else {
  n_trends <- 2
}

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
dat.dfa <- read_csv(file.path(dir.data, "DFA_trends_recruit_data.csv"))
head(dat.dfa)

# Sablefish Recruitment Data (response)
dat.rec <- read_csv(file.path(dir.data, "pred_rec_std_aligned.csv"))
head(dat.rec)

# Join Data together
dat.comb <- dat.rec %>% left_join(dat.dfa) %>% dplyr::select(-ln_rec)
head(dat.comb)

#  2) Craft Stan Input Objects =================================================

# Calculate CV (normal space)
dat.comb <- dat.comb %>% mutate("cv_norm"=std.dev/pred_rec,
                                "sd_ln"=sqrt(log(cv_norm^2 + 1)),
                                "rec_ln"=log(pred_rec)) %>% 
                         dplyr::filter(!is.na(rec_ln))

head(dat.comb)
tail(dat.comb)

# Checkup
# cv_norm_test <- sqrt(exp(dat.comb$sd_ln^2)-1) - Passed

# Create trends objects

if(n_trends==1) {
  trends <- dat.comb$model1_state1
  trends_se <- dat.comb$model1_state1_SE
}else {
  trends <- cbind(dat.comb$model2_state1, dat.comb$model2_state2)
  trends_se <- cbind(dat.comb$model2_state1_SE, dat.comb$model2_state2_SE)
}
# Convert to numeriuc
trends <- as.matrix(trends)
trends_se <- as.matrix(trends_se)



# Create Stan Input file
stan.data <- list(
  "n_year"=nrow(dat.comb),
  "rec_ln"=dat.comb$rec_ln,
  "sd_ln"=dat.comb$sd_ln,
  
  "n_trends"=n_trends,
  "trends"=trends,
  "trends_se"=trends_se
)

#  3) Function to Generate Initial Values ======================================
# Initialization Function
init_fn <- function(chain_id=1) {
  list( 
    "incpt"=rnorm(n=1, mean=mean(dat.comb$rec_ln), sd=sd(dat.comb$rec_ln)),
    "slp"=rnorm(n_trends, 0, 1)
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
                   data=stan.data,
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
