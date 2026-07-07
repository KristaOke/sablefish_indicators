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
require(yardstick)

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
run_n <- 1

# Model Name
mod <- c("model1", "model2")[1]  #Choose number of DFA trends to include
run <- c("1977to2019", "1996to2019")[run_n]
mod.name <- paste0(mod, "_", "v", version, "_", run)

# Determine number of DFA trends to fit
if(mod=="model1") {
  n_trends <- 1  
}else {
  n_trends <- 2
}

# Do we fit the model, or just load saved .rds outputs?
fit <- FALSE

fit.loo <- TRUE

# MCMC Parameters
n.chains <- 3
n.iter <- 2e4
n.thin <- 2

# Update figure and output directories
dir.figs <- file.path(dir.figs, mod.name)
dir.create(dir.figs, recursive=TRUE)

dir.output <- file.path(dir.output, mod.name)
dir.create(dir.output, recursive=TRUE)

#  1) Read in data =============================================================

# DFA Trends (predictor)
dat.dfa <- read_csv(file.path(dir.data, "DFA_trends_recruit_data_updated.csv"))
head(dat.dfa)

# Sablefish Recruitment Data (response)
dat.rec <- read_csv(file.path(dir.data, "pred_rec_std_aligned.csv"))
head(dat.rec)

# Join Data together
dat.comb <- dat.rec %>% left_join(dat.dfa) %>% dplyr::select(-ln_rec)
head(dat.comb)

#trim data 

dat.comb <- dat.comb[which(dat.comb$Year<2020),] 

#SELECT RUN HERE
run <- "1977to2019"
#run <- "1996to2019" #USE FOR SHORT TIME SERIES

# Determine number of DFA trends to fit
if(run=="1996to2019") {
  dat.comb <- dat.comb[,which(dat.comb$Year>1995)]  
}
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
  # For version 2 with forecast capacity
  trends_fcst <- 0
  trends_se_fcst <- mean(dat.comb$model1_state1_SE)
}else {
  trends <- cbind(dat.comb$model2_state1, dat.comb$model2_state2)
  trends_se <- cbind(dat.comb$model2_state1_SE, dat.comb$model2_state2_SE)
  # For version 2 with forecast capacity
  trends_fcst <- c(0,0)
  trends_se_fcst <- c(mean(dat.comb$model1_state1_SE), mean(dat.comb$model2_state2_SE))
}
# Convert to numeric
trends <- as.matrix(trends)
trends_se <- as.matrix(trends_se)

trends_fcst <- t(as.matrix(trends_fcst))
trends_se_fcst <- t(as.matrix(trends_se_fcst))



# Create Stan Input file
stan.data <- list(
  "n_year"=nrow(dat.comb),
  "rec_ln"=dat.comb$rec_ln,
  "sd_ln"=dat.comb$sd_ln,
  
  "n_trends"=n_trends,
  "trends"=trends,
  "trends_se"=trends_se,
  # v2 for forecast section
  "n_year_fcst"=1,
  "trends_fcst"=trends_fcst,
  "trends_se_fcst"=trends_se_fcst
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
                   control = list(adapt_delta = 0.99))
                   # init=init_ll)
  # Save Model Fit
  saveRDS(stan.fit, file.path(dir.output, paste0("stan.fit.rds")))

}else {
  stan.fit <- readRDS(file.path(dir.output, paste0("stan.fit.rds")))
}

# Extract parameters as a list object
pars <- rstan::extract(stan.fit)


#  5) Plot Results =============================================================


# Plot: Observed vs Predicted ==================================================

# Calculate quantiles of distribution of predicted ln rec
pred.mat <- apply(pars$pred_rec_ln, 2, quantile, probs=c(0.975, 0.75, 0.5, 0.25, 0.025))

years <- dat.comb$Year
n.years <- length(years)

y.lim <- c(min(pred.mat, dat.comb$rec_ln - 1.96*dat.comb$sd_ln),
           max(pred.mat, dat.comb$rec_ln + 1.96*dat.comb$sd_ln))

pdf(file=file.path(dir.figs, "Model Fit.pdf"), height=5, width=7)

# plot(rec_ln~Year, data=dat.comb, type='l', col="blue", ylim=y.lim,
#        ylab="Natural Log of Recruitment", main="Alaska Sablefish")
# points(rec_ln~Year, data=dat.comb, pch=21, bg="blue")

plot(rec_ln~Year, data=dat.comb, type='p', bg="blue", pch=21, ylim=y.lim,
     ylab="Natural Log of Recruitment", main="Alaska Sablefish")

# Add errorbars for 95% CI

segments(x0=years, y0=dat.comb$rec_ln, 
           x1=years, y1=dat.comb$rec_ln - 1.96*dat.comb$sd_ln,
           col='blue')
segments(x0=years, y0=dat.comb$rec_ln, 
           x1=years, y1=dat.comb$rec_ln + 1.96*dat.comb$sd_ln,
           col='blue')


# Plot Model predicted ln rec
polygon(x=c(years, rev(years)), y=c(pred.mat['97.5%',], rev(pred.mat['2.5%',])),
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
polygon(x=c(years, rev(years)), y=c(pred.mat['75%',], rev(pred.mat['25%',])),
        col=rgb(1,0,0, alpha=0.25), border=FALSE)
lines(x=years, y=pred.mat['50%',], col=rgb(1,0,0, alpha=0.5))

dev.off()

# Plot: Predicted DFA Trends vs. Observed ======================================

n_trends

pdf(file=file.path(dir.figs, "Fit to Trends.pdf"), height=ifelse(n_trends==1, 5, 7), width=6)

par(mfrow=c(n_trends,1))
# plot(rec_ln~Year, data=dat.comb, type='l', col="blue", ylim=y.lim,
#        ylab="Natural Log of Recruitment", main="Alaska Sablefish")
# points(rec_ln~Year, data=dat.comb, pch=21, bg="blue")
t <- 1
for(t in 1:n_trends) {
  # Calculate quantiles of distribution of predicted ln rec
  trends.mat <- apply(pars$pred_trends, 2, quantile, probs=c(0.975, 0.75, 0.5, 0.25, 0.025))
  
  # Determine y limit
  y.lim <- c(min(trends.mat, trends[,t] - 2*trends_se[,t]),
             max(trends.mat, trends[,t] + 2*trends_se[,t]))
  
  plot(x=years, y=trends[,t], type='p', bg="blue", pch=21, ylim=y.lim,
      ylab=paste("DFA Trend", t), main="Alaska Sablefish")

  # Add errorbars for 95% CI
  # y <- 1
  # for(y in 1:n.years) {
    segments(x0=years, y0=trends[,t], 
             x1=years, y1=trends[,t] - 2*trends_se[,t],
             col='blue')
    segments(x0=years, y0=trends[,t], 
             x1=years, y1=trends[,t] +2*trends_se[,t],
             col='blue')
  # }

  # Plot Model predicted ln rec
  polygon(x=c(years, rev(years)), y=c(trends.mat['97.5%',], rev(trends.mat['2.5%',])),
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
  polygon(x=c(years, rev(years)), y=c(trends.mat['75%',], rev(trends.mat['25%',])),
          col=rgb(1,0,0, alpha=0.25), border=FALSE)
  lines(x=years, y=trends.mat['50%',], col=rgb(1,0,0, alpha=0.5))
} # next t
dev.off()

# Plot: Intercept and Slopes
pdf(file=file.path(dir.figs, "Parameters.pdf"), height=5, width=ifelse(n_trends==1, 5, 8))

par(mfrow=c(1, n_trends+1), mar=c(5,2,2,1))

plotPost(pars$incpt, xlab="Intercept", )

t <- 1
for(t in 1:n_trends) {
  plotPost(pars$slp[,t], xlab=paste("Effect of Trend", t))
}

dev.off()

# Plot: Observed vs Predicted on X-Y Axes ======================================

# Calculate Performance metrics
??yardstick::mape

# Mean Absolute Percent Error
temp.mape <- yardstick::mape_vec(truth=dat.comb$rec_ln, estimate=pred.mat['50%',])
temp.rmse <- yardstick::rmse_vec(truth=dat.comb$rec_ln, estimate=pred.mat['50%',])

pdf(file=file.path(dir.figs, "Predicted-Observed.pdf"), height=6, width=6)
plot(x=dat.comb$rec_ln, y=pred.mat['50%',], pch=21, bg=rgb(0,0,1, alpha=0.25),
       xlab="Observed ln(recruitment)",
       ylab="Predicted ln(recruitment)",
       main=paste("Alaska Sablefish:", n_trends, "Trends"))

abline(a=0, b=1, col=rgb(1,0,0, alpha=0.5), lwd=3)

legend("topleft", legend=c(paste("Mean Abs. % Error: ", round(temp.mape, 1), "%"),
                           paste("Root Mean Squared Error: ", round(temp.rmse, 2))),
       bty='n')

dev.off()




#  6) Attempting LOOCV ===========================================================
#Fit the model
if(fit.loo==TRUE) {

  yrs <- unique(dat.comb$Year)
  
  loo.years <- yrs
  loo.obs <- vector(length=length(yrs))
  loo.pred <- vector(length=length(yrs))
  loo.postPred <- vector(length=length(yrs))
  
  # List to hold extracted parameters
  loo.pars <- vector("list", length=length(yrs))
  names(loo.pars) <- yrs
  
  i <- 1
  for(i in 1:length(dat.comb$Year)) { #update here tuesday
  
    print(paste("i:",i,"of",length(dat.comb$Year)))
    yr <- yrs[i]
    
    # Fitting data with missing year
    dat.comb.fit <- dat.comb[-i,]
    
    # Forecast data for witheld year
    dat.comb.fcst <- dat.comb[i,]
    
    # Save Truth
    loo.obs[i] <- dat.comb.fcst$rec_ln
    
    # Create trends objects
    
    if(n_trends==1) {
      trends <- dat.comb.fit$model1_state1
      trends_se <- dat.comb.fit$model1_state1_SE
      # For version 2 with forecast capacity
      trends_fcst <- dat.comb.fcst$model1_state1
      trends_se_fcst <- dat.comb.fcst$model1_state1_SE
    }else {
      trends <- cbind(dat.comb.fit$model2_state1, dat.comb.fit$model2_state2)
      trends_se <- cbind(dat.comb.fit$model2_state1_SE, dat.comb.fit$model2_state2_SE)
      # For version 2 with forecast capacity
      trends_fcst <- cbind(dat.comb.fcst$model2_state1, dat.comb.fcst$model2_state2)
      trends_se_fcst <- cbind(dat.comb.fcst$model2_state1_SE, dat.comb.fcst$model2_state2_SE)
    }
    # Convert to numeric
    trends <- as.matrix(trends)
    trends_se <- as.matrix(trends_se)
    
    trends_fcst <- as.matrix(trends_fcst)
    trends_se_fcst <- as.matrix(trends_se_fcst)
    
    # Create Stan Input file
    stan.data <- NULL
    stan.data <- list(
      "n_year"=nrow(dat.comb.fit),
      "rec_ln"=dat.comb.fit$rec_ln,
      "sd_ln"=dat.comb.fit$sd_ln,
      
      "n_trends"=n_trends,
      "trends"=trends,
      "trends_se"=trends_se,
      # v2 for forecast section
      "n_year_fcst"=1,
      "trends_fcst"=trends_fcst,
      "trends_se_fcst"=trends_se_fcst
    )
    
    stan.fit <- NULL
    stan.fit <- stan(file=file.path(dir.stan, paste0("sablefish-ss-reg-v", version, ".stan")),
                     model_name=paste0("sablefish-ss-reg-v", version),
                     data=stan.data,
                     chains=n.chains, iter=n.iter/2, thin=n.thin, # Half as many samples for testing
                     # chains=3, iter=5e3, thin=5,
                     cores=n.chains, verbose=FALSE,
                     seed=101, refresh = 0,
                     control = list(adapt_delta = 0.99))
    # init=init_ll)
    # Save Model Fit
    # saveRDS(stan.fit, file.path(dir.output, paste0("stan.fit.drop", yr, ".rds")))
    
    # Extract parameters as a list object
    pars <- rstan::extract(stan.fit)
    loo.pars[[i]] <- pars
    
    # Save prediction for ln_rec in witheld year
    loo.pred[i] <- median(pars$pred_rec_ln_fcst)
    
    # Posterior predicted values     
    loo.postPred[i] <- median(pars$pred_rec_ln_fcst_postPred)
    
  } # next i (year)
  
  # Save LOO results
  loo.results <- data.frame(loo.years, loo.obs, loo.pred, loo.postPred) %>% 
                   rename("year"="loo.years",
                          "observed"="loo.obs",
                          "predicted"="loo.pred",
                          "posteriorPred"="loo.postPred")
  
  write.csv(loo.results, file.path(dir.figs, "loo.results.csv"))
  
  # Save parameters
  saveRDS(loo.pars, file.path(dir.output, "loo.pars.rds"))

}else {
  loo.results <- read.csv(file.path(dir.figs, "loo.results.csv"))
  loo.pars <- readRDS(file.path(dir.output, "loo.pars.rds"))
}

# Plotting LOO results

yardstick::rmse_vec(truth=loo.results$observed, estimate=loo.results$predicted)
yardstick::rmse_vec(truth=loo.results$observed, estimate=loo.results$posteriorPred)



