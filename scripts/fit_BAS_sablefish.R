#==================================================================================================
#Fit BAS to training and testing datasets
#
#
#Krista, March 2023
#BORROWING HEAVILY from scripts developed by Curry Cunningham for ESP analyses
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
#TIMING:
#
#
##==================================================================================================
library(tidyverse)
library(corrplot)
library(cowplot)
library(ggplot2)
library(viridis)
library(ggthemes)
library(BAS)
library(readxl)
library(yardstick)


#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#=============================================================
#### Control Section ####

fit <- TRUE

offset <- 0

#load data===================================================

#Get data=======

# load data that is already z-scored, checked for correlations


scaled <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)


#DATA CONTROL SECTION----

#select covariates

#using ONLY the indicator subset that's been checked for collinearity
#using noncor_covars from process_data for now

#BAS CANNOT handle missing data


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


scaled_bas_dat <- scaled[,names(scaled) %in% noncor_covars3]

#leave out 2020 b/c it's a longterm mean
scaled_bas_dat <- scaled_bas_dat[which(scaled_bas_dat$Year<2020),]




#Determine Covariates

covars <- names(scaled_bas_dat)[-which(names(scaled_bas_dat) %in% c("Year", "ln_rec"))]

n.cov <- length(covars)

dat.temp <- scaled_bas_dat[-which(names(scaled_bas_dat) %in% c("Year"))]

#more parameters than data, let's try dropping the shortest time series
#Smr_temp_250m_GOA_scaled
dat.temp <- dat.temp[-which(names(dat.temp) %in% c("Smr_temp_250m_GOA_scaled"))]

# Bayesian Model Selection
bas.lm <-  bas.lm(ln_rec ~ ., data=dat.temp,
                  # prior="ZS-null",
                  modelprior=uniform(), initprobs="Uniform",
                  method='BAS', MCMC.iterations=1e5, thin=10)

summary(bas.lm)
image(bas.lm, rotate=F)

plot(bas.lm, which = 4, ask=FALSE, caption="", sub.caption="")
plot(coef(bas.lm),  ask=FALSE)
plot(bas.lm, which=4)
summary(bas.lm)


# Plot Model Predictions vs. Observed ==============================
#pdf(file.path(dir.figs,"Model Fit.pdf"), height=6, width=9)
par(oma=c(1,1,1,1), mar=c(4,4,1,1), mfrow=c(1,2))
pred.bas <- predict(bas.lm, estimator="BMA", se.fit=T)

# Omit NAs
dat.temp.na.omit <- na.omit(scaled_bas_dat[-which(names(scaled_bas_dat) %in% c("Smr_temp_250m_GOA_scaled"))])

plot(x=dat.temp.na.omit$ln_rec, y=pred.bas$Ybma,
     xlab="Observed ln(Recruitment)", ylab="Predicted ln(Recruitment)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=pred.bas$fit, y=pred.bas$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)

# Timeseries
plot(x=dat.temp.na.omit$Year, y=dat.temp.na.omit$ln_rec,
     xlab="Year", ylab="ln(Recruitment)", type='l', col=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
grid(lty=3, col='dark gray')
points(x=dat.temp.na.omit$Year, y=dat.temp.na.omit$ln_rec,
       pch=21, bg=rgb(1,0,0,alpha=0.5))
lines(x=dat.temp.na.omit$Year, y=pred.bas$Ybma, lwd=3, col=rgb(0,0,1, alpha=0.5))
points(x=dat.temp.na.omit$Year, y=pred.bas$Ybma,
       pch=21, bg=rgb(0,1,0,alpha=0.5))
#conf_int <- confint(pred.bas, parm = "pred")
#plotCI(x=dat.temp.na.omit$Year, y=pred.bas$Ybma,li=conf_int[,1], ui=conf_int[,2])

legend('topleft', legend=c("Observed","Predicted"), lty=1, col=c(rgb(1,0,0,alpha=0.5),
                                                             rgb(0,0,1, alpha=0.5)), bg="white")

dev.off()



# bas.lm.2 <-  bas.lm(ln_rec ~ ., data=dat.temp,
#                     # prior="ZS-null",
#                     modelprior=uniform(), initprobs="Uniform",
#                     method='MCMC', MCMC.iterations=1e6, thin=10)


# PLOT RESULTS ==================================================
names(summary(bas.lm))

inc.probs <- summary(bas.lm)[2:ncol(dat.temp),1]
# par(oma=c(1,1,1,1), mar=c(4,20,1,1))
# barplot(inc.probs, horiz=TRUE, xlim=c(0,1), las=2)
# abline(v=seq(from=0.2, to=0.8, by=0.2), lty=2)
# box()

bas.names <- coef(bas.lm)$namesx
inc.probs <- coef(bas.lm)$probne0
post.mean <- coef(bas.lm)$postmean
post.sd <- coef(bas.lm)$postsd
#Calcualte lower and upper 95% CI
low.95 <- post.mean - 1.96*post.sd
up.95 <- post.mean + 1.96*post.sd

# confint(coef(bas.lm), level=c(0.5))

# post.probs <- coef(bas.lm)$postprobs

cond.mean <- coef(bas.lm)$conditionalmean[,2]
cond.sd <- coef(bas.lm)$conditionalsd

names(coef(bas.lm))



#Plot it out....
par(mfrow=c(1,2), mar=c(4,1,2,1), oma=c(0,10,1,1))

plot.df <- data.frame(bas.names, inc.probs, post.mean, post.sd, low.95, up.95)
# plot.list <- melt(plot.df)

#let's make the labels pretty

plot.df$bas.names[which(plot.df$bas.names=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
plot.df$bas.names[which(plot.df$bas.names=="Smr_temp_250m_GOA_scaled")] <- "Summer 250m temperature GOA"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_biom_GOA_scaled")] <- "Spring chlorophyll A biomass GOA"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_biom_SEBS_scaled")] <- "Spring chlorophyll A biomass SEBS"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_peak_GOA_scaled")] <- "Spring chlorophyll A peak GOA"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_peak_SEBS_scaled")] <- "Spring chlorophyll A peak SEBS"
plot.df$bas.names[which(plot.df$bas.names=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
plot.df$bas.names[which(plot.df$bas.names=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
plot.df$bas.names[which(plot.df$bas.names=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
plot.df$bas.names[which(plot.df$bas.names=="smr_adult_cond_scaled")] <- "Summer adult condition"




# g <- ggplot(filter(plot.df, bas.names!='Intercept'),
#             aes(x=bas.names, post.mean, fill=bas.names)) +
#        theme_linedraw() +
#        geom_errorbar(aes(ymin=low.95, ymax=up.95), width=0.25) +
#        geom_point(pch=21) +
#        geom_hline(yintercept = 0, col='red', alpha=0.5) +
#        ylab('Effect') +
#        xlab('Covariate') +
#        coord_flip() +
#        theme(legend.position='none')
# 
# g

g <- ggplot(filter(plot.df, bas.names!='Intercept'),
            aes(x=bas.names, post.mean, fill=bas.names)) +
  theme_linedraw() +
  geom_errorbar(aes(ymin=post.mean-post.sd, ymax=post.mean+post.sd), width=0.25) +
  geom_point(pch=21) +
  geom_hline(yintercept = 0, col='red', alpha=0.5) +
  ylab('Effect') +
  xlab('Covariate') +
  coord_flip() +
  theme(legend.position='none')

g

#Inclusion prob

g2 <-  ggplot(filter(plot.df, bas.names!='Intercept'),
              aes(x=bas.names, y=inc.probs, fill=bas.names)) +
  theme_linedraw() +
  geom_bar(stat='identity', color='black') +
  ylab('Inclusion\nProbability') +
  # coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=c(0,1)) +
  theme(legend.position='none', axis.text.y = element_blank(), 
        axis.title.y=element_blank()) +
  coord_flip()
# scale_fill_continuous()
g2

# Bring Figs Together ========
g3 <- plot_grid(g,g2, nrow=1, ncol=2, rel_widths=c(3,1), align='h')
g3 #+ ggtitle('Sablefish Recruitment', subtitle=paste('Rsq:',round(,2)))
# ggsave(file=file.path(dir.figs,"BAS.png"), plot=g3, height=5, width=8, units='in',
#        dpi=500)


#PLOT OUTPUT WITHOUT RAINBOW ===========
g.b <- ggplot(filter(plot.df, bas.names!='Intercept'),
              aes(x=bas.names, post.mean, fill='blue')) +
  theme_linedraw() +
  geom_errorbar(aes(ymin=post.mean-post.sd, ymax=post.mean+post.sd), width=0.25) +
  geom_point(pch=21, fill='blue') +
  geom_hline(yintercept = 0, col='red', alpha=0.5) +
  ylab('Effect') +
  xlab('Covariate') +
  coord_flip() +
  theme(legend.position='none')

g.b

#Inclusion prob

g2.b <-  ggplot(filter(plot.df, bas.names!='Intercept'),
                aes(x=bas.names, y=inc.probs, fill=inc.probs)) +
  theme_linedraw() +
  geom_bar(stat='identity', color='black') +
  ylab('Inclusion\nProbability') +
  # coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=c(0,1)) +
  theme(legend.position='none', axis.text.y = element_blank(), 
        axis.title.y=element_blank()) +
  coord_flip() +
  scale_fill_continuous_tableau()



g2.b

# Bring Figs Together ========
g3.b <- plot_grid(g.b,g2.b, nrow=1, ncol=2, rel_widths=c(3,1), align='h')
g3.b #+ ggtitle('Sablefish Recruitment', subtitle=paste('Rsq:',round(,2)))
# ggsave(file=file.path(dir.figs,"BAS_noRainbow.png"), plot=g3.b, height=5, width=8, units='in',
#        dpi=500)

#REPEAT without short time series================================================================

#use only middleton, condition, sebs sst, adfg

#Determine Covariates

covars <- names(scaled_bas_dat)[which(names(scaled_bas_dat) %in% c("Year", "ln_rec",
                                                                   "Spr_ST_SEBS_scaled",
                                                                   "smr_adult_cond_scaled",
                                                                   "YOY_grwth_Middleton_scaled",
                                                                   "Smr_CPUE_juv_ADFG_ln_scaled"))]

n.cov <- length(covars)

dat.temp <- scaled_bas_dat[-which(names(scaled_bas_dat) %in% c("Year"))]
dat.temp <- dat.temp[which(names(dat.temp) %in% covars)]

# Bayesian Model Selection
bas.lm <-  bas.lm(ln_rec ~ ., data=dat.temp,
                  # prior="ZS-null",
                  modelprior=uniform(), initprobs="Uniform",
                  method='BAS', MCMC.iterations=1e5, thin=10)

summary(bas.lm)

plot(bas.lm, which = 4, ask=FALSE, caption="", sub.caption="")
plot(coef(bas.lm),  ask=FALSE)
plot(bas.lm, which=4)


# Plot Model Predictions vs. Observed ==============================
#pdf(file.path(dir.figs,"Model Fit.pdf"), height=6, width=9)
par(oma=c(1,1,1,1), mar=c(4,4,1,1), mfrow=c(1,2))
pred.bas <- predict(bas.lm, estimator="BMA")

# Omit NAs
#dat.temp.na.omit <- na.omit(scaled_bas_dat[-which(names(scaled_bas_dat) %in% c("Smr_temp_250m_GOA_scaled", "Smr_condition_fem_age4_GOA_scaled"))])
dat.temp.na.omit <- na.omit(scaled_bas_dat[which(names(scaled_bas_dat) %in% c("Year", "ln_rec", covars))])


plot(x=dat.temp.na.omit$ln_rec, y=pred.bas$Ybma,
     xlab="Observed ln(Recruitment)", ylab="Predicted ln(Recruitment)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=pred.bas$fit, y=pred.bas$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)

# Timeseries
plot(x=dat.temp.na.omit$Year, y=dat.temp.na.omit$ln_rec,
     xlab="Year", ylab="ln(Recruitment)", type='l', col=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
grid(lty=3, col='dark gray')
points(x=dat.temp.na.omit$Year, y=dat.temp.na.omit$ln_rec,
       pch=21, bg=rgb(1,0,0,alpha=0.5))
lines(x=dat.temp.na.omit$Year, y=pred.bas$Ybma, lwd=3, col=rgb(0,0,1, alpha=0.5))

legend('top', legend=c("Observed","Predicted"), lty=1, col=c(rgb(1,0,0,alpha=0.5),
                                                             rgb(0,0,1, alpha=0.5)), bg="white")

dev.off()



# bas.lm.2 <-  bas.lm(ln_rec ~ ., data=dat.temp,
#                     # prior="ZS-null",
#                     modelprior=uniform(), initprobs="Uniform",
#                     method='MCMC', MCMC.iterations=1e6, thin=10)


# PLOT RESULTS ==================================================
names(summary(bas.lm))

inc.probs <- summary(bas.lm)[2:ncol(dat.temp),1]
# par(oma=c(1,1,1,1), mar=c(4,20,1,1))
# barplot(inc.probs, horiz=TRUE, xlim=c(0,1), las=2)
# abline(v=seq(from=0.2, to=0.8, by=0.2), lty=2)
# box()

bas.names <- coef(bas.lm)$namesx
inc.probs <- coef(bas.lm)$probne0
post.mean <- coef(bas.lm)$postmean
post.sd <- coef(bas.lm)$postsd
#Calcualte lower and upper 95% CI
low.95 <- post.mean - 1.96*post.sd
up.95 <- post.mean + 1.96*post.sd

# confint(coef(bas.lm), level=c(0.5))

# post.probs <- coef(bas.lm)$postprobs

cond.mean <- coef(bas.lm)$conditionalmean[,2]
cond.sd <- coef(bas.lm)$conditionalsd

names(coef(bas.lm))


#Plot it out....
par(mfrow=c(1,2), mar=c(4,1,2,1), oma=c(0,10,1,1))

plot.df <- data.frame(bas.names, inc.probs, post.mean, post.sd, low.95, up.95)
# plot.list <- melt(plot.df)

plot.df$bas.names[which(plot.df$bas.names=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
plot.df$bas.names[which(plot.df$bas.names=="Smr_temp_250m_GOA_scaled")] <- "Summer 250m temperature GOA"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_biom_GOA_scaled")] <- "Spring chlorophyll A biomass GOA"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_biom_SEBS_scaled")] <- "Spring chlorophyll A biomass SEBS"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_peak_GOA_scaled")] <- "Spring chlorophyll A peak GOA"
plot.df$bas.names[which(plot.df$bas.names=="Spr_chlA_peak_SEBS_scaled")] <- "Spring chlorophyll A peak SEBS"
plot.df$bas.names[which(plot.df$bas.names=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
plot.df$bas.names[which(plot.df$bas.names=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
plot.df$bas.names[which(plot.df$bas.names=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
plot.df$bas.names[which(plot.df$bas.names=="smr_adult_cond_scaled")] <- "Summer adult condition"




# g <- ggplot(filter(plot.df, bas.names!='Intercept'),
#             aes(x=bas.names, post.mean, fill=bas.names)) +
#        theme_linedraw() +
#        geom_errorbar(aes(ymin=low.95, ymax=up.95), width=0.25) +
#        geom_point(pch=21) +
#        geom_hline(yintercept = 0, col='red', alpha=0.5) +
#        ylab('Effect') +
#        xlab('Covariate') +
#        coord_flip() +
#        theme(legend.position='none')
# 
# g

g <- ggplot(filter(plot.df, bas.names!='Intercept'),
            aes(x=bas.names, post.mean, fill=bas.names)) +
  theme_linedraw() +
  geom_errorbar(aes(ymin=post.mean-post.sd, ymax=post.mean+post.sd), width=0.25) +
  geom_point(pch=21) +
  geom_hline(yintercept = 0, col='red', alpha=0.5) +
  ylab('Effect') +
  xlab('Covariate') +
  coord_flip() +
  theme(legend.position='none')

g

#Inclusion prob

g2 <-  ggplot(filter(plot.df, bas.names!='Intercept'),
              aes(x=bas.names, y=inc.probs, fill=bas.names)) +
  theme_linedraw() +
  geom_bar(stat='identity', color='black') +
  ylab('Inclusion\nProbability') +
  # coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=c(0,1)) +
  theme(legend.position='none', axis.text.y = element_blank(), 
        axis.title.y=element_blank()) +
  coord_flip()
# scale_fill_continuous()
g2

# Bring Figs Together ========
g3 <- plot_grid(g,g2, nrow=1, ncol=2, rel_widths=c(3,1), align='h')
g3 #+ ggtitle('Sablefish Recruitment', subtitle=paste('Rsq:',round(,2)))


#PLOT OUTPUT WITHOUT RAINBOW ===========
g.b <- ggplot(filter(plot.df, bas.names!='Intercept'),
              aes(x=bas.names, post.mean, fill='blue')) +
  theme_linedraw() +
  geom_errorbar(aes(ymin=post.mean-post.sd, ymax=post.mean+post.sd), width=0.25) +
  geom_point(pch=21, fill='blue') +
  geom_hline(yintercept = 0, col='red', alpha=0.5) +
  ylab('Effect') +
  xlab('Covariate') +
  coord_flip() +
  theme(legend.position='none')

g.b

#Inclusion prob

g2.b <-  ggplot(filter(plot.df, bas.names!='Intercept'),
                aes(x=bas.names, y=inc.probs, fill=inc.probs)) +
  theme_linedraw() +
  geom_bar(stat='identity', color='black') +
  ylab('Inclusion\nProbability') +
  # coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=c(0,1)) +
  theme(legend.position='none', axis.text.y = element_blank(), 
        axis.title.y=element_blank()) +
  coord_flip() +
  scale_fill_continuous_tableau()



g2.b

# Bring Figs Together ========
g3.b <- plot_grid(g.b,g2.b, nrow=1, ncol=2, rel_widths=c(3,1), align='h')
g3.b #+ ggtitle('Sablefish Recruitment', subtitle=paste('Rsq:',round(,2)))



#LOOCV============================================================================

#LOOCV on only long time series-------



covars <- names(scaled_bas_dat)[which(names(scaled_bas_dat) %in% c("Year", "ln_rec",
                                                                   "Spr_ST_SEBS_scaled",
                                                                   "smr_adult_cond_scaled",
                                                                   "YOY_grwth_Middleton_scaled",
                                                                   "Smr_CPUE_juv_ADFG_ln_scaled"))]

n.cov <- length(covars)

#STEP 1 - Loop through training sets and fit models-------

#I am using HPM or highest probability model rather than model averaging which produces a 
#range of predictions, should revisit

scaled_loop_dat <- scaled_bas_dat

yrs <- unique(scaled_loop_dat$Year)
output_df <- data.frame(matrix(ncol=3, nrow = length(yrs)))
colnames(output_df) <- c("Year", "observed_ln_recruit", "predicted_ln_recruit")

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  
  temp_dat <- temp_dat[-which(names(temp_dat) %in% c("Year"))]
  temp_dat <- temp_dat[which(names(temp_dat) %in% covars)]
  
  dropped_yr <- scaled_loop_dat[i,]
  output_df$observed_ln_recruit[i] <- dropped_yr$ln_rec
  dropped_yr <- dropped_yr[,names(dropped_yr) %in% covars]
  dropped_yr <- dropped_yr[,!names(dropped_yr) %in% "ln_rec"]
  print(dropped_yr$Year)
  #fit model
  bas.loop <-  bas.lm(ln_rec ~ ., data=temp_dat,
                    # prior="ZS-null",
                    modelprior=uniform(), initprobs="Uniform",
                    method='BAS', MCMC.iterations=1e5, thin=10)
  
  #have model predict to missing year
  temp_predict <- predict(bas.loop, newdata=dropped_yr, estimator="BMA")
  print(temp_predict$bestmodel)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict$fit
}

#PROBLEM! Using HPM results in some loops selecting different best models!
#going back to using BMA

output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
 # geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year))+
  ylim(c(0,5)) + xlim(c(0,5)) + theme_bw()

#get MSE & MAE------

#these need to be double checked!
#BAS_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)


obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
   geom_point() + geom_smooth(method="lm")

BAS_long_rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

BAS_long_mae <- mae(output_df, truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)

write.csv(output_df, file=paste(wd,"/data/BAS_obsvpreds_long.csv", sep=""))
output_df_long <- read.csv(file=paste(wd,"/data/BAS_obsvpreds_long.csv", sep=""))


#LOOCV on short time series-------



covars <- c("Year", "ln_rec",
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

n.cov <- length(covars)

#STEP 1 - Loop through training sets and fit models-------

#I am using HPM or highest probability model rather than model averaging which produces a 
#range of predictions, should revisit

scaled_loop_dat <- scaled_bas_dat

yrs <- unique(scaled_loop_dat$Year)
output_df <- data.frame(matrix(ncol=3, nrow = length(yrs)))
colnames(output_df) <- c("Year", "observed_ln_recruit", "predicted_ln_recruit")

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  
  temp_dat <- temp_dat[-which(names(temp_dat) %in% c("Year"))]
  temp_dat <- temp_dat[which(names(temp_dat) %in% covars)]
  
  dropped_yr <- scaled_loop_dat[i,]
  output_df$observed_ln_recruit[i] <- dropped_yr$ln_rec
  dropped_yr <- dropped_yr[,names(dropped_yr) %in% covars]
  dropped_yr <- dropped_yr[,!names(dropped_yr) %in% "ln_rec"]
  print(dropped_yr$Year)
  #fit model
  bas.loop <-  bas.lm(ln_rec ~ ., data=temp_dat,
                      # prior="ZS-null",
                      modelprior=uniform(), initprobs="Uniform",
                      method='BAS', MCMC.iterations=1e5, thin=10)
  
  #have model predict to missing year
  temp_predict <- predict(bas.loop, newdata=dropped_yr, estimator="BMA")
  print(temp_predict$bestmodel)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict$fit
}

#PROBLEM! Using HPM results in some loops selecting different best models!
#going back to using BMA


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  #geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year))+
  ylim(c(0,5)) + xlim(c(0,5)) + theme_bw()

#get MSE & MAE------
library(yardstick)

#these need to be double checked!
#BAS_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)

obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

BAS_rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

BAS_mae <- mae(output_df, truth=observed_ln_recruit, 
                estimate=predicted_ln_recruit, na.rm=TRUE)


output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
  geom_point() + geom_smooth(method="lm")

write.csv(output_df, file=paste(wd,"/data/BAS_obsvpreds_reduced.csv", sep=""))
output_df_reduced <- read.csv(file=paste(wd,"/data/BAS_obsvpreds_reduced.csv", sep=""))

#repeat same plot as within-----
par(oma=c(1,1,1,1), mar=c(4,4,1,1), mfrow=c(1,2))

# Omit NAs
dat.temp <- output_df_long

plot(x=dat.temp$observed_ln_recruit, y=dat.temp$predicted_ln_recruit,
     xlab="Observed ln(Recruitment)", ylab="Predicted ln(Recruitment)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=pred.bas$fit, y=pred.bas$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)
points(x=dat.temp$observed_ln_recruit, y=output_df_reduced$predicted_ln_recruit,
       pch=21, bg=rgb(0,1,0.5,alpha=0.5))

# Timeseries
plot(x=dat.temp$Year, y=dat.temp$observed_ln_recruit,
     xlab="Year", ylab="ln(Recruitment)", type='l', col=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"), ylim=c(0,5), xlim=c(1975,2021))
grid(lty=3, col='dark gray')
# points(x=dat.temp$Year, y=dat.temp$observed_ln_recruit,
#        pch=21, bg=rgb(1,0,0,alpha=0.5))
lines(x=dat.temp$Year, y=dat.temp$predicted_ln_recruit, lwd=3, col=rgb(0,0,1, alpha=0.5))
lines(x=dat.temp$Year, y=output_df_reduced$predicted_ln_recruit, lwd=3, col=rgb(0,1,0.5, alpha=0.5))
# points(x=dat.temp$Year, y=output_df$predicted_ln_recruit,
#        pch=21, bg=rgb(0,1,0,alpha=0.5))

#maybe add the long series predictions to same plot?
