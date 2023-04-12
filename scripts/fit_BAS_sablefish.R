#==================================================================================================
#Fit BAS to training and testing datasets
#
#
#Krista, March 2023
#BORROWING HEAVILY from scripts developped by Curry Cunningham for ESP analyses
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


train1_bas_dat <- train1[,names(train1) %in% noncor_covars3]
train2_bas_dat <- train2[,names(train2) %in% noncor_covars3]
train3_bas_dat <- train3[,names(train3) %in% noncor_covars3]
train4_bas_dat <- train4[,names(train4) %in% noncor_covars3]
train5_bas_dat <- train5[,names(train5) %in% noncor_covars3]

test1_bas_dat <- testing1[,names(testing1) %in% noncor_covars3]
test2_bas_dat <- testing2[,names(testing2) %in% noncor_covars3]
test3_bas_dat <- testing3[,names(testing3) %in% noncor_covars3]
test4_bas_dat <- testing4[,names(testing4) %in% noncor_covars3]
test5_bas_dat <- testing5[,names(testing5) %in% noncor_covars3]

#remove rows missing values???

# train1_brt_dat <- train1_brt_dat[which(is.na(train1_brt_dat$ln_rec)==FALSE),]
# train2_brt_dat <- train2_brt_dat[which(is.na(train2_brt_dat$ln_rec)==FALSE),]
# train3_brt_dat <- train3_brt_dat[which(is.na(train3_brt_dat$ln_rec)==FALSE),]
# train4_brt_dat <- train4_brt_dat[which(is.na(train4_brt_dat$ln_rec)==FALSE),]
# train5_brt_dat <- train5_brt_dat[which(is.na(train5_brt_dat$ln_rec)==FALSE),]
# 
# test1_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
# test2_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
# test3_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
# test4_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]
# test5_brt_dat <- test1_brt_dat[which(is.na(test1_brt_dat$ln_rec)==FALSE),]





#Determine Covariates

covars <- names(train1_bas_dat)[-which(names(train1_bas_dat) %in% c("Year", "ln_rec"))]

n.cov <- length(covars)

dat.temp <- train1_bas_dat[-which(names(train1_bas_dat) %in% c("Year"))]

#more parameters than data, let's try dropping the shortest time series
#Smr_temp_250m_GOA_scaled
dat.temp <- dat.temp[-which(names(dat.temp) %in% c("Smr_temp_250m_GOA_scaled"))]

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
dat.temp.na.omit <- na.omit(train1_bas_dat[-which(names(train1_bas_dat) %in% c("Smr_temp_250m_GOA_scaled"))])

plot(x=dat.temp.na.omit$ln_rec, y=pred.bas$Ybma,
     xlab="Observed ln(Rec)", ylab="Predicted ln(Rec)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=pred.bas$fit, y=pred.bas$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)

# Timeseries
plot(x=dat.temp.na.omit$Year, y=dat.temp.na.omit$ln_rec,
     xlab="Year", ylab="ln(Rec)", type='l', col=rgb(1,0,0,alpha=0.5),
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
ggsave(file=file.path(dir.figs,"BAS.png"), plot=g3, height=5, width=8, units='in',
       dpi=500)


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
ggsave(file=file.path(dir.figs,"BAS_noRainbow.png"), plot=g3.b, height=5, width=8, units='in',
       dpi=500)

#REPEAT without short time series================================================================

#some covars are not important and also short, let's try removing two of them
#Smr_temp_250m_GOA_scaled starts 2005
#Smr_condition_fem_age4_GOA_scaled ends 2017


#Determine Covariates

covars <- names(train1_bas_dat)[-which(names(train1_bas_dat) %in% c("Year", "ln_rec",
                                                                    "Smr_temp_250m_GOA_scaled",
                                                                    "Smr_condition_fem_age4_GOA_scaled"))]

n.cov <- length(covars)

dat.temp <- train1_bas_dat[-which(names(train1_bas_dat) %in% c("Year"))]

#more parameters than data, let's try dropping the shortest time series
#Smr_temp_250m_GOA_scaled
dat.temp <- dat.temp[-which(names(dat.temp) %in% c("Smr_temp_250m_GOA_scaled", "Smr_condition_fem_age4_GOA_scaled"))]

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
dat.temp.na.omit <- na.omit(train1_bas_dat[-which(names(train1_bas_dat) %in% c("Smr_temp_250m_GOA_scaled", "Smr_condition_fem_age4_GOA_scaled"))])

plot(x=dat.temp.na.omit$ln_rec, y=pred.bas$Ybma,
     xlab="Observed ln(Rec)", ylab="Predicted ln(Rec)", pch=21, bg=rgb(1,0,0,alpha=0.5),
     main=paste("Sablefish"))
# plot(x=pred.bas$fit, y=pred.bas$Ybma) 
abline(a=0, b=1, col=rgb(0,0,1,alpha=0.5), lwd=3)

# Timeseries
plot(x=dat.temp.na.omit$Year, y=dat.temp.na.omit$ln_rec,
     xlab="Year", ylab="ln(Rec)", type='l', col=rgb(1,0,0,alpha=0.5),
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
ggsave(file=file.path(dir.figs,"BAS.png"), plot=g3, height=5, width=8, units='in',
       dpi=500)


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
ggsave(file=file.path(dir.figs,"BAS_noRainbow.png"), plot=g3.b, height=5, width=8, units='in',
       dpi=500)




