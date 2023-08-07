#=========================================================================================
# Fitting GAMM to sablefish indicator data from ESP
#
#=========================================================================================
#Notes:
#=========================================================================================
#packages
library(tidyverse)
library(mgcv)
library(gamm4)
library(corrplot)
library(cowplot)
library(heplots)
library(car)

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#=============================================================
#### get data ####

#first load whole dataset and look to see if any indicators look like they have nonlinear relationships
#with recruitment

scaled_dat <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)

#ONLY NONCORRELATED-----

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


subset_dat <- scaled_dat[,names(scaled_dat) %in% noncor_covars3]

plot_covars <- subset_dat %>% gather(key=type, value=value, -c(Year, ln_rec))

rec.plot <- ggplot(plot_covars, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + geom_smooth() +
  #scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
rec.plot


#plot what's left after removing NAs except 250m temp, too short

nona_covars <- na.omit(subset_dat[,-c(4)])
nona_covars <- nona_covars %>% gather(key=type, value=value, -c(Year, ln_rec))

rec.plotn <- ggplot(nona_covars, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + #geom_smooth(method="lm") +
  #scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type) +
  theme(legend.position = "NA") + ylab("ln(Recruitment)")
rec.plotn

#update with new covars if using
# labels_type <- c("Annual copepod size EGOA",
#                  "Annual copepod size WGOA",
#                  "Summer adult condition",
#                  "Summer condition age4 females",
#                  "Summer CPUE juveniles ADFG 
# large mesh survey",
# "Spring chlorophyll A biomass SEBS",
# "Spring chlorophyll A peak GOA",
# "Spring chlorophyll A peak SEBS",
# "Spring SST SEBS")

rec.plotn <- ggplot(nona_covars, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + #geom_smooth(method="lm") +
  #scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, labeller = labeller(type=labels_type)) +
  theme(legend.position = "NA") + ylab("ln(Recruitment)")
rec.plotn

#Get data=======

# load data that is already z-scored, checked for correlations


#DATA CONTROL SECTION----

#select covariates

#using ONLY the indicator subset that's been checked for collinearity
#using noncor_covars from process_data for now

#GAM cannot handle missing data, euphasiid and the every other yr GOA indicator were already removed


scaled_gam_dat <- scaled_dat[,names(scaled_dat) %in% noncor_covars3]

#leave out 2020 b/c it's a longterm mean
scaled_gam_dat <- scaled_gam_dat[which(scaled_gam_dat$Year<2020),]


#new gams-----------

#updated to run on whole dataset Jun 17, 2023

#first try
mod1_1 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
               # s(Smr_temp_250m_GOA_scaled, k=4)  +  #too short    
                s(Spr_chlA_biom_GOA_scaled, k=4)  + 
  s(Spr_chlA_biom_SEBS_scaled, k=4)  +     
  s(Spr_chlA_peak_GOA_scaled, k=4)  +     
  s(Spr_chlA_peak_SEBS_scaled, k=4)  +  
s(ann_Copepod_size_EGOA_scaled, k=4)  +     
s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
s(YOY_grwth_Middleton_scaled, k=4) +
s(smr_adult_cond_scaled, k=4), data=scaled_gam_dat) #more coefficients than data
gam.check(mod1_1) #not good
summary(mod1_1)

#next shortest are chl covars try dropping
mod1_2 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
                # s(Smr_temp_250m_GOA_scaled, k=4)  +  #too short    
               # s(Spr_chlA_biom_GOA_scaled, k=4)  + 
               # s(Spr_chlA_biom_SEBS_scaled, k=4)  +     
               # s(Spr_chlA_peak_GOA_scaled, k=4)  +     
              #  s(Spr_chlA_peak_SEBS_scaled, k=4)  +  
                s(ann_Copepod_size_EGOA_scaled, k=4)  +     
                s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
                s(YOY_grwth_Middleton_scaled, k=4) +
                s(smr_adult_cond_scaled, k=4), data=scaled_gam_dat) #runs once all chl covars removed
gam.check(mod1_2) #not good
summary(mod1_2) #nothing significantly nonlinear, only YOY_grwth has edf over 1

#change all edf < 2 to linear and try again
mod1_3 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                # s(Smr_temp_250m_GOA_scaled, k=4)  +  #too short    
                # s(Spr_chlA_biom_GOA_scaled, k=4)  + 
                # s(Spr_chlA_biom_SEBS_scaled, k=4)  +     
                # s(Spr_chlA_peak_GOA_scaled, k=4)  +     
                #  s(Spr_chlA_peak_SEBS_scaled, k=4)  +  
                ann_Copepod_size_EGOA_scaled  +     
                Smr_CPUE_juv_ADFG_ln_scaled +    
                s(YOY_grwth_Middleton_scaled, k=4) +
                smr_adult_cond_scaled, data=scaled_gam_dat) #runs once all chl covars removed
gam.check(mod1_3) #not good
summary(mod1_3) #none are significantly nonlinear

#best model
mod1_4 <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +      
                ann_Copepod_size_EGOA_scaled  +     
                Smr_CPUE_juv_ADFG_ln_scaled +    
                YOY_grwth_Middleton_scaled +
                smr_adult_cond_scaled, data=scaled_gam_dat) #
summary(mod1_4)
anova(mod1_4)
R1_4 <- residuals(mod1_4, type="pearson")
acf(R1_4)
Anova(mod1_4, type="II")

sumry1_4 <- summary(mod1_4)
efsz1_4 <- effectsize::eta_squared(car::Anova(mod1_4, type=2), partial=TRUE)

#only ADFG survey is significant!

cormod1_4 <- gls(ln_rec ~ Spr_ST_SEBS_scaled  +      
               ann_Copepod_size_EGOA_scaled  +     
               Smr_CPUE_juv_ADFG_ln_scaled +    
               YOY_grwth_Middleton_scaled +
               smr_adult_cond_scaled, 
               correlation = corCompSymm(form = ~Year),
                                           data=scaled_gam_dat, na.action = na.omit) #
AIC(cormod1_4, mod1_4)

cormod2_4 <- gls(ln_rec ~ Spr_ST_SEBS_scaled  +      
                   ann_Copepod_size_EGOA_scaled  +     
                   Smr_CPUE_juv_ADFG_ln_scaled +    
                   YOY_grwth_Middleton_scaled +
                   smr_adult_cond_scaled, 
                 correlation = corAR1(form = ~Year),
                 data=scaled_gam_dat, na.action = na.omit) #
AIC(cormod2_4, mod1_4)


#run on ONLY 4 longest time series

mod1_5 <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +        
               Smr_CPUE_juv_ADFG_ln_scaled +    
               YOY_grwth_Middleton_scaled +
               smr_adult_cond_scaled, data=scaled_gam_dat) #
summary(mod1_5)
anova(mod1_5) #still only ADFG significant

plot(mod1_5)
ggplot(scaled_gam_dat, aes(Smr_CPUE_juv_ADFG_ln_scaled, ln_rec)) + geom_point() +
  geom_smooth(method="lm")
R1_5 <- residuals(mod1_5, type="pearson")
acf(R1_5)
Anova(mod1_5, type="II")

sumry1_5 <- summary(mod1_5)
efsz1_5 <-effectsize::eta_squared(car::Anova(mod1_5, type=2), partial=TRUE)

modAL <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +      
                  Smr_temp_250m_GOA_scaled +  #    
                  Spr_chlA_biom_GOA_scaled  + 
                  Spr_chlA_biom_SEBS_scaled  +     
                  Spr_chlA_peak_GOA_scaled  +     
                  Spr_chlA_peak_SEBS_scaled  +  
                  ann_Copepod_size_EGOA_scaled  +     
                  Smr_CPUE_juv_ADFG_ln_scaled +    
                  YOY_grwth_Middleton_scaled +
                  smr_adult_cond_scaled, data=scaled_gam_dat) #
summary(modAL)
anova(modAL)

#plot effects and effect sizes------

smryplot1_4 <- as.data.frame(sumry1_4$coefficients )
smryplot1_4$var <- rownames(smryplot1_4)
smryplot1_4$stderror <- smryplot1_4$`Std. Error`

smryplot1_4$var[which(smryplot1_4$var=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
smryplot1_4$var[which(smryplot1_4$var=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
smryplot1_4$var[which(smryplot1_4$var=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
smryplot1_4$var[which(smryplot1_4$var=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
smryplot1_4$var[which(smryplot1_4$var=="smr_adult_cond_scaled")] <- "Summer adult condition"


p1 <- ggplot(smryplot1_4[which(smryplot1_4$var!="(Intercept)"),], aes(Estimate, var)) + geom_point() + ylab("Indicator") +
  xlab("Coefficient estimate") + theme_bw() + geom_errorbar(aes(xmin=Estimate-stderror,
                                                                xmax=Estimate+stderror), width = 0.1) +
  geom_vline(xintercept=0, col="red")
p1
#do same for partial eta sq

efsz1_4$Parameter[which(efsz1_4$Parameter=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
efsz1_4$Parameter[which(efsz1_4$Parameter=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
efsz1_4$Parameter[which(efsz1_4$Parameter=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
efsz1_4$Parameter[which(efsz1_4$Parameter=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
efsz1_4$Parameter[which(efsz1_4$Parameter=="smr_adult_cond_scaled")] <- "Summer adult condition"

p2 <- ggplot(efsz1_4, aes(Eta2_partial, Parameter)) + geom_point() + ylab("Indicator") +
   xlab("Partial eta squared") + theme_bw() +# geom_errorbar(aes(xmin=CI_low,
  #                                                               xmax=CI_high), width = 0.1) +
  geom_vline(xintercept=0, col="red")
p2

plot_grid(p1,p2, nrow=2, ncol=1, rel_widths=c(3,1), align='h')



#long time series

smryplot1_5 <- as.data.frame(sumry1_5$coefficients )
smryplot1_5$var <- rownames(smryplot1_5)
smryplot1_5$stderror <- smryplot1_5$`Std. Error`

smryplot1_5$var[which(smryplot1_5$var=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
smryplot1_5$var[which(smryplot1_5$var=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
smryplot1_5$var[which(smryplot1_5$var=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
smryplot1_5$var[which(smryplot1_5$var=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
smryplot1_5$var[which(smryplot1_5$var=="smr_adult_cond_scaled")] <- "Summer adult condition"


p3 <- ggplot(smryplot1_5[which(smryplot1_5$var!="(Intercept)"),], aes(Estimate, var)) + geom_point() + ylab("Indicator") +
  xlab("Coefficient estimate") + theme_bw() + geom_errorbar(aes(xmin=Estimate-stderror,
                                                                xmax=Estimate+stderror), width = 0.1) +
  geom_vline(xintercept=0, col="red")
p3
#do same for partial eta sq

efsz1_5$Parameter[which(efsz1_5$Parameter=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
efsz1_5$Parameter[which(efsz1_5$Parameter=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
efsz1_5$Parameter[which(efsz1_5$Parameter=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
efsz1_5$Parameter[which(efsz1_5$Parameter=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
efsz1_5$Parameter[which(efsz1_5$Parameter=="smr_adult_cond_scaled")] <- "Summer adult condition"

p4 <- ggplot(efsz1_5, aes(Eta2_partial, Parameter)) + geom_point() + ylab("Indicator") +
  xlab("Partial eta squared") + theme_bw() +# geom_errorbar(aes(xmin=CI_low,
  #                                                               xmax=CI_high), width = 0.1) +
  geom_vline(xintercept=0, col="red")
p4

plot_grid(p3,p4, nrow=2, ncol=1, rel_widths=c(3,1), align='h')




#within sample pred plot-----

withinpred <- predict(mod1_4, scaled_gam_dat[which(scaled_gam_dat$Year>2000),c(3:12)])
plotwithin <- scaled_gam_dat[which(scaled_gam_dat$Year>2000),]
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


#

#within sample pred plot LONG time series-----

withinpred <- predict(mod1_5, scaled_gam_dat[which(scaled_gam_dat$Year>2000),c(3:12)])
plotwithin <- scaled_gam_dat[which(scaled_gam_dat$Year>2000),]
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




#LOOCV - all time series=======

#STEP 1 - Loop through training sets and fit models-------

#
scaled_loop_dat <- scaled_gam_dat

yrs <- unique(scaled_loop_dat$Year)
output_df <- data.frame(matrix(ncol=3, nrow = length(yrs)))
colnames(output_df) <- c("Year", "observed_ln_recruit", "predicted_ln_recruit")

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  dropped_yr <- scaled_loop_dat[i,]
  output_df$observed_ln_recruit[i] <- dropped_yr$ln_rec
  dropped_yr <- dropped_yr[,!names(dropped_yr) %in% "ln_rec"]
  #fit model
  modtemp <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +      
                  ann_Copepod_size_EGOA_scaled  +     
                  Smr_CPUE_juv_ADFG_ln_scaled +    
                  YOY_grwth_Middleton_scaled +
                  smr_adult_cond_scaled, data=temp_dat)
  #have model predict to missing year
   temp_predict <- predict(modtemp, newdata=dropped_yr)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict[1]
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  #geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year))+
  ylim(c(0,5)) + xlim(c(0,5)) + theme_bw()

#STEP 2 - get MSE, MAE, and R2------

#get MSE & MAE------
library(yardstick)
#these need to be double checked!
#GAM_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)


obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
  geom_point() + geom_smooth(method="lm")

GAM_rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

GAM_mae <- mae(output_df, truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)



output_df$diff <- output_df$predicted_ln_recruit - output_df$observed_ln_recruit

ggplot(output_df, aes(Year, diff, col=as.numeric(Year))) + 
  geom_point() + geom_smooth(method="lm")

write.csv(output_df, file=paste(wd,"/data/GAM_obsvpreds_reduced.csv", sep=""))







#LOOCV - long time series=======

#STEP 1 - Loop through training sets and fit models-------

#
scaled_loop_dat <- scaled_gam_dat

yrs <- unique(scaled_loop_dat$Year)
output_df <- data.frame(matrix(ncol=3, nrow = length(yrs)))
colnames(output_df) <- c("Year", "observed_ln_recruit", "predicted_ln_recruit")

i<-1
for(i in 1:length(scaled_loop_dat$Year)){
  print(i)
  temp_dat <- scaled_loop_dat[-i,]
  dropped_yr <- scaled_loop_dat[i,]
  output_df$observed_ln_recruit[i] <- dropped_yr$ln_rec
  dropped_yr <- dropped_yr[,!names(dropped_yr) %in% "ln_rec"]
  #fit model
  modtemp <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +        
                  Smr_CPUE_juv_ADFG_ln_scaled +    
                  YOY_grwth_Middleton_scaled +
                  smr_adult_cond_scaled, data=temp_dat)
  #have model predict to missing year
  temp_predict <- predict(modtemp, newdata=dropped_yr)
  #write to output object so we can compare predicted vs obs
  output_df$Year[i] <- dropped_yr$Year
  output_df$predicted_ln_recruit[i] <- temp_predict[1]
}


output_df$predicted_ln_recruit <- as.numeric(as.character(output_df$predicted_ln_recruit))

ggplot(output_df, aes(observed_ln_recruit, predicted_ln_recruit)) + 
  #geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(observed_ln_recruit, predicted_ln_recruit, label=Year))+
  ylim(c(0,5)) + xlim(c(0,5)) + theme_bw()

#STEP 2 - get MSE, MAE, and R2------

#get MSE & MAE------

#these need to be double checked!
#GAM_MSE <- ((sum((output_df$observed_ln_recruit - output_df$predicted_ln_recruit)^2, na.rm = TRUE)))/length(output_df$observed_ln_recruit)


obs_pred_mod <- lm(predicted_ln_recruit ~ observed_ln_recruit, data=output_df)
summary(obs_pred_mod)

GAM_long_rmse <- rmse(output_df, truth=observed_ln_recruit, 
                 estimate=predicted_ln_recruit, na.rm=TRUE)

GAM_long_mae <- mae(output_df, truth=observed_ln_recruit, 
               estimate=predicted_ln_recruit, na.rm=TRUE)





