#============================================================================================================================================
# get 5-yr running avg recruitment to compare to models

#Created by Krista, 2023
#============================================================================================================================================
#Notes:
# colours:
# dark green, predicted "#1b9e77"
# dark orange, long time series predicted "#d95f02"
# purple, reduced time series predicted "#7570b3"
#============================================================================================================================================

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#Get data =======


scaled <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)

#THIS WILL VARY BY MODEL TYPE

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

scaled_5yr_dat <- scaled[,names(scaled) %in% noncor_covars3]


#loop through----

out_avgs <- as.data.frame(matrix(ncol = 2,nrow = length((1977+5):2020)))
colnames(out_avgs) <- c("Year", "avg_prev_5yr_recruit")

i<-1
counter <- 1
for(i in (1977+5):2020){
  temp_dat <- scaled_5yr_dat[which(scaled_5yr_dat$Year>(i-6) & 
                                     scaled_5yr_dat$Year<(i) ),]
  temp_mean <- mean(temp_dat$ln_rec)
  out_avgs[counter,1] <- i
    out_avgs[counter,2] <- temp_mean
    counter <- counter + 1
}

avgs_join <- left_join(scaled_5yr_dat[,c(1:2)], out_avgs)

#plot within sample predictions----


ggplot(avgs_join, aes(Year, ln_rec)) + geom_point(aes(col="red")) + 
  geom_line() + geom_point(aes(Year, avg_prev_5yr_recruit)) + geom_line(aes(Year, avg_prev_5yr_recruit)) + theme_bw()


#from BAS
# Plot Model Predictions vs. Observed ==============================
#pdf(file.path(dir.figs,"Model Fit.pdf"), height=6, width=9)
par(oma=c(1,1,1,1), mar=c(4,4,1,1), mfrow=c(1,2))

# Omit NAs
dat.temp <- avgs_join

plot(x=dat.temp$ln_rec, y=dat.temp$avg_prev_5yr_recruit,
     xlab="Observed ln(recruitment)", ylab="5yr avg ln(recruitment)", pch=21, bg="#1b9e77",
     main=paste("5-year rolling average"), ylim=c(0,5), xlim=c(0,5))
# plot(x=plotwithin$fit, y=plotwithin$Ybma) 
abline(a=0, b=1, col="dark grey", lwd=3)

# Timeseries
plot(x=dat.temp$Year, y=dat.temp$ln_rec,
     xlab="Year", ylab="ln(recruitment)", type='l', col="black",
     main=paste("5-year rolling average"))
grid(lty=3, col='dark gray')
# points(x=dat.temp$Year, y=dat.temp$ln_rec,
#        pch=21, bg="black")
lines(x=dat.temp$Year, y=dat.temp$avg_prev_5yr_recruit, lwd=3, col="#1b9e77")
#points(x=dat.temp$Year, y=plotwithin$predicted,
#      pch=21, bg=rgb(0,1,0,alpha=0.5))

legend('topleft', legend=c("Observed","5yr avg"), lty=1, col=c("black",
                                                               "#1b9e77"), bg="white")

#comparison metrics-----

obs_5_mod <- lm(ln_rec ~ avg_prev_5yr_recruit, data=dat.temp)
summary(obs_5_mod)

avg5yr_rmse <- rmse(dat.temp, truth=ln_rec, 
                 estimate=avg_prev_5yr_recruit, na.rm=TRUE)

avg5yr_mae <- mae(dat.temp, truth=ln_rec, 
               estimate=avg_prev_5yr_recruit, na.rm=TRUE)

#from 1996

obs_5_mod96 <- lm(ln_rec ~ avg_prev_5yr_recruit, data=dat.temp[which(dat.temp$Year>1995),])
summary(obs_5_mod96)

avg5yr_rmse96 <- rmse(dat.temp[which(dat.temp$Year>1995),], truth=ln_rec, 
                    estimate=avg_prev_5yr_recruit, na.rm=TRUE)

avg5yr_mae96 <- mae(dat.temp[which(dat.temp$Year>1995),], truth=ln_rec, 
                  estimate=avg_prev_5yr_recruit, na.rm=TRUE)



#from 2002

obs_5_mod02 <- lm(ln_rec ~ avg_prev_5yr_recruit, data=dat.temp[which(dat.temp$Year>2001),])
summary(obs_5_mod02)

avg5yr_rmse02 <- rmse(dat.temp[which(dat.temp$Year>2001),], truth=ln_rec, 
                    estimate=avg_prev_5yr_recruit, na.rm=TRUE)

avg5yr_mae02 <- mae(dat.temp[which(dat.temp$Year>2001),], truth=ln_rec, 
                  estimate=avg_prev_5yr_recruit, na.rm=TRUE)



