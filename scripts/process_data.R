#====================================================================================
#Data processing - run before running analyses
#
#Krista, March 2023
#====================================================================================
#Notes:
#====================================================================================



require(tidyverse)
require(ggthemes)
require(viridis)
#require(R2jags)
require(corrplot)
require(reshape2)
library(car)


#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#Get data plot data =======

#dat <- read.csv(file.path(dir.data, "sablefish_BAS_indicators_2022.csv"))
dat <- read.csv(file.path(dir.data, "sablefish_esp_2022_updated.csv"), skip=1)



#before scaling, let's look and see if any need to be normalized

covar.list <- dat %>% gather(key=type, value=value, -Year)
head(covar.list)

expore.plot <- ggplot(covar.list, aes(x=Year, y=value, fill=type)) +
  geom_point() +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
expore.plot



expore.hist <- ggplot(covar.list, aes(x=value, fill=type)) +
  theme_linedraw() +
  geom_histogram() +
  geom_density(alpha=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
expore.hist


#dat$ln_rec <- log(dat$Recruitment)
dat$ln_rec <- log(dat$Recruitment_year_class)


#log data that doesn't look normal
#Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery
#Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey


#NOT normal but not sure logging will help:
#Summer_Temperature_250m_GOA_Survey
#Annual_Arrowtooth_Biomass_GOA_Model

#euphausiid data is super sparse, will need to come out for most analyses, 
#retain for now, maybe can go into some

log.hist <- ggplot(covar.list, aes(x=log(value), fill=type)) +
  theme_linedraw() +
  geom_histogram() +
  geom_density(alpha=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
log.hist

#logging seems to help the CPUE covars but not heatwave, arrowtooth, or summer 250m GOA temp

dat$ln_Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey <- log(dat$Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey)
dat$ln_Summer_Sablefish_CPUE_Juvenile_GOA_Survey <- log(dat$Summer_Sablefish_CPUE_Juvenile_GOA_Survey)

#also incident catch in arrowtooth fishery and log + 1 of heatwave
dat$ln_Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery <- log(dat$Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery)
dat$ln_plus_Annual_Heatwave_GOA_Model <- log(dat$Annual_Heatwave_GOA_Model+1)
#heatwave gets log + 1 because it has a ton of zeros

#remove columns that have now been lagged
dat <- dat[,!names(dat) %in% c("Summer_Sablefish_CPUE_Juvenile_GOA_Survey",
                                        "Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey",
                                        "Recruitment_year_class",
                                        "Annual_Heatwave_GOA_Model",
                                        "Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery")]


covar.noR <- dat %>% gather(key=type, value=value, -c(Year, ln_rec))


ln.rec.plot <- ggplot(covar.noR, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + geom_smooth() +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
ln.rec.plot



#z-score data=======

# scaled_dat <- dat %>% #group_by(Year) %>%
#   mutate(recruit_scaled=scale(ln_rec),
#          Spr_ST_SEBS_scaled=scale(Spring_Temperature_Surface_SEBS_Satellite),
#          YOY_grwth_Middleton_scaled=scale(Annual_Sablefish_Growth_YOY_Middleton_Survey),
#          Smr_CPUE_juv_ADFG_scaled=scale(Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey),
#          spawner_mean_age_scaled=scale(Annual_Sablefish_Mean_Age_Female_Adult_Model),
#          spawner_age_evenness_scaled=scale(Annual_Sablefish_Age_Evenness_Female_Adult_Model),
#          arrowtooth_biomass_scaled=scale(Annual_Arrowtooth_Biomass_GOA_Model),
#          sablefish_bycatch_arrowtooth_fishery_scaled=scale(Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery),
#          smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Female_Adult_GOA_Survey))
# 

scaled_dat <- dat %>% #group_by(Year) %>%
  mutate(recruit_scaled=scale(ln_rec),
         ann_heatwave_GOA_scaled=scale(ln_plus_Annual_Heatwave_GOA_Model),                                      
          Spr_ST_GOA_scaled=scale(Spring_Temperature_Surface_GOA_Satellite),
         Spr_ST_SEBS_scaled=scale(Spring_Temperature_Surface_SEBS_Satellite),
         Smr_temp_250m_GOA_scaled = scale(Summer_Temperature_250m_GOA_Survey),
         Spr_chlA_biom_GOA_scaled = scale(Spring_Chlorophylla_Biomass_GOA_Satellite),                      
         Spr_chlA_biom_SEBS_scaled = scale(Spring_Chlorophylla_Biomass_SEBS_Satellite),                     
         Spr_chlA_peak_GOA_scaled = scale(Spring_Chlorophylla_Peak_GOA_Satellite),                         
         Spr_chlA_peak_SEBS_scaled = scale(Spring_Chlorophylla_Peak_SEBS_Satellite),                      
         ann_Copepod_size_EGOA_scaled = scale(Annual_Copepod_Community_Size_EGOA_Survey),                      
         ann_Copepod_size_WGOA_scaled = scale(Annual_Copepod_Community_Size_WGOA_Survey),                      
         Smr_euph_abun_Kod_scaled=scale(Summer_Euphausiid_Abundance_Kodiak_Survey),
         YOY_grwth_Middleton_scaled=scale(Annual_Sablefish_Growth_YOY_Middleton_Survey),
         Smr_CPUE_juv_ADFG_ln_scaled=scale(ln_Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey),
         Smr_CPUE_juv_GOA_ln_scaled=scale(ln_Summer_Sablefish_CPUE_Juvenile_GOA_Survey),
         spawner_mean_age_scaled=scale(Annual_Sablefish_Mean_Age_Female_Adult_Model),
         spawner_age_evenness_scaled=scale(Annual_Sablefish_Age_Evenness_Female_Adult_Model),
         Smr_condition_fem_age4_GOA_scaled=scale(Summer_Sablefish_Condition_Female_Age4_GOA_Survey),
         arrowtooth_biomass_scaled=scale(Annual_Arrowtooth_Biomass_GOA_Model),
         sablefish_bycatch_arrowtooth_fishery_scaled=scale(ln_Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery),
         smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Female_Adult_GOA_Survey))

#save
write.csv(scaled_dat, file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""))


#split data======

#NEED TO SUBSET out a training dataset

#do this five times so it can be repeated
#get sample size (n years)
datlen <- length(scaled_dat$recruit_scaled)

#randomly get seeds
seeds <- sample(1:datlen, 5)
#MY SEEDS are 31, 1, 20, 35, 3
seeds <- c(31,1,20,35,3)

set.seed(seeds[1])
train1 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),] #sample 80% of data 
train1 <- train1[order(as.numeric(row.names(train1))),] #re-order to correct order of years
testing1 <- anti_join(scaled_dat, train1) #remaining rows become testing data

set.seed(seeds[2])
train2 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train2 <- train2[order(as.numeric(row.names(train2))),]
testing2 <- anti_join(scaled_dat, train2)

set.seed(seeds[3])
train3 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train3 <- train3[order(as.numeric(row.names(train3))),]
testing3 <- anti_join(scaled_dat, train3)

set.seed(seeds[4])
train4 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train4 <- train4[order(as.numeric(row.names(train4))),]
testing4 <- anti_join(scaled_dat, train4)

set.seed(seeds[5])
train5 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train5 <- train5[order(as.numeric(row.names(train5))),]
testing5 <- anti_join(scaled_dat, train5)

#let's plot these to look at which years get included

yrs1 <- as.data.frame(train1$Year)
colnames(yrs1) <- "Year"
yrs1$set <- 1

yrs2 <- as.data.frame(train2$Year)
colnames(yrs2) <- "Year"
yrs2$set <- 2

yrs3 <- as.data.frame(train3$Year)
colnames(yrs3) <- "Year"
yrs3$set <- 3

yrs4 <- as.data.frame(train4$Year)
colnames(yrs4) <- "Year"
yrs4$set <- 4

yrs5 <- as.data.frame(train5$Year)
colnames(yrs5) <- "Year"
yrs5$set <- 5

yrs_df <- rbind(yrs1, yrs2, yrs3, yrs4, yrs5)
yrs_df$set <- as.factor(yrs_df$set)

ggplot(yrs_df, aes(Year, set)) + geom_point()

#save analysis-ready data======

write.csv(train1, file=paste(wd,"/data/dataset_training1.csv", sep=""))
write.csv(train2, file=paste(wd,"/data/dataset_training2.csv", sep=""))
write.csv(train3, file=paste(wd,"/data/dataset_training3.csv", sep=""))
write.csv(train4, file=paste(wd,"/data/dataset_training4.csv", sep=""))
write.csv(train5, file=paste(wd,"/data/dataset_training5.csv", sep=""))

write.csv(testing1, file=paste(wd,"/data/dataset_testing1.csv", sep=""))
write.csv(testing2, file=paste(wd,"/data/dataset_testing2.csv", sep=""))
write.csv(testing3, file=paste(wd,"/data/dataset_testing3.csv", sep=""))
write.csv(testing4, file=paste(wd,"/data/dataset_testing4.csv", sep=""))
write.csv(testing5, file=paste(wd,"/data/dataset_testing5.csv", sep=""))



#missing data------------------------------------------------------


#scaled only
scaled_only <- scaled_dat[,c(1,18, 24:43)]

scal_long <- pivot_longer(scaled_only, -Year, names_to = "covar", values_to = "value")

ggplot(scal_long, aes(Year, covar, size=value)) + geom_point()


#create subsets without highly correlated indicators---------------

#let's check for highly correlated indicators in steps
#first, let's look at all indicators EXCEPT Smr_euph_abun_Kod_scaled, only n=7, and 
#Smr_CPUE_juv_GOA_ln_scaled which is every 2nd yr

sub1 <- scaled_only[,-c(1,13,16)]
cov1 <- na.omit(sub1)

cov.cor1 <- cor(cov1)
corrplot(cov.cor1,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor1, upper='circle', lower='number')

#overwhelming! Let's break apart

subsub1 <- sub1[,c(1:9)]
subsub1.2 <- sub1[,c(1,10:19)]

covsub1 <- na.omit(subsub1)
covsub1.2 <- na.omit(subsub1.2)

cov.cor1.1 <- cor(covsub1)
corrplot(cov.cor1.1,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor1.1, upper='circle', lower='number')

#problem pairs
#"ann_heatwave_GOA_scaled"  + "Spr_ST_SEBS_scaled"    
#"Spr_ST_GOA_scaled"   + "Spr_ST_SEBS_scaled"         
# "ann_heatwave_GOA_scaled"  + "Spr_ST_GOA_scaled"


cov.cor1.2 <- cor(covsub1.2)
corrplot(cov.cor1.2,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor1.2, upper='circle', lower='number')

#problem pairs
#sablefish_bycatch_arrowtooth_fishery_scaled + arrowtooth_biomass_scaled
#sablefish_bycatch_arrowtooth_fishery_scaled + Smr_CPUE_juv_ADFG_ln_scaled

#ann_Copepod_size_EGOA_scaled + Spr_ST_SEBS_scaled
#ann_Copepod_size_EGOA_scaled + ann_heatwave_GOA_scaled
#ann_Copepod_size_EGOA_scaled + Spr_ST_GOA_scaled

#YOY_grwth_Middleton_scaled + Spr_ST_SEBS_scaled
#YOY_grwth_Middleton_scaled + Spr_ST_GOA_scaled

#sablefish_bycatch_arrowtooth_fishery_scaled + arrowtooth_biomass_scaled
#sablefish_bycatch_arrowtooth_fishery_scaled + Smr_CPUE_juv_ADFG_ln_scaled
#Smr_CPUE_juv_ADFG_ln_scaled + arrowtooth_biomass_scaled

#remove model derived problem indicators
sub3 <- sub1[,-c(14:15)]
cov3 <- na.omit(sub3)

cov.cor3 <- cor(cov3)
corrplot(cov.cor3,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor3, upper='circle', lower='number')


#drop those selected for removal

sub2 <- scaled_only[,-c(1,4:5,11,13,16,17:18,20:21)]
cov2 <- na.omit(sub2)

cov.cor2 <- cor(cov2)
corrplot(cov.cor2,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor2, upper='circle', lower='number')

#add juv GOA CPUE back in

sub3 <- scaled_only[,-c(1,4:5,11,13,17:18,20:21)]
cov3 <- na.omit(sub3)

cov.cor3 <- cor(cov3)
corrplot(cov.cor3,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor3, upper='circle', lower='number')
#BAD

#try only good ones
noncor_only <- scaled_only[,-c(1,4:5,11,13,16,17:18,20:21)]
cov4 <- na.omit(noncor_only)

cov.cor4 <- cor(cov4)
corrplot(cov.cor4,order='AOE',  type = 'lower', method = 'number')
corrplot.mixed(cov.cor4, upper='circle', lower='number')

#VIFs-------------------------

#let's do this without heatwave (creates problems)

vmod <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
            Spr_ST_SEBS_scaled +
             Smr_temp_250m_GOA_scaled  + 
            Spr_chlA_biom_GOA_scaled +  Spr_chlA_biom_SEBS_scaled    +  
           Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
            ann_Copepod_size_EGOA_scaled  +  
             ann_Copepod_size_WGOA_scaled   +    
           # Smr_euph_abun_Kod_scaled  +  
             YOY_grwth_Middleton_scaled  +     
            Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
            #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
            #Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
           # sablefish_bycatch_arrowtooth_fishery_scaled + 
             smr_adult_cond_scaled , data=scaled_only)


vifs <- car::vif(vmod) #not terrible

#drop worst first, that's currently ann_Copepod_size_WGOA_scaled

vmod2 <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
              Spr_ST_SEBS_scaled +
             Smr_temp_250m_GOA_scaled  + 
             Spr_chlA_biom_GOA_scaled +  
              Spr_chlA_biom_SEBS_scaled    +  
             Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
              ann_Copepod_size_EGOA_scaled  +  
             #ann_Copepod_size_WGOA_scaled   +    
             # Smr_euph_abun_Kod_scaled  +  
             YOY_grwth_Middleton_scaled  +     
             Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
             #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
            # Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
             # sablefish_bycatch_arrowtooth_fishery_scaled + 
             smr_adult_cond_scaled , data=scaled_only)


vifs2 <- car::vif(vmod2) #all fine!


#UPDATE HERE
#what are the indicators that remain?
noncor_covars <- c("Year", "ln_rec",
                   "ann_heatwave_GOA_scaled"  ,
  "Smr_temp_250m_GOA_scaled"  , 
  "Spr_chlA_biom_SEBS_scaled"  ,  
  "Spr_chlA_peak_GOA_scaled" ,
 "Spr_chlA_peak_SEBS_scaled" ,
  "ann_Copepod_size_WGOA_scaled" ,
  "YOY_grwth_Middleton_scaled" ,     
  "Smr_CPUE_juv_ADFG_ln_scaled"  ,      
  "Smr_condition_fem_age4_GOA_scaled" ,
  "smr_adult_cond_scaled")

#heatwave index causes issues, let's make a second set without it
noncor_covars2 <- c("Year", "ln_rec",
                   "Spr_ST_SEBS_scaled"  ,
                   "Smr_temp_250m_GOA_scaled"  , 
                   "Spr_chlA_biom_SEBS_scaled"  ,  
                   "Spr_chlA_peak_GOA_scaled" ,
                   "Spr_chlA_peak_SEBS_scaled" ,
                   "ann_Copepod_size_EGOA_scaled" ,
                   "ann_Copepod_size_WGOA_scaled" ,
                   "Smr_CPUE_juv_ADFG_ln_scaled"  ,      
                   "Smr_condition_fem_age4_GOA_scaled" ,
                   "smr_adult_cond_scaled")

#heatwave index and age 4 condition removed
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


noncor3_dat <- scaled_dat[,names(scaled_dat) %in% noncor_covars3]
noncor3_dat$complete <- complete.cases(noncor3_dat)

noncor3_long <- noncor3_dat %>% gather(key=type, value=value, -c(Year, ln_rec, complete))

noncor3_long$type[which(noncor3_long$type=="Spr_ST_SEBS_scaled")] <- "Spring SST SEBS"
noncor3_long$type[which(noncor3_long$type=="Smr_temp_250m_GOA_scaled")] <- "Summer 250m temperature GOA"
noncor3_long$type[which(noncor3_long$type=="Spr_chlA_biom_GOA_scaled")] <- "Spring chlorophyll A biomass GOA"
noncor3_long$type[which(noncor3_long$type=="Spr_chlA_biom_SEBS_scaled")] <- "Spring chlorophyll A biomass SEBS"
noncor3_long$type[which(noncor3_long$type=="Spr_chlA_peak_GOA_scaled")] <- "Spring chlorophyll A peak GOA"
noncor3_long$type[which(noncor3_long$type=="Spr_chlA_peak_SEBS_scaled")] <- "Spring chlorophyll A peak SEBS"
noncor3_long$type[which(noncor3_long$type=="ann_Copepod_size_EGOA_scaled")] <- "Annual copepod community size EGOA"
noncor3_long$type[which(noncor3_long$type=="YOY_grwth_Middleton_scaled")] <- "YOY growth Middleton Is. seabirds"
noncor3_long$type[which(noncor3_long$type=="Smr_CPUE_juv_ADFG_ln_scaled")] <- "Summer juvenile CPUE ADFG survey"
noncor3_long$type[which(noncor3_long$type=="smr_adult_cond_scaled")] <- "Summer adult condition"




#plot years
  ggplot(noncor3_long, aes(Year, type, size=value, col=value)) + geom_point() + theme_bw() + theme(legend.position = "none") +
    xlab("Year class") + ylab("")
  
  #plot data
  
  ggplot(noncor3_long, aes(value, ln_rec)) + geom_point() + facet_wrap(~type, nrow=4) +
    xlab("Scaled value") + ylab("ln(Recruitment)")
  
  ggplot(noncor3_long, aes(value, ln_rec, col=complete)) + geom_point() + facet_wrap(~type, nrow=4) +
    xlab("Scaled value") + ylab("ln(Recruitment)") + scale_colour_manual(values=c("red", "black"))
  

#OLD======

#now automate the process

response <-  "ln_rec"
explanatory_vars <- colnames(noncor_only[,!names(noncor_only) %in% c("Year", "ln_rec")])

form.vars <- paste(explanatory_vars, collapse=" + ")
form <- paste(response, "~",form.vars)

vif_mod <- lm(formula = form, data=noncor_only)
vif1 <- car::vif(vif_mod)
vif1df <- as.data.frame(vif1)

vif1ordered <- order(vif1df[,1]) #none of this is working!!!
highest <- rownames(vif1df[which(vif1df[,1]==max(vif1df)),] )

#uncentered vifs?-------

vif.lm <- function(object, ...) {
  V <- summary(object)$cov.unscaled
  Vi <- crossprod(model.matrix(object))
  nam <- names(coef(object))
  k <- match("(Intercept)", nam,
             nomatch = FALSE)
  v1 <- diag(V)
  v2 <- diag(Vi)
  uc.struct <- structure(v1 * v2, names = nam)
  if(k) {
    v1 <- diag(V)[-k]
    v2 <- diag(Vi)[-k] - Vi[k, -k]^2 / Vi[k, k]
    nam <- nam[-k]
    c.struct <- structure(v1 * v2, names = nam)
  
           return(c(Centered.VIF = c.struct, Uncentered.VIF = uc.struct))
  }
  else{
    warning(paste("No intercept term",
                  "detected. Uncentered VIFs computed."))
    return(Uncentered.VIF = uc.struct)
  }
}

vif.lm(vmod) #also picks Spr_chlA_biom_GOA_scaled
#after dropping that one
vif.lm(vmod2) #all under 10
