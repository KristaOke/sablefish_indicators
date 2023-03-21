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

dat$ln_Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery <- log(dat$Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery)
dat$ln_plus_Annual_Heatwave_GOA_Model <- log(dat$Annual_Heatwave_GOA_Model)

#also incident catch in arrowtooth fishery and log + 1 of heatwave

dat <- dat[,-c("Summer_Sablefish_CPUE_Juvenile_GOA_Survey",
               "Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey",
               "Recruitment_year_class")]

#not working let's brute force it for now, 
#RETURN TO THIS
dat <- dat[,c(1, 4:14, 17:20, 22:27)]


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



#split data======

#NEED TO SUBSET out a training dataset
#do this five times so it can be repeated
datlen <- length(scaled_dat$recruit_scaled)
train1 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train1 <- train1[order(row.names(train1)),]
testing1 <- anti_join(scaled_dat, train1)

train2 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train2 <- train2[order(row.names(train2)),]
testing2 <- anti_join(scaled_dat, train2)

train3 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train3 <- train3[order(row.names(train3)),]
testing3 <- anti_join(scaled_dat, train3)

train4 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train4 <- train4[order(row.names(train4)),]
testing4 <- anti_join(scaled_dat, train4)

train5 <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train5 <- train5[order(row.names(train5)),]
testing5 <- anti_join(scaled_dat, train5)

#save this training set?

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






















