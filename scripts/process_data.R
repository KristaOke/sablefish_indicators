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

dat$ln_Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery <- log(dat$Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery)
dat$ln_plus_Annual_Heatwave_GOA_Model <- log(dat$Annual_Heatwave_GOA_Model+1)
#heatwave gets log + 1 because it has a ton of zeros

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

vmod <- lm(ln_rec ~ ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
           # Spr_ST_SEBS_scaled +
             Smr_temp_250m_GOA_scaled  + 
            Spr_chlA_biom_GOA_scaled +  Spr_chlA_biom_SEBS_scaled    +  
           Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
           # ann_Copepod_size_EGOA_scaled  +  
             ann_Copepod_size_WGOA_scaled   +    
           # Smr_euph_abun_Kod_scaled  +  
             YOY_grwth_Middleton_scaled  +     
            Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
            #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
            Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
           # sablefish_bycatch_arrowtooth_fishery_scaled + 
             smr_adult_cond_scaled , data=scaled_only)


vifs <- car::vif(vmod) #ALL BAD?!

#drop worst first, that's Spr_chlA_biom_GOA_scaled  at 31

vmod2 <- lm(ln_rec ~ ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
             # Spr_ST_SEBS_scaled +
             Smr_temp_250m_GOA_scaled  + 
            # Spr_chlA_biom_GOA_scaled +  
              Spr_chlA_biom_SEBS_scaled    +  
             Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
             # ann_Copepod_size_EGOA_scaled  +  
             ann_Copepod_size_WGOA_scaled   +    
             # Smr_euph_abun_Kod_scaled  +  
             YOY_grwth_Middleton_scaled  +     
             Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
             #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
             Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
             # sablefish_bycatch_arrowtooth_fishery_scaled + 
             smr_adult_cond_scaled , data=scaled_only)


vifs2 <- car::vif(vmod2) #better! Still some above 5

#would it be better to drop SEBS?

vmod3 <- lm(ln_rec ~ ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
              # Spr_ST_SEBS_scaled +
              Smr_temp_250m_GOA_scaled  + 
               Spr_chlA_biom_GOA_scaled +  
              #Spr_chlA_biom_SEBS_scaled    +  
              Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
              # ann_Copepod_size_EGOA_scaled  +  
              ann_Copepod_size_WGOA_scaled   +    
              # Smr_euph_abun_Kod_scaled  +  
              YOY_grwth_Middleton_scaled  +     
              Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
              #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
              Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
              # sablefish_bycatch_arrowtooth_fishery_scaled + 
              smr_adult_cond_scaled , data=scaled_only)


vifs3 <- car::vif(vmod3) #higher than when GOA chlA was dropped

#if we follow up dropping GOA chlA next would be Spr_chlA_peak_GOA_scaled

vmod4 <- lm(ln_rec ~ ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
              # Spr_ST_SEBS_scaled +
              Smr_temp_250m_GOA_scaled  + 
              # Spr_chlA_biom_GOA_scaled +  
              Spr_chlA_biom_SEBS_scaled    +  
              #Spr_chlA_peak_GOA_scaled  + 
              Spr_chlA_peak_SEBS_scaled    +  
              # ann_Copepod_size_EGOA_scaled  +  
              ann_Copepod_size_WGOA_scaled   +    
              # Smr_euph_abun_Kod_scaled  +  
              YOY_grwth_Middleton_scaled  +     
              Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
              #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
              Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
              # sablefish_bycatch_arrowtooth_fishery_scaled + 
              smr_adult_cond_scaled , data=scaled_only)


vifs4 <- car::vif(vmod4) #now they're all good

#but what if we dropped the SEBS chlA-biom indicator instead since juvs should mostly be in GOA?
#next to drop would have been Spr_chlA_biom_GOA_scaled

vmod5 <- lm(ln_rec ~ ann_heatwave_GOA_scaled  + #Spr_ST_GOA_scaled +                
              # Spr_ST_SEBS_scaled +
              Smr_temp_250m_GOA_scaled  + 
               Spr_chlA_biom_GOA_scaled +  
              #Spr_chlA_biom_SEBS_scaled    +  
             # Spr_chlA_peak_GOA_scaled  + 
              Spr_chlA_peak_SEBS_scaled    +  
              # ann_Copepod_size_EGOA_scaled  +  
              ann_Copepod_size_WGOA_scaled   +    
              # Smr_euph_abun_Kod_scaled  +  
              YOY_grwth_Middleton_scaled  +     
              Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
              #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
              Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
              # sablefish_bycatch_arrowtooth_fishery_scaled + 
              smr_adult_cond_scaled , data=scaled_only)


vifs5 <- car::vif(vmod5) #Spr_chlA_biom_GOA_scaled still bad



#REPEAT with SEBS SST instead of heatwave since we know the latter might be an issue


vsebsmod <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + 
                 Spr_ST_GOA_scaled +                
             # Spr_ST_SEBS_scaled +
             Smr_temp_250m_GOA_scaled  + 
             Spr_chlA_biom_GOA_scaled +  Spr_chlA_biom_SEBS_scaled    +  
             Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
             # ann_Copepod_size_EGOA_scaled  +  
             ann_Copepod_size_WGOA_scaled   +    
             # Smr_euph_abun_Kod_scaled  +  
             YOY_grwth_Middleton_scaled  +     
             Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
             #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
             Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
             # sablefish_bycatch_arrowtooth_fishery_scaled + 
             smr_adult_cond_scaled , data=scaled_only)


vifsebs <- car::vif(vsebsmod) #ALL BAD?!

#drop worst which is Spr_chlA_biom_GOA_scaled at 21


vsebsmod2 <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + 
                 Spr_ST_GOA_scaled +                
                 # Spr_ST_SEBS_scaled +
                 Smr_temp_250m_GOA_scaled  + 
                 #Spr_chlA_biom_GOA_scaled +  
                  Spr_chlA_biom_SEBS_scaled    +  
                 Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
                 # ann_Copepod_size_EGOA_scaled  +  
                 ann_Copepod_size_WGOA_scaled   +    
                 # Smr_euph_abun_Kod_scaled  +  
                 YOY_grwth_Middleton_scaled  +     
                 Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
                 #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
                 Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
                 # sablefish_bycatch_arrowtooth_fishery_scaled + 
                 smr_adult_cond_scaled , data=scaled_only)


vifsebs2 <- car::vif(vsebsmod2)

#drop next worst, Spr_chlA_peak_GOA_scaled

vsebsmod3 <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + 
                  Spr_ST_GOA_scaled +                
                  # Spr_ST_SEBS_scaled +
                  Smr_temp_250m_GOA_scaled  + 
                  #Spr_chlA_biom_GOA_scaled +  
                  Spr_chlA_biom_SEBS_scaled    +  
                 # Spr_chlA_peak_GOA_scaled  + 
                  Spr_chlA_peak_SEBS_scaled    +  
                  # ann_Copepod_size_EGOA_scaled  +  
                  ann_Copepod_size_WGOA_scaled   +    
                  # Smr_euph_abun_Kod_scaled  +  
                  YOY_grwth_Middleton_scaled  +     
                  Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
                  #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
                  Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
                  # sablefish_bycatch_arrowtooth_fishery_scaled + 
                  smr_adult_cond_scaled , data=scaled_only)


vifsebs3 <- car::vif(vsebsmod3) #close enough! Some round down to 5 from 5.3

#what if we removed the sebs chlrA instead
vsebsmod4 <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + 
                  Spr_ST_GOA_scaled +                
                  # Spr_ST_SEBS_scaled +
                  Smr_temp_250m_GOA_scaled  + 
                  Spr_chlA_biom_GOA_scaled +  
                  #Spr_chlA_biom_SEBS_scaled    +  
                  Spr_chlA_peak_GOA_scaled  + Spr_chlA_peak_SEBS_scaled    +  
                  # ann_Copepod_size_EGOA_scaled  +  
                  ann_Copepod_size_WGOA_scaled   +    
                  # Smr_euph_abun_Kod_scaled  +  
                  YOY_grwth_Middleton_scaled  +     
                  Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
                  #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
                  Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
                  # sablefish_bycatch_arrowtooth_fishery_scaled + 
                  smr_adult_cond_scaled , data=scaled_only)


vifsebs4 <- car::vif(vsebsmod4)

#next would be Spr_chlA_peak_GOA_scaled at 9.3

vsebsmod5 <- lm(ln_rec ~ #ann_heatwave_GOA_scaled  + 
                  Spr_ST_GOA_scaled +                
                  # Spr_ST_SEBS_scaled +
                  Smr_temp_250m_GOA_scaled  + 
                  Spr_chlA_biom_GOA_scaled +  
                  #Spr_chlA_biom_SEBS_scaled    +  
                  #Spr_chlA_peak_GOA_scaled  + 
                  Spr_chlA_peak_SEBS_scaled    +  
                  # ann_Copepod_size_EGOA_scaled  +  
                  ann_Copepod_size_WGOA_scaled   +    
                  # Smr_euph_abun_Kod_scaled  +  
                  YOY_grwth_Middleton_scaled  +     
                  Smr_CPUE_juv_ADFG_ln_scaled   + #Smr_CPUE_juv_GOA_ln_scaled   +      
                  #spawner_mean_age_scaled +  spawner_age_evenness_scaled +       
                  Smr_condition_fem_age4_GOA_scaled + #arrowtooth_biomass_scaled  +        
                  # sablefish_bycatch_arrowtooth_fishery_scaled + 
                  smr_adult_cond_scaled , data=scaled_only)


vifsebs5 <- car::vif(vsebsmod5) #also fine


#now automate the process





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
