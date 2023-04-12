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

#GAM cannot handle missing data, euphasiid and the every other yr GOA indicator were already removed


train1_gam_dat <- train1[,names(train1) %in% noncor_covars3]
train2_gam_dat <- train2[,names(train2) %in% noncor_covars3]
train3_gam_dat <- train3[,names(train3) %in% noncor_covars3]
train4_gam_dat <- train4[,names(train4) %in% noncor_covars3]
train5_gam_dat <- train5[,names(train5) %in% noncor_covars3]



#new gams-----------

mod_train1_1 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
#Smr_temp_250m_GOA_scaled  +      
#Spr_chlA_biom_SEBS_scaled  +     
#Spr_chlA_peak_GOA_scaled  +     
#Spr_chlA_peak_SEBS_scaled  +  
#ann_Copepod_size_WGOA_scaled +     
s(ann_Copepod_size_EGOA_scaled, k=4)  +     
s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
s(Smr_condition_fem_age4_GOA_scaled, k=4) +
s(smr_adult_cond_scaled, k=4), data=train1_gam_dat)
gam.check(mod_train1_1) #not good
summary(mod_train1_1)

#drop lease nonlinear and try again
mod_train1_2 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
                      #Smr_temp_250m_GOA_scaled  +      
                      #Spr_chlA_biom_SEBS_scaled  +     
                      #Spr_chlA_peak_GOA_scaled  +     
                      #Spr_chlA_peak_SEBS_scaled  +  
                      #ann_Copepod_size_WGOA_scaled +     
                      s(ann_Copepod_size_EGOA_scaled, k=4)  +     
                      s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
                      s(Smr_condition_fem_age4_GOA_scaled, k=4) +
                      smr_adult_cond_scaled, data=train1_gam_dat)
gam.check(mod_train1_2) #not good
summary(mod_train1_2)

mod_train1_3 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
                      #Smr_temp_250m_GOA_scaled  +      
                      #Spr_chlA_biom_SEBS_scaled  +     
                      #Spr_chlA_peak_GOA_scaled  +     
                      #Spr_chlA_peak_SEBS_scaled  +  
                      #ann_Copepod_size_WGOA_scaled +     
                      s(ann_Copepod_size_EGOA_scaled, k=4)  +     
                      s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
                      Smr_condition_fem_age4_GOA_scaled +
                      smr_adult_cond_scaled, data=train1_gam_dat)
gam.check(mod_train1_3) #not good
summary(mod_train1_3)

mod_train1_3 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                      #Smr_temp_250m_GOA_scaled  +      
                      #Spr_chlA_biom_SEBS_scaled  +     
                      #Spr_chlA_peak_GOA_scaled  +     
                      #Spr_chlA_peak_SEBS_scaled  +  
                      #ann_Copepod_size_WGOA_scaled +     
                      s(ann_Copepod_size_EGOA_scaled, k=4)  +     
                      s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
                      Smr_condition_fem_age4_GOA_scaled +
                      smr_adult_cond_scaled, data=train1_gam_dat)
gam.check(mod_train1_3) #not good
summary(mod_train1_3)


mod_train1_4 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                      #Smr_temp_250m_GOA_scaled  +      
                      #Spr_chlA_biom_SEBS_scaled  +     
                      #Spr_chlA_peak_GOA_scaled  +     
                      #Spr_chlA_peak_SEBS_scaled  +  
                      #ann_Copepod_size_WGOA_scaled +     
                      ann_Copepod_size_EGOA_scaled  +     
                      s(Smr_CPUE_juv_ADFG_ln_scaled, k=4 ) +    
                      Smr_condition_fem_age4_GOA_scaled +
                      smr_adult_cond_scaled, data=train1_gam_dat)
gam.check(mod_train1_4) #not good
summary(mod_train1_4)

#None are nonlinear!

#add back in the covars I removed, limit k since still too few data
mod_train1_5 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                    #  s(Smr_temp_250m_GOA_scaled, k=3)  +      
                      s(Spr_chlA_biom_SEBS_scaled, k=3)  +     
                      s(Spr_chlA_peak_GOA_scaled, k=3)  +     
                      s(Spr_chlA_peak_SEBS_scaled, k=3)  +  
                      s(ann_Copepod_size_WGOA_scaled, k=3) +     
                      ann_Copepod_size_EGOA_scaled  +     
                      Smr_CPUE_juv_ADFG_ln_scaled +    
                      Smr_condition_fem_age4_GOA_scaled +
                      smr_adult_cond_scaled, data=train1_gam_dat)
gam.check(mod_train1_5) #not good
summary(mod_train1_5)

#all linear? Try and then look at resids and only try covars that look nonlin as nonlinear
mod_trainlin <- lm(ln_rec ~ Spr_ST_SEBS_scaled  +      
                      #Smr_temp_250m_GOA_scaled +       
                     Spr_chlA_biom_GOA_scaled  +     
                      Spr_chlA_biom_SEBS_scaled  +     
                      Spr_chlA_peak_GOA_scaled  +     
                      Spr_chlA_peak_SEBS_scaled  +  
                      YOY_grwth_Middleton_scaled +     
                      ann_Copepod_size_EGOA_scaled  +     
                      Smr_CPUE_juv_ADFG_ln_scaled +    
                      smr_adult_cond_scaled, data=train1_gam_dat)
summary(mod_trainlin)

mod_trainlin2 <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
                     #Smr_temp_250m_GOA_scaled +       
                     Spr_chlA_biom_GOA_scaled  +     
                     Spr_chlA_biom_SEBS_scaled  +     
                     Spr_chlA_peak_GOA_scaled  +     
                     Spr_chlA_peak_SEBS_scaled  +  
                     YOY_grwth_Middleton_scaled +     
                     ann_Copepod_size_EGOA_scaled  +     
                     Smr_CPUE_juv_ADFG_ln_scaled +    
                     smr_adult_cond_scaled, data=train1_gam_dat)
summary(mod_trainlin2)


mod_trainlin3 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                       #Smr_temp_250m_GOA_scaled +       
                       Spr_chlA_biom_GOA_scaled  +     
                       Spr_chlA_biom_SEBS_scaled  +     
                       Spr_chlA_peak_GOA_scaled  +     
                       Spr_chlA_peak_SEBS_scaled  +  
                       YOY_grwth_Middleton_scaled +     
                       ann_Copepod_size_EGOA_scaled  +     
                       s(Smr_CPUE_juv_ADFG_ln_scaled, k=4) +    
                       smr_adult_cond_scaled, data=train1_gam_dat)
summary(mod_trainlin3)

#drop all satelitte covars, they are quite short
mod_trainshort <- gam(ln_rec ~ s(Spr_ST_SEBS_scaled, k=4)  +      
                       #Smr_temp_250m_GOA_scaled +       
                       #Spr_chlA_biom_GOA_scaled  +     
                       #Spr_chlA_biom_SEBS_scaled  +     
                       #Spr_chlA_peak_GOA_scaled  +     
                       #Spr_chlA_peak_SEBS_scaled  +  
                       YOY_grwth_Middleton_scaled +     
                       #ann_Copepod_size_EGOA_scaled  +     
                       s(Smr_CPUE_juv_ADFG_ln_scaled, k=4) +    
                       smr_adult_cond_scaled, data=train1_gam_dat)
summary(mod_trainshort)

mod_trainshort2 <- gam(ln_rec ~ Spr_ST_SEBS_scaled  +      
                        #Smr_temp_250m_GOA_scaled +       
                        #Spr_chlA_biom_GOA_scaled  +     
                        #Spr_chlA_biom_SEBS_scaled  +     
                        #Spr_chlA_peak_GOA_scaled  +     
                        #Spr_chlA_peak_SEBS_scaled  +  
                        YOY_grwth_Middleton_scaled +     
                        #ann_Copepod_size_EGOA_scaled  +     
                        s(Smr_CPUE_juv_ADFG_ln_scaled, k=4), data=train1_gam_dat)
summary(mod_trainshort2)
plot(mod_trainshort2)

library(gratia)
draw(mod_trainshort2) + theme_bw()

library(itsadug)
plot_smooth(mod_trainshort2, view="Smr_CPUE_juv_ADFG_ln_scaled", ylab="ln(Recruitment)",
            xlab="Summer juvenile CPUE ADFG survey")

#OLD===============================================================================================

#gamms==============


mod1 <- gam(recruit_scaled ~ s(Spr_ST_SEBS_scaled, k=4) +
               s(ann_Copepod_size_EGOA_scaled, k=4) +
               s(Smr_CPUE_juv_ADFG_scaled, k=4) +
               s(spawner_age_evenness_scaled, k=4) +
               s(smr_adult_cond_scaled, k=4),
               data=scaled_dat)
gam.check(mod1) #NOT GOOD
summary(mod1)
plot(mod1)

#reduce to only sig nonlinear

mod2 <- gam(recruit_scaled ~ s(Spr_ST_SEBS_scaled, k=4) +
              ann_Copepod_size_EGOA_scaled +
              Smr_CPUE_juv_ADFG_scaled +
              s(spawner_age_evenness_scaled, k=4) +
              smr_adult_cond_scaled,
            data=scaled_dat)
gam.check(mod2) #NOT GOOD
summary(mod2)
plot(mod2)

#what about linear?

modL <- glm(recruit_scaled ~ Spr_ST_SEBS_scaled +
      ann_Copepod_size_EGOA_scaled +
      Smr_CPUE_juv_ADFG_scaled +
      spawner_age_evenness_scaled +
      smr_adult_cond_scaled,
    data=scaled_dat)
summary(modL)

library(lsr)
etaSquared(modL, type=2, anova=TRUE)



#Repeat WITH ALL covars------
#even those that have high correlations
#danger!

#OK with all covars has more coefficients than data

all1 <- gam(recruit_scaled ~ s(Spr_ST_SEBS_scaled, k=4) +
              s(ann_Copepod_size_EGOA_scaled, k=4) +
              s(Smr_CPUE_juv_ADFG_scaled, k=4) +
              s(spawner_mean_age_scaled, k=4) +
              s(spawner_age_evenness_scaled, k=4) +
             # s(arrowtooth_biomass_scaled, k=4) +
              s(sablefish_bycatch_arrowtooth_fishery_scaled, k=4) +
              s(smr_adult_cond_scaled, k=4),
            data=scaled_dat)
gam.check(all1) #NOT GOOD
summary(all1)
plot(all1)

#reduce to only sig nonlinear

all2 <- gam(recruit_scaled ~ s(Spr_ST_SEBS_scaled, k=4) +
              s(ann_Copepod_size_EGOA_scaled, k=4) +
              Smr_CPUE_juv_ADFG_scaled +
              spawner_mean_age_scaled +
              s(spawner_age_evenness_scaled, k=4) +
              # s(arrowtooth_biomass_scaled, k=4) +
              s(sablefish_bycatch_arrowtooth_fishery_scaled, k=4) +
              smr_adult_cond_scaled,
            data=scaled_dat)
gam.check(all2) #NOT GOOD
summary(all2)
plot(all2)

#some no longer significantly nonlinear reduce again

all3 <- gam(recruit_scaled ~ Spr_ST_SEBS_scaled +
              ann_Copepod_size_EGOA_scaled +
              Smr_CPUE_juv_ADFG_scaled +
              spawner_mean_age_scaled +
              spawner_age_evenness_scaled +
              # s(arrowtooth_biomass_scaled, k=4) +
              s(sablefish_bycatch_arrowtooth_fishery_scaled, k=4) +
              smr_adult_cond_scaled,
            data=scaled_dat)
gam.check(all3) #less awful
summary(all3)
plot(all3)

#what about linear?

allL <- glm(recruit_scaled ~ Spr_ST_SEBS_scaled +
              ann_Copepod_size_EGOA_scaled +
              Smr_CPUE_juv_ADFG_scaled +
              #spawner_mean_age_scaled +    #removed b/c too few df
              spawner_age_evenness_scaled +
              sablefish_bycatch_arrowtooth_fishery_scaled +
              smr_adult_cond_scaled,
            data=scaled_dat)
summary(allL)


etaSquared(allL, type=2, anova=TRUE)

AIC(allL, all1, all2, all3)

#interactions???
int1 <- glm(recruit_scaled ~ Spr_ST_SEBS_scaled*YOY_grwth_Middleton_scaled +
              Spr_ST_SEBS_scaled*smr_adult_cond_scaled +
              Spr_ST_SEBS_scaled*Smr_CPUE_juv_ADFG_scaled +
              #spawner_mean_age_scaled +    #removed b/c too few df
              spawner_age_evenness_scaled +
              sablefish_bycatch_arrowtooth_fishery_scaled,
            data=scaled_dat)
summary(int1)


ggplot(scaled_dat, aes(Smr_CPUE_juv_ADFG_scaled, recruit_scaled,  col=Spr_ST_SEBS_scaled))+
  geom_point()

ggplot(scaled_dat, aes(Year, Smr_CPUE_juv_ADFG_scaled, col=Spr_ST_SEBS_scaled))+
  geom_point()


