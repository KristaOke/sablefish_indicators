#============================================================================================================================================
# let's try SEM!

#Created by Krista, May 2023
#
#============================================================================================================================================
#Notes:
#============================================================================================================================================
library(ggplot2)
library(reshape2)
library(dplyr)
library(pracma)
library(MARSS)
library(tidyr)
# devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(cowplot)
library(gridExtra)
library(broom)
library(lemon)
library(MuMIn)
library(lmtest)
library(cowplot)
library(AKesp)
library(corrplot)

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#Get data=======

# load data

scaled <- read.csv(file=paste(wd,"/data/whole_dataset_scaled.csv", sep=""), row.names = 1)

#select the z-scored columns, no recruitment
scaled_sem_dat <- scaled[,c(1,24:43)]

scaled_sem_dat <- scaled_sem_dat[,!names(scaled_sem_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled",
                                                                "Smr_condition_fem_age4_GOA_scaled")]

#check to make sure years are in order!
scaled_sem_dat$Year==scaled_sem_dat$Year[order(scaled_sem_dat$Year)] #SHOULD BE ALL TRUE



#just want to look at bycatch over arrowtooth

ggplot(scaled, aes(Year, sablefish_bycatch_arrowtooth_fishery_scaled/arrowtooth_biomass_scaled)) +
  geom_point()

ggplot(scaled, aes(Year, arrowtooth_biomass_scaled)) +
  geom_point()

ggplot(scaled, aes(Year, sablefish_bycatch_arrowtooth_fishery_scaled)) +
  geom_point()

ggplot(scaled, aes(Year, ln_Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery/Annual_Arrowtooth_Biomass_GOA_Model)) +
  geom_point()

ggplot(scaled, aes(Year, Annual_Arrowtooth_Biomass_GOA_Model)) +
  geom_point()

ggplot(scaled, aes(Year, ln_Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery)) +
  geom_point()

ggplot(scaled, aes(Year, ln_rec)) +
  geom_point()














