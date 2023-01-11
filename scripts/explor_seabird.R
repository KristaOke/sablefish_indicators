#=============================================================================
#Explore USGS seabird data provided by Yumi

#=============================================================================
#Notes: Data from:
#Arimitsu, M.L., Hatch, S., 2022, Age-0 Sablefish size and growth indices from 
#seabird diets at Middleton Island, Alaska: U.S. Geological Survey data release, 
#https://doi.org/10.5066/P94KVH9X</othercit>
#Metadata, citation info in xml file in data folder
#=============================================================================

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
dir.figs <- file.path(wd,"figs")

#=============================================================
#### get data ####

birddat <- read.csv(file.path(dir.data, "MDO_age0sablefish_indices.csv"))

zbirddat <- birddat %>% #group_by(Year) %>%
  mutate(propmass_zscor=scale(propmass),
         CPUE_zscor=scale(CPUE),
         FreqOcr_zscor=scale(FO),
         growth_index_zscor=scale(growth_index),
         growth_anom_zscor=scale(growth_anomaly),
         pred_len_zscor=scale(pred_len))




















