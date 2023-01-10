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

dir.data <- file.path(wd,"data","Sablefish")
dir.output <- file.path(wd,"output","Sablefish")
dir.figs <- file.path(wd,"figs","Sablefish")

#=============================================================
#### get data ####

dat <- read.csv(file.path(dir.data, "sablefish_BAS_indicators_2021.csv"))

scaled_dat <- dat %>% #group_by(Year) %>%
  mutate(recruit_scaled=scale(Recruitment),
         Spr_ST_SEBS_scaled=scale(Spring_Surface_Temperature_SEBS),
         YOY_grwth_Middleton_scaled=scale(Sablefish_Growth_YOY_Middleton_Auklets),
         Smr_CPUE_juv_ADFG_scaled=scale(Summer_Sablefish_CPUE_Juvenile_WGOA_EAI_ADFG),
         spawner_mean_age_scaled=scale(Sablefish_Spawner_Mean_Age),
         spawner_age_evenness_scaled=scale(Sablefish_Spawner_Age_Evenness),
         arrowtooth_biomass_scaled=scale(Arrowtooth_Biomass_Assessment),
         sablefish_bycatch_arrowtooth_fishery_scaled=scale(Sablefish_Catch_Arrowtooth_Fishery),
         smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Adult_GOA_LL_Survey))

#look at covars
ggplot(scaled_dat, aes(Spr_ST_SEBS_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()

ggplot(scaled_dat, aes(YOY_grwth_Middleton_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()

ggplot(scaled_dat, aes(Smr_CPUE_juv_ADFG_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()

ggplot(scaled_dat, aes(spawner_mean_age_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()

ggplot(scaled_dat, aes(spawner_age_evenness_scaled, recruit_scaled)) + geom_point() +
  geom_smooth() #looks nonlinear

ggplot(scaled_dat, aes(arrowtooth_biomass_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()

ggplot(scaled_dat, aes(sablefish_bycatch_arrowtooth_fishery_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()

ggplot(scaled_dat, aes(smr_adult_cond_scaled, recruit_scaled)) + geom_point() +
  geom_smooth()


p1 <- ggplot(scaled_dat, aes(Year, recruit_scaled)) + geom_point() +
  geom_line()
p2 <- ggplot(scaled_dat, aes(Year, Spr_ST_SEBS_scaled)) + geom_point() +
  geom_line()

plot_grid(p1, p2, ncol=1)

#check cors==================

covar.mtx <- scaled_dat %>% select(-c(1:10))
corr.mtx <- cor(covar.mtx, use="na.or.complete")
corr.mtx
corrplot::corrplot.mixed(corr.mtx, tl.pos="lt", 
                         diag="n")
#drop arrowtooth biomass
#LOTS OF CORRELATIONS > 0.6
#drop bycatch in arrowtooth for now, highly corr w ADFG survey, which is more cor w recruitment?

covar.mtx2 <- scaled_dat %>% select(-c(1:11,17,18))
corr.mtx2 <- cor(covar.mtx2, use="na.or.complete")
corr.mtx2
corrplot::corrplot.mixed(corr.mtx2, tl.pos="lt", 
                         diag="n")
#drop mean age scaled too, could this be kind of too colinear w recruitment?
#and its corr w ADFG survey

#gamms==============


mod1 <- gam(recruit_scaled ~ s(Spr_ST_SEBS_scaled, k=4) +
               s(YOY_grwth_Middleton_scaled, k=4) +
               s(Smr_CPUE_juv_ADFG_scaled, k=4) +
               s(spawner_age_evenness_scaled, k=4) +
               s(smr_adult_cond_scaled, k=4),
               data=scaled_dat)
gam.check(mod1) #NOT GOOD
summary(mod1)
plot(mod1)

#reduce to only sig nonlinear

mod2 <- gam(recruit_scaled ~ s(Spr_ST_SEBS_scaled, k=4) +
              YOY_grwth_Middleton_scaled +
              s(Smr_CPUE_juv_ADFG_scaled, k=4) +
              spawner_age_evenness_scaled +
              s(smr_adult_cond_scaled, k=4),
            data=scaled_dat)
gam.check(mod2) #NOT GOOD
summary(mod2)
plot(mod2)
















