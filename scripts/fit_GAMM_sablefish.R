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

dat <- read.csv(file.path(dir.data, "sablefish_BAS_indicators_2022.csv"))

scaled_dat <- dat %>% #group_by(Year) %>%
  mutate(recruit_scaled=scale(Recruitment),
         Spr_ST_SEBS_scaled=scale(Spring_Temperature_Surface_SEBS_Satellite),
         YOY_grwth_Middleton_scaled=scale(Annual_Sablefish_Growth_YOY_Middleton_Survey),
         Smr_CPUE_juv_ADFG_scaled=scale(Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey),
         spawner_mean_age_scaled=scale(Annual_Sablefish_Mean_Age_Female_Adult_Model),
         spawner_age_evenness_scaled=scale(Annual_Sablefish_Age_Evenness_Female_Adult_Model),
         arrowtooth_biomass_scaled=scale(Annual_Arrowtooth_Biomass_GOA_Model),
         sablefish_bycatch_arrowtooth_fishery_scaled=scale(Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery),
         smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Female_Adult_GOA_Survey))

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
#look into possible even-odd pattern?

p1 <- ggplot(scaled_dat, aes(Year, recruit_scaled)) + geom_point() +
  geom_line()
p2 <- ggplot(scaled_dat, aes(Year, Spr_ST_SEBS_scaled)) + geom_point() +
  geom_line()

plot_grid(p1, p2, ncol=1)

p3 <- ggplot(scaled_dat, aes(Year, arrowtooth_biomass_scaled)) + geom_point() +
  geom_smooth()

p4 <- ggplot(scaled_dat, aes(Year, sablefish_bycatch_arrowtooth_fishery_scaled)) + geom_point() +
  geom_smooth()

p5 <- ggplot(scaled_dat, aes(Year, Smr_CPUE_juv_ADFG_scaled)) + geom_point() +
  geom_smooth()

plot_grid(p1, p2, p3, p4, p5, ncol=1)

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
              Smr_CPUE_juv_ADFG_scaled +
              s(spawner_age_evenness_scaled, k=4) +
              smr_adult_cond_scaled,
            data=scaled_dat)
gam.check(mod2) #NOT GOOD
summary(mod2)
plot(mod2)

#what about linear?

modL <- glm(recruit_scaled ~ Spr_ST_SEBS_scaled +
      YOY_grwth_Middleton_scaled +
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
              s(YOY_grwth_Middleton_scaled, k=4) +
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
              s(YOY_grwth_Middleton_scaled, k=4) +
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
              YOY_grwth_Middleton_scaled +
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
              YOY_grwth_Middleton_scaled +
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


