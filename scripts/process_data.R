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

dat <- read.csv(file.path(dir.data, "sablefish_BAS_indicators_2022.csv"))

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


dat$ln_rec <- log(dat$Recruitment)

covar.noR <- dat %>% gather(key=type, value=value, -c(Year, Recruitment, ln_rec))


rec.plot <- ggplot(covar.noR, aes(y=Recruitment, x=value, fill=type)) +
  geom_point() + geom_smooth() +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
rec.plot

ln.rec.plot <- ggplot(covar.noR, aes(y=ln_rec, x=value, fill=type)) +
  geom_point() + geom_smooth() +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~type, scales='free') +
  theme(legend.position = "NA")
ln.rec.plot

#z-score data=======

scaled_dat <- dat %>% #group_by(Year) %>%
  mutate(recruit_scaled=scale(ln_rec),
         Spr_ST_SEBS_scaled=scale(Spring_Temperature_Surface_SEBS_Satellite),
         YOY_grwth_Middleton_scaled=scale(Annual_Sablefish_Growth_YOY_Middleton_Survey),
         Smr_CPUE_juv_ADFG_scaled=scale(Summer_Sablefish_CPUE_Juvenile_Nearshore_GOAAI_Survey),
         spawner_mean_age_scaled=scale(Annual_Sablefish_Mean_Age_Female_Adult_Model),
         spawner_age_evenness_scaled=scale(Annual_Sablefish_Age_Evenness_Female_Adult_Model),
         arrowtooth_biomass_scaled=scale(Annual_Arrowtooth_Biomass_GOA_Model),
         sablefish_bycatch_arrowtooth_fishery_scaled=scale(Annual_Sablefish_Incidental_Catch_Arrowtooth_Target_GOA_Fishery),
         smr_adult_cond_scaled=scale(Summer_Sablefish_Condition_Female_Adult_GOA_Survey))


#split data======

#NEED TO SUBSET out a training dataset
datlen <- length(scaled_dat$recruit_scaled)
train <- scaled_dat[sample(nrow(scaled_dat),(round(datlen*0.8))),]
train <- train[order(row.names(train)),]

testing <- anti_join(scaled_dat, train)

#save this training set?

#save analysis-ready data======























