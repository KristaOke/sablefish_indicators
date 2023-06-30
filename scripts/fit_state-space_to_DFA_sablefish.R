#============================================================================================================================================
# Running state space regression on DFA outputs

# Coded by Curry Cunningham June 2023
#============================================================================================================================================
#Notes: 
#============================================================================================================================================

#=============================================================
#### Define Directory Structure ####
wd <- getwd()

dir.data <- file.path(wd,"data")
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")

#Get data =======

trends_rec_dat <- read.csv(file = paste(wd,"/data/DFA_trends_recruit_data.csv", sep=""))

#this will load trends from two different DFAs

#model 1 has a single trend onto which climate (mostly temperature) and juvenile
#abundance metrics load

#model 2 has two trends, temperature loads onto trend 1, juveniles load onto trend 2

#ideally we would run separate state space models regressing the trends from
#each model against ln_recruitment
#e.g.
# ln_rec ~ model1_state1 
# AND
# ln_rec ~ model2_state1 + model2_state2


#=============================================================







