#============================================================================================================================================
# DFA on training and testing datasets

#Created by Krista, Mar 2023
#copied to edit from older scripts
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

#DATA CONTROL SECTION----

#select the z-scored columns, no recruitment

scaled_dfa_dat <- scaled[,c(1,24:43)]

#remove problem covariates=====

#DFA can handle missing data but there is still too much missing in euphasiid data
#remove the mean age and evenness indicators
#also remove arrowtooth biomass b/c it's not converging even in best model
#same reason dropping sebs ST


scaled_dfa_dat <- scaled_dfa_dat[,!names(scaled_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled",
                                                                "Smr_condition_fem_age4_GOA_scaled")]



#manupulate data for DFA========

#check to make sure years are in order!

scaled_dfa_dat$Year==scaled_dfa_dat$Year[order(scaled_dfa_dat$Year)] #SHOULD BE ALL TRUE

#then drop Year column

#scaled_dfa_dat <- scaled_dfa_dat[,!names(scaled_dfa_dat) %in% c("Year")]


z.mat1 <- t(as.matrix(scaled_dfa_dat))
colnames(z.mat1) <- z.mat1[1,]
z.mat1 <- z.mat1[-1,]


#fit model========

# now fit DFA models with 1-3 trends and different error structures and compare

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=100000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")

#Run each and save output separately

#training set 1 DFA-----
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1:6) {  # allowing up to 3 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(z.mat1, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=FALSE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  conv=kemz$convergence,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

model.data1 <- model.data

# calculate delta-AICc scores, sort in descending order, and compare
model.data1$dAICc <- model.data1$AICc-min(model.data1$AICc)
model.data1 <- model.data1 %>%
  arrange(dAICc)
model.data1

#write_csv(model.data1, file = paste(wd,"/scripts/model_comparison_DFA_wholedataset_zscor-fixed.csv", sep=""))
write_csv(model.data1, file = paste(wd,"/scripts/model_comparison_DFA_wholedataset_zscor-fixed_conv.csv", sep=""))
model.data1 <- read_csv(file = paste(wd,"/scripts/model_comparison_DFA_wholedataset_zscor-fixed.csv", sep=""))


#---

#try different algorithm see if convergence happens faster

cntl.listB = list(minit=200, maxit=100000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             #"equalvarcov",
             "unconstrained")

#Run each and save output separately

#training set 1 DFA-----
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1:6) {  # allowing up to 3 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(z.mat1, model=dfa.model, #control=cntl.listB,
                 form="dfa", z.score=FALSE, method="BFGS")
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  conv=kemz$convergence,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

model.dataB <- model.data

# calculate delta-AICc scores, sort in descending order, and compare
model.dataB$dAICc <- model.dataB$AICc-min(model.dataB$AICc)
model.dataB <- model.dataB %>%
  arrange(dAICc)
model.dataB

#write_csv(model.dataB, file = paste(wd,"/scripts/model_comparison_DFA_wholedataset_BFGS.csv", sep=""))
write_csv(model.dataB, file = paste(wd,"/scripts/model_comparison_DFA_wholedataset_BFGS_conv.csv", sep=""))
model.dataB <- read_csv(file = paste(wd,"/scripts/model_comparison_DFA_wholedataset_BFGS.csv", sep=""))




#rotate=======

# now fit best model

model.list.1 = list(A="zero", m=1, R="diagonal and unequal") # best model by a little
cntl.list1 = list(minit=200, maxit=40000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.1 = MARSS(z.mat1, model=model.list.1, z.score=FALSE, form="dfa", control=cntl.list1)
#
autoplot(model.1)

#extract states and SE

model.1_trends <- as.data.frame(colnames(z.mat1))
colnames(model.1_trends) <- "Year"

model.1_trends$states <- as.vector(model.1$states)
model.1_trends$states_se <- as.vector(model.1$states.se)



# and rotate the loadings
Z.est = coef(model.1, type="matrix")$Z
H_inv = varimax(coef(model.1, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.1$states #

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

#trying for now IF 1 TREND
Z.rot <- Z.est
proc_rot =  model.1$states 

#Z.rot$names <- rownames(z.mat1)
#Z.rot <- arrange(Z.rot, V1)
#Z.rot <- gather(Z.rot[,c(1,2)])
#Z.rot$names <- rownames(z.mat1)
#Z.rot$plot.names <- reorder(Z.rot$names, 1:14)



#plot=========
# get CI and plot loadings...
modCI <- MARSSparamCIs(model.1)

plot.CI <- data.frame(mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                      lowCI=modCI$par.lowCI$Z)

plot.CI <- arrange(plot.CI, mean)
#plot.CI$names.order <- reorder(plot.CI$names, plot.CI$mean)
dodge <- position_dodge(width=0.9)

rec.plot <- ggplot(Z.rot, aes(names, value, fill=key)) + geom_bar(stat="identity", position="dodge") #+
# theme_bw() + ylab("Loading") + xlab("") + 
# scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
# theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
# theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) + ylim(-0.6, 0.8)

#based on nwfsc-timeseries.github.io

yr_frst <- 1977

## get number of time series
N_ts <- dim(z.mat1)[1]
## get length of time series
TT <- dim(z.mat1)[2]

## get the estimated ZZ
Z_est <- coef(model.1, type = "matrix")$Z
## get the inverse of the rotation matrix
#H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
#Z_rot = Z_est %*% H_inv
#IF ONE TREND
Z_rot <- Z_est

## rotate processes
#proc_rot = solve(H_inv) %*% model.1$states

#if 1 trend don't rotate?????
proc_rot =  model.1$states

mm <- 1 #4 processes

rec_names <- rownames(z.mat1)
ylbl <- rec_names
w_ts <- seq(dim(z.mat1)[2])
layout(matrix(c(1, 2, 3, 4), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
# jpeg("figs/ugly_DFA_trends_loadings.jpg")
par(mfcol=c(mm,2), mar = c(1.3,1.3,1.3,1.3), omi = c(0.1, 0.1, 0.1, 0.1))
## plot the processes
i<-1
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", #ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  #axis(1, 12 * (0:dim(all.clim.dat)[2]) + 1, yr_frst + 0:dim(all.clim.dat)[2])
  axis(1, 1:47, yr_frst + 0:dim(z.mat1)[2])
}
## plot the loadings
clr <- c("brown", 
         "blue", 
         "darkgreen", 
         "darkred", 
         "purple", 
         "darkorange")
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col=clr)
  for (j in 1:N_ts) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, col=clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, col=clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}
#dev.off()

par(mai = c(0.9, 0.9, 0.1, 0.1))
ccf(proc_rot[1, ], proc_rot[2, ], lag.max = 12, main = "")


#plot fits
# get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
#   ## empty list for results
#   fits <- list()
#   ## extra stuff for var() calcs
#   Ey <- MARSS:::MARSShatyt(MLEobj)
#   ## model params
#   ZZ <- coef(MLEobj, type = "matrix")$Z
#   ## number of obs ts
#   nn <- dim(Ey$ytT)[1]
#   ## number of time steps
#   TT <- dim(Ey$ytT)[2]
#   ## get the inverse of the rotation matrix
#   H_inv <- varimax(ZZ)$rotmat
#   ## check for covars
#   if (!is.null(dd)) {
#     DD <- coef(MLEobj, type = "matrix")$D
#     ## model expectation
#     fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
#   } else {
#     ## model expectation
#     fits$ex <- ZZ %*% H_inv %*% MLEobj$states
#   }
#   ## Var in model fits
#   VtT <- MARSSkfss(MLEobj)$VtT
#   VV <- NULL
#   for (tt in 1:TT) {
#     RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
#                                                          , tt] %*% t(ZZ)
#     SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
#       t(MLEobj$states[, tt, drop = FALSE])
#     VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
#   }
#   SE <- sqrt(VV)
#   ## upper & lower (1-alpha)% CI
#   fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
#   fits$lo <- qnorm(alpha/2) * SE + fits$ex
#   return(fits)
# }
# 
# ## get model fits & CI's
# mod_fit <- get_DFA_fits(model.1) #the function above not set up to work for one trend!
# ## plot the fits
# ylbl <- row.names(z.mat1)
# # par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
# #                                                              0, 0, 0))
# for (i in 1:N_ts) {
#   up <- mod_fit$up[i, ]
#   mn <- mod_fit$ex[i, ]
#   lo <- mod_fit$lo[i, ]
#   plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
#        cex.lab = 1.2, ylim = c(min(lo), max(up)))
#   axis(1,  (0:dim(z.mat1)[2]) + 1, 
#        yr_frst + 0:dim(z.mat1)[2])
#   points(w_ts, z.mat1[i, ], pch = 16)
#   lines(w_ts, up, col = "darkgray")
#   lines(w_ts, mn, col = "black", lwd = 2)
#   lines(w_ts, lo, col = "darkgray")
# }
# 
# 



#plot second best model=========


# now fit second best model

model.list.2 = list(A="zero", m=2, R="diagonal and unequal") # second best model
model.2 = MARSS(z.mat1, model=model.list.2, method="BFGS", z.score=FALSE, form="dfa")#, control=cntl.list2)
#
cntl.list2 = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001, method="BFGS" )

autoplot(model.2) #says check MARSSresiduals.tt1 warning message
MARSSresiduals(model.2, type="tt1") #seems like error b/c of NAs creating noninvertable matrix in residuals

#extract states and SE

model.2_trends <- as.data.frame(colnames(z.mat1))
colnames(model.2_trends) <- "Year"

model.2_trends$state1 <- as.vector(model.2$states[1,])
model.2_trends$state2 <- as.vector(model.2$states[2,])
model.2_trends$state1_se <- as.vector(model.2$states.se[1,])
model.2_trends$state2_se <- as.vector(model.2$states.se[2,])




# and rotate the loadings
Z.est = coef(model.2, type="matrix")$Z
H_inv = varimax(coef(model.2, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.2$states #doesn't work

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(z.mat1)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(z.mat1)
#Z.rot$plot.names <- reorder(Z.rot$names, 1:14)



# get CI and plot loadings...
modCI <- MARSSparamCIs(model.2)

plot.CI <- data.frame(mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                      lowCI=modCI$par.lowCI$Z)

plot.CI <- arrange(plot.CI, mean)
#plot.CI$names.order <- reorder(plot.CI$names, plot.CI$mean)
dodge <- position_dodge(width=0.9)

rec.plot <- ggplot(Z.rot, aes(names, value, fill=key)) + geom_bar(stat="identity", position="dodge") #+
# theme_bw() + ylab("Loading") + xlab("") + 
# scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
# theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
# theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) + ylim(-0.6, 0.8)

#based on nwfsc-timeseries.github.io

yr_frst <- 1977

## get number of time series
N_ts <- dim(z.mat1)[1]
## get length of time series
TT <- dim(z.mat1)[2]

## get the estimated ZZ
Z_est <- coef(model.2, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% model.2$states

mm <- 2 #processes

rec_names <- rownames(z.mat1)
ylbl <- rec_names
w_ts <- seq(dim(z.mat1)[2])
layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
# jpeg("figs/ugly_DFA_trends_loadings.jpg")
par(mfcol=c(mm,2), mar = c(1,1,1,1), omi = c(0, 0, 0, 0))
## plot the processes
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", #ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  #axis(1, 12 * (0:dim(all.clim.dat)[2]) + 1, yr_frst + 0:dim(all.clim.dat)[2])
  axis(1, 1:47, yr_frst + 0:dim(z.mat1)[2])
}
## plot the loadings
clr <- c("brown", 
         "blue", 
         "darkgreen", 
         "darkred", 
         "purple", 
         "darkorange")
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col=clr)
  for (j in 1:N_ts) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, col=clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, col=clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}
#dev.off()

par(mai = c(0.9, 0.9, 0.1, 0.1))
ccf(proc_rot[1, ], proc_rot[2, ], lag.max = 12, main = "")

#plot fits
get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
  ## empty list for results
  fits <- list()
  ## extra stuff for var() calcs
  Ey <- MARSS:::MARSShatyt(MLEobj)
  ## model params
  ZZ <- coef(MLEobj, type = "matrix")$Z
  ## number of obs ts
  nn <- dim(Ey$ytT)[1]
  ## number of time steps
  TT <- dim(Ey$ytT)[2]
  ## get the inverse of the rotation matrix
  H_inv <- varimax(ZZ)$rotmat
  ## check for covars
  if (!is.null(dd)) {
    DD <- coef(MLEobj, type = "matrix")$D
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
  } else {
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states
  }
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for (tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                         , tt] %*% t(ZZ)
    SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
      t(MLEobj$states[, tt, drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
  fits$lo <- qnorm(alpha/2) * SE + fits$ex
  return(fits)
}

## get model fits & CI's
mod_fit <- get_DFA_fits(model.2)
## plot the fits
ylbl <- row.names(z.mat1)
# par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
#                                                              0, 0, 0))
for (i in 1:N_ts) {
  up <- mod_fit$up[i, ]
  mn <- mod_fit$ex[i, ]
  lo <- mod_fit$lo[i, ]
  plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
       cex.lab = 1.2, ylim = c(min(lo), max(up)))
  axis(1,  (0:dim(z.mat1)[2]) + 1, 
       yr_frst + 0:dim(z.mat1)[2])
  points(w_ts, z.mat1[i, ], pch = 16)
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}


#plot third best model=========


# now fit best model

model.list.3 = list(A="zero", m=1, R="diagonal and unequal") # 
model.3 = MARSS(z.mat1, model=model.list.3, z.score=FALSE, form="dfa", method="BFGS")
#
cntl.list3 = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

autoplot(model.3)

# and rotate the loadings
Z.est = coef(model.3, type="matrix")$Z
H_inv = varimax(coef(model.3, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.3$states #doesn't work

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(z.ind.mat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(z.ind.mat)
#Z.rot$plot.names <- reorder(Z.rot$names, 1:14)



# get CI and plot loadings...
modCI <- MARSSparamCIs(model.3)

plot.CI <- data.frame(mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                      lowCI=modCI$par.lowCI$Z)

plot.CI <- arrange(plot.CI, mean)
#plot.CI$names.order <- reorder(plot.CI$names, plot.CI$mean)
dodge <- position_dodge(width=0.9)

rec.plot <- ggplot(Z.rot, aes(names, value, fill=key)) + geom_bar(stat="identity", position="dodge") #+
# theme_bw() + ylab("Loading") + xlab("") + 
# scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
# theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
# theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) + ylim(-0.6, 0.8)

#based on nwfsc-timeseries.github.io

yr_frst <- 1977

## get number of time series
N_ts <- dim(z.ind.mat)[1]
## get length of time series
TT <- dim(z.ind.mat)[2]

## get the estimated ZZ
Z_est <- coef(model.3, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% model.3$states

mm <- 1 #processes

rec_names <- rownames(z.ind.mat)
ylbl <- rec_names
w_ts <- seq(dim(z.ind.mat)[2])
layout(matrix(c(1, 2), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
# jpeg("figs/ugly_DFA_trends_loadings.jpg")
par(mfcol=c(mm,2), mar = c(1,1,1,1), omi = c(0, 0, 0, 0))
## plot the processes
i<-1
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", #ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  #axis(1, 12 * (0:dim(all.clim.dat)[2]) + 1, yr_frst + 0:dim(all.clim.dat)[2])
  axis(1, 1:47, yr_frst + 0:dim(z.ind.mat)[2])
}
## plot the loadings
clr <- c("brown", 
         "blue", 
         "darkgreen", 
         "darkred", 
         "purple", 
         "darkorange")
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col=clr)
  for (j in 1:N_ts) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, col=clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, col=clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}
#dev.off()

par(mai = c(0.9, 0.9, 0.1, 0.1))
ccf(proc_rot[1, ], proc_rot[2, ], lag.max = 12, main = "")



#plot 3 trend model=========


# fitting first model with 3 trends (7th best)

model.list.4 = list(A="zero", m=3, R="diagonal and unequal") # 
model.4 = MARSS(z.mat1, model=model.list.4, z.score=TRUE, form="dfa", method="BFGS")
#SEVERAL PARAMETERS DO NOT CONVERGE
cntl.list4 = list(minit=200, maxit=300000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

autoplot(model.4)

# and rotate the loadings
Z.est = coef(model.4, type="matrix")$Z
H_inv = varimax(coef(model.4, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.4$states #doesn't work

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(z.mat1)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(z.mat1)
#Z.rot$plot.names <- reorder(Z.rot$names, 1:14)



# get CI and plot loadings...
modCI <- MARSSparamCIs(model.4)

plot.CI <- data.frame(mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                      lowCI=modCI$par.lowCI$Z)

plot.CI <- arrange(plot.CI, mean)
#plot.CI$names.order <- reorder(plot.CI$names, plot.CI$mean)
dodge <- position_dodge(width=0.9)

rec.plot <- ggplot(Z.rot, aes(names, value, fill=key)) + geom_bar(stat="identity", position="dodge") #+
# theme_bw() + ylab("Loading") + xlab("") + 
# scale_fill_manual(values=c("Trend 1" = cb[2], "Trend 2" = cb[3])) +
# theme(legend.position = c(0.8,0.2), legend.title=element_blank()) + geom_hline(yintercept = 0) +
# theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) + ylim(-0.6, 0.8)

#based on nwfsc-timeseries.github.io

yr_frst <- 1977

## get number of time series
N_ts <- dim(z.mat1)[1]
## get length of time series
TT <- dim(z.mat1)[2]

## get the estimated ZZ
Z_est <- coef(model.4, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% model.4$states

mm <- 3 #processes

rec_names <- rownames(z.mat1)
ylbl <- rec_names
w_ts <- seq(dim(z.mat1)[2])
layout(matrix(c(1, 2), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
# jpeg("figs/ugly_DFA_trends_loadings.jpg")
par(mfcol=c(mm,2), mar = c(1,1,1,1), omi = c(0, 0, 0, 0))
## plot the processes
i<-1
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", #ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  #axis(1, 12 * (0:dim(all.clim.dat)[2]) + 1, yr_frst + 0:dim(all.clim.dat)[2])
  axis(1, 1:47, yr_frst + 0:dim(z.mat1)[2])
}
## plot the loadings
clr <- c("brown", 
         "blue", 
         "darkgreen", 
         "darkred", 
         "purple", 
         "darkorange")
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col=clr)
  for (j in 1:N_ts) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, col=clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, col=clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}
#dev.off()

par(mai = c(0.9, 0.9, 0.1, 0.1))
ccf(proc_rot[1, ], proc_rot[2, ], lag.max = 12, main = "")




#prep trends for regression-------
#join trends data to data with recruitment

#model.1
model.1_trends$Year <- as.numeric(model.1_trends$Year)
model1_reg_dat <- left_join(model.1_trends, scaled)

quickmod <- lm(ln_rec ~ states, data=model1_reg_dat)
anova(quickmod)
ggplot(model1_reg_dat, aes(states, ln_rec)) + geom_point()+ geom_smooth()

#model.2
model.2_trends$Year <- as.numeric(model.2_trends$Year)
model2_reg_dat <- left_join(model.2_trends, scaled)

quickmod2 <- lm(ln_rec ~ state1*state2, data=model2_reg_dat)
anova(quickmod2)
ggplot(model2_reg_dat, aes(state1, ln_rec)) + geom_point()+ geom_smooth()+
  geom_errorbar(aes(xmin=state1-state1_se, xmax=state1+state1_se))
ggplot(model2_reg_dat, aes(state2, ln_rec)) + geom_point() + geom_smooth()+
  geom_errorbar(aes(xmin=state2-state2_se, xmax=state2+state2_se))

quickmod2.2 <- lm(ln_rec ~ state2, data=model2_reg_dat)
anova(quickmod2.2)
summary(quickmod2.2)

#can I join to a dataset with non ln rec?
ggplot(model2_reg_dat, aes(state2, Recruitment)) + geom_point() + geom_smooth()

#compared to just cpue adfg?
quickmod2_adfg <- lm(ln_rec ~ Smr_CPUE_juv_ADFG_ln_scaled, data=model2_reg_dat)
anova(quickmod2_adfg)

ggplot(model2_reg_dat, aes(Smr_CPUE_juv_ADFG_ln_scaled, ln_rec)) + geom_point() + geom_smooth()

ggplot(model2_reg_dat, aes(sablefish_bycatch_arrowtooth_fishery_scaled, ln_rec)) + geom_point() + geom_smooth()


#join trends from diff models and save for state space====


names(model.1_trends)
names(model.2_trends)

join1 <- model.1_trends
join2 <- model.2_trends

colnames(join1) <- c("Year", "model1_state1", "model1_state1_SE")
colnames(join2) <- c("Year", "model2_state1", "model2_state2", 
                     "model2_state1_SE", "model2_state2_SE")

trends_join <- left_join(join1, join2)

trends_join$Year <- as.integer(trends_join$Year)

trends_rec_dat <- left_join(trends_join, scaled[,names(scaled) %in% c("ln_rec", "Year")])

write_csv(trends_rec_dat, file = paste(wd,"/data/DFA_trends_recruit_data.csv", sep=""))

