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


#Get data=======

# load data


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

#select the z-scored columns, no recruitment
train1_dfa_dat <- train1[,c(1,24:43)]
train2_dfa_dat <- train2[,c(1,24:43)]
train3_dfa_dat <- train3[,c(1,24:43)]
train4_dfa_dat <- train4[,c(1,24:43)]
train5_dfa_dat <- train5[,c(1,24:43)]


#remove problem covariates=====

#DFA can handle missing data but there is still too much missing in euphasiid data
#remove the mean age and evenness indicators
#also remove arrowtooth biomass b/c it's not converging even in best model

train1_dfa_dat <- train1_dfa_dat[,!names(train1_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                     "spawner_mean_age_scaled",
                                                     "spawner_age_evenness_scaled",
                                                     "arrowtooth_biomass_scaled")]
train2_dfa_dat <- train2_dfa_dat[,!names(train2_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled")]
train3_dfa_dat <- train3_dfa_dat[,!names(train3_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled")]
train4_dfa_dat <- train4_dfa_dat[,!names(train4_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled")]
train5_dfa_dat <- train5_dfa_dat[,!names(train5_dfa_dat) %in% c("Smr_euph_abun_Kod_scaled",
                                                                "spawner_mean_age_scaled",
                                                                "spawner_age_evenness_scaled",
                                                                "arrowtooth_biomass_scaled")]


#manupulate data for DFA========

#check to make sure years are in order!
train1_dfa_dat$Year==train1_dfa_dat$Year[order(train1_dfa_dat$Year)] #SHOULD BE ALL TRUE
train2_dfa_dat$Year==train2_dfa_dat$Year[order(train2_dfa_dat$Year)] #SHOULD BE ALL TRUE
train3_dfa_dat$Year==train3_dfa_dat$Year[order(train3_dfa_dat$Year)] #SHOULD BE ALL TRUE
train4_dfa_dat$Year==train4_dfa_dat$Year[order(train4_dfa_dat$Year)] #SHOULD BE ALL TRUE
train5_dfa_dat$Year==train5_dfa_dat$Year[order(train5_dfa_dat$Year)] #SHOULD BE ALL TRUE

#then drop Year column
train1_dfa_dat <- train1_dfa_dat[,!names(train1_dfa_dat) %in% c("Year")]
train2_dfa_dat <- train2_dfa_dat[,!names(train2_dfa_dat) %in% c("Year")]
train3_dfa_dat <- train3_dfa_dat[,!names(train3_dfa_dat) %in% c("Year")]
train4_dfa_dat <- train4_dfa_dat[,!names(train4_dfa_dat) %in% c("Year")]
train5_dfa_dat <- train5_dfa_dat[,!names(train5_dfa_dat) %in% c("Year")]


z.mat1 <- t(as.matrix(train1_dfa_dat))
colnames(z.mat1) <- z.mat1[1,]
z.mat1 <- z.mat1[-1,]


#fit model========

# now fit DFA models with 1-3 trends and different error structures and compare

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

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
  for(m in 1:4) {  # allowing up to 3 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(z.mat1, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
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

#---

#rotate=======

# now fit best model

model.list.1 = list(A="zero", m=1, R="diagonal and unequal") # best model by a little
cntl.list1 = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.1 = MARSS(z.mat1, model=model.list.1, z.score=TRUE, form="dfa", control=cntl.list1)
#not working yet, convergence issues


# and rotate the loadings
Z.est = coef(model.1, type="matrix")$Z
H_inv = varimax(coef(model.1, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.1$states #

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(z.mat1)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(z.mat1)
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
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% model.1$states

mm <- 2 #4 processes

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
  axis(1, 1:38, yr_frst + 0:dim(z.mat1)[2])
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







#plot second best model=========


# now fit second best model

model.list.2 = list(A="zero", m=3, R="diagonal and unequal") # second best model
model.2 = MARSS(z.ind.mat, model=model.list.2, z.score=TRUE, form="dfa", control=cntl.list2)
#DOES NOT CONVERGE bump up to 40K iter
cntl.list2 = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)


# and rotate the loadings
Z.est = coef(model.2, type="matrix")$Z
H_inv = varimax(coef(model.2, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.2$states #doesn't work

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(z.ind.mat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(z.ind.mat)
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
N_ts <- dim(z.ind.mat)[1]
## get length of time series
TT <- dim(z.ind.mat)[2]

## get the estimated ZZ
Z_est <- coef(model.2, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% model.2$states

mm <- 3 #processes

rec_names <- rownames(z.ind.mat)
ylbl <- rec_names
w_ts <- seq(dim(z.ind.mat)[2])
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





#plot third best model=========


# now fit best model

model.list.3 = list(A="zero", m=1, R="diagonal and equal") # 
model.3 = MARSS(z.ind.mat, model=model.list.3, z.score=TRUE, form="dfa", control=cntl.list)
#DOES NOT CONVERGE bump up to 60K iter
cntl.list3 = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)


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





#plot obs vs fitted========

#NOT UPDATED BELOW HERE JAN 2023!!!!!!!!!!!!!!--------------------------------------------

#from online course "nwfsc-timeseries.github.io"

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

#demean data
# y_bar <- apply(all.clim.dat, 1, mean, na.rm = TRUE)
# dat <- all.clim.dat - y_bar
# rownames(dat) <- rownames(all.clim.dat)

# dat <- scale(log.rec.mat)
# 
# head(log.rec.mat)
# 
# log.rec.mat.std <- log.rec.mat

# i <- 1
# for(i in 1:nrow(z.ind.mat)) {
#   z.ind.mat.std[i,] <- (z.ind.mat[i,]-mean(z.ind.mat[i,], na.rm=TRUE))/sd(z.ind.mat[i,], na.rm=TRUE)  
# }
#Double checking
# apply(z.ind.mat.std, 1, mean, na.rm=TRUE)
# apply(z.ind.mat.std, 1, sd, na.rm=TRUE)
dat <- z.ind.mat
#plot demeaned data

driv <- rownames(z.ind.mat)
clr <- c("brown", "blue", "darkgreen", "darkred", "purple",
         "brown", "blue", "darkgreen", "darkred", "purple",
         "brown", "blue", "darkgreen", "darkred", "purple",
         "brown", "blue", "darkgreen", "darkred")
cnt <- 1
# par(mfrow = c(N_ts, 2), mar = c(1, 1,1.5,1), omi = c(0.1, 
#                                                      0.1, 0.1, 0.1))
par(mfrow = c(N_ts, 1), mar = c(1, 1,1.5,1), omi = c(0.1, 
                                                     0.1, 0.1, 0.1))
for (i in driv) {
  plot(dat[i, ], xlab = "", ylab = "", bty = "L", 
       xaxt = "n", pch = 16, col = clr[cnt], type = "b")
  axis(1,  (0:dim(z.ind.mat)[2]) + 1, yr_frst + 0:dim(z.ind.mat)[2])
  title(i)
  cnt <- cnt + 1
}

#edit below here (and above!)

## get model fits & CI's
mod_fit <- get_DFA_fits(model.1)
## plot the fits
par(mfrow = c(5, 4), mar = c(1, 1, 1, 1), omi = c(0, 
                                                  0, 0, 0))
for (i in 1:N_ts) {
  up <- mod_fit$up[i, ]
  mn <- mod_fit$ex[i, ]
  lo <- mod_fit$lo[i, ]
  plot(w_ts, mn, xlab = "", xaxt = "n", type = "n", 
       cex.lab = 1.2, ylim = c(min(lo), max(up)))
  axis(1,  (0:dim(z.ind.mat)[2]) + 1, yr_frst + 0:dim(z.ind.mat)[2])
  points(w_ts, dat[i, ], pch = 16, col = clr[i])
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}

