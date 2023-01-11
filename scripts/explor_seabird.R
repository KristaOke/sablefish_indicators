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
library(MARSS)

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


longbird <- pivot_longer(zbirddat[,c(1:3, 12:17)], -c(year, n_samples, n_sablefish), names_to='index',
                         values_to = 'value')

#plot thru time
ggplot(longbird, aes(year, value)) + geom_point() + geom_line() + facet_wrap(~index)


ggplot(longbird, aes(year, value, col=index)) + geom_point() + geom_line()

#plot against each other

ggplot(zbirddat, aes(propmass_zscor, growth_index_zscor)) + geom_point() +
  geom_smooth(method="lm")

ggplot(zbirddat, aes(propmass_zscor, pred_len_zscor)) + geom_point() +
  geom_smooth(method="lm")


ggplot(zbirddat, aes( CPUE_zscor, growth_index_zscor)) + geom_point() +
  geom_smooth(method="lm")

ggplot(zbirddat, aes( CPUE_zscor, pred_len_zscor)) + geom_point() +
  geom_smooth(method="lm")


ggplot(zbirddat, aes(FreqOcr_zscor, growth_index_zscor)) + geom_point() +
  geom_smooth(method="lm")

ggplot(zbirddat, aes(FreqOcr_zscor, pred_len_zscor)) + geom_point() +
  geom_smooth(method="lm")




#DFA=========

birdsub <- zbirddat[,c(1,12:17)] 

#manupulate data for DFA========

seabird.mat <- t(as.matrix(birdsub))
colnames(seabird.mat) <- seabird.mat[1,]
seabird.mat <- seabird.mat[-1,]




#fit model========

# now fit DFA models with 1-3 trends and different error structures and compare

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1:4) {  # allowing up to 3 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(seabird.mat, model=dfa.model, control=cntl.list,
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

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data



#rotate=======

# now fit best model

model.list.1 = list(A="zero", m=4, R="diagonal and equal") # best model
model.1 = MARSS(seabird.mat, model=model.list.1, z.score=TRUE, form="dfa", control=cntl.list)


# and rotate the loadings
Z.est = coef(model.1, type="matrix")$Z
H_inv = varimax(coef(model.1, type="matrix")$Z)$rotmat
Z.rot = as.data.frame(Z.est %*% H_inv)

proc_rot = solve(H_inv) %*% model.1$states #

# reverse trend 2 to plot
Z.rot[,2] <- -Z.rot[,2]

Z.rot$names <- rownames(seabird.mat)
Z.rot <- arrange(Z.rot, V1)
Z.rot <- gather(Z.rot[,c(1,2)])
Z.rot$names <- rownames(seabird.mat)
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
N_ts <- dim(seabird.mat)[1]
## get length of time series
TT <- dim(seabird.mat)[2]

## get the estimated ZZ
Z_est <- coef(model.1, type = "matrix")$Z
## get the inverse of the rotation matrix
H_inv <- varimax(Z_est)$rotmat

## rotate factor loadings
Z_rot = Z_est %*% H_inv
## rotate processes
proc_rot = solve(H_inv) %*% model.1$states

mm <- 4 #4 processes

rec_names <- rownames(seabird.mat)
ylbl <- rec_names
w_ts <- seq(dim(seabird.mat)[2])
layout(matrix(c(1, 2, 3, 4, 5, 6,7,8), mm, 2), widths = c(2, 1))
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
  axis(1, 1:47, yr_frst + 0:dim(seabird.mat)[2])
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











