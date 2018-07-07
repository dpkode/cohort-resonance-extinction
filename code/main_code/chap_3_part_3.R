# Phase 3 of Chapter 3 analyses: 
# Spectrum of population simulations under increasingly reddened noise ####

# make sure we start with a clean slate: clear all objects from current environment
rm(list = ls()) 
# functions for conducting chap 3 analysis
setwd("/Users/patrickkilduff/projects/cohort_resonance_risk")
source("./code/functions.R")

# Load packages
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(biwavelet)
library(mail)
library(grid)
library(TSA)

# If script run from BATCH method, then the following will record session related info for 
# reproducibility
system("hostname")  # record name of the machine
date() # record the date
sessionInfo() # documents the version of R and any included packages, 
# thereby enabling restoration of the environment used to create the results

## Set random seed for reproducibility ####
r_seed <- 64 #123
set.seed(r_seed) # ensure the reproducibility of results

todays_date <-  Sys.Date()
diss_dir_name <- file.path(".", "output", "sim_results", paste0("dissertation_files_", todays_date))
if(!exists(diss_dir_name)) dir.create(diss_dir_name, recursive = TRUE)

support_dir_name <- file.path(".", "output", "sim_results", paste0("support_files_", todays_date))
if (!exists(support_dir_name)) dir.create(support_dir_name, recursive = TRUE)


reps <- 1000
N <- 1424

# Simulate time series white noise random variables 
# using the random sine wave approach

# WHITE NOISE ####

# create "reps" number of white noise time series of length N using random sine wave approach
# with selected frequency contents
white_n <- customFR2ts(N = N, # number of time steps
                       reps = reps,
                       r_seed = r_seed,
                       amp = mk_white(N)) # mean = 0, variance/sd = 1

surv1 <- 0.02 # first year ocean survival
surv2 <- surv3 <- 0.8 # ocean survival in later years
wanted_frac <- 0.5 # parameter to determine the fraction spawning at age-3 and 4 (equilibrium).

# First, the paper as is does not show the spectral sensitivity of the population.  
# I think you should insert the frequency response of the population so that readers can see 
# how the inputs are related to that.  I would do that for several different levels of survival, 
# so that they can see the different degrees of cohort resonance.  Also, they will see that CR 
# involves greater sensitivity at low frequencies, not just at the generational frequencies.

# Frequency response of "population" at three levels of survival for white noise
alphaMult <- 4
EQsp <- 7500

# BASE
meanPS <- c(0.275, 0.5, 0.8)
sigPSmult <- seq(0.1, 0.5, by = 0.1)

# Generate increasingly reddened noise series ####

# Make matrices to hold AR1 noise with AR(1) coefficients
# from 0.1 to 0.9, by 0.1
ar_0 <- ar_1 <- ar_2 <- ar_3 <- ar_4 <- ar_5 <- ar_6 <- ar_7 <- ar_8 <- ar_9 <- matrix(NA, nrow = N, ncol = reps)
# put them in a list so they can be looped over
noiseList <- list(ar_0, ar_1, ar_2, ar_3, ar_4, ar_5, ar_6, ar_7, ar_8, ar_9)
# declare the AR(1) autocorrelation values - phis
phis <- seq(0, 0.9, len = length(noiseList))

# Redden white noise by increasing autocorrelation
# All ts have same underlying random variation and variance, differing only in AR(1) coefficient
for (i in 1:length(noiseList)) {
  for (j in 1:ncol(white_n)) {
    noiseList[[i]][,j] <- ar1_redden(white_n[,j], 0, 1, phi = phis[i]) 
  }
}

# Examine the spectral sensitivity to low frequencies ####

# Set up appropriate data.table structures for holding data

meanPS_r <- rep(meanPS, each = length(sigPSmult)*length(alphaMult)*length(EQsp)*reps*N)
sigPSmult_r <- rep(sigPSmult, each = length(alphaMult)*length(EQsp)*reps*N, times = length(meanPS))
alphaMult_r <- rep(alphaMult, each = length(EQsp)*reps*N, times = length(meanPS)*length(sigPSmult))
EQsp_r <- rep(EQsp, each = reps*N, times = length(meanPS)*length(sigPSmult)*length(alphaMult))
reps_r <- rep(1:reps, each = N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp))
n_r <- rep(1:N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp)*reps)

storageP <- data.table(meanPS_c = meanPS_r,
                       sigPSmult_c = sigPSmult_r,
                       alphaMult_c = alphaMult_r,
                       EQsp_c = EQsp_r,
                       reps_c = reps_r,
                       N = n_r,
                       ar_0 = 0,
                       ar_1 = 0,
                       ar_2 = 0,
                       ar_3 = 0,
                       ar_4 = 0,
                       ar_5 = 0,
                       ar_6 = 0,
                       ar_7 = 0,
                       ar_8 = 0,
                       ar_9 = 0) 

# setting the key for "fast" indexing of data.table storage
setkey(storageP, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)

# split storageP by meanPS so each component of the
# then I'll put all the lists back together at the end of the simulation
french <- vector("list", length = length(meanPS)) # random name, another version I used italian named in honor of LWB

for (i in 1:length(meanPS)) {
  french[[i]] <- storageP[i = meanPS_c == meanPS[i]]
}

# pre-compile functions to improve run time.
library(compiler)
popSimPSvaryCmp <- cmpfun(popSimPSvary)

parSimCmp <- function(dt, simLen = 1000, survMean) {
  for (i in 1:length(noiseList)) {
    for (k in sigPSmult) {
      # NORMALLY DISTRIBUTED NOISE
      surv <- matrix(NA, nrow = nrow(noiseList[[i]]), ncol = ncol(noiseList[[i]]))
      surv[1:(nrow(surv)-simLen), ] <- noiseList[[i]][1:(nrow(surv)-simLen), ] * 0.01 + survMean
      surv[(nrow(surv)-(simLen - 1)):nrow(surv), ] <- noiseList[[i]][(nrow(surv)-(simLen - 1)):nrow(surv), ] * k + survMean
      surv[surv > 1] <- 1
      surv[surv < 0] <- 0
      for (l in alphaMult) {
        for (m in EQsp) {
          
          dt[i = (sigPSmult_c == k & alphaMult_c == l & EQsp_c == m), 
             j = i+6 := list(melt(popSimPSvaryCmp(rand_surv = surv, 
                                                  surv1 = surv1, 
                                                  surv2 = surv2, 
                                                  surv3 = surv3, 
                                                  EQsp = m, 
                                                  wanted_frac = wanted_frac,
                                                  alpha_scale = l)[[2]])[,3]), with = FALSE]
        }
      }
    }
  }
  return(dt)
} # end of parSimCmp()

# run simulations using foreach framework to send jobs to multiple cores
# for each chunk of storageP - this should cut simulation times down to ~ 68/9 = 7.5 hours

# load foreach package for parallel processing
library(foreach)
# set up backend to do parallel processing
library(doParallel)
# detectCores()
registerDoParallel() # defaults to half of the cores 


system.time(
  storage <- foreach(h = 1:length(meanPS),.packages="reshape2") %dopar% {
    parSimCmp(french[[h]], simLen = 1024, survMean = meanPS[h])
  }
)

storage <- rbindlist(storage)

rm(storageP, french)


# plot spectral frequency responses for a given survival as a function of AR(1)
# Use white_n as a base case

# Main: Plot frequency response at 3 survival levels
pdf(file.path(diss_dir_name, "freqResp3SurvReddened.pdf"), width = 10, height = 10)
old <- par(mfrow = c(4,1), mar = c(1.7,5, 1, 1))
plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1 ], N = 1024, surv = meanPS[1], scale = "CV", AR_col = "ar_0", yaxis_lim = c(0, 3))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1 ], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV", AR_col = "ar_0")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1 ], N = 1024, surv = meanPS[3], line_color = "grey50", scale = "CV", AR_col = "ar_0")
text(0.48, 1, expression(phi == 0), cex = 2)
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey50"), lwd = 3, cex = 1.5)

plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], scale = "CV", AR_col = "ar_3", yaxis_lim = c(0, 3))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV", AR_col = "ar_3")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[3], line_color = "grey50", scale = "CV", AR_col = "ar_3")
text(0.48, 1, expression(phi == 0.3), cex = 2)
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey50"), lwd = 3, cex = 1.5)

plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], scale = "CV", AR_col = "ar_6", yaxis_lim = c(0, 3))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV", AR_col = "ar_6")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[3], line_color = "grey50", scale = "CV", AR_col = "ar_6")
text(0.48, 1, expression(phi == 0.6), cex = 2)
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey50"), lwd = 3, cex = 1.5)

plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], scale = "CV", AR_col = "ar_9", yaxis_lim = c(0, 3))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV", AR_col = "ar_9")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[3], line_color = "grey50", scale = "CV", AR_col = "ar_9")
text(0.48, 1, expression(phi == 0.6), cex = 2)
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey50"), lwd = 3, cex = 1.5)
par(old)
dev.off()

       
pdf(file.path(diss_dir_name, "freqRespReddenedLowSurvCV.pdf"), width = 10, height = 10)       
plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1 ], N = 1024, surv = meanPS[1], scale = "CV", AR_col = "ar_0", yaxis_lim = c(0, 3))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], line_color = "black", scale = "CV", AR_col = "ar_5")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], line_color = "grey30", scale = "CV", AR_col = "ar_9")
legend("topright", legend = c("white", "0.5", "0.9"), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)
title(paste("CV low summer survival:", meanPS[1]))
dev.off()

pdf(file.path(support_dir_name, "freqRespReddenedLowSurvNORM.pdf"), width = 10, height = 10)
plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], scale = "sd", AR_col = "ar_0", yaxis_lim = c(0, 5))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], line_color = "black", scale = "sd", AR_col = "ar_5")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], line_color = "grey30", scale = "sd", AR_col = "ar_9")
legend("topright", legend = c("white", "0.5", "0.9"), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)
title(paste("Normalized low summer survival:", meanPS[1]))
dev.off()

pdf(file.path(diss_dir_name, "freqRespReddenedModerateSurvCV.pdf"), width = 10, height = 10)  
plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1 ], N = 1024, surv = meanPS[2], scale = "CV", AR_col = "ar_0", yaxis_lim = c(0, 0.5))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV", AR_col = "ar_5")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "grey30", scale = "CV", AR_col = "ar_9")
legend("topright", legend = c("white", "0.5", "0.9"), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)
title(paste("CV mid summer survival:", meanPS[2]))
dev.off()

pdf(file.path(support_dir_name, "freqRespReddenedModerateSurvNORM.pdf"), width = 10, height = 10)  
plotMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1 ], N = 1024, surv = meanPS[2], scale = "sd", AR_col = "ar_0", yaxis_lim = c(0, 6))
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "black", scale = "sd", AR_col = "ar_5")
linesMeanFR_DTmanyAR(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "grey30", scale = "sd", AR_col = "ar_9")
legend("topright", legend = c("white", "0.5", "0.9"), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)
title(paste("Normalized mid summer survival:", meanPS[2]))
dev.off()
# plot(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[1] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[1] & reps_c == 1, j = white]),
#      type = "l", col = "black", lty = 2, lwd = 1.5, 
#      ylab = "CV")
# plot(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[2] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[2] & reps_c == 1, j = white]),
#      type = "l", col = "black", lty = 1, lwd = 2, 
#      ylab = "CV")
# plot(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[3] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[3] & reps_c == 1, j = white]),
#      type = "l", col = "grey30", lty = 1, lwd = 2, 
#      ylab = "CV")
# par(old)

fileName <- file.path(diss_dir_name, "chap_3_part_3.RData")

save.image(file = fileName)
sendmail("dpkilduff@ucdavis.edu", "chinook simulation", "It (chap_3_part_3.R) be done or it crashed!")

