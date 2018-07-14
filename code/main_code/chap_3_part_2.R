# Phase 2 of Chapter 3 analyses: 
# Spectral response of population model under white noise at multiple levels of survival with white noise forcing

# make sure we start with a clean slate: clear all objects from current environment
rm(list = ls()) 
# functions for conducting chap 3 analysis
source("./code_ms/functions.R")

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

# Store simulated noise sets in a list for use in simualtions (easier indexing)
noiseList <- list(noise_white = white_n)
freqCont <- "white"
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
# ALT
# meanPS <- seq(0.3, 1, by = 0.1)
# sigPSmult <- seq(0.1, 0.3, by = 0.2)

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
                       white = 0) 

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

parSimCmp <- function(dt, simLen = 1024, survMean) {
  for (i in 1:length(freqCont)) {
    for (k in sigPSmult) {
      # NORMAL DISTRIBUTED NOISE
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
                                                  alpha_scale = l)[[2]])[,3])]
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

# Plot mean frequency response for each survival level.
pdf(file.path(support_dir_name, "white_noise_popFreqResp_CV.pdf"), width = 8, height = 8)
old <- par(mfrow = c(3,1), mar = c(1,4,1,1))
plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1]) 
plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2]) 
plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[3]) 
par(old)
dev.off()

# plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c <= 0.31 & sigPSmult_c >= 0.29], N = 1024, surv = meanPS[1]) 
# plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c <= 0.31 & sigPSmult_c >= 0.29], N = 1024, surv = meanPS[2]) 
# plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c <= 0.31 & sigPSmult_c >= 0.29], N = 1024, surv = meanPS[3]) 
# 
# plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c == 0.5], N = 1024, surv = meanPS[1]) 
# plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c == 0.5], N = 1024, surv = meanPS[2]) 
# plotMeanFreqR_DT(storage[ i = N > 400 & sigPSmult_c == 0.5], N = 1024, surv = meanPS[3]) 

# Main: Plot frequency response at 3 survival levels with corresponding time series
pdf(file.path(diss_dir_name, "white_noise_popFreqResp_with_TimeSeries_CV.pdf"), width = 8, height = 8)
old <- par(mfrow = c(4,1), mar = c(1,5, 1, 1))
plotMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[1], scale = "CV", yaxis_lim = c(0,4))
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = meanPS[3], line_color = "grey30", scale = "CV")
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)


plot(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[1] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[1] & reps_c == 1, j = white]),
     type = "l", col = "black", lty = 2, lwd = 1.5, 
     ylab = "Deviations", ylim = c(0,4))
plot(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[2] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[2] & reps_c == 1, j = white]),
     type = "l", col = "black", lty = 1, lwd = 2, 
     ylab = "Deviations", ylim = c(0,4))
plot(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[3] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[3] & reps_c == 1, j = white]),
     type = "l", col = "grey30", lty = 1, lwd = 2, 
     ylab = "Deviations", ylim = c(0,4))
par(old)
dev.off()

# Main: Plot frequency response at 3 survival levels higher variability
pdf(file.path(support_dir_name, "white_noise_popFreqResp_CV_hiVar.pdf"), width = 8, height = 8)
old <- par(mfrow = c(4,1), mar = c(1,5, 1, 1))
plotMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c >= 0.3 - 0.01 & sigPSmult_c <= 0.3 + 0.01 ], N = 1024, surv = meanPS[1], scale = "CV", yaxis_lim = c(0,150))
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c >= 0.3 - 0.01 & sigPSmult_c <= 0.3 + 0.01 ], N = 1024, surv = meanPS[2], line_color = "black", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c >= 0.3 - 0.01 & sigPSmult_c <= 0.3 + 0.01 ], N = 1024, surv = meanPS[3], line_color = "grey30", scale = "CV")
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)


plot(storage[i = N > 400 & sigPSmult_c >= 0.3 - 0.01 & sigPSmult_c <= 0.3 + 0.01 & meanPS_c == meanPS[1] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[1] & reps_c == 1, j = white]),
     type = "l", col = "black", lty = 2, lwd = 1.5, 
     ylab = "Deviations", ylim = c(0,10))
plot(storage[i = N > 400 & sigPSmult_c >= 0.3 - 0.01 & sigPSmult_c <= 0.3 + 0.01 & meanPS_c == meanPS[2] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[2] & reps_c == 1, j = white]),
     type = "l", col = "black", lty = 1, lwd = 2, 
     ylab = "Deviations", ylim = c(0,10))
plot(storage[i = N > 400 & sigPSmult_c >= 0.3 - 0.01 & sigPSmult_c <= 0.3 + 0.01 & meanPS_c == meanPS[3] & reps_c == 1, j = white]/mean(storage[i = N > 400 & sigPSmult_c == 0.1 & meanPS_c == meanPS[3] & reps_c == 1, j = white]),
     type = "l", col = "grey30", lty = 1, lwd = 2, 
     ylab = "Deviations", ylim = c(0,10))
par(old)
dev.off()
# ALTERNATE: PLOTS of FREQUENCY RESPONSE for meanPS from 0.3 to 1, by 0.1

# ALT
meanPS <- seq(0.3, 1, by = 0.1)
sigPSmult <- seq(0.1, 0.3, by = 0.2)

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
                       white = 0) 

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
  for (i in 1:length(freqCont)) {
    for (k in sigPSmult) {
      # NORMAL DISTRIBUTED NOISE
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

# Plots of frequency responses for many survivals for many survivals when plots scaled by mean or normalized (hi and lo variability)

pdf(file.path(diss_dir_name, "white_noise_popFreqResp_CV_lowVar_manySurvivals.pdf"), width = 8, height = 8)
plotMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.3, scale = "CV", yaxis_lim = c(0, 2))
title("Lower variabiliity (0.1) CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.4, line_color = "red", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.5, line_color = "orange", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.6, line_color = "goldenrod", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.7, line_color = "green", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.8, line_color = "blue", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.9, line_color = "purple", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 1.0, line_color = "black", scale = "CV")
legend("topright", legend = seq(0.3, 1, by = 0.1), lty = c(2,1,1,1,1,1,1,1), 
       col = c("black", "red", "orange", "goldenrod", "green", "blue", "purple", "black"), lwd = 3)
dev.off()

pdf(file.path(support_dir_name, "white_noise_popFreqResp_NORM_lowVar_manySurvivals.pdf"), width = 8, height = 8)
plotMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.3, scale = "sd", yaxis_lim = c(0, 8))
title("Lower variabiliity (0.1) normalized")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.4, line_color = "red", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.5, line_color = "orange", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.6, line_color = "goldenrod", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.7, line_color = "green", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.8, line_color = "blue", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 0.9, line_color = "purple", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.1], N = 1024, surv = 1.0, line_color = "black", scale = "norm")
legend("topright", legend = seq(0.3, 1, by = 0.1), lty = c(2,1,1,1,1,1,1,1), 
       col = c("black", "red", "orange", "goldenrod", "green", "blue", "purple", "black"), lwd = 3)
dev.off()

pdf(file.path(support_dir_name, "white_noise_popFreqResp_CV_highVar_manySurvivals.pdf"), width = 8, height = 8)
plotMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.3, scale = "CV", yaxis_lim = c(0, 60))
title("Higher variabiliity (0.3) CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.4, line_color = "red", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.5, line_color = "orange", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.6, line_color = "goldenrod", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.7, line_color = "green", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.8, line_color = "blue", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.9, line_color = "purple", scale = "CV")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 1.0, line_color = "black", scale = "CV")
legend("topright", legend = seq(0.3, 1, by = 0.1), lty = c(2,1,1,1,1,1,1,1), 
       col = c("black", "red", "orange", "goldenrod", "green", "blue", "purple", "black"), lwd = 3)
dev.off()

pdf(file.path(support_dir_name, "white_noise_popFreqResp_NORM_highVar_manySurvivals.pdf"), width = 8, height = 8)
plotMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.3, scale = "sd", yaxis_lim = c(0, 6))
title("Higher variabiliity (0.3) normalized")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.4, line_color = "red", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.5, line_color = "orange", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.6, line_color = "goldenrod", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.7, line_color = "green", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.8, line_color = "blue", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 0.9, line_color = "purple", scale = "norm")
linesMeanFR_DTmany(storage[ i = N > 400 & sigPSmult_c == 0.3], N = 1024, surv = 1.0, line_color = "black", scale = "norm")
legend("topright", legend = seq(0.3, 1, by = 0.1), lty = c(2,1,1,1,1,1,1,1), 
       col = c("black", "red", "orange", "goldenrod", "green", "blue", "purple", "black"), lwd = 3)
dev.off()

## END OF SCRIPT

fileName <- file.path(diss_dir_name, "chap_3_part_2.RData")

save.image(file = fileName)
sendmail("dpkilduff@ucdavis.edu", "chinook simulation", "It (chap_3_part_2.R) be done or it crashed!")


# # Generate increasingly reddened noise series ####
# 
# # Make matrices to hold AR1 noise with AR(1) coefficients
# # from 0.1 to 0.9, by 0.1
# ar_0 <- ar_1 <- ar_2 <- ar_3 <- ar_4 <- ar_5 <- ar_6 <- ar_7 <- ar_8 <- ar_9 <- matrix(NA, nrow = N, ncol = reps)
# # put them in a list so they can be looped over
# noiseList <- list(ar_0, ar_1, ar_2, ar_3, ar_4, ar_5, ar_6, ar_7, ar_8, ar_9)
# # declare the AR(1) autocorrelation values - phis
# phis <- seq(0, 0.9, len = length(ar_n))
# 
# # Redden white noise by increasing autocorrelation
# # All ts have same underlying random variation and variance, differing only in AR(1) coefficient
# for (i in 1:length(noiseList)) {
#   for (j in 1:ncol(white_n)) {
#     ar_n[[i]][,j] <- ar1_redden(white_n[,j], 0, 1, phi = phis[i]) 
#   }
# }
# 
# # Examine the spectral sensitivity to low frequencies ####
# 
# # Set up appropriate data.table structures for holding data
# storageP <- data.table(meanPS_c = meanPS_r,
#                        sigPSmult_c = sigPSmult_r,
#                        alphaMult_c = alphaMult_r,
#                        EQsp_c = EQsp_r,
#                        reps_c = reps_r,
#                        N = n_r,
#                        ar_0 = 0,
#                        ar_1 = 0,
#                        ar_2 = 0,
#                        ar_3 = 0,
#                        ar_4 = 0,
#                        ar_5 = 0,
#                        ar_6 = 0,
#                        ar_7 = 0,
#                        ar_8 = 0,
#                        ar_9 = 0) 
# 
# # setting the key for "fast" indexing of data.table storage
# setkey(storageP, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)
# 
# # split storageP by meanPS so each component of the
# # then I'll put all the lists back together at the end of the simulation
# french <- vector("list", length = length(meanPS)) # random name, another version I used italian named in honor of LWB
# 
# for (i in 1:length(meanPS)) {
#   french[[i]] <- storageP[i = meanPS_c == meanPS[i]]
# }
# 
# parSimCmp <- function(dt, simLen = 1000, survMean) {
#   for (i in 1:length(noiseList)) {
#     for (k in sigPSmult) {
#       # NORMALLY DISTRIBUTED NOISE
#       surv <- matrix(NA, nrow = nrow(noiseList[[i]]), ncol = ncol(noiseList[[i]]))
#       surv[1:(nrow(surv)-simLen), ] <- noiseList[[i]][1:(nrow(surv)-simLen), ] * 0.01 + survMean
#       surv[(nrow(surv)-(simLen - 1)):nrow(surv), ] <- noiseList[[i]][(nrow(surv)-(simLen - 1)):nrow(surv), ] * k + survMean
#       surv[surv > 1] <- 1
#       surv[surv < 0] <- 0
#       for (l in alphaMult) {
#         for (m in EQsp) {
#           
#           dt[i = (sigPSmult_c == k & alphaMult_c == l & EQsp_c == m), 
#              j = i+6 := list(melt(popSimPSvaryCmp(rand_surv = surv, 
#                                                   surv1 = surv1, 
#                                                   surv2 = surv2, 
#                                                   surv3 = surv3, 
#                                                   EQsp = m, 
#                                                   wanted_frac = wanted_frac,
#                                                   alpha_scale = l)[[2]])[,3]), with = FALSE]
#         }
#       }
#     }
#   }
#   return(dt)
# } # end of parSimCmp()
# 
# system.time(
#   storage <- foreach(h = 1:length(meanPS),.packages="reshape2") %dopar% {
#     parSimCmp(french[[h]], simLen = 1024, survMean = meanPS[h])
#   }
# )
# 
# storage <- rbindlist(storage)
# 
# rm(storageP, french)
#