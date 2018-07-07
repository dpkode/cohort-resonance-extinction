# Phase 4 of Chapter 3 analyses: 
# Extinction dynamics with increasing redness ####

# For each level of the 3 levels of CR you could look at pQE with white noise, 
# then various levels of reddening obtained by various levels of combining sequential values in 
# the noise, e.g., .9 and .1, .8 and .2, .7 and .3, etc.  Results could be shown as a plot of pQE 
# vs. redness, with a line for each level of CR.

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


# First generated the normal random variables ####

reps <- 1000
N <- 500

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
sigPSmult <- seq(0.1, 0.5, by = 0.025)

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

parSimCmp <- function(dt, simLen = 100, survMean) {
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
    parSimCmp(french[[h]], simLen = 100, survMean = meanPS[h])
  }
)

storage <- rbindlist(storage)

rm(storageP, french)

## PLOT TIME SERIES OF SPAWNER ABUNDANCE to make sure they look okay! ####
# he should have 
# length(unique(storage$N))
# length(unique(storage$reps))
spPlot <- copy(storage[ i = N > 350 & reps_c <= 10 & sigPSmult_c %in% c(0.1, 0.2, 0.3, 0.4, 0.5)])
spPlot_m <- melt(spPlot, id = 1:6)
pdf(file.path(support_dir_name, "spawnerTSplot_AR_noise.pdf"), width = 8, height = 6)
ggplot(data = spPlot_m, aes(x = N, y = value, color = factor(meanPS_c))) + 
  facet_grid(sigPSmult_c ~ variable) +
  geom_point(alpha = 0.08, ) +
  geom_line(alpha = 0.08) +
  scale_colour_hue() + #manual(values = c("red", "pink", "black")) +
  ylab("Spawner abundance") +
  xlab("Time (Years)") +
  guides(colour = guide_legend(title = "Survival")) +
  theme_bw()
dev.off()

## ANALYSIS PLOTS: POPULATION VARIABILTY, QE TIME DISTRIBUTION, CV (mean and sd), PQE ####

# The first 400 years of the simulation are run with low noise so population runs ~ equlibrium before
# running the model for 100 years. Thus drop the first 400 years

storage_sub <- copy(storage[ i = (N > 400)]) # Only the last 100 years, the first years are "burn-in"

# CALCULATE QE YEAR ####

# A function wrapper format for data.tables
calcQEyr <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)]
}

# qeLev <- 100
# run_length <- 4
# Call the function
sb_qeyr <- calcQEyr(dt = storage_sub, expr = list( as.integer(JA_consec(ar_0, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_1, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_2, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_3, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_4, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_5, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_6, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_7, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_8, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(ar_9, run_length = 4, qeLev = 100))) )

# need to name the new columns
setnames(sb_qeyr, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
         c("ar_0", "ar_1", "ar_2", "ar_3", "ar_4", "ar_5", "ar_6", "ar_7", "ar_8", "ar_9"))


## PLOT QE TIME DISTRIBUTIONS

pdf(file.path(diss_dir_name, "AR_Noise_QE_time_Distribution.pdf"), width = 10, height = 6)

#for (j in 1:length(sigPSmult)) {
  qet_tmp_sb_m <- melt(copy(sb_qeyr[ i = (sigPSmult_c == sigPSmult[13] & meanPS_c %in% c(meanPS[1]))]), id = c(1:5))
  #droplevels(subset(qet_tmp_sb_m, !is.na(value) )
  x <- ggplot(qet_tmp_sb_m[ which(qet_tmp_sb_m$variable %in% c("ar_0", "ar_1", "ar_2", "ar_3", "ar_4", "ar_5")), ], aes(x = value)) + 
    geom_histogram() + 
    facet_grid(meanPS_c ~ variable) +
    #ggtitle(paste("vertical survival \n noise sd =", sigPSmult[j], "\n", "SR slope mult = 4") ) +
    ylim(c(0,350)) +
    xlab("Time (Year)") +
    ylab("Count") +
    theme_bw()
  print(x)
#}
dev.off()

# CALCULATE PQE ####

calc_pQE <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)]
}

sb_pQE <- calc_pQE(sb_qeyr, list(ar0_qe = length(which(!is.na(ar_0)))/1000,
                                 ar1_qe = length(which(!is.na(ar_1)))/1000,
                                 ar2_qe = length(which(!is.na(ar_2)))/1000,
                                 ar3_qe = length(which(!is.na(ar_3)))/1000,
                                 ar4_qe = length(which(!is.na(ar_4)))/1000,
                                 ar5_qe = length(which(!is.na(ar_5)))/1000,
                                 ar6_qe = length(which(!is.na(ar_6)))/1000,
                                 ar7_qe = length(which(!is.na(ar_7)))/1000,
                                 ar8_qe = length(which(!is.na(ar_8)))/1000,
                                 ar9_qe = length(which(!is.na(ar_9)))/1000)  )

pdf(file.path(diss_dir_name, "AR_Noise_pQE.pdf"), width = 10, height = 10)
for (i in 1:length(alphaMult)) {
  pqe_tmp_m <- melt(copy(sb_pQE[ i =  alphaMult_c == alphaMult[i]]), id = c(1:4))
  x <- ggplot(pqe_tmp_m[ which(qet_tmp_sb_m$variable %in% c("ar_0", "ar_1", "ar_2", "ar_3", "ar_4", "ar_5")), ], aes(x = sigPSmult_c, y = value)) + 
    geom_line(colour = "slateblue", size = 1.2) + 
    facet_grid(meanPS_c ~ variable) +
    ggtitle(paste("p(Quasi-Extinction) \n rows: survival \n columns: AR(1) coeff") ) +
    ylab("Probability of Quasi-Extinction") +
    xlab(expression(paste(sigma)))+
    theme(axis.text = element_text(size = 8, colour = "black"))+
    theme_bw() 
  print(x)
}
dev.off()

pdf(file.path(diss_dir_name, "AR_Noise_pQE1plot.pdf"), width = 10, height = 6)

  pqe_tmp_m <- melt(copy(sb_pQE[ i =  alphaMult_c == alphaMult & meanPS_c == meanPS[1]]), id = c(1:4))  
  x <- ggplot(pqe_tmp_m[which(pqe_tmp_m$variable %in% c("ar0_qe", "ar1_qe", "ar2_qe", "ar3_qe", "ar4_qe", "ar5_qe")),], aes(x = sigPSmult_c, y = value)) + 
    geom_line(colour = "slateblue", size = 1.2) + 
    facet_grid(meanPS_c ~ variable) +
    #ggtitle(paste("p(Quasi-Extinction) \n rows: survival \n columns: AR(1) coeff") ) +
    ylab("Probability of Quasi-Extinction") +
    xlab(expression(paste(sigma)))+
    theme(axis.text = element_text(size = 8, colour = "black")) +
  theme_bw()
  print(x)

dev.off()

fileName <- file.path(diss_dir_name, "chap_3_part_4.RData")

save.image(file = fileName)
sendmail("dpkilduff@ucdavis.edu", "chinook simulation", "It (chap_3_part_4.R) be done or it crashed!")

