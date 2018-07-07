# Script to produce plots for working chapter 3


# make sure we start with a clean slate: clear all objects from current environment
rm(list = ls()) 
# functions for conducting chap 3 analysis
source("./code/functions.R")

# Load packages
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(biwavelet)
library(mail)
library(grid)

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


## DETERMINE realized mean and sd oflognormal survival rates based on specified mean and sd ####

# Examples of truncated log-normal and normal distributions of white noise

# The goal is to examine log-normal and normal survivals between 0 and 1 generated
# that have representative mean and sd based on "observations" from 
# the projections from Thompson et al. 2012

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

# COLORED NOISE ####
# Simulate different 'colors' of noise [really just specified bandwidths]

# Bandpass period 3-4 (frequencies 1/4 to 1/3) in white noise signals.

rsin_34_n <- customFR2ts(N = N, # number of time steps
                         reps = reps,
                         r_seed = r_seed,
                         amp = mk_rsin(N, highF=1/3, lowF=1/4)) # mean = 0, variance/sd = 1

# Bandreject period 3-4 (frequencies 1/4 to 1/3) in white noise signals.

rsin_34_reject_n <- customFR2ts(N = N,
                                reps = reps,
                                r_seed = r_seed,
                                amp = mk_rsin_reject(N, lowF=1/4, highF=1/3) )
# Band-pass greater than period 4 (frequencies lower than 0.25)
rsin_gt4_n <- customFR2ts(N = N,
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_rsin(N, lowF=0, highF=1/4) )

# Band-pass less than period 3 (frequencies higher than 0.33)
rsin_lt3_n <- customFR2ts(N = N,
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_rsin(N, lowF=1/3, highF=1/2) )

# Band-pass greater than period 3 (frequencies lower than 0.33)
rsin_gt3_n <- customFR2ts(N = N,
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_rsin(N, lowF=0, highF=1/3) )

# Band-pass less than period 4 (frequencies higher than 0.25)
rsin_lt4_n <- customFR2ts(N = N,
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_rsin(N, lowF=1/4, highF=1/2) )

# Band-pass greater than period 10 (frequencies lower than 0.1)
rsin_gt10_n <- customFR2ts(N = N,
                           reps = reps,
                           r_seed = r_seed,
                           amp = mk_rsin(N, lowF=0, highF=1/10) )  

# Band-pass greater than period 10 and period 3-4 (frequencies lower than 0.1 plus period 3-4)
# By adding the period 10 and greater noise to the period 3-4 only noise (divide by 2) 
rsin_34_gt10_n <- matrix(NA, nrow = N, ncol = reps)
for (i in 1:reps) {
  rsin_34_gt10_n[,i] <- (rsin_gt10_n[,i] + rsin_34_n[,i])/2
}
rsin_34_gt10_n <- apply(rsin_34_gt10_n, 2, scale)

# Store simulated noise sets in a list for use in simualtions (easier indexing)
noiseList <- list(noise_white = white_n,
                  noise_34 = rsin_34_n,
                  nose_34_reject = rsin_34_reject_n,
                  noise_gt4 = rsin_gt4_n,
                  noise_lt3 = rsin_lt3_n,
                  noise_gt3 = rsin_gt3_n,
                  noise_lt4 = rsin_lt4_n,
                  noise_gt10 = rsin_gt10_n,
                  noise_34gt10 = rsin_34_gt10_n)

# rm individual sets of noise
rm(white_n, rsin_34_n, rsin_34_reject_n, rsin_gt4_n, rsin_lt3_n, rsin_gt3_n, rsin_lt4_n, rsin_gt10_n, rsin_34_gt10_n)

## ANALYSIS OF TRUNCATED SURVIVAL DISTRIBUTIONS #####

# Set the parameters for the survival values (mean and variability)

# "fixed" parameters for subsequent noise generation and population modeling
freqCont <- c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10")
surv1 <- 0.02 # first year ocean survival
surv2 <- surv3 <- 0.8 # ocean survival in later years
wanted_frac <- 0.5 # parameter to determine the fraction spawning at age-3 and 4 (equilibrium).
# Before generating the random variables of noise to force survival
# determine the mean survival level based on the (low) survival rate
# the produces stock collapse [1/SPR == alpha] and then set survival
# rates at incrementally higher levels.

# Run simulations that approach the persistence threshold for each alpha level
# by dropping meanPS from some arbitrary distance above 1/SPR collapse 

# calculate the fraction spawning early
delta_e <- calc_de(wanted_frac, surv3)
# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = 1)  
alphaMult <- 4 # multiplier of replacement line slope to set alpha [from 2-10]
alpha <- alphaMult * 1/spr
sps_crash <- 1/{alpha*surv1*surv2*{delta_e + surv3*{1 - delta_e}}}

#sps_adds <- seq(0.05, 0.3, by = 0.05) # 0.05 increments in survival
sps_mults <- seq(1.1, 1.5, by = 0.1) #seq(1.1, 1.5, by = 0.2)

# SET REMAINING PARAMETERS FOR EVALUATING TRUNCATED NORMAL AND LOGNORMAL SURVIVAL RATES 0 to 1
# Set up matrix of survival rates based on persistence threshold conditions
# for each alphaMult

meanPS <- outer(sps_mults, sps_crash, '*') 
sigPSmult <- seq(.1, 1.5, by = 0.2)


# HISTOGRAMS with specified mean and sd and realized of all replicates #####

# LOGNORMAL ####
for (i in 1:length(freqCont)) {
  f_path <- paste0(support_dir_name, "/trunc_lognormal_", freqCont[i], ".pdf", sep = "")
  pdf(f_path)
  for (j in 1:length(meanPS)) {
    for (k in sigPSmult)  {
      lsurv <- exp(log(meanPS[j]) + noiseList[[i]] * k )
      lsurv[lsurv > 1] <- 1
      lsurv[lsurv < 0] <- 0
      old <- par(mar = c(2,2,5,1))
      hist(lsurv, 
           col = "goldenrod", main = "", 
           xlim = c(0,1))
      title(main = paste("Truncated lognormal survival rates \n Frequency content:", 
                         freqCont[i], 
                         "\n prespawn surv =", meanPS[j],
                         "\n sigma_env = ", k,
                         "\n obs SD:", round(sd(lsurv),2),
                         "\n obsMean:", round(mean(lsurv),2) ),
            cex.main = 0.8)
      par(old)
    }
  }
  dev.off()
}

## NORMAL #####

for (i in 1:length(freqCont)) {
  f_path <- paste0(support_dir_name, "/trunc_normal_", freqCont[i], ".pdf", sep = "")
  pdf(f_path)
  for (j in 1:length(meanPS)) {
    for (k in sigPSmult) {
      surv <- noiseList[[i]] * k + meanPS[j]
      surv[surv > 1] <- 1
      surv[surv < 0] <- 0
      old <- par(mar = c(2,2,5,1))
      hist(surv, 
           col = "goldenrod", main = "", 
           xlim = c(0,1))
      title(main = paste("Truncated normal survival rates \n Frequency content:", 
                         freqCont[i],  
                         "\n prespawn surv =", meanPS[j],
                         "\n sigma_env = ", k,
                         "\n obs SD:", round(sd(surv),2),
                         "\n obsMean:", round(mean(surv),2) ),
            cex.main = 0.8)
      par(old)
    }
  }
  dev.off()
}

# data.frame of SPECIFIED and REALIZED TRUNCATED SURVIVAL RATES for each FREQUENCY BAND ####

freqCont_r <- rep(freqCont, each = length(meanPS)*length(sigPSmult))
meanPS_r <-rep(meanPS, each = length(sigPSmult), times = length(freqCont))
sigPSmult_r <- rep(sigPSmult, times = length(freqCont)*length(meanPS))

survCompDF <- survCompDFnorm <- data.frame(freqCont = freqCont_r,
                                           specMean = meanPS_r,
                                           specSD = sigPSmult_r,
                                           realMean = 0,
                                           realSD = 0)
# LOGNORMAL ####
for (i in 1:length(freqCont)) {
  for (j in 1:length(meanPS)) {
    for(k in sigPSmult ) {
      lsurv <- exp(log(meanPS[j]) + noiseList[[i]] * k )
      lsurv[lsurv > 1] <- 1
      lsurv[lsurv < 0] <- 0
      
      ind <- which(survCompDF$freqCont == freqCont[i] & survCompDF$specMean == meanPS[j]  & survCompDF$specSD == k )
      survCompDF[ind, 4] <- mean(colMeans(lsurv))
      survCompDF[ind, 5]<- mean(apply(lsurv, 2, sd))
    }
  }
}

# NORMAL ####
for (i in 1:length(freqCont)) {
  for (j in 1:length(meanPS)) {
    for (k in sigPSmult) {
      surv <- noiseList[[i]] * k + meanPS[j]
      surv[surv > 1] <- 1
      surv[surv < 0] <- 0
      
      ind <- which(survCompDF$freqCont == freqCont[i] & survCompDF$specMean == meanPS[j]  & survCompDF$specSD == k )
      survCompDFnorm[ind, 4] <- mean(colMeans(surv))
      survCompDFnorm[ind, 5]<- mean(apply(surv, 2, sd))
    }
  }
}



# PLOT of specified mean vs. realized mean ####
pdf(paste0(support_dir_name, "/specMeanRealMeanLOGNORM.pdf", sep = ""), width = 8, height = 8)
ggplot(data = survCompDF, aes(x=specMean, y = realMean, colour = specSD)) +
  geom_point() +
  scale_colour_gradientn(colours = rainbow(7)) +
  facet_wrap( ~ freqCont) +
  geom_abline(intercept = 0, slope = 1, colour = "red", size = 2) +
  xlab("Specified mean survival") +
  ylab("Actual mean survival") +
  ggtitle("Lognormal") +
  theme_bw()
dev.off()

pdf(paste0(support_dir_name, "/specMeanRealMeanNORM.pdf", sep = ""), width = 8, height = 8)
ggplot(data = survCompDFnorm, aes(x=specMean, y = realMean, colour = specSD)) +
  geom_point() +
  scale_colour_gradientn(colours = rainbow(7)) +
  facet_wrap( ~ freqCont) +
  geom_abline(intercept = 0, slope = 1, colour = "red", size = 2) +
  xlab("Specified mean survival") +
  ylab("Actual mean survival") +
  ggtitle("Normal") +
  theme_bw()
dev.off()

# PLOT of specified sd vs. realized sd ####
pdf(paste0(support_dir_name, "/specSD_RealSD_LOGNORM.pdf", sep = ""), width = 8, height = 8)
ggplot(data = survCompDF, aes(x=specSD, y = realSD, colour = specMean)) +
  geom_point() +
  scale_colour_gradientn(colours = rainbow(7)) +
  facet_wrap( ~ freqCont) +
  geom_abline(intercept = 0, slope = 1, colour = "red", size = 2) +
  xlab("Specified mean survival") +
  ylab("Actual mean survival") +
  ggtitle("Lognormal") +
  theme_bw()
dev.off()

pdf(paste0(support_dir_name, "/specSD_RealSD_NORM.pdf", sep = ""), width = 8, height = 8)
ggplot(data = survCompDFnorm, aes(x=specSD, y = realSD, colour = specMean)) +
  geom_point() +
  scale_colour_gradientn(colours = rainbow(7)) +
  facet_wrap( ~ freqCont) +
  geom_abline(intercept = 0, slope = 1, colour = "red", size = 2) +
  xlab("Specified mean survival") +
  ylab("Actual mean survival") +
  ggtitle("Normal") +
  theme_bw()
dev.off()

## LOOMING pQE QUESTION ####
# Why are their bumps in the quasi-extinction where pQE jumps for low to moderate variability (of specified sd)
# is it an error or a real thing
# Look at 3 cases from initial analyses

## SIMULATIONS ####

# Set up the vectors describing the columns of the data.table (~ data.frame - just faster) 
# called "storage"

# since we have separate survival rates to monitor for each alpha, 
# there needs to be separate data.table for each alphaMult and its 
# corresponding meanPS levels
# set-up "empty" data.tables for each alphaMult with corresponding meanPS --> persistence threshold
# -- a list [storageP] of data.tables for each alphaMult (rather than for meanPS - since that now varies with alphaMult)
sigPSmult <-  seq(.1, 0.5, by = 0.1) # seq(.1, 0.5, by = 0.2)  seq(.1, 0.6, by = 0.025)
alphaMult <- 4
meanPS <- c(0.35, 0.40)
EQsp <- 7500

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
                       white = 0,
                       p34 = 0,
                       pNo34 = 0,
                       pgt4 = 0,
                       plt3 = 0,
                       pgt3 = 0,
                       plt4 = 0,
                       pgt10 = 0,
                       p34gt10 = 0) 

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
pdf(paste0(support_dir_name, "/spawnerTSplot.pdf", sep = ""), width = 12, height = 8)
ggplot(data = spPlot_m, aes(x = N, y = value, color = factor(meanPS_c))) + 
  facet_grid(sigPSmult_c ~ variable) +
  geom_point(alpha = 0.08) +
  geom_line(alpha = 0.08) +
  scale_colour_manual(values = c("red","blue")) +
  ylab("Spawner abundance") +
  xlab("Time (Years)") +
  theme(axis.text = element_text(size = 8, colour = "black")) +
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
sb_qeyr <- calcQEyr(dt = storage_sub, expr = list( as.integer(JA_consec(white, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(p34, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(pNo34, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(pgt4, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(plt3, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(pgt3, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(plt4, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(pgt10, run_length = 4, qeLev = 100)),
                                                   as.integer(JA_consec(p34gt10, run_length = 4, qeLev = 100)) ) )

# need to name the new columns
setnames(sb_qeyr, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))


## PLOT QE TIME DISTRIBUTIONS

pdf(paste0(support_dir_name, "/QE_BUMP_QE_time_Distribution.pdf", sep = ""), width = 12, height = 8)

for (j in 1:length(sigPSmult)) {
  qet_tmp_sb_m <- melt(copy(sb_qeyr[ i = (sigPSmult_c == sigPSmult[3])]), id = c(1:5))
  #droplevels(subset(qet_tmp_sb_m, !is.na(value) )
  x <- ggplot(qet_tmp_sb_m, aes(x = value)) + 
    geom_histogram() + 
    facet_grid(meanPS_c ~ variable) +
    ggtitle(paste("vertical survival \n noise sd =", sigPSmult[j], "\n", "SR slope mult = 4") ) +
    ylim(c(0,350)) +
    xlab("Time (Year)") +
    ylab("Count") +
    theme(axis.text = element_text(size = 8, colour = "darkgrey")) +
    theme_bw()
  print(x)
}

dev.off()

# CALCULATE PQE ####

calc_pQE <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)]
}

sb_pQE <- calc_pQE(sb_qeyr, list(white_qe = length(which(!is.na(white)))/1000,
                                 p34_qe = length(which(!is.na(p34)))/1000,
                                 pNo34_qe = length(which(!is.na(pNo34)))/1000,
                                 pgt4_qe = length(which(!is.na(pgt4)))/1000,
                                 plt3_qe = length(which(!is.na(plt3)))/1000,
                                 pgt3_qe = length(which(!is.na(pgt3)))/1000,
                                 plt4_qe = length(which(!is.na(plt4)))/1000,
                                 pgt10_qe = length(which(!is.na(pgt10)))/1000,
                                 p34gt10_qe = length(which(!is.na(p34gt10)))/1000 )  )

pdf(paste0(support_dir_name, "/QE_BUMP_pQE.pdf", sep = ""), width = 12, height = 8)
for (i in 1:length(alphaMult)) {
  pqe_tmp_m <- melt(copy(sb_pQE[ i =  alphaMult_c == alphaMult[i]]), id = c(1:4))
  x <- ggplot(pqe_tmp_m, aes(x = sigPSmult_c, y = value)) + 
    geom_line(colour = "slateblue", size = 1.2) + 
    facet_grid(meanPS_c ~ variable) +
    ggtitle(paste("p(Quasi-Extinction) \n rows: survival \n columns: frequency bands") ) +
    ylab("Probability of Quasi-Extinction") +
    xlab(expression(paste(sigma)))+
    theme_bw()
  print(x)
}
dev.off()

# CV, SD, MEAN -- > story to loo #### 

calcCVsp <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)]
}

sb_CVsp <- calcCVsp(dt = storage_sub, expr = list( as.double(sd(white, na.rm = TRUE)/mean(white, na.rm = TRUE)),
                                                   as.double(sd(p34, na.rm = TRUE)/mean(p34, na.rm = TRUE)),
                                                   as.double(sd(pNo34, na.rm = TRUE)/mean(pNo34, na.rm = TRUE)),
                                                   as.double(sd(pgt4, na.rm = TRUE)/mean(pgt4, na.rm = TRUE)),
                                                   as.double(sd(plt3, na.rm = TRUE)/mean(plt3, na.rm = TRUE)),
                                                   as.double(sd(pgt3, na.rm = TRUE)/mean(pgt3, na.rm = TRUE)),
                                                   as.double(sd(plt4, na.rm = TRUE)/mean(plt4, na.rm = TRUE)),
                                                   as.double(sd(pgt10, na.rm = TRUE)/mean(pgt10, na.rm = TRUE)),
                                                   as.double(sd(p34gt10, na.rm = TRUE)/mean(p34gt10, na.rm = TRUE)) ) ) 

setnames(sb_CVsp, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

pdf(paste0(support_dir_name, "/QE_BUMP_CVsp.pdf", sep = ""), width = 12, height = 8)
for (i in 1:length(alphaMult)) {
  cv_tmp_m <- melt(copy(sb_CVsp[ i =   alphaMult_c == alphaMult[i]]), id = c(1:5))
  x <- ggplot(cv_tmp_m, aes(x = factor(sigPSmult_c), y = value)) + 
    geom_violin(fill = "slateblue") + 
    facet_grid(meanPS_c ~ variable) +
    ggtitle(paste("Coefficient of variation spawner abunance \n rows: survival \n columns: frequency bands") ) +
    ylab("CV spaawner abundance") +
    xlab(expression(paste(sigma))) +
    theme(axis.text.x = element_text(size = 7.5)) +
    scale_x_discrete(breaks=seq(0.1,0.5, by = 0.1)) +
    theme_bw()
  print(x)
}
dev.off()

# MEAN spawning abundance ####

calcMeansp <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)]
}

sb_Meansp <- calcMeansp(dt = storage_sub, expr = list( as.double(mean(white, na.rm = TRUE)),
                                                       as.double(mean(p34, na.rm = TRUE)),
                                                       as.double(mean(pNo34, na.rm = TRUE)),
                                                       as.double(mean(pgt4, na.rm = TRUE)),
                                                       as.double(mean(plt3, na.rm = TRUE)),
                                                       as.double(mean(pgt3, na.rm = TRUE)),
                                                       as.double(mean(plt4, na.rm = TRUE)),
                                                       as.double(mean(pgt10, na.rm = TRUE)),
                                                       as.double(mean(p34gt10, na.rm = TRUE)) ) ) 

setnames(sb_Meansp, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

pdf(paste0(support_dir_name, "/QE_BUMP_mean.pdf", sep = ""), width = 12, height = 8)
for (i in 1:length(alphaMult)) {
  cv_tmp_m <- melt(copy(sb_Meansp[ i =   alphaMult_c == alphaMult[i]]), id = c(1:5))
  x <- ggplot(cv_tmp_m, aes(x = factor(sigPSmult_c), y = value)) + 
    geom_violin(fill = "slateblue") + 
    facet_grid(meanPS_c ~ variable) +
    ggtitle(paste("Mean spawner abundance \n rows: survival \n columns: frequency bands") ) +
    ylab("Mean spaawner abundance") +
    xlab(expression(paste(sigma))) +
    theme(axis.text.x = element_text(size = 7.5)) +
    scale_x_discrete(breaks=seq(0.1,0.5, by = 0.1)) + 
    theme_bw()
  print(x)
}
dev.off()

# STANDARD DEVIATION spawning abundance ####

calcSDsp <- function(dt, expr) {
  e <- substitute(expr) 
  dt[,eval(e), by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)]
}

sb_SDsp <- calcSDsp(dt = storage_sub, expr = list( as.double(sd(white, na.rm = TRUE)),
                                                   as.double(sd(p34, na.rm = TRUE)),
                                                   as.double(sd(pNo34, na.rm = TRUE)),
                                                   as.double(sd(pgt4, na.rm = TRUE)),
                                                   as.double(sd(plt3, na.rm = TRUE)),
                                                   as.double(sd(pgt3, na.rm = TRUE)),
                                                   as.double(sd(plt4, na.rm = TRUE)),
                                                   as.double(sd(pgt10, na.rm = TRUE)),
                                                   as.double(sd(p34gt10, na.rm = TRUE)) ) ) 

setnames(sb_SDsp, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

pdf(paste0(support_dir_name, "/QE_BUMP_SD.pdf", sep = ""), width = 12, height = 8)
for (i in 1:length(alphaMult)) {
  cv_tmp_m <- melt(copy(sb_SDsp[ i =  alphaMult_c == alphaMult[i]]), id = c(1:5))
  x <- ggplot(cv_tmp_m, aes(x = factor(sigPSmult_c), y = value)) + 
    geom_violin(fill = "slateblue") + 
    facet_grid(meanPS_c ~ variable) +
    ggtitle(paste("SD spawner abundance \n rows: survival \n columns: frequency bands") ) +
    ylab("SD spaawner abundance") +
    xlab(expression(paste(sigma))) +
    theme(axis.text.x = element_text(size = 7.5)) +
    scale_x_discrete(breaks=seq(0.1,0.5, by = 0.1)) + 
    theme_bw()
  print(x)
}
dev.off()

## PLOT OF SPAWNER ABUNDANCES - HIGH AND LOW SPECIFIED STANDARD DEVIATION
rm(spPlot, spPlot_m)
spPlot <- copy(storage[ i = N > 400 & sigPSmult_c %in% c(0.1, 0.2) & reps_c <= 20])
spPlot_m <- melt(spPlot, id = 1:6)

pdf(paste0(support_dir_name, "/spawnerTSplot_BUMP.pdf", sep = ""), width = 8, height = 6)
for (i in 1:length(freqCont)) {
  sub_sp_m <- droplevels(subset(spPlot_m, variable == freqCont[i]))
  x<- ggplot(data = sub_sp_m, aes(x = N, y = value, color = factor(sigPSmult_c))) + 
    facet_grid(meanPS_c ~ sigPSmult_c ) +
    geom_point(alpha = 0.05) +
    geom_line(alpha = 0.1) +
    geom_abline(intercept = 100, slope = 0) +
    scale_colour_manual(values = c("blue", "red")) +
    ylab("Spawner abundance") +
    xlab("Time (Years)") +
    guides(colour = guide_legend(title = "Specified \nSD"))+
    ggtitle(paste0("Spawner abundance: ", freqCont[i], "\n rows: survival \n columns: SD") ) +
    theme_bw()
  print(x)
}
dev.off()

## PLOT OF ONLY SUB-QE ABUNDANCES 
rm(spPlot, spPlot_m)

spPlot <- copy(storage[ i = N > 400 & sigPSmult_c %in% c(0.1, 0.2)])
spPlot_m <- as.data.table(melt(spPlot, id = 1:6))
spPlot_m <- spPlot_m[ i = value <= 100]
pdf(paste0(support_dir_name, "/subQEspawnerTimeSeries.pdf", sep = ""), width = 8, height = 8)
for (i in 1:length(freqCont)) {
  sub_sp_m <- droplevels(subset(spPlot_m, variable == freqCont[i]))
  x<- ggplot(data = sub_sp_m, aes(x = N, y = value, color = factor(sigPSmult_c))) + 
    facet_grid(meanPS_c ~ sigPSmult_c ) +
    geom_point(alpha = 0.75) +
    geom_abline(intercept = 100, slope = 0) +
    scale_colour_manual(values = c("blue", "red")) +
    ylab("Spawner abundance") +
    xlab("Time (Years)") +
    guides(colour = guide_legend(title = "Specified \nSD"))+
    ggtitle(paste0("Spawner abundance: ", freqCont[i], "\n rows: survival \n columns: SD") ) +
    theme_bw()
  print(x)
}
dev.off()
                 
## COMPARE a few HAPHAZARD selected COUNTS BY HAND TO JA_CONSECUTIVE RESULTS FOR SD 0.1 and 0.2 SD for plt4 at mean survival of 0.35
## WHY: there is a higher pQE at 0.1 sd than at 0.2 sd., but plots of spawner abundance show very few dips below 100 when
## surv sd == 0.1 but several at 0.2.
## PICKING OUTPUT 
## Either it's Jaime's function or the data.table wrapper that's not working
## First I'll test Jaime's JA_consec() on example series that supposedly went extinct
## to see if that's the problem, if not I have to test the data.table approach.

print(sb_qeyr[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1) & !is.na(plt4), 
               j = c(1:5, 12), with = FALSE], nrow = Inf)

# sb_qeyr "says" that the rep 3 goes extinct in year 41. 
sb_qeyr[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1, 0.2) &  reps_c == 3, 
         j = c(1:5, 12), with = FALSE]

plt4_3 <- storage[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1, 0.2) & reps_c == 3 & N > 400, 
                   j = c(1:6, 13), with = FALSE]
JA_consec(plt4_3[ i = sigPSmult_c == 0.1,
                  j = plt4])
JA_consec(plt4_3[ i = sigPSmult_c == 0.2,
                  j = plt4])

# sb_qeyr "says" that the rep 37 goes extinct in year 30 
sb_qeyr[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1, 0.2) &  reps_c == 37, 
         j = c(1:5, 12), with = FALSE]
plt4_37 <- storage[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1, 0.2) & reps_c == 37 & N > 400, 
                    j = c(1:6, 13), with = FALSE]
JA_consec(plt4_37[ i = sigPSmult_c == 0.1,
                   j = plt4])
JA_consec(plt4_37[ i = sigPSmult_c == 0.2,
                   j = plt4])


# sb_qeyr "says" that the rep 787 goes extinct in year 81
sb_qeyr[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1, 0.2) &  reps_c == 787, 
         j = c(1:5, 12), with = FALSE]        
plt4_787 <- storage[ i = meanPS_c == 0.35 & sigPSmult_c %in% c(0.1, 0.2) & reps_c == 787 & N > 400, 
                     j = c(1:6, 13), with = FALSE]
JA_consec(plt4_787[ i = sigPSmult_c == 0.1,
                    j = plt4])
JA_consec(plt4_787[ i = sigPSmult_c == 0.2,
                    j = plt4])


## THE SOLUTION--The problem was that QE events were triggered if there were only 1-3 dips below the QE threshold, 
# rather than 4 consecutive dips. When there were no dips below the QE threshold or >= 4 dips below the QE threshold 
# there were no problems. When the sd of environmental noise was about 0.2, that was the sweet spot which had the 
# erroneously high number of QE events - the bumps that you correctly saw no explanation.
## This fixed the pQE BUMP!
#### OLD --------------------------------- ####
## NONE OF THESE ACTUALL GO EXTINCT AND APPLYING JA_consec() by hand gives the correct answer
## CONCLUSION: test data.table wrapper approach. 
#### OLD --------------------------------- ####


# MOTIVATION analysis and plots


# Plot the stock recruitment curve. Declining survival story.

spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = 1)
alpha <- alphaMult * 1/spr
beta <- calc_beta(alpha, EQsp, spr)
spawners <- seq(0, max(EQsp)*1.5, len = 200)

recruits <- matrix(0, nrow = length(spawners), ncol = length(alpha))
alpha_names <- list()
for (i in 1:length(spawners)) {
  for (j in 1:length(alpha)) {
    recruits[i,j] <- alpha[j]*spawners[i]/(1 + beta[j]*spawners[i])
    alpha_names[j] <- paste("alpha_",toString(signif(alpha[j], 2)), sep = "")
  }
}
# make and melt df
sr <- data.frame(spawners, recruits)
names(sr)[2:(length(alpha)+1)] <- alpha_names
srm <- melt(sr, id.vars = "spawners")
srm$spEsc <- with(srm, value*surv1*surv2*delta_e + value*surv1*surv2*surv3*(1-delta_e))
names(srm)[2:3] <- c("Alpha", "recruits")

# SR Plot
pdf(paste0(diss_dir_name, "/stock_recruit.pdf", sep = ""), width = 8, height = 6)
ggplot(data = srm, aes(spawners, recruits) ) +
  geom_line(colour = "slateblue", size = 2) +
  theme_bw() +
  xlab("Spawners") +
  ylab("Recruits") +
  theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
  theme(legend.position = "right") +
  geom_abline(intercept = 0, slope = 1/spr, colour = "forestgreen", linetype = 2, size = 1.5) +
  geom_abline(intercept = 0, slope = ((1/spr)+alpha)/3, colour = "goldenrod", linetype = 2, size = 1.5) +
  geom_abline(intercept = 0, slope = ((1/spr)+alpha)/2, colour = "darkorange", linetype = 2, size = 1.5) +
  geom_abline(intercept = 0, slope = alpha, colour = "red", linetype = 2, size = 1.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme_bw()
dev.off()

## SUMMARY PLOT OF NOISE GENERATION AND FREQUENCY CONTENT noise, survival, and spawning female abundance

## FREQ CONTENT | TIME SERIES | WAVELET POWER SPECTRUM | GLOBAL WAVELET POWER SPECTRUM

## SUMMARY PLOT: noise
pdf(paste0(diss_dir_name, "/summaryFreqContTS_Noise.pdf", sep = ""), width = 8, height = 12)
summaryCh3tsPlotNoise(noise = noiseList, 
                      n = 1, 
                      J1 = trunc((log(32/(2 * 1))/log(2))/0.01))
dev.off()
## SUMMARY PLOT: truncated survival -- not stored currently
## SUMMARY PLOT: spawning female abundance

pdf(paste0(diss_dir_name, "/summaryFreqContTS_SpawningFemales.pdf", sep = ""), width = 8, height = 12)
summaryCh3tsPlotSpawners(spawners = storage, 
                         meanSurv = 0.35,
                         sigma = 0.1, 
                         n = 1, 
                         J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) 
dev.off()

## POPULATION (SPAWNER) VARIABILITY AS A FUNCTION OF FREQUENCY CONTENT

## SD 
SDsp <- melt(sb_SDsp[ i = sigPSmult_c %in% c(0.2, 0.4)], id = 1:5)

pdf(paste0(support_dir_name, "/SD_SpawningFemalesQEBUMP.pdf", sep = ""), width = 8, height = 8)
ggplot(data = SDsp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Standard Deviation Spawning Abundance") +
  theme_bw()
dev.off()

## COEFF VARIATION

CVsp <- melt(sb_CVsp[ i = sigPSmult_c %in% c(0.2, 0.4)], id = 1:5)

pdf(paste0(support_dir_name, "/CoeffVar_SpawningFemales_QEBUMP.pdf", sep = ""), width = 8, height = 8)
ggplot(data = CVsp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Coefficient of Variation of Spawning Abundance") +
  theme_bw()
dev.off()
## MEAN 

Meansp <- melt(sb_Meansp[ i = sigPSmult_c %in% c(0.2, 0.4)], id = 1:5)

pdf(paste0(support_dir_name, "/Mean_SpawningFemalesQEBUMP.pdf", sep = ""), width = 8, height = 8)
ggplot(data = Meansp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Spawning Abundance") +
  theme_bw()
dev.off() 


## FREQUENCY CONTENT on DISTRIBUTION OF QE TIME ####
# Different distributions of QE times as you move down to the equilbrium abundance_collapse (meanPS decreasing)
# for the survival rates with different frequency contents
# check out how initial abundance is specified for logNormal runs
# Run for 100 years at different specifed levels of sigPSmult

# ONE SLOPE 4

# Three SURVIVALS seq(0.275, 0.5, 0.8)

# Multiple Noise

# calculate the fraction spawning early
delta_e <- calc_de(wanted_frac, surv3)
# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = 1)  
alphaMult <- 4 # multiplier of replacement line slope to set alpha [from 2-10]
alpha <- alphaMult * 1/spr
# sps_crash <- 1/{alpha*surv1*surv2*{delta_e + surv3*{1 - delta_e}}}

meanPS <- c(0.275, 0.5, 0.8)
sigPSmult <- seq(0.1, 0.5, by = 0.05) # c(0.2, 0.4) 
EQsp <- 7500

meanPS_r <- rep(meanPS, each = length(sigPSmult)*length(alphaMult)*length(EQsp)*reps*N)
sigPSmult_r <- rep(sigPSmult, each = length(alphaMult)*length(EQsp)*reps*N, times = length(meanPS))
alphaMult_r <- rep(alphaMult, each = length(EQsp)*reps*N, times = length(meanPS)*length(sigPSmult))
EQsp_r <- rep(EQsp, each = reps*N, times = length(meanPS)*length(sigPSmult)*length(alphaMult))
reps_r <- rep(1:reps, each = N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp))
n_r <- rep(1:N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp)*reps)
rm(storage, storageP, french)
storageP <- data.table(meanPS_c = meanPS_r,
                       sigPSmult_c = sigPSmult_r,
                       alphaMult_c = alphaMult_r,
                       EQsp_c = EQsp_r,
                       reps_c = reps_r,
                       N = n_r,
                       white = 0,
                       p34 = 0,
                       pNo34 = 0,
                       pgt4 = 0,
                       plt3 = 0,
                       pgt3 = 0,
                       plt4 = 0,
                       pgt10 = 0,
                       p34gt10 = 0) 

# setting the key for "fast" indexing of data.table storage
setkey(storageP, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)

# split storageP by meanPS so each component of the
# then I'll put all the lists back together at the end of the simulation
french <- vector("list", length = length(meanPS)) # random name, another version I used italian named in honor of LWB

for (i in 1:length(meanPS)) {
  french[[i]] <- storageP[i = meanPS_c == meanPS[i]]
}

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

# A function wrapper format for data.tables

# qeLev <- 100
# run_length <- 4
# Call the function
storage_sub <- copy(storage[ i = N > 400])

sb_qeyr <- calcQEyr(dt = storage_sub, expr = list( as.integer(JA_consec(white, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(p34, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pNo34, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pgt4, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(plt3, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pgt3, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(plt4, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pgt10, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(p34gt10, run_length = 4, qeLev = 100)) ) )

# need to name the new columns
setnames(sb_qeyr, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))


## PLOT QE TIME DISTRIBUTIONS

# SURVIVAL RATE
pdf(paste0(diss_dir_name, "/Surv_Freq_QE_time_Dist_allSurv.pdf", sep = ""), width = 12, height =12)
rm(qet_tmp_sb_m)
qet_tmp_sb_m <- as.data.table(melt(copy(sb_qeyr[ i = sigPSmult_c == 0.4 & N > 400]), id = c(1:5)))
x <- ggplot(qet_tmp_sb_m, aes(x = value)) + 
  geom_histogram() + 
  facet_grid(variable ~ meanPS_c) +
  xlab("Time (Year)") +
  ylab("Count") +
  theme_bw()
print(x)

dev.off()

pdf(paste0(diss_dir_name, "/Surv_Freq_QE_time_Dist_lowSurv.pdf", sep = ""), width = 5, height = 8)
rm(qet_tmp_sb_m)
qet_tmp_sb_m <- as.data.table(melt(copy(sb_qeyr[ i = sigPSmult_c == 0.4 & N > 400 & meanPS_c == 0.275]), id = c(1:5)))
x <- ggplot(qet_tmp_sb_m, aes(x = value)) + 
  geom_histogram() + 
  facet_grid(variable ~ meanPS_c) +
  xlab("Time (Year)") +
  ylab("Count") +
  theme_bw()
print(x)

dev.off()


# VARIABILTY 


pdf(paste0(support_dir_name, "/Sigma_Freq_QE_time_Dist.pdf", sep = ""), width = 12, height = 12)
rm(qet_tmp_sb_m)
qet_tmp_sb_m <- as.data.table(melt(copy(sb_qeyr[ i = meanPS_c == 0.275 & N > 400]), id = c(1:5)))
x <- ggplot(qet_tmp_sb_m, aes(x = value)) + 
  geom_histogram() + 
  facet_grid(variable ~ sigPSmult_c) +
  xlab("Time (Year)") +
  ylab("Count") +
  theme_bw()

print(x)

dev.off()

## pQE,  ####

sb_pQE <- calc_pQE(sb_qeyr, list(white_qe = length(which(!is.na(white)))/1000,
                                 p34_qe = length(which(!is.na(p34)))/1000,
                                 pNo34_qe = length(which(!is.na(pNo34)))/1000,
                                 pgt4_qe = length(which(!is.na(pgt4)))/1000,
                                 plt3_qe = length(which(!is.na(plt3)))/1000,
                                 pgt3_qe = length(which(!is.na(pgt3)))/1000,
                                 plt4_qe = length(which(!is.na(plt4)))/1000,
                                 pgt10_qe = length(which(!is.na(pgt10)))/1000,
                                 p34gt10_qe = length(which(!is.na(p34gt10)))/1000 )  )

pdf(paste0(diss_dir_name, "/sigma_vs_pQE.pdf", sep = ""), width = 12, height = 12)
for (i in 1:length(alphaMult)) {
  pqe_tmp_m <- melt(copy(sb_pQE[ i = alphaMult_c == alphaMult[i]]), id = c(1:4))
  x <- ggplot(pqe_tmp_m, aes(x = sigPSmult_c, y = value)) + 
    geom_line(colour = "slateblue", size = 1.2) + 
    facet_grid(meanPS_c ~ variable) +
    theme(axis.text = element_text(size = 8, colour = "black")) +
    ggtitle(paste("p(Quasi-Extinction) \n rows: survival \n columns: frequency bands") ) +
    ylab("Probability of Quasi-Extinction") +
    xlab(expression(paste(sigma))) +
    theme_bw()
  print(x)
}
dev.off()

pdf(paste0(diss_dir_name, "/sigma_vs_pQE2row.pdf", sep = ""), width = 12, height = 9)
for (i in 1:length(alphaMult)) {
  pqe_tmp_m <- melt(copy(sb_pQE[ i = alphaMult_c == alphaMult[i] & meanPS_c %in% c(meanPS[1], meanPS[2])]), id = c(1:4))
  x <- ggplot(pqe_tmp_m, aes(x = sigPSmult_c, y = value)) + 
    geom_line(colour = "slateblue", size = 1.2) + 
    facet_grid(meanPS_c ~ variable) +
    theme(axis.text = element_text(size = 8, colour = "black")) +
    # ggtitle(paste("p(Quasi-Extinction) \n rows: survival \n columns: frequency bands") ) +
    ylab("Probability of Quasi-Extinction") +
    xlab(expression(paste(sigma))) +
    theme_bw()
  print(x)
}
dev.off()


## SD 

sb_SDsp <- calcSDsp(dt = storage_sub, expr = list( as.double(sd(white, na.rm = TRUE)),
                                                   as.double(sd(p34, na.rm = TRUE)),
                                                   as.double(sd(pNo34, na.rm = TRUE)),
                                                   as.double(sd(pgt4, na.rm = TRUE)),
                                                   as.double(sd(plt3, na.rm = TRUE)),
                                                   as.double(sd(pgt3, na.rm = TRUE)),
                                                   as.double(sd(plt4, na.rm = TRUE)),
                                                   as.double(sd(pgt10, na.rm = TRUE)),
                                                   as.double(sd(p34gt10, na.rm = TRUE)) ) ) 

setnames(sb_SDsp, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

SDsp <- melt(sb_SDsp[ i = meanPS_c %in% c(0.275, 0.5) & sigPSmult_c %in% c(0.2, 0.4)], id = 1:5)

pdf(paste0(diss_dir_name, "/SD_SpawningFemales.pdf", sep = ""), width = 8, height = 8)
ggplot(data = SDsp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Standard Deviation Spawning Abundance") +
  theme_bw()
dev.off()

SDsp <- melt(sb_SDsp[ i = meanPS_c == 0.275 & sigPSmult_c == 0.4], id = 1:5)

pdf(paste0(diss_dir_name, "/SD_SpawningFemales1plot.pdf", sep = ""), width = 8, height = 8)
ggplot(data = SDsp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Standard Deviation Spawning Abundance") +
  theme_bw()
dev.off()


## COEFF VARIATION

sb_CVsp <- calcCVsp(dt = storage_sub, expr = list( as.double(sd(white, na.rm = TRUE)/mean(white, na.rm = TRUE)),
                                                   as.double(sd(p34, na.rm = TRUE)/mean(p34, na.rm = TRUE)),
                                                   as.double(sd(pNo34, na.rm = TRUE)/mean(pNo34, na.rm = TRUE)),
                                                   as.double(sd(pgt4, na.rm = TRUE)/mean(pgt4, na.rm = TRUE)),
                                                   as.double(sd(plt3, na.rm = TRUE)/mean(plt3, na.rm = TRUE)),
                                                   as.double(sd(pgt3, na.rm = TRUE)/mean(pgt3, na.rm = TRUE)),
                                                   as.double(sd(plt4, na.rm = TRUE)/mean(plt4, na.rm = TRUE)),
                                                   as.double(sd(pgt10, na.rm = TRUE)/mean(pgt10, na.rm = TRUE)),
                                                   as.double(sd(p34gt10, na.rm = TRUE)/mean(p34gt10, na.rm = TRUE)) ) ) 

setnames(sb_CVsp, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

CVsp <- melt(sb_CVsp[ i = meanPS_c %in% c(0.275, 0.5) & sigPSmult_c %in% c(0.2, 0.4)], id = 1:5)

pdf(paste0(diss_dir_name, "/CoeffVar_SpawningFemales.pdf", sep = ""), width = 8, height = 8)
ggplot(data = CVsp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Coefficient of Variation of Spawning Abundance") +
  theme_bw()
dev.off()

CVsp <- melt(sb_CVsp[ i = meanPS_c == 0.275 & sigPSmult_c == 0.4], id = 1:5)

pdf(paste0(diss_dir_name, "/CoeffVar_SpawningFemales1plot.pdf", sep = ""), width = 8, height = 8)
ggplot(data = CVsp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Coefficient of Variation of Spawning Abundance") +
  theme_bw()
dev.off()



## MEAN 

sb_Meansp <- calcMeansp(dt = storage_sub, expr = list( as.double(mean(white, na.rm = TRUE)),
                                                       as.double(mean(p34, na.rm = TRUE)),
                                                       as.double(mean(pNo34, na.rm = TRUE)),
                                                       as.double(mean(pgt4, na.rm = TRUE)),
                                                       as.double(mean(plt3, na.rm = TRUE)),
                                                       as.double(mean(pgt3, na.rm = TRUE)),
                                                       as.double(mean(plt4, na.rm = TRUE)),
                                                       as.double(mean(pgt10, na.rm = TRUE)),
                                                       as.double(mean(p34gt10, na.rm = TRUE)) ) ) 

setnames(sb_Meansp, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

Meansp <- melt(sb_Meansp[ i = meanPS_c %in% c(0.275, 0.5) & sigPSmult_c %in% c(0.2, 0.4)], id = 1:5)

pdf(paste0(diss_dir_name, "/Mean_SpawningFemales.pdf", sep = ""), width = 8, height = 8)
ggplot(data = Meansp, aes(x = variable, y = value)) +
  facet_grid( meanPS_c ~ sigPSmult_c  ) +
  geom_boxplot() +
  xlab("") +
  ylab("Mean Spawning Abundance") +
  theme_bw()
dev.off() 



## POPULATION VARIABILTY PLOTTED AS A FUNCTION OF SLOPES

# DO SECOND SIMULATION WITH alphaMult of 2
alphaMult <- 2 # multiplier of replacement line slope to set alpha 
alpha <- alphaMult * 1/spr
sps_crash <- 1/{alpha*surv1*surv2*{delta_e + surv3*{1 - delta_e}}}
sigPSmult <- c(0.2, 0.4)
meanPS2 <- seq(0.55, 0.675, 0.025)

meanPS_r <- rep(meanPS2, each = length(sigPSmult)*length(alphaMult)*length(EQsp)*reps*N)
sigPSmult_r <- rep(sigPSmult, each = length(alphaMult)*length(EQsp)*reps*N, times = length(meanPS2))
alphaMult_r <- rep(alphaMult, each = length(EQsp)*reps*N, times = length(meanPS2)*length(sigPSmult))
EQsp_r <- rep(EQsp, each = reps*N, times = length(meanPS2)*length(sigPSmult)*length(alphaMult))
reps_r <- rep(1:reps, each = N, times = length(meanPS2)*length(sigPSmult)*length(alphaMult)*length(EQsp))
n_r <- rep(1:N, times = length(meanPS2)*length(sigPSmult)*length(alphaMult)*length(EQsp)*reps)
rm(storageP, french)
storageQ <- data.table(meanPS_c = meanPS_r,
                       sigPSmult_c = sigPSmult_r,
                       alphaMult_c = alphaMult_r,
                       EQsp_c = EQsp_r,
                       reps_c = reps_r,
                       N = n_r,
                       white = 0,
                       p34 = 0,
                       pNo34 = 0,
                       pgt4 = 0,
                       plt3 = 0,
                       pgt3 = 0,
                       plt4 = 0,
                       pgt10 = 0,
                       p34gt10 = 0) 

# setting the key for "fast" indexing of data.table storage
setkey(storageQ, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)

# split storageQ by meanPS so each component of the
# then I'll put all the lists back together at the end of the simulation
italian <- vector("list", length = length(meanPS2)) # random name, another version I used italian named in honor of LWB

for (i in 1:length(meanPS2)) {
  italian[[i]] <- storageQ[i = meanPS_c == meanPS2[i]]
}

# run simulations using foreach framework to send jobs to multiple cores
# for each chunk of storageP - this should cut simulation times down to ~ 68/9 = 7.5 hours

# load foreach package for parallel processing
library(foreach)
# set up backend to do parallel processing
library(doParallel)
# detectCores()
registerDoParallel() # defaults to half of the cores 

system.time(
  storage2 <- foreach(h = 1:length(meanPS2),.packages="reshape2") %dopar% {
    parSimCmp(italian[[h]], simLen = 100, survMean = meanPS2[h])
  }
)

storage2 <- rbindlist(storage2)

rm(storageQ, italian)

## FOR EACH SLOPE (2 & 4 x 1/spr) plot SD, CV, and pQE 
## PQE

## FIRST TALLY QE EVENTS 

storage2_sub <- copy(storage2[ i = N > 400])
sb_qeyr2 <- calcQEyr(dt = storage2_sub, expr = list( as.integer(JA_consec(white, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(p34, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pNo34, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pgt4, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(plt3, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pgt3, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(plt4, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(pgt10, run_length = 4, qeLev = 100)),
                                               as.integer(JA_consec(p34gt10, run_length = 4, qeLev = 100)) ) )
# need to name the new columns
setnames(sb_qeyr2, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         c("white", "p34", "pNo34", "pgt4", "plt3", "pgt3", "plt4", "pgt10", "p34gt10"))

## PLOT QE TIME AS A FUNCTION OF SLOPE

## SINCE EACH SLOPE HAS DIFFERENT SURVIVAL RATES,
## MANUALLY SPLIT PAGE FOR PLOTTING TOP
# vp4 <- viewport(x = 0.5, y = 0.75, width = 1, height = 0.5)
# vp2 <- viewport(x = 0.5, y = 0.25, width = 1, height = 0.5)
# 
# tmp4 <- as.data.table(copy(sb_qeyr[ i = sigPSmult_c == 0.4]))
# tmp2 <- as.data.table(copy(sb_qeyr2[ i = sigPSmult_c == 0.4]))
# 
# pdf(paste0(support_dir_name, "/QE_time_plots_Slope_Freq_cont.pdf", sep = ""), width = 8, height = 8) 
# for ( i in 1:length(freqCont)) {
#   dt4 <- copy(tmp4[ j = c(1:5, 5 + i), with = FALSE])
#   names(dt4)[6] <- "bandwidth"
#   
#   dt2 <- copy(tmp2[ j = c(1:5, 5 + i), with = FALSE])
#   names(dt2)[6] <- "bandwidth"
#   
#   p4 <- ggplot(dt4, aes(x = bandwidth)) +
#     geom_histogram() +
#     facet_grid(meanPS_c ~ .) + 
#     xlab("Year") +
#     ylab("Count")
#     
#   p2 <- ggplot(dt2, aes(x = bandwidth)) +
#     geom_histogram() +
#     facet_grid(meanPS_c ~ .) + 
#     xlab("Year") +
#     ylab("Count")
#   
#   print(p4, vp = vp4)
#   print(p2, vp = vp2)
# }
# dev.off()
# ## The other approach is to plot all freqContent for each survival level on each page
# 
# ## 
# ## melt tmp4 and tmp2
# # meanPS <- seq(0.275, 0.4, 0.025)
# # meanPS2 <- seq(0.5, 0.625, 0.025)
#   tmp4 <- melt(sb_qeyr[i = sigPSmult_c == 0.4 & meanPS_c %in% c(meanPS[1], meanPS[length(meanPS)]) ], id = c(1:5))
#   tmp2 <- melt(sb_qeyr2[i = sigPSmult_c == 0.4 & meanPS_c %in% c(meanPS2[1], meanPS2[length(meanPS2)]) ], id = c(1:5))
# ## 
# pdf(paste0(support_dir_name, "/QE_time_plots_SLOPE_Freq_pg.pdf", sep = ""), width = 4.5, height = 8) 
#   q4 <- ggplot(tmp4, aes(x = value)) +
#   geom_histogram(fill = "black") +
#   facet_grid(variable ~ meanPS_c) + 
#   xlab("") + 
#   ylab("Count") +
#   ggtitle(paste0("SR slope: ", unique(tmp4$alphaMult))) + 
#   theme(strip.text.y = element_text(angle = 0, size = 10))
#   
#   q2 <- ggplot(tmp2, aes(x = value)) +
#   geom_histogram(fill = "black") +
#   facet_grid(variable ~ meanPS_c) + 
#   xlab("") + 
#   ylab("Count") +
#   ggtitle(paste0("SR slope: ", unique(tmp2$alphaMult))) + 
#   theme(strip.text.y = element_text(angle = 0, size = 10))
#   
# 
#   print(q4, vp = vp4)
#   print(q2, vp = vp2)
#   dev.off()

## THEN CALC PQE
sb_pQE2 <- calc_pQE(sb_qeyr2, list(white_qe = length(which(!is.na(white)))/1000,
                                 p34_qe = length(which(!is.na(p34)))/1000,
                                 pNo34_qe = length(which(!is.na(pNo34)))/1000,
                                 pgt4_qe = length(which(!is.na(pgt4)))/1000,
                                 plt3_qe = length(which(!is.na(plt3)))/1000,
                                 pgt3_qe = length(which(!is.na(pgt3)))/1000,
                                 plt4_qe = length(which(!is.na(plt4)))/1000,
                                 pgt10_qe = length(which(!is.na(pgt10)))/1000,
                                 p34gt10_qe = length(which(!is.na(p34gt10)))/1000 )  )

pQE4_m <- melt(sb_pQE[ i = sigPSmult_c == 0.2], id = 1:4)
pQE2_m <- melt(sb_pQE2[ i = sigPSmult_c == 0.2], id = 1:4)

vp4 <- viewport(x = 0.5, y = 0.75, width = 1, height = 0.5)
vp2 <- viewport(x = 0.5, y = 0.25, width = 1, height = 0.5)

pdf(paste0(support_dir_name, "/pQE_w_slope_by_Freq.pdf", sep = ""), width = 4, height = 8) 

p4 <- ggplot(data = pQE4_m, aes(x = meanPS_c, y = value)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ . ) + 
  xlab("Pre-spawn survival") + 
  ylab("Probability of Quasi-Extinction") + 
  ggtitle(paste0("SR slope: ", unique(pQE4_m$alphaMult))) + 
  theme(strip.text.y = element_text(angle = 0, size = 10)) +
  theme_bw()

p2 <- ggplot(data = pQE2_m, aes(x = meanPS_c, y = value)) +
  geom_line() +
  geom_point() +
  facet_grid(variable ~ . ) + 
  xlab("Pre-spawn survival") + 
  ylab("Probability of Quasi-Extinction") + 
  ggtitle(paste0("SR slope: ", unique(pQE2_m$alphaMult))) + 
  theme(strip.text.y = element_text(angle = 0, size = 10)) +
  theme_bw()

  print(p4, vp = vp4)
  print(p2, vp = vp2)
  dev.off()

## END OF SCRIPT

fileName <- paste(diss_dir_name, "/chap_3_part_1.RData", sep = "")

save.image(file = fileName)
sendmail("dpkilduff@yahoo.com", "chinook simulation", "It (chap_3_part_1.R) be done or it crashed!")