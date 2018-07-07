# main.R
# Main file to call all code for chapter 3 of dissertation.


# Load packages and project functions
library(reshape2)
library(ggplot2)
library(plyr)
library(xtable)
library(ggthemes)
library(biwavelet)
library(data.table)
library(RSEIS)

source("./code/functions.R")

mgmtScenarios <- c("BAU", "NoDiversion", "ColdWater") # BAU, mgmt option 1, and mgmt option 2;
                                                      #missing management combo - same as NoDiversion
                #c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion",
                #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
climateScenarios <- c("A2","B1")
gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")

# Extract and prepare SALMOD data ####
SALMODextract()
SALMODsumSpawnJuv()
mergeSALMODextract()
spawnFemSurvSALMOD()
# bootSRsurvSALMOD()
egg2fry2outSurv()

SALMODfecundity()
SALMODrelMort()

# Exploratory plots and analyses ####
plotEggsPerSpawner()
plotMortEggsFryDD()
plotMortEggsFryTemp()
plotPrespawnMort() # plot pre-spawn mortality for each water management and climate scenario

# Population modeling ####

# Age-structured model projections, plots, and extinction risk calculations
# assuming constant early ocean survival 

## EQsp 7500 and 3000 going to zero.

# Analysis outline ####

# Compare deterministic population model on time to extinction
  # QE levels: 0, 20, 100
  # Influence of alpha
  # Influence of spawning capacity (i.e., alpha/beta) and selected QE threshold
  # So what?: Simplest way to show that population dynamics influences time to extinction calculations

# Stochastic model frequency response
  # High and low survival, SR slope at origin (alpha)
  # Spawners and recruits
  # Influence of white, and other colors of noise on Frequency responses
    # So what?: If we know what the frequency response of the population model then we understand how 
    # noise will increase variability of output (spawners - the metric of interest)

# Characterize ocean environmental variability (Wavelet - PDO, NPGO, ENSO) -- use Global power spectrum to generate TS (unit variability)
  # - periods of greatest variability? (Scaled WPS)
  # 

# Characterize FW variability of survival (Wavelet) -- Use Global power spectrum to generate TS (unit variability)
  # - periods of greatest variability? (Scaled WPS) - consistency within/across model output (cross wavelet spectra)

    # call functions for mass production of time series diagnostic and analysis plots
    # for each wtrMgmt, emmission, climate model and FW life stage (prespawn, egg fry)

salmod_surv_dp() # makes several of the diagnostic plots in Zuur 2010 for each 
                 # wtrMgmt, emmission, climate model and FW life stage (prespawn, egg fry)

wavelet_SALMOD_series() # Plots ts, wavelet power spectrum, plus the frequency-averaged Global WPS and
                        # the time-averaged Scaled WPS - need to check on how to treat/scale magnitudes 

# Deterministic simulations

# Case1: uses baseline parameters for survival rates and the fraction of early spawners [results in 50% of 
# age3 and age4 spawners] to investigate the influence of varying the slope of the BevHolt SR at the origin,
# alpha, for given number of spawning females (7,500)

# The BAU_A1 
source("./code/functions.R")
case1A2 <- popProjConstMarine(surv1 = 0.02, 
                            surv2 = 0.8, 
                            surv3 = 0.8,
                            survPS = "vary", # 1 or "vary"
                            wanted_frac = 0.5,
                            alpha_vec = c(1.5, 3, 5, 10),
                            EQsp = c(7500),
                            wtrMgmt = "BAU", # can be "" for test case of survPS = 1
                            climate = "A2", # can be "" for test case of survPS = 1
                            yrs = 2010:2099,
                            make_plot = FALSE) 

#  Time series plots - influence of changing alpha on time to extinciton
plotPopProjConst(droplevels(subset(case1A2[[1]], age >= 3)), case1A2[[2]], case1A2[[3]], vary = "alpha")

# Calculate time to Quasi extinction (last of four years)

out1A2qe0 <- calcQEtimeConstPers(case1A2, mgmt = "BAU", climate = "A2",QEthr = 0)
out1A2qe20 <- calcQEtimeConstPers(case1A2, mgmt = "BAU", climate = "A2",QEthr = 20)
out1A2qe100 <- calcQEtimeConstPers(case1A2, mgmt = "BAU", climate = "A2",QEthr = 100)

# operations to get qe output in table form
# xtable(out1A2qe[which(out1A2qe$SRinfo == 1), 2:4])
# colMeans(out1A2qe[which(out1A2qe$SRinfo == 1), 3:4])
# colMeans(dcast(out1A2qe, gcm ~ SRinfo, value.var = "yrs2ext")[,-1])

# This is gener
case1B1 <- popProjConstMarine(surv1 = 0.02, 
                              surv2 = 0.8, 
                              surv3 = 0.8,
                              survPS = "vary", # 1 or "vary"
                              wanted_frac = 0.5,
                              alpha_vec = c(1.5, 3, 5, 10),
                              EQsp = c(7500),
                              wtrMgmt = "NoD", # can be "" for test case of survPS = 1
                              climate = "B1", # can be "" for test case of survPS = 1
                              yrs = 2010:2099,
                              make_plot = FALSE) 
# calc_de(surv3=0.8, wanted_frac=0.5, plot_wanted = "yes")
plotPopProjConst(droplevels(subset(case1B1[[1]], age >= 3)), case1B1[[2]], case1B1[[3]], vary = "alpha")

# Calculate time to Quasi extinction (last of four years)
out1B1qe0 <- calcQEtimeConstPers(case1B1, mgmt = "NoDiversion", climate = "B1",QEthr = 0)
out1B1qe20 <- calcQEtimeConstPers(case1B1, mgmt = "NoDiversion", climate = "B1",QEthr = 20)
out1B1qe100 <- calcQEtimeConstPers(case1B1, mgmt = "NoDiversion", climate = "B1",QEthr = 100)


# Case2: uses baseline parameters for survival rates and the fraction of early spawners [results in 50% of 
# age3 and age4 spawners] to investigate the influence of varying the carry capacity of the system EQsp
# which is alpha/beta scaled by survival and maturation rates to produce a deterministic number of spawners 
case2A2 <- popProjConstMarine(surv1 = 0.2, 
                            surv2 = 0.8, 
                            surv3 = 0.8,
                            survPS = "vary", # 1 or "vary"
                            wanted_frac = 0.5,
                            alpha_vec = c(5),
                            EQsp = c(2500, 7500, 30000),
                            wtrMgmt = "BAU", # can be "" for test case of survPS = 1
                            climate = "A2", # can be "" for test case of survPS = 1
                            yrs = 2010:2099,
                            make_plot = FALSE) # not used in this run

plotPopProjConst(droplevels(subset(case2A2[[1]], age >= 3)), case2A2[[2]], case2A2[[3]], vary = "EQsp")

# Calculate time to Quasi extinction (last of four years)
out2A2qe0 <- calcQEtimeConstPers(case2A2, mgmt = "BAU", climate = "A2",QEthr = 0)
out2A2qe20 <- calcQEtimeConstPers(case2A2, mgmt = "BAU", climate = "A2",QEthr = 20)
out2A2qe100 <- calcQEtimeConstPers(case2A2, mgmt = "BAU", climate = "A2",QEthr = 100)


case2B1 <- popProjConstMarine(surv1 = 0.2, 
                              surv2 = 0.8, 
                              surv3 = 0.8,
                              survPS = "vary", # 1 or "vary"
                              wanted_frac = 0.5,
                              alpha_vec = c(5),
                              EQsp = c(2500, 7500, 30000),
                              wtrMgmt = "NoD", # can be "" for test case of survPS = 1
                              climate = "B1", # can be "" for test case of survPS = 1
                              yrs = 2010:2099,
                              make_plot = FALSE) # not used in this run

plotPopProjConst(droplevels(subset(case2B1[[1]], age >= 3)), case2B1[[2]], case2B1[[3]], vary = "EQsp")

# Calculate time to Quasi extinction (last of four years)
out2B1qe0 <- calcQEtimeConstPers(case2B1, mgmt = "BAU", climate = "A2", QEthr = 0)
out2B1qe20 <- calcQEtimeConstPers(case2B1, mgmt = "BAU", climate = "A2", QEthr = 20)
out2B1qe100 <- calcQEtimeConstPers(case2B1, mgmt = "BAU", climate = "A2", QEthr = 100)

# Calculate statistics for write-up


# Stochastic simulations #### 

# for each type of noise
# calculate the probability of extinction (mean time to extinction) as a function of variability 
# for alpha ranging from 2-10 times 1/SPR
# for EQsp ranging from 1000 to 50000 
# for mean survivals 0.1 to 0.8
# for QEthr 20 and 100 at 50 and 1000. (lowest count 14)
# 

r_seed <- 123

# Psuedo-code for stochastic simulations

# Model parameters
set.seed(r_seed)
freqCont <- c("white", "rsw34", "rsw34no", "rswgt3", "rswlt4", "rswgt4", "rswlt3", "rswgt10", "rsw34gt10")
meanPS <- seq(0.1, 0.9, by = 0.1)
# Notes on "observed" meanPS from the SALMOD output files
# Showing these values to justify ranges used in simulations
# Early = first 33 years, Mid = Middle 33 years, and Late = final 34 years
# BAU and A2: mean: Early max = 0.6, min = 0.34  
#                   Mid max = 0.51, min = 0.18
#                   Late max = 0.23, min = 0.01
# BAU and A2: sd: Early max = 0.34, min = 0.25
#                   Mid max = 0.34, min = 0.2
#                   Late max = 0.24., min = 0.01
# NoD and B1: mean: Early max = 0.92, min = 0.59 
#                   Mid max = 0.81, min = 0.32
#                   Late max = 0.65, min = 0.19
# NoD and B1: sd: Early max = 0.38, min = 0.16
#                   Mid max = 0.42, min = 0.31
#                   Late max = 0.41, min = 0.27
sigPSmult <- seq(.1, 0.5, by = 0.1) # seq(.1, 0.5, by = 0.05)
alphaMult <- 2:10 # changing meanPS also changes SPR, so don't need separate SPR values
EQsp <-  c(7500) # seq(2000, 20000, by = 2000)
# QEthr <- c(20, 100) # don't use
reps <- 50 # 1000
N <- 100 
surv1 <- 0.02
surv2 <- surv3 <- 0.8
wanted_frac <- 0.5

# Simulate survival time series

# create "reps" number of white noise time series of length N using random sine wave approach
white_n <- customFR2ts(N = N, # number of time steps
                       reps = reps,
                       r_seed = r_seed,
                       amp = mk_white(N)) # mean = 0, variance/sd = 1

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
# plotMeanFreqR(rsin_34_gt10_n, N=N)
# plot(1:N, rsin_gt10_n[,1], type = "l")
# lines(1:N, rsin_34_n[,1], type = "l", col = "red")
# lines(1:N, (rsin_gt10_n[,1] + rsin_34_n[,1])/2, type = "l", col = "slateblue", lwd = 2)

# For test runs
# noiseList <- list(noise_white = white_n,
#                   noise_34 = rsin_34_n)
# For Full run
noiseList <- list(noise_white = white_n,
                  noise_34 = rsin_34_n,
                  nose_34_reject = rsin_34_reject_n,
                  noise_gt4 = rsin_gt4_n,
                  noise_lt3 = rsin_lt3_n,
                  noise_gt3 = rsin_gt3_n,
                  noise_lt4 = rsin_lt4_n,
                  noise_gt10 = rsin_gt10_n,
                  noise_34gt10 = rsin_34_gt10_n)
rm(white_n, rsin_34_n, rsin_34_reject_n, rsin_gt4_n, rsin_lt3_n, rsin_gt3_n, rsin_lt4_n, rsin_gt10_n, rsin_34_gt10_n)
# array for storing spawner data
# spawnOutArray <- array(NA, c(length(freqCont), length(meanPS), length(sigPSmult), length(alphaMult),
#                              length(EQsp), reps, N))
# dimnames(spawnOutArray) <- list(freqCont, paste("meanPS_", meanPS, sep = ""), 
#                                 paste("sigPSmult_", sigPSmult, sep = ""),
#                                 paste("alphaMult_", alphaMult, sep = ""),
#                                 paste("EQsp_", EQsp, sep = ""),
#                                 paste("r_", 1:reps, sep = ""),
#                                 paste("yr_", 1:N, sep = ""))

# array for storing abundance data at each age
# abundOutArray <- array(NA, c(length(freqCont), length(meanPS), length(sigPSmult), length(alphaMult),
#                             length(EQsp), reps, N, 4))

# try a data.table solution

# testing params
meanPS <-  c(0.1, 0.3, 0.5) # seq(0.1, 0.9, by = 0.1)
M <- -log(meanPS)
sigPSmult <- c(0.1, 0.5) # seq(0.1, 0.5, by = 0.1) #c(0.1, 0.5)
alphaMult <- c(2, 5) # 2:10 # c(2, 5)
EQsp <- 7500
reps <- 50
N <- 100
# i <- j <- k <- l <- m <- 1
meanPS_r <- rep(meanPS, each = length(sigPSmult)*length(alphaMult)*length(EQsp)*reps*N)
sigPSmult_r <- rep(sigPSmult, each = length(alphaMult)*length(EQsp)*reps*N, times = length(meanPS))
alphaMult_r <- rep(alphaMult, each = length(EQsp)*reps*N, times = length(meanPS)*length(sigPSmult))
EQsp_r <- rep(EQsp, each = reps*N, times = length(meanPS)*length(sigPSmult)*length(alphaMult))
reps_r <- rep(1:reps, each = N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp))
n_r <- rep(1:N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp)*reps)

# Since the simulation output is large, I opted to use data.tables  
# in hopes of speeding up simulation time

storage <- data.table(meanPS_c = meanPS_r,
                       sigPSmult_c = sigPSmult_r,
                       alphaMult_c = alphaMult_r,
                       EQsp_c = EQsp_r,
                       reps_c = reps_r,
                       N = n_r,
                      white = 0,
                      p34 = 0,
                      pNo34 = 0,
                      pgt3 = 0,
                      plt4 = 0,
                      pgt4 = 0,
                      plt3 = 0,
                      pgt10 = 0,
                      p34gt10 = 0) 

# setting the key for "fast" indexing of data.table storage
setkey(storage, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)

# colsRefs <- c("white", "p34", "pNo34", "pgt3", ""
# plt4 = 0,
# pgt4 = 0,
# plt3 = 0,
# pgt10 = 0,
# p34gt10 = 0)


for (i in 1:length(freqCont)) {
  for (j in 1:length(meanPS)) {
    for (k in 1:length(sigPSmult)) {
      for (l in 1:length(alphaMult)) {
        for (m in 1:length(EQsp)) {
          # noise to survival: surv <- noise * sigPSmult + meanPS 
          surv <- noiseList[[i]] * sigPSmult[k] + meanPS[j]
          # surv <- exp(-(noiseList[[i]]* sigPSmult[k]  + M[j])))
          surv[surv > 1] <- 1
          surv[surv < 0] <- 0
          # truncate survival: surv[surv > 1] <- 1; surv[surv < 0] <- 0
          # improved survival truncation (control variance - 0's and 1's)
          # run population model
          out <- melt(popSimPSvary(rand_surv = surv, 
                                 surv1 = surv1, 
                                 surv2 = surv2, 
                                 surv3 = surv3, 
                                 EQsp = EQsp[m], 
                                 wanted_frac = wanted_frac,
                                 alpha_scale = alphaMult[l])[[2]])
          #storage[meanPS[j], sigPSmult[k], alphaMult[l], EQsp[m], white := out]
          #colnames(modOut[[2]]) <- paste("r",1:reps,sep="")
          storage[i = (i = (meanPS_c == meanPS[j] & sigPSmult_c == sigPSmult[k] & alphaMult_c == alphaMult[l] & EQsp_c == EQsp[m])), 
                  j = i+6 := list(out[,3]), with = FALSE]
          rm(out, surv)
        }
      }
    }
  }
} # End of simulation code


# Munge and analyze output of "storage" data.table
# Do this on chinook, too!

spawnerDT <- copy(storage)

wrap_QE_counterDT <- function(spawnerDT, qeLev =  100) {
  # This function tallies the time (year) a population either went extinct or did not 
  # for each type of noise in the spawner abundance data.table for each replicate simulation
  # (i.e., by meanPS, sigPSmult, alphaMult, and EQsp)
  
  # create a data.table in which to store output
  
  meanPS <- unique(spawnerDT$meanPS_c)
  sigPSmult <- unique(spawnerDT$sigPSmult_c)
  alphaMult <- unique(spawnerDT$alphaMult_c)
  EQsp <- unique(spawnerDT$EQsp)
  reps <- unique(spawnerDT$reps_c)
  
  meanPS_r <- rep(meanPS, each = length(sigPSmult)*length(alphaMult)*length(EQsp)*length(reps))
  sigPSmult_r <- rep(sigPSmult, each = length(alphaMult)*length(EQsp)*length(reps), times = length(meanPS))
  alphaMult_r <- rep(alphaMult, each = length(EQsp)*length(reps), times = length(meanPS)*length(sigPSmult))
  EQsp_r <- rep(EQsp, each = length(reps), times = length(meanPS)*length(sigPSmult)*length(alphaMult))
  reps_r <- rep(reps, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp))
  #n_r <- rep(1:N, times = length(meanPS)*length(sigPSmult)*length(alphaMult)*length(EQsp)*reps)
  storage_ext <- data.table(meanPS_c = meanPS_r,
                        sigPSmult_c = sigPSmult_r,
                        alphaMult_c = alphaMult_r,
                        EQsp_c = EQsp_r,
                        reps_c = reps_r,
                        white = 0,
                        p34 = 0,
                        pNo34 = 0,
                        pgt3 = 0,
                        plt4 = 0,
                        pgt4 = 0,
                        plt3 = 0,
                        pgt10 = 0,
                        p34gt10 = 0) 

  setkey(storage_ext, meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c, reps_c)




}


  

QE_counter_DT <- function(spawnerDT, qeLev = 100) {
  #works on long(ish) rather than wide data

  # year of extinction will be in the 'named' noise columns
  
  # output will be used to calculate a pQE, but it will also be interesting to see how 
  # the histograms of the noise types depict the timing of the QE crossings
  # really dark noise (C&Y - black), showed that there where initial extinctions but
  # few afterwards. That is, if the slowly changing, random trajectory was going down 
  # to start, populations are/were in trouble, but otherwise tended to persist for a long time.
  
  # Step 1: construct a for loop that works through each subset* to tally the first year
  # when the qeLev has been crossed for 4 consecutive years for a certain set of noises.
  # 
  # *subeset by meanPS, sigPSmult, alphaMult, EQsp, and rep
  # Extinctions tallied by extinction year
  # Non-extinctions are tallied as 101 lends to pQE psuedo-code: length(x < 100)/reps
  
  di_len <- length(unique(spawnerDT$meanPS_c))
  j_len <- length(unique(spawnerDT$sigPSmult_c))
  k_len <- length(unique(spawnerDT$alphaMult_c))
  l_len <- length(unique(spawnerDT$EQsp_c))
  m_len <- length(unique(spawnerDT$reps_c))
 
  for (h in 1:length(freqCont))  {
    for (i in 1:i_len) {
      for (j in 1:j_len) {
        for (k in 1:k_len) {
          for (l in 1:l_len) {
            for (m in 1:m_len) {
              #storage_ext[] <-QE_counter_DT
              #print(h , i) #, j, k, l, m)
              spVec <- spawnerDT[i = (meanPS_c == meanPS[i] & sigPSmult_c == sigPSmult[j] & alphaMult_c == alphaMult[k] & EQsp_c == EQsp[l] & reps_c == m ),
                               j = h+6, with = FALSE ]
              qeInd <- which(spVec <= qeLev )  # note since there are "m_len" number of reps"licates" 
                                                                 # that are from 1:reps just set reps_c = the index "m".
               
               extRuns <- vector("integer", 0)
               kk <- 1
               if (length(qeInd) >= 4) {
                 
                 for (n in 1:(length(qeInd)-3)) {
                   
                   if (qeInd[n] == (qeInd[n+1] - 1) & qeInd[n] == (qeInd[n+2] - 2) & 
                         
                         qeInd[n] == (qeInd[n+3] - 3)  ) {
                     extRuns[kk] <- (qeInd[n] + 3)
                     kk <- kk + 1
                   
                   }  
                 
                 }
                 
               } else {
                 
                 extRuns[kk] <- (N + 1) 
                 kk <- kk+1
                 
               }
              if(length(extRuns > 0)) extInd <- min(extRuns) else extInd <- N + 1
                  
            storage_ext[i = (meanPS_c == meanPS[i] & sigPSmult_c == sigPSmult[j] & alphaMult_c == alphaMult[k] & EQsp_c == EQsp[l] & reps_c == m ),
                      j = h+5 := extInd, with = FALSE ] 
            }
          }
        }
      }
    }
  }
  
  
  ggplot(storage_ext, aes(x=white, fill = factor(alphaMult_c))) + 
    geom_density(alpha = 0.2) + 
    facet_grid(meanPS_c ~ sigPSmult_c)+ 
    theme_minimal() +
    xlim(0,100)
  
#   ggplot(storage_ext, aes(x=white, fill = factor(alphaMult_c))) + 
#     geom_histogram(alpha = 0.9, aes(y = 2 * ..density..)) + 
#     facet_grid(meanPS_c ~ sigPSmult_c)+ 
#     theme_minimal()
  
  
# determine pQE from storage_ext
  
  # for each parameter combination, count the number of extinction times < 100 and divide by 
  # the number of reps
  reps <- length(unique(storage_ext$reps_c))

storage_pQE <- copy(storage_ext[, list(count = .N, 
                        white_qe = length(which(white < 100))/reps,
                        p34_qe= length(which(p34 < 100))/reps,
                        pNo34_qe= length(which(pNo34 < 100))/reps,
                        pgt3_qe= length(which(pgt3 < 100))/reps,
                        plt4_qe= length(which(plt4 < 100))/reps,
                        pgt4_qe= length(which(pgt4 < 100))/reps,
                        plt3_qe= length(which(plt3 < 100))/reps,
                        pgt10_qe= length(which(pgt10 < 100))/reps,
                        p34gt10_qe= length(which(p34gt10 < 100))/reps),
                        by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)])

  
  
for (i in 1:9) {
  #i = 9
  qed <- copy(storage_pQE[,c(1:5,i+5), with = FALSE])
  names(qed)[6] <- "pQE"
  
  x <- ggplot(data = qed, aes(x = sigPSmult_c, y = pQE)) + 
    geom_line(colour = "slateblue", size = 1) +
    geom_point() +
    facet_grid(meanPS_c ~ alphaMult_c) +
    theme_minimal() +
    ggtitle(paste(toupper(names(storage_pQE)[i+5]))) +
    xlab(expression(paste(sigma))) +
    ylab("Probability of quasi-extinction (T = 100 years)")
  rm(qed)     
  print(x)
}

}


pQE_DT <- function() {
  # uses results of QE_counter_DT to calculate a pQE
}
# Results of stochastic simulations

#source functions.R (if need be)
source("~/NEP_salmon/chap_3/code/functions.R")
#load data
load(...)
#set parameters
qeLev <- 100




# Plot results

# pQE vs. variability in survPS



# mean and quantiles of time to QE for each freqContent




##-------------------OLD AND NOT USED-------------------------##
# surv1mid <- 0.025 # higher than cwt: natural and some fw mortality in cwt incorporated into SR??
# surv1low <- 0.01 #
# surv1high <- 0.05 #
# surv1vhigh <- 0.10
# surv2 <- 0.8
# surv3 <- 0.8
# survPS <- 1
# 
# 
# # The delta_e calcs are based on longterm/equilibrium assumptions.
# # CADFG data are based on recovery BC SRCS CWT recoveries. Might not be definitive, but provide 
# # a range to test
# de_50 <- calc_de(surv3, wanted_frac = 0.5) # Base case: from CADFG reports this is about the norm
# de_20 <- calc_de(surv3, wanted_frac = 0.2) # but sometimes many more come back at age 4
# de_80 <- calc_de(surv3, wanted_frac = 0.8) # and other times many return at age 3
# 
# # Calculate SPR, use survPS = 1 (no oversummer prespawning mortality) to get maximum SPR for assumptions
# # about surv1 and delta_e (assume surv2 and surv3 fixed)
# # 50% age-3 spawners
# sprS1mid <- SPR_srcs(surv1 = surv1mid, surv2, surv3, delta_e = de_50, survPS = 1)
# sprS1low <- SPR_srcs(surv1 = surv1low, surv2, surv3, delta_e = de_50, survPS = 1)
# sprS1high <- SPR_srcs(surv1 = surv1high, surv2, surv3, delta_e = de_50, survPS = 1)
# sprS1vhigh <- SPR_srcs(surv1 = surv1vhigh, surv2, surv3, delta_e = de_50, survPS = 1)
# 
# # 20% age-3 spawners
# sprS1mid20de <- SPR_srcs(surv1 = surv1mid, surv2, surv3, delta_e = de_20, survPS = 1)
# sprS1low20de <- SPR_srcs(surv1 = surv1low, surv2, surv3, delta_e = de_20, survPS = 1)
# sprS1high20de <- SPR_srcs(surv1 = surv1high, surv2, surv3, delta_e = de_20, survPS = 1)
# sprS1vhigh20de <- SPR_srcs(surv1 = surv1vhigh, surv2, surv3, delta_e = de_20, survPS = 1)
# 
# # 80% age-3 spawners
# sprS1mid80de <- SPR_srcs(surv1 = surv1mid, surv2, surv3, delta_e = de_80, survPS = 1)
# sprS1low80de <- SPR_srcs(surv1 = surv1low, surv2, surv3, delta_e = de_80, survPS = 1)
# sprS1high80de <- SPR_srcs(surv1 = surv1high, surv2, surv3, delta_e = de_80, survPS = 1)
# sprS1vhigh80de <- SPR_srcs(surv1 = surv1vhigh, surv2, surv3, delta_e = de_80, survPS = 1)
# 
# # Plot the relationship between survival (early ocean) and delta_e with SPR [for above parameters]
# 
# sprVec <- c(sprS1low20de, sprS1mid20de, sprS1high20de, sprS1vhigh20de,
#             sprS1low, sprS1mid, sprS1high, sprS1vhigh,
#             sprS1low80de, sprS1mid80de, sprS1high80de, sprS1vhigh80de)
# 
# sprDF <- data.frame(SPR = sprVec, delta_e = rep(c(0.2, 0.5, 0.8), each = 4), surv1 = rep(c(0.01, 0.025, 0.05, 0.10), times = 3))
# 
# plotSPR_de_surv1(sprDF)
# 
# vec <- c(1.5, 3, 5, 10)  # Assumptions about the productivity of the population at
# # low abundances relatiove to the slope of the replacement
# # line 2 = indicates alpha = 2*(1/SPR)
# 
# 
# 
# 
# # calc betas for SPR, alphas, EQsp
# # EQsp = 2500
# EQsp <- 2500
# # de = .2
# # S1 = low
# # alpha = spr
# alpha <- (1/spr)*vec # Hypotheses about alpha
# betaSPRS1lowde20EQsp2500alpha1_5 <- calc_beta(alpha = alpha[1], EQsp = EQsp, spr = sprVec[1])
# betaSPRde_EQspalpha3 <- calc_beta(alpha = alpha[1], EQsp = EQsp, spr = sprVec[1])
# betaSPRde_EQspalpha5 <- calc_beta(alpha = alpha[1], EQsp = EQsp, spr = sprVec[1])
# betaSPRde_EQspalpha110 <- calc_beta(alpha = alpha[1], EQsp = EQsp, spr = sprVec[1])
# 
# # EQsp = 7500
# EQsp <- 7500
# 
# 
# # EQsp = 15000
# EQsp <- 15000
# 
# # EQsp = 30000
# EQsp <- 30000
# 
# # EQsp inf (beta == 0)
# # The no density dependence case
# EQsp <- "infinity"
# beta0 <- 0
# 
# 
# 
# 
# # Load survPS data 
# 
# 
# # Call age-structured model [function] that holds early ocean survival constant [low, mid, high, very high]
# # Many scenarios
# 
# 
# 
# 
# 
# 
# # Call age-structured model [function] that allows time-varying ocean survival 
# # Many more scenarios
# 
# 
# 
# 
# # for (i in 1:length(mgmtScenarios)) {
# #   for (j in 1:length(climateScenarios)) {
# #     baseDetLM(mgmtScenarios = mgmtScenarios[i],
# #               climateScenarios = climateScenarios[j],
# #               n0 = c(10000, 5000, 5000, 5000),
# #               f3 = 5530/2,
# #               f4 = 5530/2,
# #               s_eo = 0.02,
# #               s2 = 0.5,
# #               s3 = 0.8,
# #               d3 = 0.6) 
# #   }
# # }
# UNUSED STOCHASTIC SIMULATION CODE ####
# ggplot(storage[i = alphaMult_c == 5], aes(x=N, y=white, colour = factor(reps_c))) + geom_line() + facet_grid(meanPS_c ~ sigPSmult_c, scale = "free_y") 
# 
# 
# save(spawnOutArray, "~/NEP_salmon/chap_3/data/popModelOut/spawnOutArray.Rdata")
# # Calculate statistics
# 
# whiteSp <- spawnOutArray["white",,,,,,,]
# rsw34Sp <- spawnOutArray["rsw34",,,,,,,]
# rsw34noSp <- spawnOutArray["rsw34no",,,,,,,]
# rswgt3Sp <- spawnOutArray["rswgt3",,,,,,,]
# rswlt4 <- spawnOutArray["rswlt4",,,,,,,]
# rswgt4 <- spawnOutArray["rswgt4",,,,,,,]
# rswlt3 <- spawnOutArray["rswlt3",,,,,,,]
# rswgt10 <- spawnOutArray["rswgt10",,,,,,,]
# rsw34gt10 <- spawnOutArray["rsw34gt10",,,,,,,]
# 
# rm(spawnOutArray)
# 
# # Calculate statistics
# wrap_QE_counter <- function(spawnerData, qeLev) {
# 
#   i_len <- dim(spawnerData)[1]
#   j_len <- dim(spawnerData)[2]
#   k_len <- dim(spawnerData)[3]
# 
#   QE_out <- array(NA, c(i_len, j_len, k_len, reps) )
#   
# for (i in 1:i_len) {
#   for (j in 1:j_len) {
#     for (k in 1:k_len) {
#         QE_out[i,j,k,] <- QE_counter(data = spawnerData[i,j,k,"QEthr_100",,], QE = qeLev)
#     }
#   } 
# }
#   return(QE_out)
# } # end of wrap_QE_counter function
# QElev <- 100
# QE_white_100 <- wrap_QE_counter(spawnerData = whiteSp, qeLev = QElev)
# system.time(QE_34_100 <- wrap_QE_counter(spawnerData = rsw34Sp, qeLev = QElev))
# QE_no34_100 <- wrap_QE_counter(spawnerData = rsw34noSp, qeLev = QElev)
# QE_lt4_100 <- wrap_QE_counter(spawnerData = rswlt4, qeLev = QElev)
# QE_gt3_100 <- wrap_QE_counter(spawnerData = rswgt4, qeLev = QElev)
# QE_lt3_100 <- wrap_QE_counter(spawnerData = rswlt3, qeLev = QElev)
# QE_gt10_100 <- wrap_QE_counter(spawnerData = rswgt10 , qeLev = QElev)
# QE_34gt10_100 <- wrap_QE_counter(spawnerData = rsw34gt10, qeLev = QElev)

#   for (i in 1:i_len) {
#     for (j in 1:j_len) {
#       for (k in 1:k_len) {
#         QE_out[i,j,k,] <- QE_counter(data = spawnerData[i,j,k,"QEthr_100",,], QE = qeLev)
#       }
#     } 
#   }

#   storage_ext[(i = (meanPS_c == meanPS[j] & sigPSmult_c == sigPSmult[k] & alphaMult_c == alphaMult[l] & EQsp_c == EQsp[m])), 
#          j = i+5 := list(out[,3]), with = FALSE] <- QE_counter_DT(data = spawnerDT[i = (i = (meanPS_c == meanPS[j] & sigPSmult_c == sigPSmult[k] & alphaMult_c == alphaMult[l] & EQsp_c == EQsp[m])), 
#                                                                                 j = i+5 := list(out[,3]), with = FALSE])
#  
#  return(QE_outDT)
  # bob <- copy(storage_ext)
  #  bob[, names(bob)[6:14] := lapply(.SD, function(x) (length(x < qeLev)/reps)), .SDcols=6:14]#,
  #      , by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c)]
  #  ddply(bob, .(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c), function(x) length(x[6:14] < 100))
  #DT[,lapply(.SD,sum),by=grp,.SDcols=c("x","y")]
  # 
  # bob[,names(bob)[6:14] := lapply(.SD, function(x) sum(x<100)) ,  # /reps
  #       by = list(meanPS_c, sigPSmult_c, alphaMult_c, EQsp_c), .SDcols=6:14]
  # bob
  