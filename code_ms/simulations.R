# simulations.R

# Overview
## Fig. 2: Frequency response of Spring Run Chinook salmon to white 
## noise at different levels of mean prespawning survival 
## (0.275 - just above collapse, 0.5, and 0.8).

## basic session info & logging

source("./code_ms/functions.R")

# libraries
library(ggplot2)
library(reshape2)
# library(plyr)
library(data.table)
library(TSA)
library(biwavelet)
# library(grid)
library(compiler)
# load foreach package for parallel processing
library(foreach)
# set up backend to do parallel processing
library(doParallel)
# detectCores()
registerDoParallel() # defaults to half of the cores 

system("hostname")  # record name of the machine
date() # record the date
sessionInfo() # documents the version of R and any included packages, 
# to reproduce the environment that yielded the results

# set random seed
r_seed <- 64
set.seed(r_seed)

# prelims

# todays_date <-  Sys.Date()
diss_dir_name <- file.path(".", "output_ms", "sim_results", "dissertation_files")
if(!exists(diss_dir_name)) dir.create(diss_dir_name, recursive = TRUE)

support_dir_name <- file.path(".", "output_ms", "sim_results", "support_files")
if (!exists(support_dir_name)) dir.create(support_dir_name, recursive = TRUE)

##  params

# "fixed" parameters for subsequent noise generation and population modeling
reps <- 10
N <- 1424
freqCont <- c("white", "p34", "pgt10", "p34gt10", "one_over_f")
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


# create white noise

# create "reps" number of white noise time series of length N using random sine wave approach
# with selected frequency contents
white_n <- customFR2ts(N = N, # number of time steps
                       reps = reps,
                       r_seed = r_seed,
                       amp = mk_white(N)) # mean = 0, variance/sd = 1

# Filtered/selected bandwidths

# Bandpass period 3-4 (frequencies 1/4 to 1/3) in white noise signals.

rsin_34_n <- customFR2ts(N = N, # number of time steps
                         reps = reps,
                         r_seed = r_seed,
                         amp = mk_rsin(N, highF=1/3, lowF=1/4)) 

rsin_gt10_n <- customFR2ts(N = N,
                           reps = reps,
                           r_seed = r_seed,
                           amp = mk_rsin(N, lowF=0, highF=1/10) ) 

red_beta_1 <- customFR2ts(N = N, # number of time steps
                          reps = reps,
                          r_seed = r_seed,
                          amp = mk_1_over_f_beta(N, beta = .5)) 

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
                  noise_gt10 = rsin_gt10_n,
                  noise_34gt10 = rsin_34_gt10_n,
                  noise_red_beta_1 = red_beta_1)

# rm individual sets of noise
# rm(white_n, rsin_34_n, rsin_gt10_n, rsin_34_gt10_n, red_beta_1)



# Frequency response of "population" at three levels of survival for white noise
alphaMult <- 4
EQsp <- 7500

# BASE
meanPS <- as.character(c(0.3, 0.4, 0.5))
survRange <- as.character(seq(0.2, 0.8, by = 0.2))


meanPS_r <- rep(meanPS, each = length(survRange)*length(alphaMult)*length(EQsp)*reps*N)
survRange_r <- rep(survRange, each = length(alphaMult)*length(EQsp)*reps*N, times = length(meanPS))
alphaMult_r <- rep(alphaMult, each = length(EQsp)*reps*N, times = length(meanPS)*length(survRange))
EQsp_r <- rep(EQsp, each = reps*N, times = length(meanPS)*length(survRange)*length(alphaMult))
reps_r <- rep(1:reps, each = N, times = length(meanPS)*length(survRange)*length(alphaMult)*length(EQsp))
n_r <- rep(1:N, times = length(meanPS)*length(survRange)*length(alphaMult)*length(EQsp)*reps)


storageP <- data.table::data.table(meanPS_c = meanPS_r,
                                   survRange_c = survRange_r,
                                   alphaMult_c = alphaMult_r,
                                   EQsp_c = EQsp_r,
                                   reps_c = reps_r,
                                   N = n_r,
                                   white = 0) 

data.table::setkey(storageP, meanPS_c, survRange_c, alphaMult_c, EQsp_c, reps_c)

# split storageP by meanPS so each component of the
# then I'll put all the lists back together at the end of the simulation
french <- vector("list", length = length(meanPS)) # random name, another version I used italian named in honor of LWB

for (i in 1:length(meanPS)) {
  french[[i]] <- storageP[i = meanPS_c == meanPS[i]]
}

# pre-compile functions to improve run time.
popSimPSvaryCmp <- cmpfun(popSimPSvary)

parSimCmp <- function(dt, 
                      noise_list,
                      freq_cont = c("white", "p34", "pgt10", "p34gt10", "one_over_f"),  
                      sim_len = 1024, 
                      surv_mean, 
                      surv_range, 
                      alpha_mult = 4,
                      EQ_sp = EQsp,
                      frac_wanted = 0.5) {
  for (i in 1:length(freq_cont)) {
    # for (j in surv_mean) {
      for (k in surv_range) {
        surv <- make_surv_mat(noise_dat = noise_list[[i]], 
                              mean_surv = as.numeric(surv_mean), 
                              sim_len,
                              range_surv = as.numeric(k))
        for (l in alpha_mult) {
          for (m in EQ_sp) {
            dt[i = (survRange_c == k & alphaMult_c == l & EQsp_c == m), 
               j = freq_cont[i] := list(melt(popSimPSvaryCmp(rand_surv = surv, 
                                                             surv1 = surv1, 
                                                             surv2 = surv2, 
                                                             surv3 = surv3, 
                                                             EQsp = m, 
                                                             wanted_frac = frac_wanted,
                                                             alpha_scale = l)[[2]])[,3])]
          }
        }
      } 
    # }
  }
  return(dt)
} # end of parSimCmp()


# run simulations using foreach framework to send jobs to multiple cores
# for each chunk of storageP - this should cut simulation times down to ~ 68/9 = 7.5 hours

system.time(
  storage <- foreach(h = 1:length(meanPS),.packages="reshape2") %dopar% {
    parSimCmp(dt = french[[h]], 
              noise_list = noiseList,
              freq_cont = c("white", "p34", "pgt10", "p34gt10", "one_over_f"),  
              sim_len = 1024, 
              surv_mean = meanPS[h], 
              surv_range = survRange)
  }
) 

storage <- rbindlist(storage)



# Fig 2: Plot frequency response at 3 survival levels 
pdf(file.path(".", "output_ms", "Fig_2_white_noise_popFreqResp_with_TimeSeries_CV.pdf"), width = 8, height = 6)
# old <- par(mfrow = c(4,1), mar = c(1,5, 1, 1))
old <- par(mar = c(1,5,1,1))
plot_dat <- storage[ i = N > 400 & survRange_c == 0.4]
plotMeanFR_DTmany(plot_dat, N = 1024, surv = meanPS[1], scale = "CV", yaxis_lim = c(0,1))
title(xlab = "Frequency")
linesMeanFR_DTmany(plot_dat, N = 1024, surv = meanPS[2], line_color = "black", scale = "CV")
linesMeanFR_DTmany(plot_dat, N = 1024, surv = meanPS[3], line_color = "grey30", scale = "CV")
legend("topright", legend = c(meanPS[1], meanPS[2], meanPS[3]), lty = c(2,1,1), col = c("black", "black", "grey30"), lwd = 3)
par(old)
dev.off()


## Fig 3 Summary plots of noise signals
pdf(file.path(".", "output_ms", "/Fig_3_summaryFreqContTS_Noise.pdf"), width = 8, height = 6)
plot_gen_freq_wvlt(noise = noiseList, 
                      n = 1, 
                      J1 = trunc((log(32/(2 * 1))/log(2))/0.01))
dev.off()

## Fig 4. Summary plot of spawning female abundance

pdf(file.path(".", "output_ms", "/Fig_4_summaryFreqContTS_SpawningFemales.pdf"), width = 8, height = 6)
plot_surv_spawn_ts(spawners = storage,
                   noise = noiseList,
                   meanSurv = 0.3,
                   rangeSurv = 0.4,
                   n = 1,
                   J1 = trunc((log(32/(2 * 1))/log(2))/0.01))
dev.off()


