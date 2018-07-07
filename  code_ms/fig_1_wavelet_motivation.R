# Figure 1. Motivation: Summer survival rates of Butte Creek spring-run Chinook 
# salmon are projected to decline as summer stream temperatures increase and flows 
# decrease from climate change. (a) Example time series of summer adult survival 
# rates from Thompson et al. (2012). (b) Periods of environmental variability equal 
# to the mean generation time (approx. 3-4 years) of these salmon are evident in 
# the wavelet power spectrum from the middle 2060s to the late 2080s. The white 
# dashed line denotes the cone of influence that indicates the region where which 
# edge effects distort the spectrum. Thick black contour lines denote regions 
# where variance is significantly greater (α = 0.05) than a red-noise process with 
# the same lag-1 autocorrelation as the original time series (Torrence and Compo 1998).

## basic session info & logging

system("hostname")  # record name of the machine
date() # record the date
sessionInfo() # documents the version of R and any included packages, 
# to reproduce the environment that yielded the results

# set random seed
r_seed <- 64
set.seed(r_seed)

## libraries
library(biwavelet)
library(RSEIS)
## data

load_data_SALMOD_series <- function(management = c("BAU", "NoDiversion", "ColdWater", "ForeCast", "ForecastError", 
                                                   "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade"),
                                    emissions = c("A2","B1"),
                                    global_mods = c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")) {
  
  return()
}

## parameters



## wrangle

get_simulation_data_wavelet_plots <- function(management, # c("BAU", "NoDiversion")
                                              emissions, # c("A2","B1")
                                              global_mods) # c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1") 
{
  # get Thompson et al. prespawn survival SALMOD simulation data to generate
  # the motivating wavelet plot
  # params to generate this plot
  # management = "NoDiversion"
  # emissions = "B1"
  # global_mods = "mpiecham5"
  sim_path <- file.path(".", "data_ms", "simulation_results")
  
  mgmtScenarios <- management 
  climateScenarios <- emissions
  
  salmod_out <- read.delim(file.path(sim_path, mgmtScenarios, 
                                     paste0(climateScenarios, "_", global_mods), 
                                     "SALMODsumOutMerge.txt"), sep = ",")
  
  mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
  mortFW$global_mods <- global_mods
  years <- 2010:2099
  prespSurv <-  with(mortFW, data.frame(year = years, gcms = global_mods, PreSpawn = allMortSF/AFem))
  eggSurv <- with(mortFW, data.frame(year = years, gcms = global_mods, Egg = FryGrad/Eggs))
  frySurv <- with(mortFW, data.frame(year = years, gcms = global_mods, Fry = FryExit/FryGrad))
  eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
  frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
  
  return(list(prespSurv, eggSurv, frySurv))
}

wavelet_plot <- function(df, 
                         management,
                         emissions, 
                         global_mods) {
  mgmtScenarios <- management 
  climateScenarios <- emissions
  gcms <- global_mods 
  gcmScenarios <- numeric(length(climateScenarios)*length(global_mods))
  
  # number of scales minus - 1
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) 
  
  # Calculate and visualize wavelet spectra
  
  out_wavelet <- file.path(".", "output_ms")
  if (!dir.exists(out_wavelet)) { dir.create(out_wavelet, recursive = TRUE) }
  
  path4file <- file.path(out_wavelet, paste0("Fig_1_wavelet_", mgmtScenarios, "_", climateScenarios, "_", global_mods, ".pdf"))
  pdf(file = path4file)
  
  # Time series
  par(fig = c(0, 1, 0.5, 1))
  old <- par(mar = c(3, 3, 3, 1) )
  plot(df$year, df$PreSpawn, type = "l", lwd = 2, col = "darkgrey",
       axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
  axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
       lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
       las = 1, cex.axis = 1, tck = 0.02)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       las = 2, cex.axis = 1, tck = 0.02)
  box(lwd = 2)
  # title(paste0("GCM: ", gcms, " Prespawn survival emissions scenario: ", 
  #              climateScenarios, " \n Water managament alternative: ", 
  #              mgmtScenarios))
  par(old)
  
  # Wavelet Power Spectrum
  par(fig= c(0, 1, 0, 0.5), new = TRUE)
  old <- par(mar = c(3, 3, 1, 1) )
  test_cwt <- wt(cbind(df$year, detrend(as.numeric(scale(df$PreSpawn)))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
  axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
       lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
       las = 1, cex.axis = 1, tck = 0.02)
  par(old)
  
  dev.off()
}

dat <- get_simulation_data_wavelet_plots(management = "NoDiversion",
                                         emissions = "B1", 
                                         global_mods = "mpiecham5")
ps_dat <- dat[[1]]

wavelet_plot(df = ps_dat, 
             management = "NoDiversion",
             emissions = "B1", 
             global_mods = "mpiecham5") 

