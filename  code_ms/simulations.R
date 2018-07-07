# simulations.R

## basic session info & logging

source("./code_ms/functions.R")

system("hostname")  # record name of the machine
date() # record the date
sessionInfo() # documents the version of R and any included packages, 
# to reproduce the environment that yielded the results

# set random seed
r_seed <- 64
set.seed(r_seed)


# libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(biwavelet)
library(grid)

# prelims

todays_date <-  Sys.Date()
diss_dir_name <- file.path(".", "output_ms", "sim_results", paste0("dissertation_files_", todays_date))
if(!exists(diss_dir_name)) dir.create(diss_dir_name, recursive = TRUE)

support_dir_name <- file.path(".", "output_ms", "sim_results", paste0("support_files_", todays_date))
if (!exists(support_dir_name)) dir.create(support_dir_name, recursive = TRUE)

# params

reps <- 1000
N <- 500

