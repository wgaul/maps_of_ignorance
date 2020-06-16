################################
## This script organises the analysis for the millipede maps of ignorance
## 
## TODO:
##  - put all environmental data in this project directory
##  - add soil data
##  - add geology data
##  
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 16 June 2020
#################################
warning("TODO: 21 May 2020: Re-set the number of spatial undersamplingn block configurations to make to a higher number (e.g. 2000 or something much more than the number of models to test and fit).")

rm(list = ls())
dbg <- F
calc_1k_distances <- F # run distances for 1km grid (might take a long time)
seed <- 23012020  # 23 Jan 2020
set.seed(seed) 

run_rf <- T
make_sampling_plots <- F # map sampling coverage in env. and geographic space
make_spatial_blocks <- T # takes a few minutes. Set to T for final run
get_partial_dependence <- T # calculate partial dependence (time consuming)

analysis_resolution <- 1000 # analysis resolution (10000 or 1000 m rid squares)
n_folds <- 3 # number of cross-validation folds to use
n_cv_trials <- 3 # number of different cross-validation fold layouts to use
cv_block_sizes <- c("random", 30000) # sizes of CV spatial blocks (in meters)
n_subsamp_block_draws <- 10 # number of spatial subsampling block configurations to make
block_range_spat_undersamp <- 30000 # spatial undersampling grid block size (m)

library(wgutil)
library(Hmisc)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(fields)
library(gstat)
library(raster)
library(parallel)
library(sf) # can't install sf on sonic as of 8 Jan 2020
library(fasterize)
library(rgdal)
library(dismo)
#library(rgeos)
library(lubridate)
library(tidyverse)
library(randomForest)

source("functions_maps_of_ignorance.R")

n_cores <- 1

sp_to_fit <- list("Ommatoiulus sabulosus")

# sp_to_fit <- list("Julus scandinavius", "Tachypodoiulus niger",
#                   "Ommatoiulus sabulosus", "Cylindroiulus punctatus",
#                   "Boreoiulus tenuis", "Proteroiulus fuscus",
#                   "Blaniulus guttulatus", "Cylindroiulus latestriatus",
#                   "Proteroiulus fuscus", "Ophiodesmus albonanus",
#                   "Glomeris marginata", "Macrosternodesmus palicola")
#  
names(sp_to_fit) <- sp_to_fit

source("prepare_data.R")
source("prepare_objects_for_SDM.R")
mod_names <- c("day_ll_rf", "spat_ll_rf", "env_ll_rf", "env_spat_ll_rf") # 
mods_for_pd_plots <- c("spat_ll_rf", "env_spat_ll_rf")

if(run_rf) source("fit_rf.R")

## evaluate models
evals <- data.frame() # data frame to hold evaluation results
# Specify the model(s) using the string that was used as the beginning of the 
# file names that hold the fitted model objects.  e.g. here use "day_ll_rf" 
# for the file day_ll_rf_noSubSamp_fits_Julus_scaninavius.rds

for(mod_name in mod_names) {
  source("evaluate_models.R")
}

## make plots
source("plots_main.R")
source("plots_supplementary.R")

if(make_sampling_plots) source("sampling_coverage_maps.R")



save.image("millipede_maps_sonic.RData")
