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
## last modified: 24 Jan 2020
#################################
dbg <- F
calc_1k_distances <- F # run distances for 1km grid (might take a long time)
seed <- 23012020  # 23 Jan 2020
set.seed(seed) 

fit_brt <- F
fit_rf <- T
make_sampling_plots <- F
n_folds <- 3

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

source("functions_maps_of_ignorance.R")

n_cores <- 2

source("prepare_data.R")
source("prepare_objects_for_SDM.R")

if(fit_rf) source("fit_rf.R")
if(fit_brt) source("fit_brt.R")

if(make_sampling_plots) source("sampling_coverage_maps.R")



save.image("millipede_maps_sonic.RData")
