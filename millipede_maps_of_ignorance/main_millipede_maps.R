################################
## This script organises the analysis for the millipede maps of ignorance
## 
## TODO:
##  - put all environmental data in this project directory
##  - add soil data
##  - add geology data
##  - get all hectads in CORINE
##  
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 25 Oct 2019
#################################
dbg <- F

library(wgutil)
library(Hmisc)
library(rgdal)
library(gdalUtils)
library(ggplot2)
library(fields)
library(gstat)
library(raster)
library(parallel)
library(sf)
library(fasterize)
library(rgdal)
library(rgeos)
library(lubridate)
library(tidyverse)

source("../functions_maps_of_ignorance.R")

source("prepare_data.R")

## make query points ---------------------------------------------------------
query_points <- cbind(hec_names,  
                      data.frame(raster::extract(pred_brick, hec_names_spat, 
                                                 df = TRUE)))
query_points$year_csc <- mill$year_csc[mill$year == max(mill$year)][1]
query_points$eastings_csc <- query_points$eastings - east_mean # centre
query_points$eastings_csc <- query_points$eastings_csc/spat_sd # scale
query_points$northings_csc <- query_points$northings - north_mean # centre
query_points$northings_csc <- query_points$northings_csc/spat_sd # scale

## measure spatial distance ---------------------------------------------------
dist_sp <- dist_to_nearest_record(mill, 
                                  query_points = query_points, 
                                  coords = c("eastings_csc", 
                                             "northings_csc"), 
                                  parallel = T, ncores = 3, 
                                  chunk.size = 1621)
ggplot(data = dist_sp, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  geom_point(data = mill, aes(x = eastings, y = northings),
             size = 1, color = "orange") +
  # geom_point(data = mill[tp, ], aes(x = eastings, y = northings),
  #            color = "red") + 
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial distance")

## measure spatial & environmental distance -----------------------------------
dist_sp_env <- dist_to_nearest_record(mill, 
                                      query_points = query_points, 
                                      coords = c("eastings_csc", 
                                                 "northings_csc", 
                                                 "mean_tn", "mean_tx", 
                                                 "mean_rr", 
                                                 "artificial_surfaces", 
                                                 "forest_seminatural_l1", 
                                                 "wetlands_l1", "pasture_l2", 
                                                 "arable_l2", "coast_dist", 
                                                 "elev"), 
                                      parallel = T, ncores = 3, 
                                      chunk.size = 1621)
ggplot(data = dist_sp_env, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  geom_point(data = mill, aes(x = eastings, y = northings),
             size = 1, color = "orange") +
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial & environmental distance")

## measure spatial, environmental, & year distance -----------------------------


## measure spatial, environmental, year, & day of year distance ---------------
# make day of year be July 1
query_points$doy_csc <- mill$doy_csc[mill$day_of_year == yday("2017-07-01")][1]


