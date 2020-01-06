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
## last modified: 6 Jan 2020
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

t_size <- 19
line_size <- 1.5
c_f = 1 # coord fixed for plots

## What is spatial resolution of records?
table(mill$Precision)

## What is temporal resolution of records?
table(as.numeric(mill$temp_resolution))

temp_res_df <- data.frame(temp_resolution = c("one_day", "one_month", "one_year"), 
                          nrec = NA)
temp_res_df$nrec[temp_res_df$temp_resolution == "one_day"] <- length(which(
  as.numeric(mill$temp_resolution) < 2))
temp_res_df$nrec[temp_res_df$temp_resolution == "one_month"] <- length(which(
  as.numeric(mill$temp_resolution) >= 2 & as.numeric(mill$temp_resolution) < 31))
temp_res_df$nrec[temp_res_df$temp_resolution == "one_year"] <- length(which(
  as.numeric(mill$temp_resolution) > 31 & as.numeric(mill$temp_resolution) < 366))

ggplot(data = temp_res_df, aes(x = factor(temp_resolution, 
                                          levels = c("one_day", "one_month", 
                                                     "one_year"), 
                                          labels = c("One Day", "One Month", 
                                                     "One Year")), 
                               y = nrec)) + 
  geom_bar(stat = "identity") +
  xlab("Temporal Resolution of Records") + 
  ylab("Number of Records") + 
  theme_bw()


## map spatial bias - data density
# count number of checklists per hectad
nlist_hec <- select(mill, hectad, checklist_ID) %>%
  group_by(hectad) %>%
  summarise(nlist = n_distinct(checklist_ID)) %>%
  full_join(., hec_names, by = "hectad") %>%
  replace_na(list(nlist = 0))

# specify locations of annotations for compass arrow and scale bar
annot <- data.frame(x1 = c(265000, 310000, 60000, 60000), 
                    x2 = c(365000, 310000, 60000, 60000), 
                    y1 = c(60000, 40000, 400000, 380000), 
                    y2 = c(60000, 40000, 455000, 380000),
                    label = c(NA, "100 km", NA, "N"), 
                    bias = "e)")

ggplot() + 
  geom_raster(data = nlist_hec, 
              aes(x = eastings, y = northings, fill = nlist)) + 
  coord_fixed(c_f) + 
  scale_fill_gradient(name = "Number\nof\nChecklists\n", 
                      low = "light grey", high = "black") + 
  geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label)) + 
  geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
               arrow = arrow(length = unit(0.1, "npc"))) + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(hjust = -0.01))


## make query points ---------------------------------------------------------
query_points <- cbind(hec_names,  
                      data.frame(raster::extract(pred_brick, hec_names_spat, 
                                                 df = TRUE, 
                                                 method = "simple", 
                                                 cellnumbers = TRUE)))
query_points$year_csc <- mill$year_csc[mill$year == max(mill$year)][1]
query_points$eastings_csc <- query_points$eastings - east_mean # centre
query_points$eastings_csc <- query_points$eastings_csc/spat_sd # scale
query_points$northings_csc <- query_points$northings - north_mean # centre
query_points$northings_csc <- query_points$northings_csc/spat_sd # scale

## measure spatial distance ---------------------------------------------------
tp <- sample(1:nrow(mill), size = 600) # just for testing
dist_sp <- dist_to_nearest_record(mill[tp, ], 
                                  query_points = query_points, 
                                  coords = c("eastings_csc", 
                                             "northings_csc"), 
                                  parallel = T, ncores = 3, 
                                  chunk.size = 200) # 1621
ggplot(data = dist_sp, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  # geom_point(data = mill, aes(x = eastings, y = northings),
  #            size = 1, color = "orange") +
  geom_point(data = mill[tp, ], aes(x = eastings, y = northings),
             color = "red") +
  # scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial distance")

## measure environmental distance -----------------------------------
dist_sp_env <- dist_to_nearest_record(mill[tp, ], 
                                      query_points = query_points, 
                                      coords = c("mean_tn", "mean_tx", 
                                                 "mean_rr", 
                                                 "artificial_surfaces", 
                                                 "forest_seminatural_l1", 
                                                 "wetlands_l1", "pasture_l2", 
                                                 "arable_l2", "coast_dist", 
                                                 "elev"), 
                                      parallel = T, ncores = 3, 
                                      chunk.size = 200) #1621
ggplot(data = dist_sp_env, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 7) + 
  geom_point(data = mill[tp, ], 
             aes(x = eastings, y = northings),
             size = 1, color = "orange") +
  # scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("environmental distance")



## measure spatial & environmental distance -----------------------------------


## measure spatial & year distance -------------------------------------------


## measure spatial, environmental, & year distance -----------------------------


## measure spatial, environmental, year, & day of year distance ---------------
# make day of year be July 1
query_points$doy_csc <- mill$doy_csc[mill$day_of_year == yday("2017-07-01")][1]


## histogram of years ---------------------------------------------------------
list_year_summary <- select(mill, checklist_ID, StartDate, year) %>%
  distinct()

ggplot(data = list_year_summary, aes(x = year)) +
  geom_histogram(bins = max(mill$year) - min(mill$year)) + 
  xlab("Year") + 
  ylab("Number of Checklists") + 
  theme_bw()


## histogram of day of year ----------------------------------------------------
# ggplot(data = mill, aes(x = day_of_year)) + 
#   geom_histogram(breaks = seq(0, 365, by = 14)) +
#   coord_polar(start = 0) +
#   scale_x_continuous("", breaks = c(0, 31, 59, 91, 120, 151, 181, 212, 243, 
#                                     273, 304, 334), 
#                      labels = c("Jan", "Feb", "March", "April", "May", "June", 
#                                 "July", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
#   theme_bw()



list_week_summary <- select(mill, checklist_ID, StartDate, day_of_year) %>%
  distinct()

ggplot(data = list_week_summary, aes(x = day_of_year)) + 
  geom_histogram(breaks = seq(0, 365, by = 14)) +
  scale_x_continuous("", breaks = c(0, 31, 59, 91, 120, 151, 181, 212, 243, 
                                    273, 304, 334) + 15, 
                     labels = c("Jan", "Feb", "March", "April", "May", "June", 
                                "July", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ylab("Number of Checklists") + 
  theme_bw()

# mill$week <- lubridate::week(mill$StartDate)
# week_nrec_summary <- group_by(mill, week) %>%
#   summarise(nrecs_week = n())
# 
# ggplot(data = week_nrec_summary, aes(x = week, y = nrecs_week)) + 
#   geom_line() + 
#   coord_polar()


## measure spatial distance by month -----------------------------------------
test_dates <- c("2017-01-15", "2017-03-15", "2017-06-15", "2017-09-15") #
test_dates <- yday(test_dates)            

query_points$day_lag_csc <- (0 - 90)/sd(-90:90) #lag of zero "centred" and scaled in same way millipede day lag will be below
mill$day_lag_csc <- NA


make_plot <- function(test_dates, mill, query_points) {
  mill$day_lag_csc <- abs(mill$day_of_year - test_dates)
  mill$day_lag_csc[abs(mill$day_lag_csc) > 365/2] <- 365 - 
    mill$day_lag_csc[mill$day_lag_csc > 365/2]
  mill$day_lag_csc <- mill$day_lag_csc - 90 # roughly centre day lag
  mill$day_lag_csc <- mill$day_lag_csc/sd(-90:90) # scale by the same amount every time
  
  dist_sp_doy <- dist_to_nearest_record(mill, 
                                        query_points = query_points, 
                                        coords = c("eastings_csc", 
                                                   "northings_csc", 
                                                   "day_lag_csc"), 
                                        parallel = T, ncores = 3, 
                                        chunk.size = 334) #1621
  ggplot(data = dist_sp_doy, aes(x = eastings, y = northings)) + 
    geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
    scale_color_gradient(limits = c(0, 1), name = "Distance\nto nearest\nrecord") + 
    ggtitle(paste0("spatial/seasonal distance at doy ", test_dates))
}

warning("Spatial/temporal distances will take a long time to calculate.")
dist_by_month_plots <- lapply(test_dates, FUN = make_plot, 
                              mill = mill, query_points = query_points)

multiplot(plotlist = dist_by_month_plots, cols = 3)
