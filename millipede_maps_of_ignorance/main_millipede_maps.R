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
## last modified: 8 Jan 2020
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
#library(sf) # can't install sf on sonic as of 8 Jan 2020
library(fasterize)
library(rgdal)
#library(rgeos)
library(lubridate)
library(tidyverse)

source("functions_maps_of_ignorance.R")

n_cores <- 22

source("prepare_data.R")

# tp <- sample(1:nrow(mill), size = 600) # test points - subset of records for testing
tp <- 1:nrow(mill)
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
dim(nlist_hec)
head(nlist_hec)

# specify locations of annotations for compass arrow and scale bar
annot <- data.frame(x1 = c(265000, 310000, 60000, 60000), 
                    x2 = c(365000, 310000, 60000, 60000), 
                    y1 = c(60000, 40000, 400000, 380000), 
                    y2 = c(60000, 40000, 455000, 380000),
                    label = c(NA, "100 km", NA, "N"), 
                    bias = "e)")

p_data_density <- tryCatch({ggplot() + 
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
          strip.text = element_text(hjust = -0.01))}, 
    error = function(x) NA)

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
dist_sp <- dist_to_nearest_record(mill[tp, ], 
                                  query_points = query_points, 
                                  coords = c("eastings_csc", 
                                             "northings_csc"), 
                                  parallel = T, ncores = n_cores, 
                                  chunk.size = floor(nrow(mill)/n_cores)) # 1621
p_spat_dist <- tryCatch({ggplot(data = dist_sp, aes(x = eastings, y = northings)) + 
    geom_raster(aes(fill = dist_to_nearest_rec)) + 
    # geom_point(data = mill[tp, ], aes(x = eastings, y = northings),
    #            color = "orange") +
    scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                        trans = "reverse") + 
    coord_fixed(c_f) +
    ggtitle("spatial distance") + 
    theme(text = element_text(size = t_size*1.2), 
          legend.key.width = unit(1.8*t_size, "points"))}, 
    error = function(x) NA)

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
                                      parallel = T, ncores = n_cores, 
                                      chunk.size = floor(nrow(mill[tp, ])/n_cores)) #1621
p_dist_env <- tryCatch({ggplot(data = dist_sp_env, aes(x = eastings, y = northings)) + 
    geom_raster(aes(fill = dist_to_nearest_rec)) + 
    # geom_point(data = mill[tp, ], 
    #            aes(x = eastings, y = northings),
    #            size = 1, color = "orange") +
    scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                        trans = "reverse") + 
    coord_fixed(c_f) + 
    ggtitle("environmental distance") + 
    theme_bw() + 
    theme(text = element_text(size = t_size*1.2), 
          legend.key.width = unit(1.8*t_size, "points"))}, 
    error = function(x) NA)



## measure spatial & environmental distance -----------------------------------


## measure spatial & year distance -------------------------------------------


## measure spatial, environmental, & year distance -----------------------------


## measure spatial, environmental, year, & day of year distance ---------------
# make day of year be July 1
query_points$doy_csc <- mill$doy_csc[mill$day_of_year == yday("2017-07-01")][1]


## histogram of years ---------------------------------------------------------
list_year_summary <- select(mill, checklist_ID, StartDate, year) %>%
  distinct()

p_year_hist <- tryCatch({ggplot(data = list_year_summary, aes(x = year)) +
    geom_histogram(bins = max(mill$year) - min(mill$year)) + 
    xlab("Year") + 
    ylab("Number of Checklists") + 
    theme_bw()}, 
    error = function(x) NA)


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

p_week_hist <- tryCatch({ggplot(data = list_week_summary, aes(x = day_of_year)) + 
    geom_histogram(breaks = seq(0, 365, by = 14)) +
    scale_x_continuous("", breaks = c(0, 31, 59, 91, 120, 151, 181, 212, 243, 
                                      273, 304, 334) + 15, 
                       labels = c("Jan", "Feb", "March", "April", "May", "June", 
                                  "July", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
    ylab("Number of Checklists") + 
    theme_bw()}, 
    error = function(x) NA)

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


calc_month_dists <- function(test_dates, mill, query_points) {
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
                                        parallel = T, ncores = n_cores, 
                                        chunk.size = floor(nrow(mill)/n_cores)) 
}

warning("Spatial/temporal distances will take a long time to calculate.")
dists_by_month <- tryCatch({lapply(test_dates, FUN = calc_month_dists, 
                                   mill = mill[tp, ], 
                                   query_points = query_points)}, 
                           error = function(x) NA)

# Distance at doy 15
p_season_winter <- ggplot(data = dists_by_month[[1]], 
                          aes(x = eastings, y = northings)) + 
  geom_raster(aes(fill = dist_to_nearest_rec)) + 
  scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                      trans = "reverse") + 
  coord_fixed(c_f) + 
  ggtitle("spatial/seasonal distance\nWinter") + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(t_size, "points"))

# Distance at doy 74
p_season_spring <- ggplot(data = dists_by_month[[2]], 
                          aes(x = eastings, y = northings)) + 
  geom_raster(aes(fill = dist_to_nearest_rec)) + 
  scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                      trans = "reverse") + 
  coord_fixed(c_f) + 
  ggtitle("spatial/seasonal distance\nSpring") + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(t_size, "points"))

# Distance at doy 166
p_season_summer <- ggplot(data = dists_by_month[[3]], 
                          aes(x = eastings, y = northings)) + 
  geom_raster(aes(fill = dist_to_nearest_rec)) + 
  scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                      trans = "reverse") + 
  coord_fixed(c_f) + 
  ggtitle("spatial/seasonal distance\nSummer") + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(t_size, "points"))

# Distance at doy 258
p_season_autumn <- ggplot(data = dists_by_month[[4]], 
                          aes(x = eastings, y = northings)) + 
  geom_raster(aes(fill = dist_to_nearest_rec)) + 
  scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                      trans = "reverse") + 
  coord_fixed(c_f) + 
  ggtitle("spatial/seasonal distance\nAutumn") + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(t_size, "points"))
### end distances by month ----------------------------------------------------



### Spatial & Environmental distance in 3 main recording periods --------------
# 1970 - 1985
dist_sp_hec_1970_85 <- dist_to_nearest_record(
  mill[mill$year >= 1970 & 
         mill$year <= 1985, ], 
  query_points = query_points, 
  coords = c("eastings_csc", 
             "northings_csc"), 
  parallel = T, ncores = n_cores, 
  chunk.size = floor(nrow(mill[mill$year >= 1970 & mill$year <= 1985, ])/n_cores))

p_sp_1970_85 <- tryCatch({ggplot(data = dist_sp_hec_1970_85, 
                       aes(x = eastings, y = northings)) + 
    geom_raster(aes(fill = dist_to_nearest_rec)) + 
    # geom_point(data = mill[mill$year >= 1970 & mill$year <= 1985, ], 
    #            aes(x = eastings, y = northings), color = "orange") + 
    scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                        trans = "reverse") + 
    coord_fixed(c_f) +
    ggtitle("spatial distance\n1970 - 1985") + 
    theme(text = element_text(size = t_size), 
          legend.key.width = unit(t_size, "points"))}, 
    error = function(x) NA)

# 1986 - 2005
dist_sp_hec_1986_2005 <- dist_to_nearest_record(
  mill[mill$year >= 1986 & 
         mill$year <= 2005, ], 
  query_points = query_points, 
  coords = c("eastings_csc", 
             "northings_csc"), 
  parallel = T, ncores = n_cores, 
  chunk.size = floor(nrow(mill[mill$year >= 1986 & 
                                 mill$year <= 2005, ])/n_cores))

p_sp_1986_05 <- tryCatch({ggplot(data = dist_sp_hec_1986_2005, 
                                 aes(x = eastings, y = northings)) + 
    geom_raster(aes(fill = dist_to_nearest_rec)) + 
    # geom_point(data = mill[mill$year >= 1986 & mill$year <= 2005, ], 
    #           aes(x = eastings, y = northings), color = "orange") + 
    scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                        trans = "reverse") + 
    coord_fixed(c_f) +
    ggtitle("spatial distance\n1986 - 2005") + 
    theme(text = element_text(size = t_size), 
          legend.key.width = unit(t_size, "points"))}, 
    error = function(x) NA)

# 2006 - 2019
dist_sp_hec_2006_2019 <- dist_to_nearest_record(
  mill[mill$year >= 2006 & 
         mill$year <= 2019, ], 
  query_points = query_points, 
  coords = c("eastings_csc", 
             "northings_csc"), 
  parallel = T, ncores = n_cores, 
  chunk.size = floor(nrow(mill[mill$year >= 2006 & 
                                 mill$year <= 2019, ])/n_cores))

p_sp_2006_19 <- tryCatch({ggplot(data = dist_sp_hec_2006_2019, 
                                 aes(x = eastings, y = northings)) + 
    geom_raster(aes(fill = dist_to_nearest_rec)) + 
    # geom_point(data = mill[mill$year >= 2006 & mill$year <= 2019, ],
    #            aes(x = eastings, y = northings), color = "orange") +
    scale_fill_gradient(name = "Euclidean\ndistance\nto nearest\nrecord", 
                        trans = "reverse") + 
    coord_fixed(c_f) +
    ggtitle("spatial distance\n2006 - 2019") + 
    theme(text = element_text(size = t_size), 
          legend.key.width = unit(t_size, "points"))}, 
    error = function(x) NA)

### end distance in 3 main recording periods ----------------------------------

pdf("millipede_maps.pdf") 
print(p_data_density)
print(p_spat_dist)
print(p_dist_env)
print(p_year_hist)
print(p_week_hist)
print(multiplot(plotlist = dist_by_month_plots, cols = 2))
print(multiplot(p_sp_1970_85, p_sp_1986_05, p_sp_2006_19, 
                layout = matrix(c(1,2,3,3), nrow = 2, byrow = TRUE)))
print(multiplot(
  p_sp_1970_85 + geom_point(
    data = mill[mill$year >= 1970 & mill$year <= 1985, ], 
    aes(x = eastings, y = northings), color = "orange"), 
  p_sp_1986_05 + geom_point(
    data = mill[mill$year >= 1986 & mill$year <= 2005, ], 
    aes(x = eastings, y = northings), color = "orange"), 
  p_sp_2006_19 + geom_point(
    data = mill[mill$year >= 2006 & mill$year <= 2019, ],
    aes(x = eastings, y = northings), color = "orange"), 
  layout = matrix(c(1,2,3,3), nrow = 2, byrow = TRUE)))
dev.off()

## save jpg versons of plots
try(ggsave("seasonal_spatial_distance.jpg", 
       multiplot(
         p_season_winter, p_season_spring, p_season_summer, p_season_autumn, 
         layout = matrix(c(1,2,3,4), nrow = 2, byrow = TRUE)), 
       width = 25, height = 25, units = "cm", 
       device = "jpeg"))

try(ggsave("seasonal_spatial_distance_2.jpg", 
           print(multiplot(
             p_season_winter, p_season_spring, p_season_summer, p_season_autumn, 
             layout = matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))), 
           width = 25, height = 25, units = "cm", 
           device = "jpeg"))

try(ggsave("survey_period_spatial_dists_no_points.jpg", 
       multiplot(p_sp_1970_85, p_sp_1986_05, p_sp_2006_19, 
                layout = matrix(c(1,2,3), nrow = 1, byrow = TRUE)), 
       width = 45, height = 15, units = "cm", 
       device = "jpeg"))

try(ggsave("survey_period_spatial_dists.jpg", multiplot(
  p_sp_1970_85 + geom_point(
    data = mill[mill$year >= 1970 & mill$year <= 1985, ], 
    aes(x = eastings, y = northings), color = "orange"), 
  p_sp_1986_05 + geom_point(
    data = mill[mill$year >= 1986 & mill$year <= 2005, ], 
    aes(x = eastings, y = northings), color = "orange"), 
  p_sp_2006_19 + geom_point(
    data = mill[mill$year >= 2006 & mill$year <= 2019, ],
    aes(x = eastings, y = northings), color = "orange"), 
  layout = matrix(c(1,2,3), nrow = 1, byrow = TRUE)), 
  width = 45, height = 15, units = "cm", 
  device = "jpeg"))

## alternately, save jpgs individually
try(ggsave("winter_spatial_distance.jpg", p_season_winter, 
           width = 25, height = 25, units = "cm", 
           device = "jpg"))
try(ggsave("spring_spatial_distance.jpg", p_season_spring, 
           width = 25, height = 25, units = "cm", 
           device = "jpg"))
try(ggsave("summer_spatial_distance.jpg", p_season_summer, 
           width = 25, height = 25, units = "cm", 
           device = "jpg"))
try(ggsave("fall_spatial_distance.jpg", p_season_autumn, 
           width = 25, height = 25, units = "cm", 
           device = "jpg"))

save.image("millipede_maps_sonic_test_8Jan2020.RData")
