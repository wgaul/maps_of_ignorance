####################################
## Maps of Ignorance - spatial coverage
## Distance to nearest record in space
## 
## author: Willson Gaul
## created: 22 Oct 2019 
## last modified: 24 Oct 2019
#####################################

setwd("~/Documents/Data_Analysis/UCD/maps_of_ignorance/")
library(wgutil)
library(Hmisc)
library(raster)
library(spatstat)
library(GGally)
library(parallel)
library(lubridate)
library(tidyverse)
source("~/Documents/Data_Analysis/R_templates/biorecording-master/os2eastnorth.R")

### load data --------------------------------------------------------------
# bryophyte training and test data
load("~/Documents/Data_Analysis/UCD/bryophytes/data/bry_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/bryophytes/data/test_vault/bry_test_data.Rdata")

# NBDC training and test data
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/NBDC_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

hec_names <- read_csv("~/Documents/Data_Analysis/UCD/mapping/data/Irish_land_hectads.csv")
### end load data --------------------------------------------------------------

### format data ------------------------------------------------------------
bry_test_data$split <- rep("test", nrow(bry_test_data))
bry_training_data$split <- rep("training", nrow(bry_training_data))

# combine test and training data for all groups
bry <- bind_rows(bry_training_data, bry_test_data)
nbdc <- test_data
for (i in 1:length(names(nbdc))) {
  taxon = names(nbdc)[i]
  nbdc[[i]] <- bind_rows(training_data[which(names(training_data) == taxon)], 
                         test_data[which(names(test_data) == taxon)])
}
# remove NBDC bryophyte groups, since I have the BBS bry data
nbdc <- nbdc[-c(which(names(nbdc) %in% c("liverwort", "hornwort", "moss")))]

# convert NBDC gridRefs to hectads (many currently 1km blocks)
for (i in 1:length(nbdc)) {
  nbdc[[i]]$hectad <- gridref_to_hec(nbdc[[i]]$GridReference)
}

# remove british hecs from bryophyte data
bry <- bry[which(nchar(bry$OSGR_10km) == 3), ]

# make a version of bird data without Bird Atlas 2007 - 2011
nbdc$bird_noAtlas <- nbdc$bird[nbdc$bird$DataSetHash != 
                                 "AB385E4D6E73BFA8C42CD208BD00B1E0", ]
# make a version of bird data with only Bird Atlas 2007 - 2011
nbdc$bird_Atlas <- nbdc$bird[nbdc$bird$DataSetHash == 
                               "AB385E4D6E73BFA8C42CD208BD00B1E0", ]

# make unique checklist identifier
nbdc <- lapply(nbdc, FUN = function(x) {x$checklist_ID <- paste(
  x$StartDate, x$GridReference, x$SiteName, x$RecorderHash, sep = "_"); x})
bry$checklist_ID <- paste(bry$Year, bry$Month, bry$Date, bry$OSGR,
                          bry$Locality, bry$observer_1,
                          bry$observer_2, bry$observer_3, bry$observer_4,
                          bry$observer_5, sep = "_")

rm(test_data, training_data, bry_test_data, bry_training_data) # clean workspace
### end format data ------------------------------------------------------------


euc_dist <- function(x, y) {
  # compute the euclidean distance between points x and y 
  # ARGS: x, y - vectors of matching lengths giving the coordinates of each
  #         point x and y in any number of dimensions
  # VALUE:  The euclidean distance between x and y
  if(length(x) != length(y)) stop("x and y must be vectors of the same length.")
  ndim <- length(x) # get the number of dimensions to use
  sides <- mapply(FUN = function(a, b) {a-b}, a = x, b = y)
  sqrt(sum(sides^2)) # return euclidean distance
}

dist_to_nearest_record <- function(df, query_points, coords, parallel = FALSE, 
                                   ncores = NULL, ...) {
  # calculate the euclidean distance to the nearest record in df from the 
  # points in query_points.
  # All dimensions should be scaled before being sent into this function.
  # ARGS: df - data frame of observations containing columns for coordinates
  #       query_points - data frame giving the points from which distances to 
  #           nearest records are desired.  This should have two columns for 
  #           coordinates with the same names as those in df and coords
  #       coords - character vector of any length giving the names of the columns
  #           containing the coordinates of each record
  # VALUE: a data frame with all the points in query_points and the distances 
  #         from those points to the nearest records in df
  list2env(list(...), envir = environment())
  # Check that dimensions have been scaled by making sure SDs are within 2
  # (all SDs should be close to 1 after scaling, so differences between SDs should
  # be small, though SDs might not be exactly 1 if some sites are removed after
  # scaling variables)
  df_sds <- lapply(df[, coords], function(x) sd(x, na.rm = T))
  df_sds <- as.matrix(dist(df_sds))
  if(any(df_sds > 2)) stop("Scale dimensions of records before using this function.")
  q_sds <- lapply(query_points[, coords], function(x) sd(x, na.rm = T))
  q_sds <- as.matrix(dist(q_sds))
  if(any(q_sds > 2)) stop("Scale dimensions of query points before using this function.")
  
  dists <- query_points # make df to hold resulting minimum distances
  dists$dist_to_nearest_rec <- NA
 
  # subset columns to only coords and coerce to data frame (in case it is a tibble)
  df <- data.frame(df[, coords]) 
  query_points <- data.frame(query_points[, coords])
  # coerce all columns to numeric
  for (i in 1:ncol(df)) {df[, i] <- as.numeric(df[, i])}
  for (i in 1:ncol(query_points)) {
    query_points[, i] <- as.numeric(query_points[, i])}

  if(parallel) {
    cl <- makeCluster(ncores)
    for (i in 1:nrow(query_points)) {
      # find smallest distance between query point i and a record in df
      ds <- parApply(cl = cl, X = df, MARGIN = 1, FUN = euc_dist, 
                     y = query_points[i, ], chunk.size = chunk.size)
      if(any(!is.na(ds))) {
        dists$dist_to_nearest_rec[i] <- min(ds, na.rm = TRUE)
      } else dists$dist_to_nearest_rec[i] <- NA
    }
    stopCluster(cl)
  }
  if(!parallel) {
    for (i in 1:nrow(query_points)) {
      # find smallest distance between query point i and a record in df
      dists$dist_to_nearest_rec[i] <- min(apply(df, MARGIN = 1, FUN = euc_dist, 
                                                y = query_points[i, ]))
    }
  }
  
  # make distances relative so most distant grid query point has distance = 1
  dists$dist_to_nearest_rec <- dists$dist_to_nearest_rec/max(
    dists$dist_to_nearest_rec, na.rm = T)
  dists # return df with distance from each query point to nearest record
}

# test 
od <- nbdc$`insect - dragonfly (Odonata)`
#rm(nbdc)
query_points <- hec_names
query_points$StartDate <- as.Date("2017-10-23")

# scale dimensions to 0 to 1 range
# eastings and northings should be scaled in the same way so space isn't more
# important in one direction
spat_range <- max(c(max(od$eastings) - min(od$eastings), 
                max(od$northings) - min(od$northings)))
od$eastings_sc <- od$eastings/spat_range
od$northings_sc <- od$northings/spat_range

date_range <- max(as.numeric(od$StartDate)) - min(as.numeric(od$StartDate))
od$StartDate_sc <- as.numeric(od$StartDate)/date_range

# a hectad (10 km) is about as much "distance" as a year in this dataset if 
# dimensions are scaled to 0 to 1:
spat_range/10000 # spatial range is about 44 hectads, so 0 to 1 will span 44 hectads
year(max(od$StartDate)) - year(min(od$StartDate)) # 0 to 1 will span 47 years

# center and scale query point variables
query_points$eastings_sc <- query_points$eastings/spat_range
query_points$northings_sc <- query_points$northings/spat_range

query_points$StartDate_sc <- as.numeric(query_points$StartDate)/date_range

## measure spatial/temporal distance
test_dist_spTime <- dist_to_nearest_record(od, 
                                           query_points = query_points, 
                                           coords = c("eastings_sc", 
                                                      "northings_sc", 
                                                      "StartDate_sc"), 
                                           parallel = T, ncores = 3, 
                                           chunk.size = 1000)
ggplot(data = test_dist_spTime, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  # geom_point(data = od[1:1000, ], aes(x = eastings, y = northings,
  #                                     size = StartDate),
  #            color = "orange") +
  # geom_point(data = od[tp, ], aes(x = eastings, y = northings),
  #            color = "red") + 
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial / temporal distance")

## measure spatial distance only
test_dist_space <- dist_to_nearest_record(od[1:300, ], 
                                           query_points = query_points, 
                                           coords = c("eastings_sc", 
                                                      "northings_sc"), 
                                          parallel = T, ncores = 3, 
                                          chunk.size = 100)
ggplot(data = test_dist_space, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  geom_point(data = od[1:300, ], aes(x = eastings, y = northings),
             color = "orange") +
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial distance")

## measure spatial & environmental distance
# join environmental variables to data
load("~/Documents/Data_Analysis/UCD/predictor_variables/ETOPO1/elevation_hec_ETOPO1.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/annual_precip_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/summer_tx_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/winter_tn_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/CORINE/corine_label_1_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/CORINE/corine_label_2_hectad.RData")
load("~/Documents/Data_Analysis/UCD/bird_SDM_IE/coast_dist_hec.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/CCM/ccm21/river_length_hec.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/WFDLakeWaterbodies_19022016/lake_area_hec.RData")
## combine all predictors into a single df ------------------------------------
# make a data frame of predictors
pred_rast_brick <- brick(list(
  mean_tn = resample(krg_mean_tn_rast, krg_mean_rr_rast), 
  mean_tx = resample(krg_mean_tx_rast, krg_mean_rr_rast), 
  mean_rr = krg_mean_rr_rast, 
  artificial_surfaces = resample(artificial_surfaces_l1_rast, krg_mean_rr_rast), 
  forest_seminatural_l1 = resample(forest_seminatural_l1_rast, krg_mean_rr_rast),
  wetlands_l1 = resample(wetlands_l1_rast, krg_mean_rr_rast), 
  pasture_l2 = resample(pasture_l2_rast, krg_mean_rr_rast), 
  arable_l2 = resample(arable_land_l2_rast, krg_mean_rr_rast), 
  lake_area = resample(lk_area_hec_rast, krg_mean_rr_rast), 
  river_length = resample(rv_length_hec_rast, krg_mean_rr_rast), 
  coast_dist = resample(coast_dist_hec_rast, krg_mean_rr_rast), 
  elev = resample(elev_hec, krg_mean_rr_rast)))
# scale predictor variables
pred_rast_brick_csc <- raster::scale(pred_rast_brick, center = T, scale = T)

# make odonata data spatial
od_spat <- SpatialPointsDataFrame(coords = od[, c("eastings", "northings")], 
                                    data = od, 
                                    proj4string = CRS("+init=epsg:29903"))
od_spat <- spTransform(od_spat, raster::projection(pred_rast_brick_csc))

# match predictor values to od observations
raw_pred_df <- cbind(od[, colnames(od) == "hectad"], 
                     data.frame(raster::extract(pred_rast_brick_csc, od_spat, 
                                                df = TRUE)))
# using sp for the previous transformation, I must remove an "ID" column that 
# gets added somehow
raw_pred_df <- select(raw_pred_df, -ID) %>%
  distinct()
colnames(raw_pred_df)[1] <- "hectad"
raw_pred_df <- full_join(raw_pred_df, hec_names, by = "hectad")
raw_pred_df <- raw_pred_df[raw_pred_df$hectad %in% hec_names$hectad, ]

# join predictors to query points
query_points <- left_join(query_points, raw_pred_df)
# join predictors to odonata data
od <- left_join(od, raw_pred_df[, which(colnames(raw_pred_df) %nin% 
                                          c("eastings", "northings"))], 
                by = "hectad")

dist_sp_env <- dist_to_nearest_record(od[1:300, ], 
                                      query_points = query_points, 
                                      coords = c("eastings_sc", "northings_sc", 
                                                 "mean_tn", "mean_tx", "mean_rr", 
                                                 "artificial_surfaces", 
                                                 "forest_seminatural_l1", 
                                                 "wetlands_l1", "pasture_l2", 
                                                 "arable_l2", "lake_area", 
                                                 "river_length", "coast_dist", 
                                                 "elev"), 
                                      parallel = T, ncores = 3, 
                                      chunk.size = 100)

ggplot(data = dist_sp_env, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  geom_point(data = od[1:300, ], aes(x = eastings, y = northings),
             color = "orange") +
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial & environmental distance")

# measure spatial, temporal, & environmental distance
dist_sp_tm_env <- dist_to_nearest_record(
  od[1:20000, ], 
  query_points = query_points, 
  coords = c("eastings_sc", "northings_sc", 
             "StartDate_sc", 
             "mean_tn", "mean_tx", "mean_rr", 
             "artificial_surfaces", 
             "forest_seminatural_l1", 
             "wetlands_l1", "pasture_l2", 
             "arable_l2", "lake_area", 
             "river_length", "coast_dist", 
             "elev"), 
  parallel = T, ncores = 3, 
  chunk.size = 1000)

ggplot(data = dist_sp_tm_env, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = log(dist_to_nearest_rec)), 
             size = 5) + 
  # geom_point(data = od[1:300, ], aes(x = eastings, y = northings),
  #            color = "orange") +
  # scale_color_gradient(limits = c(0, 1), name = "") + 
  ggtitle("spatial, temporal, & environmental distance")

# measure evenness of distances
simpson_even(dist_sp_tm_env$dist_to_nearest_rec)

