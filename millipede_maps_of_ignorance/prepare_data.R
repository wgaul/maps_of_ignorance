#############################
## prepare data for millipede maps of ignorance project
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 25 Oct 2019
##############################

### load millipede data
# NBDC training and test data
load("./data/6_Oct_download/NBDC_training_data.Rdata")
load("./data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

hec_names <- read_csv("./data/Irish_land_hectads.csv")

# combine test and training data for all groups
test_data <- test_data$millipede
training_data <- training_data$millipede
mill <- bind_rows(test_data, training_data)
rm(test_data, training_data)

### prepare and load environmental data ---------------------------------------
if(!all(file.exists("annual_precip_hectad.RData") & 
        file.exists("summer_tx_hectad.RData") & 
        file.exists("winter_tn_hectad.RData") & 
        file.exists("elevation_hec_ETOPO1.rds") & 
        file.exists("corine_label_1_hectad.RData") & 
        file.exists("corine_label_2_hectad.RData") & 
        file.exists("coast_dist_hec.RData") & 
        file.exists("coast_dist_1km.RData"))) {
  source("./eobs_annual_precipitation_hectad.R")
  source("./eobs_max_summer_temp_hec.R")
  source("./eobs_min_winter_temp_hec.R")
  source("./interp_elev_hec_etopo1.R")
  source("./prep_corine.R")
  source("./distance_to_coast_calc.R")
}
load("annual_precip_hectad.RData")
load("summer_tx_hectad.RData")
load("winter_tn_hectad.RData")
elev <- read_rds("elevation_hec_ETOPO1.rds")
load("corine_label_1_hectad.RData")
load("corine_label_2_hectad.RData")
load("coast_dist_hec.RData")

pred_brick <- brick(list(
  mean_tn = resample(krg_mean_tn_rast, krg_mean_rr_rast), 
  mean_tx = resample(krg_mean_tx_rast, krg_mean_rr_rast), 
  mean_rr = krg_mean_rr_rast, 
  artificial_surfaces = resample(artificial_surfaces_l1_rast, krg_mean_rr_rast), 
  forest_seminatural_l1 = resample(forest_seminatural_l1_rast, krg_mean_rr_rast),
  wetlands_l1 = resample(wetlands_l1_rast, krg_mean_rr_rast), 
  pasture_l2 = resample(pasture_l2_rast, krg_mean_rr_rast), 
  arable_l2 = resample(arable_land_l2_rast, krg_mean_rr_rast), 
  coast_dist = resample(coast_dist_hec_rast, krg_mean_rr_rast), 
  elev = resample(elev, krg_mean_rr_rast)))
# mask pred brick by one of the CORINE layers to get only Irish land cells
pred_brick <- mask(pred_brick, pred_brick$artificial_surfaces)
# scale and centre environmental predictors over study extent
pred_brick <- scale(pred_brick)

# make hec_names spatial 
hec_names_spat <- SpatialPointsDataFrame(coords = hec_names[, c("eastings", "northings")], 
                                         data = hec_names, 
                                         proj4string = CRS("+init=epsg:29903"))
# make sure millipede data is in same projection as predictor data
hec_names_spat <- spTransform(hec_names_spat, raster::projection(pred_brick))
### end prepare predictor variables -------------------------------------------

### prepare millipede data ----------------------------------------------------
## make Julian day and year variables
mill$year <- year(mill$StartDate)
mill$day_of_year <- yday(mill$StartDate)
mill$temp_resolution <- mill$EndDate - mill$StartDate

# join predictor variables to millipede data
# make millipede data spatial 
mill_spat <- SpatialPointsDataFrame(coords = mill[, c("eastings", "northings")], 
                                    data = mill, 
                                    proj4string = CRS("+init=epsg:29903"))
# make sure millipede data is in same projection as predictor data
mill_spat <- spTransform(mill_spat, raster::projection(pred_brick))

pred_df <- data.frame(raster::extract(pred_brick, mill_spat, 
                                      df = TRUE))
mill <- cbind(mill, pred_df)
### end prepare millipede data ------------------------------------------------

## scale years, day of year, and spatial coordinate dimensions
# eastings and northings should be scaled in the same way so space isn't more
# important in one direction
east_mean <- mean(mill$eastings)
north_mean <- mean(mill$northings)
mill$eastings_csc <- mill$eastings - east_mean
mill$northings_csc <- mill$northings - north_mean
spat_sd <- sd(c(mill$eastings_csc, mill$northings_csc))
mill$eastings_csc <- mill$eastings_csc/spat_sd
mill$northings_csc <- mill$northings_csc/spat_sd

mill$year_csc <- mill$year - mean(mill$year)
mill$year_csc <- mill$year_csc/sd(mill$year_csc)
mill$doy_csc <- mill$day_of_year - mean(mill$day_of_year)
mill$doy_csc <- mill$doy_csc/sd(mill$doy_csc)

# I think 10 km in space is about as big as a year and as big as about 8 julian
# days when all variables are scaled:
max(mill$year) - min(mill$year)
(max(mill$northings) - min(mill$northings))/10000
365/43
