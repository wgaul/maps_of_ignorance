#############################
## prepare data for millipede maps of ignorance project
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 25 Oct 2019
## last modified: 6 Jan 2020
##############################

### load millipede data
# NBDC training and test data
load("./data/6_Oct_download/NBDC_training_data.Rdata")
load("./data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

hec_names <- read_csv("./data/Irish_land_hectads.csv")
# hopefully make hectads line up with grid cells of predictor raster brick
# hec_names$eastings = hec_names$eastings - 4900 
# hec_names$northings = hec_names$northings - 4900

# combine test and training data for all groups
test_data <- test_data$millipede
training_data <- training_data$millipede
mill <- bind_rows(test_data, training_data)
rm(test_data, training_data)

### prepare and load environmental data ---------------------------------------
if(!all(file.exists("annual_precip_hectad.rds") & 
        file.exists("annual_precip_1km.rds") &
        file.exists("summer_tx_hectad.rds") & 
        file.exists("summer_tx_1km.rds") & 
        file.exists("winter_tn_hectad.rds") &
        file.exists("winter_tn_1km.rds") & 
        file.exists("elevation_hec_ETOPO1.rds") & 
        file.exists("elevation_1km_ETOPO1.rds") & 
        file.exists("corine_label_1_hectad.RData") & 
        file.exists("corine_label_1_1km.RData") & 
        file.exists("corine_label_2_hectad.RData") & 
        file.exists("corine_label_2_1km.RData") &
        file.exists("coast_dist_hec.RData") & 
        file.exists("coast_dist_1km.RData") & 
        file.exists("soil_IFS_10km_brick.grd") & 
        file.exists("soil_IFS_1km_brick.grd") & 
        file.exists("soil_drainage_10km_brick.grd") & 
        file.exists("soil_drainage_1km_brick.grd"))) {
  source("./eobs_annual_precipitation_hectad.R")
  source("./eobs_max_summer_temp_hec.R")
  source("./eobs_min_winter_temp_hec.R")
  source("./interp_elev_hec_etopo1.R")
  source("./prep_corine.R")
  source("./distance_to_coast_calc.R")
  source("./soil_SIS_prep.R")
}
load("annual_precip_hectad.RData")
rr_rast_1k <- read_rds("annual_precip_1km.rds")
load("summer_tx_hectad.RData")
mean_tx_rast_1k <- read_rds("summer_tx_1km.rds")
load("winter_tn_hectad.RData")
mean_tn_rast_1k <- read_rds("winter_tn_1km.rds")
elev <- read_rds("elevation_hec_ETOPO1.rds")
elev_1k <- read_rds("elevation_1km_ETOPO1.rds")
load("corine_label_1_hectad.RData")
load("corine_label_1_1km.RData")
load("corine_label_2_hectad.RData")
load("corine_label_2_1km.RData")
load("coast_dist_hec.RData")
load("coast_dist_1km.RData")

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
# pred_brick <- mask(pred_brick, pred_brick$artificial_surfaces)
# scale and centre environmental predictors over study extent
pred_brick <- scale(pred_brick)

pred_brick_1k <- brick(list(
  mean_tn = mean_tn_rast_1k, 
  mean_tx = mean_tx_rast_1k, 
  mean_rr = rr_rast_1k, 
  artificial_surfaces = artificial_surfaces_l1_rast_1km, 
  forest_seminatural_l1 = forest_seminatural_l1_rast_1km, 
  wetlands_l1 = wetlands_l1_rast_1km, 
  pasture_l2 = pasture_l2_rast_1km,
  arable_l2 = arable_land_l2_rast_1km,
  coast_dist = coast_dist_1km_rast, 
  elev = elev_1k
))
pred_brick_1k <- scale(pred_brick_1k)

# make hec_names spatial 
hec_names_spat <- SpatialPointsDataFrame(coords = hec_names[, c("eastings", "northings")], 
                                         data = hec_names, 
                                         proj4string = CRS("+init=epsg:29903"))
# make sure millipede data is in same projection as predictor data
hec_names_spat <- spTransform(hec_names_spat, raster::projection(pred_brick))
### end prepare predictor variables -------------------------------------------

### prepare millipede data ----------------------------------------------------
## make checlist ID variable
mill$checklist_ID <- paste0(mill$RecorderHash, mill$StartDate, mill$EndDate, 
                            mill$SiteName, mill$Latitude, mill$Longitude)
## calculate list lengths for each checklist
mill$list_length <- NA
for(i in 1:length(unique(mill$checklist_ID))) {
  id <- unique(mill$checklist_ID)[i]
  ll <- length(which(mill$checklist_ID == id))
  mill$list_length[mill$checklist_ID == id] <- ll
}

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

pred_df <- data.frame(raster::extract(pred_brick, mill_spat, method = "simple",
                                      df = TRUE, cellnumbers = TRUE))
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
mill$doy_csc <- mill$day_of_year - 182.5 # center day of year
mill$doy_csc <- mill$doy_csc/sd(mill$doy_csc)

# I think 10 km in space is about as big as a year and as big as about 8 julian
# days when all variables are scaled:
max(mill$year) - min(mill$year)
(max(mill$northings) - min(mill$northings))/10000
365/43

