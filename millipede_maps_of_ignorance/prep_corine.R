####################################
## Prepare Ireland CORINE landcover data for use at hectad scale
## 
## author: Willson Gaul wgaul@hotmail.com (based on a script by Hannah White
##          vasc_fres_all.R shared with wg via GitHub in April 2018)
## created: 2 May 2018 
## last modified: 25 Oct 2019 - change directory to millipede maps of ignorance project
####################################

corine <- raster("./data/CORINE_IE.grd")
c_legend <- read_csv("./data/clc_legend.csv")
hecs <- read_csv("./data/Irish_land_hectads.csv")

corine_sp <- rasterToPoints(corine, fun = function(x){x!=0}, spatial = T)
corine_sp <- spTransform(corine_sp, CRS("+init=epsg:29903")) # transform to IE grid
corine_df <- as.data.frame(corine_sp)

# add land cover text descriptions
corine_df <- left_join(corine_df, c_legend, 
                       by = c("g100_clc12_V18_5" = "GRID_CODE"))
# convert land cover text to valid names (to be used for column headers later)
corine_df$LABEL1 <- make.names(corine_df$LABEL1, unique = FALSE)
corine_df$LABEL2 <- make.names(corine_df$LABEL2, unique = FALSE)
#corine_df$LABEL3 <- make.names(corine_df$LABEL3, unique = FALSE)

# make array to hold the number of Label1-level land use classes in each hectad
clc_l1_counts_hecs <- array(0, dim = c(nrow(hecs), 
                                       length(unique(c_legend$LABEL1))))
colnames(clc_l1_counts_hecs) <- make.names(unique(c_legend$LABEL1), unique = T)
# make array to hold proportion of each Label1-level land use class in hectads
clc_l1_props_hecs <- array(0, dim = c(nrow(hecs), 
                                      length(unique(c_legend$LABEL1)))) 
colnames(clc_l1_props_hecs) <- make.names(unique(c_legend$LABEL1), unique = T)

# make array to hold the number of LABEL2-level land use classes in each hectad
clc_l2_counts_hecs <- array(0, dim = c(nrow(hecs), 
                                       length(unique(c_legend$LABEL2))))
colnames(clc_l2_counts_hecs) <- make.names(unique(c_legend$LABEL2), unique = T)
# make array to hold proportion of each LABEL2-level land use class in hectads
clc_l2_props_hecs <- array(0, dim = c(nrow(hecs), 
                                      length(unique(c_legend$LABEL2)))) 
colnames(clc_l2_props_hecs) <- make.names(unique(c_legend$LABEL2), unique = T)

for (i in 1:nrow(hecs)) {
  # find centers of corine raster cells that lie within each hectad
  # (some cells may cross hectad boundary but this is ignored)
  corine_sub = corine_df[which(corine_df$x >= hecs$eastings[i] & 
                                 corine_df$x < hecs$eastings[i] + 10^4 & 
                                 corine_df$y >= hecs$northings[i] & 
                                 corine_df$y < hecs$northings[i] + 10^4), ]
  
  ## LABEL1 classes
  # count number of corine LABEL1 classes within this hectad
  clc_l1_counts_hecs[i, ] <- sapply(colnames(clc_l1_counts_hecs), 
                                 FUN = function(x) {
                                   sum(x == corine_sub$LABEL1)})
  # calculate proportions of corine LABEL1 classes
  clc_l1_props_hecs[i, ] <- clc_l1_counts_hecs[i, ] / sum(
    clc_l1_counts_hecs[i, ])
  
  ## LABEL2 classes
  # count number of corine LABEL2 classes within this hectad
  clc_l2_counts_hecs[i, ] <- sapply(colnames(clc_l2_counts_hecs), 
                                    FUN = function(x) {sum(
                                      x == corine_sub$LABEL2)})
  # calculate proportions of corine LABEL1 classes
  clc_l2_props_hecs[i, ] <- clc_l2_counts_hecs[i, ] / sum(
    clc_l2_counts_hecs[i, ])
}

clc_l1_counts_hecs <- cbind(hecs, clc_l1_counts_hecs)
clc_l1_props_hecs <- cbind(hecs, clc_l1_props_hecs)
clc_l2_counts_hecs <- cbind(hecs, clc_l2_counts_hecs)
clc_l2_props_hecs <- cbind(hecs, clc_l2_props_hecs)

## promote results to spatial objects -----------------------------------------
coordinates(clc_l1_props_hecs) <- as.matrix(
  clc_l1_props_hecs[, c("eastings", "northings")])
coordinates(clc_l2_props_hecs) <- as.matrix(
  clc_l2_props_hecs[, c("eastings", "northings")])

proj4string(clc_l1_props_hecs) <- CRS("+init=epsg:29903")
proj4string(clc_l2_props_hecs) <- CRS("+init=epsg:29903")

gridded(clc_l1_props_hecs) <- TRUE # promote to SpatialPixelsDataFrame
gridded(clc_l2_props_hecs) <- TRUE 

clc_l1_props_hecs <- as(clc_l1_props_hecs, "SpatialGridDataFrame")
clc_l2_props_hecs <- as(clc_l2_props_hecs, "SpatialGridDataFrame")

## Make rasters of each variable  --------------------------------------------
## make a raster for each land cover class giving the proportion of the hectad
# (raster cell) that is that class.
# These 5 (or 16) rasters will be the predictor variables.

# LABEL 1 level
artificial_surfaces_l1_rast <- raster::raster(
  clc_l1_props_hecs["Artificial.surfaces"])
agricultural_l1_rast <- raster::raster(
  clc_l1_props_hecs["Agricultural.areas"])
forest_seminatural_l1_rast <- raster::raster(
  clc_l1_props_hecs["Forest.and.semi.natural.areas"])
wetlands_l1_rast <- raster::raster(
  clc_l1_props_hecs["Wetlands"])
water_l1_rast <- raster::raster(
  clc_l1_props_hecs["Water.bodies"])
# unclassified_l1_rast <- raster::raster(
#   clc_l1_props_hecs["UNCLASSIFIED"])

# LABEL2 level
urban_fabric_l2_rast <- raster::raster(
  clc_l2_props_hecs["Urban.fabric"])
industrial_commercial_transport_l2_rast <- raster::raster(
  clc_l2_props_hecs["Industrial..commercial.and.transport.units"])
mine_dump_construction_l2_rast <- raster::raster(
  clc_l2_props_hecs["Mine..dump.and.construction.sites"])
artificial_non_ag_vegetated_l2_rast <- raster::raster(
  clc_l2_props_hecs["Artificial..non.agricultural.vegetated.areas"])
arable_land_l2_rast <- raster::raster(
  clc_l2_props_hecs["Arable.land"])
permanent_crops_l2_rast <- raster::raster(
  clc_l2_props_hecs["Permanent.crops"])
pasture_l2_rast <- raster::raster(
  clc_l2_props_hecs["Pastures"])
heterogeneous_ag_l2_rast <- raster::raster(
  clc_l2_props_hecs["Heterogeneous.agricultural.areas"])
forest_l2_rast <- raster::raster(
  clc_l2_props_hecs["Forests"])
scrub_herbaceous_l2_rast <- raster::raster(
  clc_l2_props_hecs["Scrub.and.or.herbaceous.vegetation.associations"])
open_space_no_veg_l2_rast <- raster::raster(
  clc_l2_props_hecs["Open.spaces.with.little.or.no.vegetation"])
inland_wetlands_l2_rast <- raster::raster(
  clc_l2_props_hecs["Inland.wetlands"])
maritime_wetlands_l2_rast <- raster::raster(
  clc_l2_props_hecs["Maritime.wetlands"])
inland_water_l2_rast <- raster::raster(
  clc_l2_props_hecs["Inland.waters"])
marine_water_l2_rast <- raster::raster(
  clc_l2_props_hecs["Marine.waters"])

## end make rasters -----------------------------------------------------------

## clean workspace except for rasters and results dfs
rm(list = c("c_legend", "corine_df", "corine_sub", "hecs", "corine", 
            "corine_sp", "i"))

## save outputs --------------------------------------------------------------
## write out RData with rasters and results dfs
save(clc_l1_props_hecs, artificial_surfaces_l1_rast, agricultural_l1_rast, 
     forest_seminatural_l1_rast, wetlands_l1_rast, water_l1_rast,
     file = "corine_label_1_hectad.RData")

save(urban_fabric_l2_rast, industrial_commercial_transport_l2_rast, 
     mine_dump_construction_l2_rast, artificial_non_ag_vegetated_l2_rast, 
     arable_land_l2_rast, permanent_crops_l2_rast, pasture_l2_rast, 
     heterogeneous_ag_l2_rast, forest_l2_rast, scrub_herbaceous_l2_rast, 
     open_space_no_veg_l2_rast, inland_wetlands_l2_rast, 
     maritime_wetlands_l2_rast, inland_water_l2_rast, marine_water_l2_rast, 
     file = "corine_label_2_hectad.RData")