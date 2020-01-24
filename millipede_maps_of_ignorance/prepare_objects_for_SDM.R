#############################
## Prepare observation data for millipede SDMs
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_data.R'.  This puts the millipede data 
## into a format for use on sonit by fit_SDM.R. This script also makes the
## spatial block cross-validation folds using spatial packages that are not
## available on sonic.  This script saves the required objects as .rds files
## which can be loaded on sonic for fitting SDMs
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 24 Jan 2020
## last modified: 24 Jan 2020
##############################
library(blockCV)
library(sf)

# For now, just fit a single environmental SDM and a single coordinates + survey
# effort SDM to Julius scandinavius to see how it looks
# with these predictors (which do not yet include soil), vs. a spatial only model.

# spread mill data frame to have a column for each species
mill_fewer_vars <- select(mill, -DataSetHash, -SurveyHash, -RecorderHash, 
                          -CommonName, -TaxonGroupName, -StartMonth, -EndMonth, 
                          -Genus, -Species, -TaxonName, 
                          -(Name_matched_TNRS:Accepted_name_lsid_TNRS), 
                          -ID, -cells) 

mill_wide <- mill_fewer_vars %>% 
  mutate(observed = 1) %>%  
  spread(Genus_species, value = observed, fill = 0) %>%
  select(-RecordId) %>%
  group_by(checklist_ID) %>%
  summarise_at(vars("Adenomeris gibbosa":"Thalassisobates littoralis"), sum) %>%
  left_join(mill_fewer_vars[, colnames(mill_fewer_vars) %nin% 
                                            c("RecordId", "StartDate", "EndDate",
                                              "SiteName", "GridReference", 
                                              "Latitude", "Longitude", 
                                              "Precision")], 
            by = "checklist_ID")


### Make spatial block CV folds ------------------------------------------------
# Find best block size using only a subset of all data (for memory efficiency)
spatialAutoRange(pred_brick, sampleNumber = 1000, 
                 nCores = 3, showPlots = TRUE, plotVariograms = TRUE)

# make spatial blocks 
block_mill_10k <- spatialBlock(mill_spat, 
                               theRange = 100000, 
                               k = 5, 
                               selection = "random", 
                               iteration = 5, 
                               showBlocks = TRUE, 
                               rasterLayer = pred_brick$pasture_l2, 
                               biomod2Format = FALSE)

mill_wide <- SpatialPointsDataFrame(
  coords = mill_wide[, c("eastings", "northings")], 
  data = mill_wide, proj4string = CRS("+init=epsg:29903"))
mill_wide <- st_join(st_as_sf(mill_wide), st_as_sf(block_mill_10k$blocks))
mill_wide <- data.frame(mill_wide)
mill_wide <- select(mill_wide, -Easting, -Northing, -geometry, -layer)
### end make spatial blocks ---------------------------------------------------

# make new data with standardized recording effort
newdata <- cbind(hec_names,  
                 data.frame(raster::extract(pred_brick, hec_names_spat, 
                                            df = TRUE, 
                                            method = "simple", 
                                            cellnumbers = TRUE)))
newdata <- SpatialPointsDataFrame(coords = newdata[, c("eastings", "northings")], 
                                  data = newdata, 
                                  proj4string = CRS("+init=epsg:29903"))
newdata <- st_join(st_as_sf(newdata), st_as_sf(block_mill_10k$blocks))

newdata_temp <- data.frame(newdata)
newdata_temp <- select(newdata_temp, -Easting, -Northing, -geometry, -layer)

newdata <- data.frame()

for(i in 1:182) {
  dt <- newdata_temp
  dt$day_of_year <- i*2
  newdata <- bind_rows(newdata, dt)
}
newdata$list_length <- 4

### save objects for use on sonic ---------------------------------------------
saveRDS(mill_wide, "mill_wide.rds")
saveRDS(newdata, "newdata.rds")
