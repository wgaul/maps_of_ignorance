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
  summarise_at(vars("Adenomeris gibbosa":"Thalassisobates littoralis"), sum)

# get dataset of just Julus scandinavius detection/non-detection and predictor
# variables
J_scand <- select(mill_wide, -("Adenomeris gibbosa":"Glomeris marginata"), 
                  -c("Leptoiulus belgicus":"Thalassisobates littoralis")) %>%
  left_join(mill_fewer_vars[, colnames(mill_fewer_vars) %nin% 
                              c("Genus_species", "RecordId")], 
            by = "checklist_ID") %>%
  distinct()

# make J. scandinavius data spatial
J_scand_spat <- SpatialPointsDataFrame(
  coords = J_scand[, c("eastings", "northings")], 
  data = J_scand, 
  proj4string = CRS("+init=epsg:29903"))
# make sure millipede data is in same projection as predictor data
J_scand_spat <- spTransform(J_scand_spat, raster::projection(pred_brick))


### Make spatial block CV folds ------------------------------------------------
# Find best block size using only a subset of all data (for memory efficiency)
spatialAutoRange(pred_brick, sampleNumber = 1000, 
                 nCores = 3, showPlots = TRUE, plotVariograms = TRUE)

# make spatial blocks 
block_mill_10k <- spatialBlock(J_scand_spat, 
                               theRange = 90000, 
                               k = 5, 
                               selection = "random", 
                               iteration = 5, 
                               showBlocks = TRUE, 
                               rasterLayer = pred_brick$pasture_l2, 
                               biomod2Format = FALSE)
fold_index_10k <- block_mill_10k$foldID
### end make spatial blocks ---------------------------------------------------

# make new data with standardized recording effort
newdata <- cbind(hec_names,  
                 data.frame(raster::extract(pred_brick, hec_names_spat, 
                                            df = TRUE, 
                                            method = "simple", 
                                            cellnumbers = TRUE)))
newdata$list_length <- 10

### save objects for use on sonic ---------------------------------------------
saveRDS(J_scand, "J_scand.rds")
saveRDS(fold_index_10k, "fold_index_10k.rds")
saveRDS(newdata, "newdata.rds")
