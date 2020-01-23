#############################
## Fit SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_data.R'
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 23 Jan 2020
## last modified: 23 Jan 2020
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

### end make spatial blocks ---------------------------------------------------


### fit boosted regression tree -----------------------------------------------
fit_brt <- function(test_fold, fold_index, sp_df, pred_names) {
  # ARGS: test_fold - integer giving the CV fold to be used as test data
  #       fold_index - integer vector saying which CV fold each row of sp_df 
  #                     belongs in
  #       sp_df - data frame with species observations
  #       pred_names - character vector giving the variables to use as predictors
  mod <- tryCatch(gbm.step(data = sp_df[fold_index != test_fold], 
                           gbm.x = which(colnames(sp_df) %in% pred_names), 
                           gbm.y = "Julus scandinavius", 
                           tree.complexity = 5, learning.rate = 0.005, 
                           n.trees = 500, step.size = 100,
                           family = "bernoulli", max.trees = 10000, plot.main = F), 
                  error = function(x) NA)
  
  # make predictions for measuring AUC (predict only to data subset used for 
  # overall model trainign)
  f_pred <- tryCatch(predict(f_m, newdata = sp_df[fold_index == test_fold, ], 
                             type = "response",
                             n.trees = f_m$gbm.call$best.trees), 
                     error = function(x) NA)
  f_auc <- tryCatch(
    roc(response = factor(sp_df[fold_index == test_fold, "plantago_presence"], 
                          levels = c("0", "1")),
        predictor = f_pred,
        auc = T)$auc[1], 
    error = function(x) NA)
  
  # make dataframe of predictions with standardized recording effort
  newdata <- newdata[fold_index == test_fold, ]
  # make predictions to all grid cells in test_df
  f_pred <- tryCatch(predict(f_m, newdata = newdata, 
                             type = "response",
                             n.trees = f_m$gbm.call$best.trees), 
                     error = function(x) NA)
  
  newdata$pred <- tryCatch(f_pred, error = function(x) NA)
  # return fitted model, AUC value, and predictions for this model
  tryCatch(list(m = f_m, auc = f_auc, predictions = newdata), 
           error = function(x) "No list exported from fit_brt.")
}

# make new data with standardized recording effort
newdata <- cbind(hec_names,  
                 data.frame(raster::extract(pred_brick, hec_names_spat, 
                                            df = TRUE, 
                                            method = "simple", 
                                            cellnumbers = TRUE)))
newdata$list_length <- 10

J_scand_brt_fits <- mclapply(1:5, FUN = fit_brt, 
                             fold_index = block_mill_10k$foldID, 
                             sp_df = J_scand, 
                             pred_names = c("eastings", "northings", 
                                            "list_length"), 
                             mc.cores = n_cores)

### end fit brt ---------------------------------------------------------------