#############################
## Fit SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_objects_for_SDM.R'
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 23 Jan 2020
## last modified: 28 Jan 2020
##############################
library(dismo)
library(pROC)
library(parallel)
library(dplyr)
library(tibble)
set.seed(01242020) # Jan 24 2020
n_cores <- 1
on_sonic <- T

# load objects that are needed to fit SDMs
mill_wide <- readRDS("mill_wide.rds")
newdata <- readRDS("newdata.rds")

### fit boosted regression tree -----------------------------------------------
fit_brt <- function(test_fold, sp_name, sp_df, pred_names) {
  # ARGS: test_fold - integer giving the CV fold to be used as test data
  #       sp_name - character string giving the name of the column with 
  #                 detection/non-detection data to use
  #       sp_df - data frame with species observations and a column called 
  #               "folds" giving the cross-validation folds each record is in
  #       pred_names - character vector giving the variables to use as predictors

  if(is_tibble(sp_df)) {
    sp_df <- data.frame(sp_df)
  }
  sp_name <- gsub(" ", ".", sp_name)
  colnames(sp_df) <- gsub(" ", ".", colnames(sp_df))
  
  sp_df[sp_df[, sp_name] > 1, sp_name] <- 1 
  
  # fit boosted regression tree
  mod <- tryCatch(gbm.step(data = sp_df[sp_df$folds != test_fold, ], 
                           gbm.x = which(colnames(sp_df) %in% pred_names), 
                           gbm.y = sp_name, 
                           tree.complexity = 3, learning.rate = 0.005, 
                           n.trees = 500, step.size = 100,
                           family = "bernoulli", max.trees = 1000,
                           plot.main = F), 
                  error = function(x) NA)

  # make predictions for measuring AUC (predict only to data subset used for 
  # overall model training)
  f_pred <- tryCatch(predict(mod, newdata = sp_df[sp_df$fold == test_fold, ], 
                             type = "response",
                             n.trees = mod$gbm.call$best.trees), 
                     error = function(x) NA)
  f_auc <- tryCatch(
    roc(response = factor(sp_df[sp_df$folds == test_fold, sp_name], 
                          levels = c("0", "1")),
        predictor = f_pred, auc = T), 
    error = function(x) NA)
  
  # make dataframe of predictions with standardized recording effort
  newdata$pred[newdata$folds == test_fold] <- tryCatch({
    predict(mod, newdata = newdata[newdata$folds == test_fold, ], 
            type = "response", n.trees = mod$gbm.call$best.trees)}, 
    error = function(x) NA)
  
  # return fitted model, AUC value, and predictions for this model
  tryCatch(list(m = mod, auc = f_auc$auc[1], predictions = newdata, 
                roc_object = f_auc), 
           error = function(x) "No list exported from fit_brt.")
}

call_fit_brt <- function(sp_name, test_fold, sp_df, pred_names, ...) {
  # Function to call fit_rf on each species in a list of species dfs
  lapply(1:5, FUN = fit_brt, sp_name = sp_name, 
         sp_df = sp_df, pred_names = pred_names)
}



## fit spatial models
sp_to_fit <- list("Glomeris marginata")
# "Ophyiulus pilosus", "Blaniulus guttulatus", 
# "Tachypodoiulus niger", "Julus scandinavius", 

names(sp_to_fit) <- sp_to_fit
spatial_brt_fits <- mclapply(sp_to_fit, 
                            FUN = call_fit_brt, 
                            sp_df = mill_wide, 
                            pred_names = c("eastings", "northings", 
                                           "day_of_year", 
                                           "list_length"), 
                            mc.cores = n_cores)
try(print(pryr::object_size(spatial_brt_fits)))
try(saveRDS(spatial_brt_fits, "spatial_brt_fits.rds"))

## fit environmental model
env_brt_fits <- mclapply(sp_to_fit, 
                         FUN = call_fit_brt, 
                         sp_df = mill_wide, 
                         pred_names = c("mean_tn", "mean_tx", 
                                        "mean_rr", "artificial_surfaces", 
                                        "forest_seminatural_l1", 
                                        "wetlands_l1", "pasture_l2", 
                                        "arable_l2", "elev", 
                                        "day_of_year", 
                                        "list_length"), 
                         mc.cores = n_cores)
try(print(pryr::object_size(env_brt_fits)))
try(saveRDS(env_brt_fits, "env_brt_fits.rds"))
### end fit brt ---------------------------------------------------------------


# view AUC (averaged over 5 folds)
lapply(spatial_brt_fits, 
       FUN = function(x) {mean(sapply(x, FUN = function(y) {
         tryCatch(y$auc, error = function(x) NA)}), 
         na.rm = T)})
lapply(env_brt_fits, 
       FUN = function(x) {mean(sapply(x, FUN = function(y) {
         tryCatch(y$auc, error = function(x) NA)}), 
         na.rm = T)})

# Get predictions with standardized survey effort
mill_predictions_spatial_brt <- lapply(
  spatial_brt_fits, 
  FUN = function(x) {lapply(x, FUN = function(x) {
    tryCatch(x$predictions, error = function(x) NA)})
  })
mill_predictions_spatial_brt <- lapply(
  mill_predictions_spatial_brt, 
  FUN = function(x) {
    bind_rows(x[!is.na(x)])})

mill_predictions_env_brt <- lapply(
  env_brt_fits, 
  FUN = function(x) {lapply(x, FUN = function(x) {
    tryCatch(x$predictions, error = function(x) NA)})})
mill_predictions_env_brt <- lapply(
  mill_predictions_env_brt, 
  FUN = function(x) {
    bind_rows(x[!is.na(x)])})

# get variable importance (averaged over 5 folds) 
mill_var_imp_spatial_brt <- lapply(
  spatial_brt_fits, 
  FUN = function(x) {
    bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
      df <- data.frame(y$m$contributions)
      df$variable <- rownames(df)
      df[, c("variable", "rel.inf")]}))}
)
mill_var_imp_env_brt <- lapply(
  env_brt_fits, 
  FUN = function(x) {
    bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
      df <- data.frame(y$m$contributions)
      df$variable <- rownames(df)
      df[, c("variable", "rel.inf")]}))}
)


# save results
try(saveRDS(mill_predictions_spatial_brt, "mill_predictions_spatial_brt.rds"))
try(saveRDS(mill_predictions_env_brt, "mill_predictions_env_brt.rds"))




### exploratory plotting - this is not for running on sonic
if(!on_sonic) {
  # plot average of predictions from all 5 folds (so 4 predictions will be to
  # training data, one prediction to test data in each grid cell)
  
  # spatial model
  mill_predictions_spatial_brt <- lapply(
    mill_predictions_spatial_brt, FUN= function(x) {
      group_by(x, hectad) %>%
        summarise(mean_prediction = mean(pred, na.rm = T), 
                  eastings = mean(eastings), northings = mean(northings))
    })
  
  for(i in 1:length(mill_predictions_spatial_brt)) {
    print(ggplot() + 
            geom_tile(data = mill_predictions_spatial_brt[[i]], 
                      aes(x = eastings, y = northings, fill = mean_prediction)) + 
            geom_point(data = mill, aes(x = eastings, y = northings),
                       color = "light grey", size = 0.5) +
            geom_point(data = mill[mill$Genus_species ==
                                     names(mill_predictions_env_brt)[i], ],
                       aes(x = eastings, y = northings), color = "orange") +
            ggtitle(paste0(names(mill_predictions_spatial_brt)[i], 
                           " - spatial model"))
    )
  }
  
  
  # environmental variables model
  mill_predictions_env_brt <- lapply(
    mill_predictions_env_brt, FUN= function(x) {
      group_by(x, hectad) %>%
        summarise(mean_prediction = mean(pred, na.rm = T), 
                  eastings = mean(eastings), northings = mean(northings))
    })
  
  for(i in 1:length(mill_predictions_env_brt)) {
    print(ggplot() + 
            geom_tile(data = mill_predictions_env_brt[[i]], 
                      aes(x = eastings, y = northings, fill = mean_prediction)) + 
            geom_point(data = mill, aes(x = eastings, y = northings),
                       color = "light grey", size = 0.5) +
            geom_point(data = mill[mill$Genus_species ==
                                     names(mill_predictions_env_brt)[i], ],
                       aes(x = eastings, y = northings), color = "orange") +
            ggtitle(paste0(names(mill_predictions_env_brt)[i], 
                           " - environmental model"))
    )
  }
  
  ## plot variable importance
  mill_var_imp_spatial_brt <- lapply(
    mill_var_imp_spatial_brt, FUN= function(x) {
      group_by(x, variable) %>%
        summarise(mean_rel.inf = mean(rel.inf))
    })
  mill_var_imp_env_brt <- lapply(
    mill_var_imp_env_brt, FUN= function(x) {
      group_by(x, variable) %>%
        summarise(mean_rel.inf = mean(rel.inf))
    })
  
  for(i in 1:length(mill_var_imp_spatial_brt)) {
    print(ggplot(data = mill_var_imp_spatial_brt[[i]], 
                 aes(x = variable, y = MeanDecreaseGini)) + 
            geom_bar(stat = "identity") + 
            coord_flip() + 
            ggtitle(paste0(names(mill_var_imp_spatial_brt)[i], 
                           " - spatial model"))
    )
  }
  
  for(i in 1:length(mill_var_imp_spatial_brt)) {
    print(ggplot(data = mill_var_imp_spatial_brt[[i]], 
                 aes(x = variable, y = mean_rel.inf)) + 
            geom_bar(stat = "identity") + 
            coord_flip() + 
            ggtitle(paste0(names(mill_var_imp_spatial_brt)[i], 
                           " - spatial model"))
    )
  }
  
  for(i in 1:length(mill_var_imp_spatial_brt)) {
    print(ggplot(data = mill_var_imp_env_brt[[i]], 
                 aes(x = variable, y = mean_rel.inf)) + 
            geom_bar(stat = "identity") + 
            coord_flip() + 
            ggtitle(paste0(names(mill_var_imp_env_brt)[i], 
                           " - environmental model"))
    )
  }
}

if(on_sonic) quit(save = "no")