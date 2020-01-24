#############################
## Fit SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_objects_for_SDM.R'
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 23 Jan 2020
## last modified: 24 Jan 2020
##############################
library(dismo)
set.seed(01242020) # Jan 24 2020

# load objects that are needed to fit SDMs
J_scand <- readRDS("J_scand.rds")
fold_index_10k <- readRDS("fold_index_10k.rds")
newdata <- readRDS("newdata.rds")

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
                           tree.complexity = 3, learning.rate = 0.01, 
                           n.trees = 500, step.size = 500,
                           family = "bernoulli", max.trees = 2000, plot.main = F), 
                  error = function(x) NA)
  
  # make predictions for measuring AUC (predict only to data subset used for 
  # overall model training)
  f_pred <- tryCatch(predict(f_m, newdata = sp_df[fold_index == test_fold, ], 
                             type = "response",
                             n.trees = f_m$gbm.call$best.trees), 
                     error = function(x) NA)
  f_auc <- tryCatch(
    roc(response = factor(sp_df[fold_index == test_fold, "Julus scandinavius"], 
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

# fit model
J_scand_brt_fits <- mclapply(1:5, FUN = fit_brt, 
                             fold_index = fold_index_10k, 
                             sp_df = J_scand, 
                             pred_names = c("eastings", "northings", 
                                            "list_length"), 
                             mc.cores = n_cores)

### end fit brt ---------------------------------------------------------------

J_scand_preds <- bind_rows(lapply(J_scand_brt_fits, 
                        FUN = function(x) {x$predictions}))