#############################
## Fit random forest SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_objects_for_SDM.R'
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 23 Jan 2020
## last modified: 24 Jan 2020
##############################
library(randomForest)
library(pROC)
set.seed(01242020) # Jan 24 2020
on_sonic <- F
n_cores <- 2

# load objects that are needed to fit SDMs
mill_wide <- readRDS("mill_wide.rds")
newdata <- readRDS("newdata.rds")

### fit random forest ---------------------------------------------------------
fit_rf <- function(test_fold, sp_name, sp_df, pred_names) {
  # ARGS: test_fold - integer giving the CV fold to be used as test data
  #       sp_name - character string giving the name of the column with 
  #                 detection/non-detection data to use
  #       sp_df - data frame with species observations and a column called 
  #               "folds" giving the cross-validation folds each record is in
  #       pred_names - character vector giving the variables to use as predictors
  browser()
  if(is_tibble(sp_df)) {
    sp_df <- data.frame(sp_df)
  }
  sp_name <- gsub(" ", ".", sp_name)
  colnames(sp_df) <- gsub(" ", ".", colnames(sp_df))
  # optimize number of variables to use at each node
  mtry_err <- c() # make list to hold errors for the mtry tests
  nvar <- length(pred_names) # number of variables
  # mtry values to test - default, half that, and twice that
  mtry_tests <- c(ceiling(sqrt(nvar)/2), ceiling(sqrt(nvar)), 
                  ceiling(sqrt(nvar)*2)) 
  
  for(k in 1:length(mtry_tests)) {
    m_k <- tryCatch(randomForest(
      x = sp_df[, colnames(sp_df) %in% pred_names],
      y = factor(sp_df[, colnames(sp_df) == sp_name]),
      ntree = 2000, 
      mtry = mtry_tests[k], 
      nodesize = 1, 
      replace = TRUE, classwt = NULL, 
      importance = FALSE, 
      keep.forest = FALSE), 
      error = function(x) NA)
    
    # if error from this model is lowest so far, keep this model
    if(class(m_k) == "randomForest" && 
       !is.na(m_k$err.rate[nrow(m_k$err.rate), "OOB"])) {
      mtry_err[k] <- m_k$err.rate[nrow(m_k$err.rate), "OOB"] # error for this mtry
    }
  }
  
  # get best mtry value.  If multiply mtry values tied for best error rate,
  # use the one that is closest to the square root of the number of variables
  if(length(which(mtry_err == min(mtry_err))) > 1) {
    # calculate the distance of each tied "best" mtry value from the default
    dist_from_default <- tryCatch({
      abs(mtry_tests[mtry_err == min(mtry_err)] - sqrt(nvar))}, 
      error = function(x) NA)
    # keep the mtry value that is closest to the default
    mtry_best <- tryCatch({
      mtry_tests[mtry_err == min(mtry_err)][dist_from_default == 
                                              min(dist_from_default)]}, 
      error = function(x) NA)
  } else mtry_best <- tryCatch(mtry_tests[mtry_err == min(mtry_err)], 
                               error = function(x) NA)
  
  # fit model with optimum mtry 
  # use 2000 trees which is hopefully high enough to get stable variable
  # importance if I want it (see manual linked in help documentation)
  mod <- tryCatch({
    randomForest(
      x = sp_df[, which(colnames(sp_df) %in% pred_names)],
      y = factor(sp_df[, which(colnames(sp_df) == sp_name)]),
      ntree = 2000, 
      mtry = mtry_best, 
      nodesize = 1, 
      replace = TRUE, classwt = NULL, 
      importance = FALSE, 
      keep.forest = TRUE)}, error = function(x) NA)
  
  # make predictions for measuring AUC (predict only to data subset not used for 
  # model training)
  f_pred <- tryCatch(predict(mod, newdata = sp_df[sp_df$folds == test_fold, ], 
                             type = "prob")[, "1"], 
                     error = function(x) NA)
  f_auc <- tryCatch(
    roc(response = factor(sp_df[sp_df$folds == test_fold, sp_name], 
                          levels = c("0", "1")),
        predictor = f_pred,
        auc = T), 
    error = function(x) NA)
  
  # make dataframe of predictions with standardized recording effort
  # make predictions to ALL grid cells (both test and training)
  newdata$pred[newdata$folds == test_fold] <- tryCatch(predict(
    mod, newdata = newdata[newdata$folds == test_fold, ], type = "prob")[, "1"], 
    error = function(x) NA)
  # return fitted model, AUC value, and predictions for this model
  tryCatch(list(m = mod, auc = f_auc$auc[1], predictions = newdata, 
                proc_object = f_auc), 
           error = function(x) "No list exported from fit_rf.")
}


call_fit_rf <- function(sp_name, test_fold, sp_df, pred_names, ...) {
  # Function to call fit_rf on each species in a list of species dfs
  lapply(1:5, FUN = fit_rf, sp_name = sp_name, 
         sp_df = sp_df, pred_names = pred_names)
}

# fit spatial models
sp_to_fit <- list("Julus scandinavius", "Glomeris marginata")
# "Ophyiulus pilosus", "Blaniulus guttulatus", 
# "Tachypodoiulus niger"

names(sp_to_fit) <- sp_to_fit
spatial_rf_fits <- mclapply(sp_to_fit, 
                            FUN = call_fit_rf, 
                            sp_df = mill_wide, 
                            pred_names = c("eastings", "northings", 
                                           "day_of_year", 
                                           "list_length"), 
                            mc.cores = n_cores)
# fit environmental model
env_rf_fits <- mclapply(sp_to_fit, 
                        FUN = call_fit_rf, 
                        sp_df = mill_wide, 
                        pred_names = c("mean_tn", "mean_tx", 
                                       "mean_rr", "artificial_surfaces", 
                                       "forest_seminatural_l1", 
                                       "wetlands_l1", "pasture_l2", 
                                       "arable_l2", "elev", 
                                       "day_of_year", 
                                       "list_length"), 
                        mc.cores = n_cores)

### end fit random forest ----------------------------------------------------

# view AUC (averaged over 5 folds)
lapply(spatial_rf_fits, 
            FUN = function(x) {mean(sapply(x, FUN = function(y) {y$auc}), 
                                    na.rm = T)})
lapply(env_rf_fits, 
       FUN = function(x) {mean(sapply(x, FUN = function(y) {y$auc}), 
                               na.rm = T)})

# Get predictions with standardized survey effort
mill_predictions_spatial_rf <- lapply(
  spatial_rf_fits, 
  FUN = function(x) {
    bind_rows(lapply(x, FUN = function(x) x$predictions))
  })
mill_predictions_env_rf <- lapply(
  env_rf_fits, 
  FUN = function(x) {
    bind_rows(lapply(x, FUN = function(x) x$predictions))
  })

# get variable importance (averaged over 5 folds) 
mill_var_imp_spatial_rf <- lapply(
  spatial_rf_fits, 
  FUN = function(x) {
    bind_rows(lapply(x, FUN = function(y) {
      df <- data.frame(y$m$importance)
      df$variable <- rownames(df)
      df[, c("variable", "MeanDecreaseGini")]}))}
)
mill_var_imp_env_rf <- lapply(
  env_rf_fits, 
  FUN = function(x) {
    bind_rows(lapply(x, FUN = function(y) {
      df <- data.frame(y$m$importance)
      df$variable <- rownames(df)
      df[, c("variable", "MeanDecreaseGini")]}))}
)

# save results
saveRDS(mill_predictions_spatial_rf, "mill_predictions_spatial_rf.rds")
saveRDS(mill_predictions_env_rf, "mill_predictions_env_rf.rds")



### exploratory plotting - this is not for running on sonic
# plot average of predictions from all 5 folds (so 4 predictions will be to
# training data, one prediction to test data in each grid cell)

# spatial model
mill_predictions_spatial_rf <- lapply(
  mill_predictions_spatial_rf, FUN= function(x) {
    group_by(x, hectad) %>%
      summarise(mean_prediction = mean(pred, na.rm = T), 
                eastings = mean(eastings), northings = mean(northings))
  })

for(i in 1:length(mill_predictions_spatial_rf)) {
  print(ggplot() + 
    geom_tile(data = mill_predictions_spatial_rf[[i]], 
              aes(x = eastings, y = northings, fill = mean_prediction)) + 
  # geom_point(data = J_scand, aes(x = eastings, y = northings, 
  #                                color = factor("Julus scandinavius"))) + 
    ggtitle(paste0(names(mill_predictions_spatial_rf)[i], 
                   " - spatial model"))
  )
}


# environmental variables model
mill_predictions_env_rf <- lapply(
  mill_predictions_env_rf, FUN= function(x) {
    group_by(x, hectad) %>%
      summarise(mean_prediction = mean(pred, na.rm = T), 
                eastings = mean(eastings), northings = mean(northings))
  })

for(i in 1:length(mill_predictions_env_rf)) {
  print(ggplot() + 
          geom_tile(data = mill_predictions_env_rf[[i]], 
                    aes(x = eastings, y = northings, fill = mean_prediction)) + 
          # geom_point(data = mill, aes(x = eastings, y = northings),
          #            color = "light grey", size = 0.5) +
          # geom_point(data = mill[mill$Genus_species ==
          #                          names(mill_predictions_env_rf)[i], ],
          #            aes(x = eastings, y = northings), color = "orange") +
          ggtitle(paste0(names(mill_predictions_env_rf)[i], 
                         " - environmental model"))
  )
}

## plot variable importance
mill_var_imp_spatial_rf <- lapply(
  mill_var_imp_spatial_rf, FUN= function(x) {
    group_by(x, variable) %>%
      summarise(MeanDecreaseGini = mean(MeanDecreaseGini))
  })
mill_var_imp_env_rf <- lapply(
  mill_var_imp_env_rf, FUN= function(x) {
    group_by(x, variable) %>%
      summarise(MeanDecreaseGini = mean(MeanDecreaseGini))
  })

for(i in 1:length(mill_var_imp_spatial_rf)) {
  print(ggplot(data = mill_var_imp_spatial_rf[[i]], 
               aes(x = variable, y = MeanDecreaseGini)) + 
          geom_bar(stat = "identity") + 
          coord_flip() + 
          ggtitle(paste0(names(mill_var_imp_spatial_rf)[i], 
                         " - spatial model"))
  )
}

for(i in 1:length(mill_var_imp_env_rf)) {
  print(ggplot(data = mill_var_imp_env_rf[[i]], 
               aes(x = variable, y = MeanDecreaseGini)) + 
          geom_bar(stat = "identity") + 
          coord_flip() + 
          ggtitle(paste0(names(mill_var_imp_env_rf)[i], 
                         " - spatial model"))
  )
}

if(on_sonic) quit(save = "no")

