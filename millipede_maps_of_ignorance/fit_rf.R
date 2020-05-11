#############################
## Fit random forest SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.  It
## relies on the outputs from 'prepare_objects_for_SDM.R'
## 
## TODO:  - make all testing be on spatially sub-sampled data (and use the
##          same spatially sub-sampled dataset to test all models)
##        - add list length model
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 23 Jan 2020
## last modified: 11 May 2020
##############################
library(randomForest)
library(pROC)
library(parallel)
library(tidyverse)
on_sonic <- F

sp_to_fit <- list("Tachypodoiulus niger", 
                  "Ommatoiulus sabulosus")
# "Ophyiulus pilosus", "Blaniulus guttulatus", "Glomeris marginata", 
# "Julus scandinavius", 
names(sp_to_fit) <- sp_to_fit

# load objects that are needed to fit SDMs
mill_wide <- readRDS("mill_wide.rds")
newdata <- readRDS("newdata.rds")

### fit random forest ---------------------------------------------------------
fit_rf <- function(test_fold, sp_name, sp_df, pred_names, newdata) {
  # ARGS: test_fold - integer giving the CV fold to be used as test data
  #       sp_name - character string giving the name of the column with 
  #                 detection/non-detection data to use
  #       sp_df - data frame with species observations and a column called 
  #               "folds" giving the cross-validation folds each record is in
  #       pred_names - character vector giving the variables to use as 
  #       predictors
  #       newdata - new data for predictions with standardized recording effort
  # browser()
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
      x = sp_df[sp_df$folds != test_fold, colnames(sp_df) %in% pred_names],
      y = factor(sp_df[sp_df$folds != test_fold, colnames(sp_df) == sp_name]),
      ntree = 1000, 
      mtry = mtry_tests[k], 
      nodesize = 3, 
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
      x = sp_df[sp_df$folds != test_fold, which(colnames(sp_df) %in% 
                                                  pred_names)],
      y = factor(sp_df[sp_df$folds != test_fold, which(colnames(sp_df) == 
                                                         sp_name)]),
      ntree = 1000, 
      mtry = mtry_best, 
      nodesize = 1, 
      replace = TRUE, classwt = NULL, 
      importance = FALSE, 
      keep.forest = TRUE)}, error = function(x) NA)
  
  # make predictions for model testing (predict only to data subset not used 
  # for model training)
  f_pred <- tryCatch(predict(mod, newdata = sp_df[sp_df$folds == test_fold, ], 
                             type = "prob")[, "1"], 
                     error = function(x) NA)
  
  # make dataframe of predictions with standardized recording effort
  newdata$pred[newdata$folds == test_fold] <- tryCatch(predict(
    mod, newdata = newdata[newdata$folds == test_fold, ], type = "prob")[, "1"], 
    error = function(x) NA)
  # return fitted model, and predictions for this model
  tryCatch(list(m = mod, predictions = newdata), 
           error = function(x) "No list exported from fit_rf.")
}


call_fit_rf <- function(fold_assignments, sp_name, test_fold, sp_df, 
                        pred_names, spatial.under.sample, ...) {
  # Function to call fit_rf on each species in a list of species dfs
  # First spatially sub-sample non-detection records to even absence records
  # and improve class balance (presuming sp. is rare)  
  # Then fit RF.
  #
  # ARGS: fold_assignments - object of class SpatialBlock with the 
  #         cross-validation fold that each observation (row) of sp_df is
  #         assigned to
  #       sp_name - character string giving the name of the column with species
  #           detection/non-detection data
  #       test_fold - integer givin the CV fold to be used for testing
  #       sp_df - data frame with species observations and predictor variables
  #       pred_names - character vector with names of columns to be used as 
  #           predictor variables
  #       spatial.under.sample - T/F indicating whether to perform spatial
  #           under sampling (Robinson et al. 2018).  If TRUE, the spatial
  #           undersampling grid cells must be provided as a column named
  #           "spat_subsamp_cell"
  browser()
  # add CV fold assignments to sp_df
  sp_df <- st_join(
    st_as_sf(sp_df), 
    st_as_sf(fold_assignments$blocks[, names(fold_assignments$blocks) ==
                                       "folds"]))
  if(is_tibble(sp_df) | "sf" %in% class(sp_df)) {
    sp_df <- data.frame(sp_df)
  }
  sp_name <- gsub(" ", ".", sp_name)
  colnames(sp_df) <- gsub(" ", ".", colnames(sp_df))
  
  if(spatial.under.sample) {
    # separate presence and absence checklists.  Keep all presence checklists.
    presences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 1, ]
    absences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 0, ]
    # spatially sub-sample absence checklists to 1 per cell
    cell_abs_tab <- table(absences$spat_subsamp_cell)
    keep_ab_rows <- c()
    for(ri in 1:length(unique(absences$spat_subsamp_cell))) {
      cell <- unique(absences$spat_subsamp_cell)[ri]
      keep_ab_rows <- c(keep_ab_rows, 
                        sample(which(absences$spat_subsamp_cell == cell), 
                               size = 1))
    }
    absences <- absences[keep_ab_rows, ] 
    
    # combine spatially sub-sampled non-detection data with all detection data
    sp_df <- bind_rows(absences, presences)
  }
 
  newdata <- st_join(st_as_sf(newdata), st_as_sf(fold_assignments$blocks))
  newdata <- data.frame(newdata)
  
  # fit RF
  lapply(1:n_folds, FUN = fit_rf, sp_name = sp_name, 
         sp_df = sp_df, pred_names = pred_names, newdata = newdata)
}


#### Fit Random Forest -------------------------------------------------------
# fit Day of Year + List Length models
day_ll_rf_fits <- sp_to_fit
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  day_ll_rf_fits[[i]] <- mclapply(fold_assignments_10k, 
                                  sp_name = sp_name, 
                                  FUN = call_fit_rf, 
                                  sp_df = mill_wide, 
                                  pred_names = c("day_of_year", "list_length"),
                                  spatial.under.sample = TRUE, 
                                  mc.cores = n_cores)
  try(print(pryr::object_size(day_ll_rf_fits)))
  try(saveRDS(day_ll_rf_fits, paste0("day_ll_rf_fits_", gsub(" ", "_", sp_name),
                                     ".rds")))
  rm(day_ll_rf_fits)
}






# fit spatial models
spatial_rf_fits <- mclapply(sp_to_fit, 
                            FUN = call_fit_rf, 
                            sp_df = mill_wide, 
                            pred_names = c("eastings", "northings", 
                                           "day_of_year", 
                                           "list_length"),
                            spatial.under.sample = TRUE, 
                            mc.cores = n_cores)
try(print(pryr::object_size(spatial_rf_fits)))
try(saveRDS(spatial_rf_fits, "spatial_rf_fits.rds"))

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
                        spatial.under.sample = TRUE, 
                        mc.cores = n_cores)
try(print(pryr::object_size(env_rf_fits)))
try(saveRDS(env_rf_fits, "env_rf_fits.rds"))
### end fit random forest ----------------------------------------------------


# view AUC (averaged over CV folds)
# rf_performance <- expand.grid(species = names(sp_to_fit), 
#                               model = c("doy_ll", "spatial", "environmental"), 
#                               stringsAsFactors = FALSE)
# rf_performance$metric <- "AUC"
# rf_performance$value <- NA
# rf_performance$value[rf_performance$model == "doy_ll"] <- sapply(
#   day_ll_rf_fits, FUN = function(x) {mean(sapply(x, FUN = function(y) {
#     tryCatch(y$auc, error = function(x) NA)}), 
#     na.rm = T)})
# rf_performance$value[rf_performance$model == "spatial"] <- sapply(
#   spatial_rf_fits, FUN = function(x) {mean(sapply(x, FUN = function(y) {
#     tryCatch(y$auc, error = function(x) NA)}), 
#     na.rm = T)})
# rf_performance$value[rf_performance$model == "environmental"] <- sapply(
#   env_rf_fits, FUN = function(x) {mean(sapply(x, FUN = function(y) {
#     tryCatch(y$auc, error = function(x) NA)}), 
#     na.rm = T)})

# get AUC for each fold & each model
rf_performance <- expand.grid(fold = 1:n_folds, 
                              species = names(sp_to_fit),
                              model = c("doy_ll", "spatial", "environmental"),
                              stringsAsFactors = FALSE)
rf_performance$metric <- "AUC"
rf_performance$value <- NA
rf_performance$value[rf_performance$model == "doy_ll"] <- sapply(
  day_ll_rf_fits, FUN = function(x) {sapply(x, FUN = function(y) {
    tryCatch(y$auc, error = function(x) NA)})})
rf_performance$value[rf_performance$model == "spatial"] <- sapply(
  spatial_rf_fits, FUN = function(x) {sapply(x, FUN = function(y) {
    tryCatch(y$auc, error = function(x) NA)})})
rf_performance$value[rf_performance$model == "environmental"] <- sapply(
  env_rf_fits, FUN = function(x) {sapply(x, FUN = function(y) {
    tryCatch(y$auc, error = function(x) NA)})})


# Get predictions with standardized survey effort
mill_predictions_spatial_rf <- lapply(
  spatial_rf_fits, 
  FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
    tryCatch(x$predictions, error = function(x) NA)})
  })
mill_predictions_spatial_rf <- lapply(
  mill_predictions_spatial_rf, 
  FUN = function(x) {
    bind_rows(x[!is.na(x)])})

mill_predictions_env_rf <- lapply(
  env_rf_fits, 
  FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
    tryCatch(x$predictions, error = function(x) NA)})})
mill_predictions_env_rf <- lapply(
  mill_predictions_env_rf, 
  FUN = function(x) {
    bind_rows(x[!is.na(x)])})

# get variable importance (averaged over 5 folds) 
mill_var_imp_spatial_rf <- lapply(
  spatial_rf_fits, 
  FUN = function(x) {
    bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
      df <- data.frame(y$m$importance)
      df$variable <- rownames(df)
      df[, c("variable", "MeanDecreaseGini")]}))}
)
mill_var_imp_env_rf <- lapply(
  env_rf_fits, 
  FUN = function(x) {
    bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
      df <- data.frame(y$m$importance)
      df$variable <- rownames(df)
      df[, c("variable", "MeanDecreaseGini")]}))}
)

# save results
try(saveRDS(mill_predictions_spatial_rf, "mill_predictions_spatial_rf.rds"))
try(saveRDS(mill_predictions_env_rf, "mill_predictions_env_rf.rds"))



### exploratory plotting - this is not for running on sonic

ggplot(data = rf_performance, aes(x = model, y = value, color = species)) + 
  geom_point() + 
  geom_boxplot() + 
  facet_wrap(~metric) + 
  theme_bw()


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
      geom_point(data = mill, aes(x = eastings, y = northings),
                 color = "light grey", size = 0.5) +
      geom_point(data = mill[mill$Genus_species ==
                               names(mill_predictions_spatial_rf)[i], ],
                 aes(x = eastings, y = northings), color = "orange") +
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
          geom_point(data = mill, aes(x = eastings, y = northings),
                     color = "light grey", size = 0.5) +
          geom_point(data = mill[mill$Genus_species ==
                                   names(mill_predictions_env_rf)[i], ],
                     aes(x = eastings, y = northings), color = "orange") +
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
                         " - environmental model"))
  )
}

if(on_sonic) quit(save = "no")
