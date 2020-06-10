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
## last modified: 10 June 2020
##############################

on_sonic <- F

# load objects that are needed to fit SDMs
mill_wide <- readRDS("mill_wide.rds")
newdata <- readRDS("newdata.rds")
fold_assignments <- readRDS("fold_assignments.rds")
block_subsamp_10k <- readRDS("block_subsamp_10k.rds")

### fit random forest ---------------------------------------------------------
fit_rf <- function(test_fold, sp_name, sp_df, pred_names, newdata, 
                   sp_df_original, mtry, block_cv_range) {
  # ARGS: test_fold - integer giving the CV fold to be used as test data
  #       sp_name - character string giving the name of the column with 
  #                 detection/non-detection data to use
  #       sp_df - data frame to be used for model fitting (possibly 
  #               undersampled), with species observations and a column called 
  #               "folds" giving the cross-validation folds each record is in
  #       pred_names - character vector giving the variables to use as 
  #       predictors
  #       newdata - new data for predictions with standardized recording effort
  #       sp_df_original - the original (not undersampled) data
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
 
  # label which observations are in the test fold
  sp_df$test_fold <- sp_df$folds == test_fold
  sp_df_original$test_fold <- sp_df_original$folds == test_fold

  # ### ---- test mtry. (if I want to do this later) ---------------------
  # ### As of 20 May, I will not do this and will just use the square root of the
  # ### number of predictor variables, following Robinson et al. (2018) and 
  # ### references in that
  # for(k in 1:length(mtry_tests)) {
  #   m_k <- tryCatch(randomForest(
  #     x = sp_df[sp_df$folds != test_fold, colnames(sp_df) %in% pred_names],
  #     y = factor(sp_df[sp_df$folds != test_fold, colnames(sp_df) == sp_name]),
  #     ntree = 1000, 
  #     mtry = mtry_tests[k], 
  #     nodesize = 3, 
  #     replace = TRUE, classwt = NULL, 
  #     importance = FALSE, 
  #     keep.forest = FALSE), 
  #     error = function(x) NA)
  #   
  #   # if error from this model is lowest so far, keep this model
  #   if(class(m_k) == "randomForest" && 
  #      !is.na(m_k$err.rate[nrow(m_k$err.rate), "OOB"])) {
  #     mtry_err[k] <- m_k$err.rate[nrow(m_k$err.rate), "OOB"] # error for this mtry
  #   }
  # }
  # 
  # # get best mtry value.  If multiply mtry values tied for best error rate,
  # # use the one that is closest to the square root of the number of variables
  # if(length(which(mtry_err == min(mtry_err))) > 1) {
  #   # calculate the distance of each tied "best" mtry value from the default
  #   dist_from_default <- tryCatch({
  #     abs(mtry_tests[mtry_err == min(mtry_err)] - sqrt(nvar))}, 
  #     error = function(x) NA)
  #   # keep the mtry value that is closest to the default
  #   mtry_best <- tryCatch({
  #     mtry_tests[mtry_err == min(mtry_err)][dist_from_default == 
  #                                             min(dist_from_default)]}, 
  #     error = function(x) NA)
  # } else mtry_best <- tryCatch(mtry_tests[mtry_err == min(mtry_err)], 
  #                              error = function(x) NA)
  
  # use 2000 trees which is hopefully high enough to get stable variable
  # importance if I want it (see manual linked in help documentation)
  mod <- tryCatch({
    randomForest(
      x = sp_df[sp_df$folds != test_fold, which(colnames(sp_df) %in% 
                                                  pred_names)],
      y = factor(sp_df[sp_df$folds != test_fold, which(colnames(sp_df) == 
                                                         sp_name)]),
      ntree = 1000, 
      mtry = mtry,   # use mtry_best if I want to fit model with optimum mtry 
      nodesize = 1, 
      replace = TRUE, classwt = NULL, 
      importance = FALSE, 
      keep.forest = TRUE)}, error = function(x) NA)

  # make predictions for model testing (predict to ALL folds, and predict to 
  # ALL locations for which there is data, not just the spatially 
  # undersampled locations)
  # For some reason, predict will fail if any columns have NAs, even if those
  # are not predictor columns.  Because in some rare cases a fold is not 
  # assigned to some cells (b/c they are on the boundary of a fold), there are
  # occasionally NAs in the "folds" column, so that column should not be passed
  # to predict.
  f_pred <- tryCatch({
    predict(mod, newdata = sp_df_original[, which(
      colnames(sp_df_original) != "folds")], 
            type = "prob")[, "1"]}, 
    error = function(x) NA)
  # select columns to keep in df of predictions
  preds <- sp_df_original[ , c("checklist_ID", "eastings", "northings", 
                               "hectad", "folds", "test_fold", sp_name)]
  preds$pred <- f_pred # add predictions to df

  # make dataframe of predictions with standardized recording effort
  # predict to ALL sites (both in training and test sets)
  newdata$pred <- tryCatch({
    predict(mod, newdata = newdata[, which(colnames(newdata) != "folds")], 
            type = "prob")[, "1"]}, 
    error = function(x) NA)
  
  # label test folds in case this is ever needed later
  if("folds" %in% colnames(newdata)) {
    newdata$test_fold <- tryCatch(newdata$folds == test_fold, 
                                  error = function(x) NA)
  } else {
    newdata$folds <- NA
    newdata$test_fold <- NA}
  
  # drop predictor variables from predictions dataframe 
  newdata <- select(newdata, hectad, eastings, northings, cells,
                    list_length, day_of_year, folds, test_fold, pred)
  
  ## calculate class balance (proportion of checklists with a detection)
  prop_dets <- tryCatch({
    length(which(sp_df[sp_df$folds != test_fold, 
                                  which(colnames(sp_df) == sp_name)] != 0)) / 
    nrow(sp_df)}, error = function (x) NA)
  
  ## calculate Simpson's evenness for training and test datasets
  # Only calculate this for fold # 1.  The same dataset is used in multiple 
  # folds, so the value will be identical for all folds using that dataset.
  if(test_fold == 1) {
    table_nobs <- table(sp_df$hectad)
    simps_train <- simpson_even(as.numeric(table_nobs))
  } else {
    simps_train <- NA
  }
  
  
  # return fitted model, and predictions for this model
  tryCatch(list(
    m = mod, preds = preds, standardized_preds = newdata, 
    train_sites = unique(newdata$hectad[newdata$folds != test_fold]), 
    test_sites = unique(newdata$hectad[newdata$folds == test_fold]), 
    block_cv_range = block_cv_range, 
    n_detections_test_fold = try(sum(sp_df[sp_df$folds == test_fold, 
                                           sp_name])),
    simpson_training = simps_train, 
    proportion_detections = prop_dets), 
    error = function(x) "No list exported from fit_rf.")
}


call_fit_rf <- function(fold_assignments, sp_name, test_fold, sp_df, 
                        pred_names, spatial.under.sample, block_subsamp, mtry, 
                        ...) {
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

  # add CV fold assignments to sp_df
  if(class(fold_assignments) == "SpatialBlock") {
    sp_df <- st_join(
      st_as_sf(sp_df), 
      st_as_sf(fold_assignments$blocks[, names(fold_assignments$blocks) ==
                                         "folds"]))
  }
  if(class(fold_assignments) == "list") { # add random CV blocks
    sp_df$folds <- fold_assignments$blocks
  }
  
  # convert to df from spatial
  if(is_tibble(sp_df) | "sf" %in% class(sp_df) | 
     "SpatialPointsDataFrame" %in% class(sp_df)) {
    sp_df <- data.frame(sp_df)
  }
  colnames(sp_df) <- gsub(" ", ".", colnames(sp_df))
  sp_name <- gsub(" ", ".", sp_name)
 
  # convert species record counts to p/a
  sp_df[colnames(sp_df) == sp_name] <- pa(sp_df[colnames(sp_df) == sp_name])
  sp_df_original <- sp_df
  
  if(spatial.under.sample) {
    # add a column of subsampling block assignments by chosing randomly from
    # the many allocatins created in "prepare_objects_for_SDM.R"
    sp_df$spat_subsamp_cell <- block_subsamp[, sample(2:(ncol(block_subsamp)-1), 
                                                         size = 1)]
    # spatially sub-sample absence checklists to 1 per cell
    # separate presence and absence checklists.  Keep all presence checklists.
    presences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 1, ]
    absences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 0, ]
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

  # if using block CV, assign newdata grid cells to CV fold
  # if using random CV, the folds pertain to checklists, not locations, so all
  # newdata locations are essentially in test set
  if(class(fold_assignments) == "SpatialBlock") {
    newdata <- st_join(
      st_as_sf(newdata), 
      st_as_sf(fold_assignments$blocks[, names(fold_assignments$blocks) == 
                                         "folds"]))
  }
  newdata <- data.frame(newdata)
  
  # fit RF
  lapply(1:n_folds, FUN = fit_rf, sp_name = sp_name, 
         sp_df = sp_df, pred_names = pred_names, newdata = newdata, 
         sp_df_original = sp_df_original, mtry = mtry, 
         block_cv_range = fold_assignments$range)
}


#### Fit Random Forest -------------------------------------------------------
### fit Day of Year + List Length models -------------------------------------
# train with raw data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  day_ll_rf_fits <- mclapply(fold_assignments, 
                             sp_name = sp_name, 
                             FUN = call_fit_rf, 
                             sp_df = mill_wide, 
                             pred_names = c("day_of_year", "list_length"),
                             block_subsamp = block_subsamp_10k, 
                             spatial.under.sample = FALSE, 
                             mtry = 1, 
                             mc.cores = n_cores)
  try(print(pryr::object_size(day_ll_rf_fits)))
  try(saveRDS(day_ll_rf_fits, paste0("day_ll_rf_noSubSamp_fits_", 
                                     gsub(" ", "_", sp_name),
                                     ".rds")))
  rm(day_ll_rf_fits)
}

# train with spatially subsampled data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  day_ll_rf_fits <- mclapply(fold_assignments, 
                             sp_name = sp_name, 
                             FUN = call_fit_rf, 
                             sp_df = mill_wide, 
                             pred_names = c("day_of_year", "list_length"),
                             block_subsamp = block_subsamp_10k, 
                             spatial.under.sample = TRUE, 
                             mtry = 1,
                             mc.cores = n_cores)
  try(print(pryr::object_size(day_ll_rf_fits)))
  try(saveRDS(day_ll_rf_fits, paste0("day_ll_rf_SubSamp_fits_", 
                                     gsub(" ", "_", sp_name),
                                     ".rds")))
  rm(day_ll_rf_fits)
}
### end Day of Year + List Length --------------------------------------------

### fit Spatial + List Length + DOY models -----------------------------------
# train with raw data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  spat_rf_fits <- mclapply(fold_assignments, 
                           sp_name = sp_name, 
                           FUN = call_fit_rf, 
                           sp_df = mill_wide, 
                           pred_names = c("eastings", "northings", 
                                          "day_of_year", "list_length"),
                           block_subsamp = block_subsamp_10k, 
                           spatial.under.sample = FALSE, 
                           mtry = 2, 
                           mc.cores = n_cores)
  try(print(pryr::object_size(spat_rf_fits)))
  try(saveRDS(spat_rf_fits, paste0("spat_ll_rf_noSubSamp_fits_", 
                                     gsub(" ", "_", sp_name),
                                     ".rds")))
  
  # retrieve predictions with standardized sampling effort
  mill_predictions <- lapply(
    spat_rf_fits, 
    FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
      tryCatch(x$standardized_preds, error = function(x) NA)})
    })
  mill_predictions <- bind_rows(lapply(mill_predictions, 
                             FUN = function(x) {bind_rows(x[!is.na(x)])}))
  
  try(saveRDS(mill_predictions, paste0("mill_predictions_spat_ll_rf_noSubSamp_", 
                                       gsub(" ", "_", sp_name),
                                       ".rds")))
  rm(spat_rf_fits, mill_predictions)
}

# train with spatially undersampled data 
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  spat_rf_fits <- mclapply(fold_assignments, 
                           sp_name = sp_name, 
                           FUN = call_fit_rf, 
                           sp_df = mill_wide, 
                           pred_names = c("eastings", "northings", 
                                          "day_of_year", "list_length"),
                           block_subsamp = block_subsamp_10k, 
                           spatial.under.sample = TRUE, 
                           mtry = 2, 
                           mc.cores = n_cores)
  try(print(pryr::object_size(spat_rf_fits)))
  try(saveRDS(spat_rf_fits, paste0("spat_ll_rf_SubSamp_fits_", 
                                   gsub(" ", "_", sp_name),
                                   ".rds")))
  
  # retrieve predictions with standardized sampling effort
  mill_predictions <- lapply(
    spat_rf_fits, 
    FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
      tryCatch(x$standardized_preds, error = function(x) NA)})
    })
  mill_predictions <- bind_rows(
    lapply(mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
  
  try(saveRDS(mill_predictions, paste0("mill_predictions_spat_ll_rf_SubSamp_", 
                                       gsub(" ", "_", sp_name),
                                       ".rds")))
  
  rm(spat_rf_fits, mill_predictions)
}
### end spatial + List Length + DOY models ------------------------------------


### fit environmental + LL + DOY model ----------------------------------------
# train with raw data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  env_rf_fits <- mclapply(fold_assignments, 
                          sp_name = sp_name, 
                          FUN = call_fit_rf, 
                          sp_df = mill_wide, 
                          pred_names = c("mean_tn", "mean_tx", 
                                         "mean_rr", "artificial_surfaces", 
                                         "forest_seminatural_l1", 
                                         "wetlands_l1", "pasture_l2", 
                                         "arable_l2", "elev", 
                                         "day_of_year", "list_length"),
                          block_subsamp = block_subsamp_10k, 
                          spatial.under.sample = FALSE, 
                          mtry = 3, 
                          mc.cores = n_cores)
  try(print(pryr::object_size(env_rf_fits)))
  try(saveRDS(env_rf_fits, paste0("env_ll_rf_noSubSamp_fits_", 
                                  gsub(" ", "_", sp_name),
                                  ".rds")))
  
  # retrieve predictions with standardized sampling effort
  mill_predictions <- lapply(
    env_rf_fits, 
    FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
      tryCatch(x$standardized_preds, error = function(x) NA)})
    })
  mill_predictions <- bind_rows(lapply(mill_predictions, 
                                       FUN = function(x) {bind_rows(x[!is.na(x)])}))
  
  try(saveRDS(mill_predictions, paste0("mill_predictions_env_ll_rf_noSubSamp_", 
                                       gsub(" ", "_", sp_name),
                                       ".rds")))
  rm(env_rf_fits, mill_predictions)
}

# train with spatially subsampled data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  env_rf_fits <- mclapply(fold_assignments, 
                          sp_name = sp_name, 
                          FUN = call_fit_rf, 
                          sp_df = mill_wide, 
                          pred_names = c("mean_tn", "mean_tx", 
                                         "mean_rr", "artificial_surfaces", 
                                         "forest_seminatural_l1", 
                                         "wetlands_l1", "pasture_l2", 
                                         "arable_l2", "elev", 
                                         "day_of_year", "list_length"),
                          block_subsamp = block_subsamp_10k, 
                          spatial.under.sample = TRUE, 
                          mtry = 3, 
                          mc.cores = n_cores)
  try(print(pryr::object_size(env_rf_fits)))
  try(saveRDS(env_rf_fits, paste0("env_ll_rf_SubSamp_fits_", 
                                  gsub(" ", "_", sp_name),
                                  ".rds")))
  
  # retrieve predictions with standardized sampling effort
  mill_predictions <- lapply(
    env_rf_fits, 
    FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
      tryCatch(x$standardized_preds, error = function(x) NA)})
    })
  mill_predictions <- bind_rows(
    lapply(mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
  
  try(saveRDS(mill_predictions, paste0("mill_predictions_env_ll_rf_SubSamp_", 
                                       gsub(" ", "_", sp_name),
                                       ".rds")))
  rm(env_rf_fits, mill_predictions)
}
### end environmental + LL + DOY model ------------------------------------

### fit environmental + lat + long + LL + DOY model ----------------------------
# train with raw data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  env_rf_fits <- mclapply(fold_assignments, 
                          sp_name = sp_name, 
                          FUN = call_fit_rf, 
                          sp_df = mill_wide, 
                          pred_names = c("mean_tn", "mean_tx", 
                                         "mean_rr", "artificial_surfaces", 
                                         "forest_seminatural_l1", 
                                         "wetlands_l1", "pasture_l2", 
                                         "arable_l2", "elev", 
                                         "eastings", "northings", 
                                         "day_of_year", "list_length"),
                          block_subsamp = block_subsamp_10k, 
                          spatial.under.sample = FALSE, 
                          mtry = 4, 
                          mc.cores = n_cores)
  try(print(pryr::object_size(env_rf_fits)))
  try(saveRDS(env_rf_fits, paste0("env_spat_ll_rf_noSubSamp_fits_", 
                                  gsub(" ", "_", sp_name),
                                  ".rds")))
  
  # retrieve predictions with standardized sampling effort
  mill_predictions <- lapply(
    env_rf_fits, 
    FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
      tryCatch(x$standardized_preds, error = function(x) NA)})
    })
  mill_predictions <- bind_rows(lapply(mill_predictions, 
                                       FUN = function(x) {bind_rows(x[!is.na(x)])}))
  
  try(saveRDS(mill_predictions, 
              paste0("mill_predictions_env_spat_ll_rf_noSubSamp_", 
                     gsub(" ", "_", sp_name),
                     ".rds")))
  rm(env_rf_fits, mill_predictions)
}

# train with spatially subsampled data
for(i in 1:length(sp_to_fit)) {
  sp_name <- names(sp_to_fit)[i]
  env_rf_fits <- mclapply(fold_assignments, 
                          sp_name = sp_name, 
                          FUN = call_fit_rf, 
                          sp_df = mill_wide, 
                          pred_names = c("mean_tn", "mean_tx", 
                                         "mean_rr", "artificial_surfaces", 
                                         "forest_seminatural_l1", 
                                         "wetlands_l1", "pasture_l2", 
                                         "arable_l2", "elev", 
                                         "eastings", "northings", 
                                         "day_of_year", "list_length"),
                          block_subsamp = block_subsamp_10k, 
                          spatial.under.sample = TRUE, 
                          mtry = 4, 
                          mc.cores = n_cores)
  try(print(pryr::object_size(env_rf_fits)))
  try(saveRDS(env_rf_fits, paste0("env_spat_ll_rf_SubSamp_fits_", 
                                  gsub(" ", "_", sp_name),
                                  ".rds")))
  
  # retrieve predictions with standardized sampling effort
  mill_predictions <- lapply(
    env_rf_fits, 
    FUN = function(x) {lapply(x[sapply(x, is.list)], FUN = function(x) {
      tryCatch(x$standardized_preds, error = function(x) NA)})
    })
  mill_predictions <- bind_rows(
    lapply(mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
  
  try(saveRDS(mill_predictions, 
              paste0("mill_predictions_env_spat_ll_rf_SubSamp_", 
                     gsub(" ", "_", sp_name), ".rds")))
  rm(env_rf_fits, mill_predictions)
}
### end environmental + LL + DOY model ------------------------------------
### end fit random forest ----------------------------------------------------


for(mod_name in mod_names) {
  for(i in 1:length(sp_to_fit)) {
    ### trained with raw data
    fits <- readRDS(paste0(mod_name, "_noSubSamp_fits_", 
                           gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
    sp_name <- names(sp_to_fit)[i]
    
    # Get predictions with standardized survey effort
    mill_predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$standardized_preds, error = function(x) NA)
        preds$cv <- as.character(x$block_cv_range)
        return(preds)})})
    mill_predictions <- bind_rows(lapply(
      mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(mill_predictions, 
                paste0("standard_predictions_", mod_name, "_noSubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), ".rds")))
    
    # get variable importance (averaged over 5 folds) 
    var_imp <- lapply(fits, FUN = function(x) {
      bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
        df <- data.frame(y$m$importance)
        df$variable <- rownames(df)
        df[, c("variable", "MeanDecreaseGini")]
        df$species <- sp_name
        df$model <- mod_name
        df$train_dat <- "raw"
        df$cv <- as.character(y$block_cv_range)
        return(df)}))}
    )
    var_imp <- bind_rows(var_imp)
    # save results
    try(saveRDS(var_imp, 
                paste0("var_importance_", mod_name, "_noSubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), ".rds")))
    rm(fits)
    
    ### trained with spatial subsampling
    fits <- readRDS(paste0(mod_name, "_SubSamp_fits_", 
                           gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
    
    # Get predictions with standardized survey effort
    mill_predictions <- lapply(fits, FUN = function(x) {
      lapply(x[sapply(x, is.list)], FUN = function(x) {
        preds <- tryCatch(x$standardized_preds, error = function(x) NA)
        preds$cv <- as.character(x$block_cv_range)
        return(preds)})})
    mill_predictions <- bind_rows(lapply(
      mill_predictions, FUN = function(x) {bind_rows(x[!is.na(x)])}))
    # save results
    try(saveRDS(mill_predictions, 
                paste0("standard_predictions_", mod_name, "_SubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), ".rds")))
    
    # get variable importance (averaged over 5 folds) 
    var_imp <- lapply(fits, FUN = function(x) {
      bind_rows(lapply(x[sapply(x, is.list)], FUN = function(y) {
        df <- data.frame(y$m$importance)
        df$variable <- rownames(df)
        df[, c("variable", "MeanDecreaseGini")]
        df$species <- sp_name
        df$model <- mod_name
        df$train_dat <- "spat_subsamp"
        df$cv <- as.character(y$block_cv_range)
        return(df)}))}
    )
    var_imp <- bind_rows(var_imp)
    # save results
    try(saveRDS(var_imp, 
                paste0("var_importance_", mod_name, "_SubSamp_", 
                       gsub(" ", "_", sp_to_fit[[i]]), ".rds")))
    rm(fits)
  }
}






if(on_sonic) quit(save = "no")
