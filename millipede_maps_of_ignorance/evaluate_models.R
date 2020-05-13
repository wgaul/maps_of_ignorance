#############################
## Evaluate SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
## 
## TODO:  - test on original data
##        - test on spatially sub-sampled data (and use the
##          same spatially sub-sampled dataset to test all models)
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 12 May 2020
## last modified: 13 May 2020
##############################
library(pROC)

evals <- data.frame() # data frame to hold evaluation results

## Get spatially subsampled test points ---------------------------------------
test_points_ss <- sapply(sp_to_fit, FUN = function(x) NULL)
names(test_points_ss) <- gsub(" ", ".", names(test_points_ss))
mill_wide_df <- data.frame(mill_wide)
for(i in 1:length(test_points_ss)) {
  mill_wide_df[, colnames(mill_wide_df) == names(test_points_ss)[i]] <- pa(
    mill_wide_df[, colnames(mill_wide_df) == names(test_points_ss)[i]])
  # separate presence and absence checklists.  Keep all presence checklists.
  presences <- mill_wide_df[mill_wide_df[colnames(mill_wide_df) == 
                                           names(test_points_ss)[i]] > 0, ]
  absences <- mill_wide_df[mill_wide_df[colnames(mill_wide_df) == 
                                          names(test_points_ss)[i]] == 0, ]
  
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
  test_points_ss[[i]]$test_points_subsampled <- bind_rows(absences, presences)
}
## end get spatially subsampled test points ------------------------------------


## evaluate Day of Year + List Length models ----------------------------------
# trained with no spatial subsampling
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("day_ll_rf_noSubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original data 
  # AUC --------------------------------------------------
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "day_ll_rf", 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev) # end AUC --------------------------------------
  
  # sensitivity --------------------------------------
  
  # end sensitivity ----------------------------------
  
  ## test against spatially subsampled data
  aucs <- lapply(fits, FUN = function(x, sp_name, test_data) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_data) {
      # get only subsampled test data that is in the spatial test blocks
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" are in different CV folds b/c cv block 
      # boundary is within the hectad and checklists are at finer spatial 
      # resolution)
      test_data <- test_data[test_data$test_fold == T, ] 
      tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name, test_data = test_data)
  }, sp_name = sp_to_fit[[i]], 
  test_data = test_points_ss[[i]]$test_points_subsampled)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "day_ll_rf", 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  
  
  
  rm(ev, fits)
}

# trained with spatial subsampling
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("day_ll_rf_SubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original data 
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "day_ll_rf", 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev)
  
  ## test against spatially subsampled data
  aucs <- lapply(fits, FUN = function(x, sp_name, test_data) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_data) {
      # get only subsampled test data that is in the spatial test blocks
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" are in different CV folds b/c cv block 
      # boundary is within the hectad and checklists are at finer spatial 
      # resolution)
      test_data <- test_data[test_data$test_fold == T, ] 
      tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                gsub(" ", ".", sp_name)], 
                           predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name, test_data = test_data)
  }, sp_name = sp_to_fit[[i]], 
  test_data = test_points_ss[[i]]$test_points_subsampled)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "day_ll_rf", 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev, fits)
}
### end DOY + LL evaluation ---------------------------------------------------

### evaluate spatial + LL model -----------------------------------------------
# trained with raw data
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("spat_rf_noSubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original data 
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "spat_ll_rf", 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev)
  
  ## test against spatially subsampled data
  aucs <- lapply(fits, FUN = function(x, sp_name, test_data) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_data) {
      # get only subsampled test data that is in the spatial test blocks
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" are in different CV folds b/c cv block 
      # boundary is within the hectad and checklists are at finer spatial 
      # resolution)
      test_data <- test_data[test_data$test_fold == T, ] 
      tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name, test_data = test_data)
  }, sp_name = sp_to_fit[[i]], 
  test_data = test_points_ss[[i]]$test_points_subsampled)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "spat_ll_rf", 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev, fits)
}

# trained with spatially subsampled data
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("spat_rf_SubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original data 
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "spat_ll_rf", 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev)
  
  ## test against spatially subsampled data
  aucs <- lapply(fits, FUN = function(x, sp_name, test_data) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_data) {
      # get only subsampled test data that is in the spatial test blocks
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" are in different CV folds b/c cv block 
      # boundary is within the hectad and checklists are at finer spatial 
      # resolution)
      test_data <- test_data[test_data$test_fold == T, ] 
      tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name, test_data = test_data)
  }, sp_name = sp_to_fit[[i]], 
  test_data = test_points_ss[[i]]$test_points_subsampled)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "spat_ll_rf", 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev, fits)
}
### end spatial + LL evaluation -----------------------------------------------

### evaluate env + LL model --------------------------------------------------
# trained with raw data
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("env_rf_noSubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original data 
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "env_ll_rf", 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev)
  
  ## test against spatially subsampled data
  aucs <- lapply(fits, FUN = function(x, sp_name, test_data) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_data) {
      # get only subsampled test data that is in the spatial test blocks
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" are in different CV folds b/c cv block 
      # boundary is within the hectad and checklists are at finer spatial 
      # resolution)
      test_data <- test_data[test_data$test_fold == T, ] 
      tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name, test_data = test_data)
  }, sp_name = sp_to_fit[[i]], 
  test_data = test_points_ss[[i]]$test_points_subsampled)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "env_ll_rf", 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev, fits)
}

# trained with spatially subsampled data
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0("env_rf_SubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original data 
  aucs <- lapply(fits, FUN = function(x, sp_name) {
    auc_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "env_ll_rf", 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev)
  
  ## test against spatially subsampled data
  aucs <- lapply(fits, FUN = function(x, sp_name, test_data) {
    auc_iter <- lapply(x, FUN = function(y, sp_name, test_data) {
      # get only subsampled test data that is in the spatial test blocks
      test_data <- test_data[test_data$hectad %in% y$test_sites |
                               test_data$hectad %in% 
                               y$preds$hectad[y$preds$test_fold == T], ]
      test_data <- left_join(test_data, y$preds) # join predictions to test data
      # drop checklists not in the test fold (sometimes checlists in hectads 
      # that are in "y$test_sites" are in different CV folds b/c cv block 
      # boundary is within the hectad and checklists are at finer spatial 
      # resolution)
      test_data <- test_data[test_data$test_fold == T, ] 
      tryCatch({
        as.numeric(pROC::roc(response = test_data[, colnames(test_data) == 
                                                    gsub(" ", ".", sp_name)], 
                             predictor = test_data$pred)$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name, test_data = test_data)
  }, sp_name = sp_to_fit[[i]], 
  test_data = test_points_ss[[i]]$test_points_subsampled)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = "env_ll_rf", 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs))
  evals <- bind_rows(evals, ev)
  rm(ev, fits)
}
### end env + LL evaluation ---------------------------------------------------




# # get AUC for each fold & each model
# rf_performance <- expand.grid(fold = 1:n_folds, 
#                               species = names(sp_to_fit),
#                               model = c("doy_ll", "spatial", "environmental"),
#                               stringsAsFactors = FALSE)
# rf_performance$metric <- "AUC"
# rf_performance$value <- NA
# rf_performance$value[rf_performance$model == "doy_ll"] <- sapply(
#   day_ll_rf_fits, FUN = function(x) {sapply(x, FUN = function(y) {
#     tryCatch(y$auc, error = function(x) NA)})})
# rf_performance$value[rf_performance$model == "spatial"] <- sapply(
#   spatial_rf_fits, FUN = function(x) {sapply(x, FUN = function(y) {
#     tryCatch(y$auc, error = function(x) NA)})})
# rf_performance$value[rf_performance$model == "environmental"] <- sapply(
#   env_rf_fits, FUN = function(x) {sapply(x, FUN = function(y) {
#     tryCatch(y$auc, error = function(x) NA)})})

