#############################
## Evaluate DOY + LL SDMs for millipedes
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
## 
## TODO:  - test on original data
##        - test on spatially sub-sampled data
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 12 May 2020
## last modified: 25 May 2020
##############################
library(pROC)
library(psych)

## Get spatially subsampled test points ---------------------------------------
test_points_ss <- sapply(sp_to_fit, FUN = function(x) NULL)
names(test_points_ss) <- gsub(" ", ".", names(test_points_ss))
mill_wide_df <- data.frame(mill_wide)
for(i in 1:length(test_points_ss)) {
  mill_wide_df[, colnames(mill_wide_df) == names(test_points_ss)[i]] <- pa(
    mill_wide_df[, colnames(mill_wide_df) == names(test_points_ss)[i]])
  
  test_points_ss[[i]] <- list()
  # make a large number of different spatially subsampled test points datasets 
  # to choose from so each test (i.e. each model in a CV fold) is tested with a 
  # different set.  Use the subsample block designations created earlier.
  for(j in 1:n_subsamp_block_draws) {
    # add a column of subsampling block assignments by chosing randomly from
    # the many allocatins created in "prepare_objects_for_SDM.R"
    mill_wide_df$spat_subsamp_cell <- block_subsamp_10k[, sample(
      2:(ncol(block_subsamp_10k)-1), size = 1)]
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
    test_points_ss[[i]][[j]] <- bind_rows(absences, presences) 
  }
}
## end get spatially subsampled test points ------------------------------------


## evaluate models ------------------------------------------------------------
# trained with raw data
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0(mod_name, "_noSubSamp_fits_", 
                         gsub(" ", "_", sp_to_fit[[i]]), ".rds"))
  sp_name <- names(sp_to_fit)[i]
  
  ## test against original (raw) data 
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
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end AUC --------------------------------------
  
  # Cohen's Kappa --------------------------------------
  kappa_calc <- function(x, resp, pred) {
    ## Function to calculate Cohen's Kappa using the best threshold
    # make a 2 column df with observed responses (0 or 1) and a T/F whether the 
    # predicted value is above threshold x.  These will be used to compute
    # Cohen's kappa for agreement between categorical response (0, 1) and 
    # categorical predictor (0, 1)
    # ARGS: x - cutoff to test
    #       resp - vector of observed responses (0 or 1)
    #       pred - numeric vector of continuous predictions
    vals <- data.frame(resp = factor(resp[!is.na(resp)]), 
                       pred = factor(as.numeric(pred[!is.na(pred)] > x)))
    psych::cohen.kappa(vals)$kappa
  }
  
  kappas <- seq(from = 0, to = 1, by = 0.1) # Thresholds to try
  names(kappas) <- kappas
  
  # get kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas) {
      tryCatch({
        k_res <- sapply(kappas, kappa_calc, 
                        resp = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == sp_name], 
                        pred = y$preds$pred[y$preds$test_fold == T])
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  
  # sensitivity --------------------------------------------------
  sens <- lapply(fits, FUN = function(x, sp_name) {
    sens_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        
        
        as.numeric(pROC::roc(
          response = y$preds[y$preds$test_fold == T, 
                             colnames(y$preds) == gsub(" ", ".", sp_name)], 
          predictor = y$preds$pred[y$preds$test_fold == T])$auc)}, 
        error = function(x) NA)
    }, sp_name = sp_name)
  }, sp_name = sp_to_fit[[i]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "sensitivity", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  
  
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
  test_data = test_points_ss[[i]][[sample(1:length(test_points_ss[[i]]), 
                                          size = 1)]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end AUC ----------------------------------------
  
  # Kappa ------------------------------------
  # get kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas, test_data) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas, test_data) {
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
        k_res <- sapply(kappas, kappa_calc, 
                        resp = test_data[, colnames(test_data) == sp_name], 
                        pred = test_data$pred)
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas, test_data = test_data)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas, 
  test_data = test_points_ss[[i]][[sample(1:length(test_points_ss[[i]]),                                           size = 1)]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  rm(fits)
}

# trained with spatial subsampling
for(i in 1:length(sp_to_fit)) {
  fits <- readRDS(paste0(mod_name, "_SubSamp_fits_", 
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
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end AUC -----------------------------------------
  
  # Kappa ----------------------------------------------
  # get Kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas) {
      tryCatch({
        k_res <- sapply(kappas, kappa_calc, 
                        resp = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == sp_name], 
                        pred = y$preds$pred[y$preds$test_fold == T])
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  
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
  test_data = test_points_ss[[i]][[sample(1:length(test_points_ss[[i]]),                                           size = 1)]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) 
  
  # Kappa ------------------------------------
  # get kappa for each fold
  kp <- lapply(fits, function(x, sp_name, kappas, test_data) {
    kp_iter <- lapply(x, FUN = function(y, sp_name, kappas, test_data) {
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
        k_res <- sapply(kappas, kappa_calc, 
                        resp = test_data[, colnames(test_data) == sp_name], 
                        pred = test_data$pred)
        # return kappa and threshold that maximised kappa
        k_res[k_res == max(k_res)][1]}, 
        error = function(x) NA)
    }, sp_name = sp_name, kappas = kappas, test_data = test_data)
  }, sp_name = gsub(" ", ".", sp_to_fit[[i]]), kappas = kappas, 
  test_data = test_points_ss[[i]][[sample(1:length(test_points_ss[[i]]),                                           size = 1)]])
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block", metric = "Kappa", 
                   value = unlist(kp), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end Kappa ----------------------------------
  rm(fits)
}
### end evaluation ---------------------------------------------------