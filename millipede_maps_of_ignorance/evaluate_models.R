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
## last modified: 26 May 2020
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
    # the many allocations created in "prepare_objects_for_SDM.R"
    mill_wide_df$spat_subsamp_cell <- block_subsamp_10k[, sample(
      2:(ncol(block_subsamp_10k)-1), size = 1)]

    # spatially sub-sample absence checklists to 1 per cell
    # for testing, spatially subsample ALL data, so take only 1 checklist from 
    # each cell (i.e. do not undersample by keeping all detections and only 
    # spatially subsampling non-detections)
    # cell_abs_tab <- table(absences$spat_subsamp_cell)
    keep_rows <- c()
    for(ri in 1:length(unique(mill_wide_df$spat_subsamp_cell))) {
      cell <- unique(mill_wide_df$spat_subsamp_cell)[ri]
      keep_rows <- c(keep_rows, 
                     sample(which(mill_wide_df$spat_subsamp_cell == cell), 
                            size = 2))
    }
    # store this test datasets
    test_points_ss[[i]][[j]] <- mill_wide_df[keep_rows, ]
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
  
  ndets <- lapply(fits, FUN = function(x, sp_name) {
    det_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        # count number of detections in test fold
        length(which(y$preds[y$preds$test_fold == TRUE, colnames(y$preds) %in% 
                               gsub(" ", ".", sp_name)] > 0))}, 
        error = function(x) NA)
    } , sp_name = sp_name)}, sp_name = sp_to_fit[[i]])
  nnondets <- lapply(fits, FUN = function(x, sp_name) {
    det_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        # count number of detections in test fold
        length(which(y$preds[y$preds$test_fold == FALSE, colnames(y$preds) %in% 
                               gsub(" ", ".", sp_name)] > 0))}, 
        error = function(x) NA)
    } , sp_name = sp_name)}, sp_name = sp_to_fit[[i]])
  
  # find number of detections and non-detections in each test fold
  n_det_nondet <- data.frame(species = sp_name, model = mod_name,
                             block_range = unlist(block_ranges), 
                             n_det = unlist(ndets), 
                             n_non_det = unlist(nnondets))
  n_dets_df <- bind_rows(n_dets_df, n_det_nondet)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges, n_det_nondet) # end AUC ---------------------------------
  
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
  # use threshold that maximised Kappa for each model
  sens <- mapply(FUN = function(x, thresh, sp_name) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  # use threshold that maximised Kappa for each model
  specif <- mapply(FUN = function(x, thresh, sp_name) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        specificity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "raw", cv = "block", metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  
  
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
  
  # sensitivity --------------------------------------------------
  sens <- mapply(FUN = function(x, thresh, sp_name) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block", 
                   metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  specif <- mapply(FUN = function(x, thresh, sp_name) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        specificity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "raw", 
                   test_data = "spat_subsamp", cv = "block",
                   metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
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
  
  ndets <- lapply(fits, FUN = function(x, sp_name) {
    det_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        # count number of detections in test fold
        length(which(y$preds[y$preds$test_fold == TRUE, colnames(y$preds) %in% 
                               gsub(" ", ".", sp_name)] > 0))}, 
        error = function(x) NA)
    } , sp_name = sp_name)}, sp_name = sp_to_fit[[i]])
  nnondets <- lapply(fits, FUN = function(x, sp_name) {
    det_iter <- lapply(x, FUN = function(y, sp_name) {
      tryCatch({
        # count number of detections in test fold
        length(which(y$preds[y$preds$test_fold == FALSE, colnames(y$preds) %in% 
                               gsub(" ", ".", sp_name)] > 0))}, 
        error = function(x) NA)
    } , sp_name = sp_name)}, sp_name = sp_to_fit[[i]])
  
  # find number of detections and non-detections in each test fold
  n_det_nondet <- data.frame(species = sp_name, model = mod_name,
                             block_range = unlist(block_ranges), 
                             n_det = unlist(ndets), 
                             n_non_det = unlist(nnondets))
  n_dets_df <- bind_rows(n_dets_df, n_det_nondet)
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", metric = "AUC", 
                   value = unlist(aucs), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges, n_det_nondet) # end AUC ---------------------------------
  
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
  
  # sensitivity --------------------------------------------------
  sens <- mapply(FUN = function(x, thresh, sp_name) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block", 
                   metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  specif <- mapply(FUN = function(x, thresh, sp_name) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        specificity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "raw", cv = "block",
                   metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  
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
  
  # sensitivity --------------------------------------------------
  sens <- mapply(FUN = function(x, thresh, sp_name) {
    sens_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        sensitivity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block", 
                   metric = "sensitivity", 
                   value = unlist(sens), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end sensitivity --------------------------------------
  
  # specificity --------------------------------------------------
  specif <- mapply(FUN = function(x, thresh, sp_name) {
    specif_iter <- mapply(FUN = function(y, thresh, sp_name) {
      tryCatch({
        specificity(threshold = thresh, 
                    response = y$preds[y$preds$test_fold == T, 
                                       colnames(y$preds) == gsub(" ", ".", 
                                                                 sp_name)], 
                    predictions = y$preds$pred[y$preds$test_fold == T])}, 
        error = function(x) NA)
    }, x, thresh, MoreArgs = list(sp_name = sp_name))
  }, fits, kp, MoreArgs = list(sp_name = sp_to_fit[[i]]), SIMPLIFY = F)
  
  block_ranges <- lapply(fits, FUN = function(x) {
    range_iter <- lapply(x, FUN = function(y) {
      tryCatch({y$block_cv_range}, error = function(x) NA)
    })})
  
  # put evaluation metrics for every fold into df
  ev <- data.frame(species = sp_name, model = mod_name, 
                   train_data = "spat_subsamp", 
                   test_data = "spat_subsamp", cv = "block",
                   metric = "specificity", 
                   value = unlist(specif), 
                   block_cv_range = unlist(block_ranges))
  evals <- bind_rows(evals, ev)
  rm(ev, block_ranges) # end specificity --------------------------------------
  rm(fits)
}
### end evaluation ---------------------------------------------------