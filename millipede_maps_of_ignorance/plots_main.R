#############################
## Plots for millipedes SDMs
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
## 
## TODO:  - test on original data
##        - test on spatially sub-sampled data (and use the
##          same spatially sub-sampled dataset to test all models)
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 13 May 2020
## last modified: 27 May 2020
##############################
t_size <- 12

## plot results using random CV -----------------------------------------------
# how many detections and non-detections in test folds?
ndet_df <- pivot_longer(evals[evals$test_data == "spat_subsamp", ], 
                        n_dets_in_test:n_nondets_in_test, 
                        names_to = "detection", values_to = "count")
ggplot(data = ndet_df, 
       aes (x = factor(as.character(block_cv_range)), y = count, 
            color = detection)) + 
  geom_boxplot() + 
  ylab("Number in test dataset") + 
  ggtitle("When testing on spatially subsampled data") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))


# plot AUC for random CV
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                                     "\nEnvironment +\nDOY +\nList Length\n", 
                                     "\nLat + Lon +\nDOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("AUC\n(block Cross-Validated)") + 
  ggtitle("Random CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

ggplot(data = evals[evals$metric == "Kappa" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year (DOY) +\nList Length\n", 
                                     "\nEnvironment + DOY +\nList Length\n", 
                                     "\nLat + Lon + DOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Kappa\n(block Cross-Validated)") + 
  ggtitle("Random CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

# plot sensitivity
ggplot(data = evals[evals$metric == "sensitivity" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year (DOY) +\nList Length\n", 
                                     "\nEnvironment + DOY +\nList Length\n", 
                                     "\nLat + Lon + DOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Sensitivity\n(block Cross-Validated)") + 
  ggtitle("Random CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))


# plot specificity
ggplot(data = evals[evals$metric == "specificity" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year (DOY) +\nList Length\n", 
                                     "\nEnvironment + DOY +\nList Length\n", 
                                     "\nLat + Lon + DOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Specificity\n(block Cross-Validated)") + 
  ggtitle("Random CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end random CV -------------------------------------------------------------


### plot results using 100km block CV -----------------------------------------
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                                     "\nEnvironment +\nDOY +\nList Length\n", 
                                     "\nLat + Lon +\nDOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("AUC\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

ggplot(data = evals[evals$metric == "Kappa" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year (DOY) +\nList Length\n", 
                                     "\nEnvironment + DOY +\nList Length\n", 
                                     "\nLat + Lon + DOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Kappa\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

# plot sensitivity
ggplot(data = evals[evals$metric == "sensitivity" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year (DOY) +\nList Length\n", 
                                     "\nEnvironment + DOY +\nList Length\n", 
                                     "\nLat + Lon + DOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Sensitivity\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))


# plot specificity
ggplot(data = evals[evals$metric == "specificity" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year (DOY) +\nList Length\n", 
                                     "\nEnvironment + DOY +\nList Length\n", 
                                     "\nLat + Lon + DOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Specificity\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end 100km block CV ---------------------------------------------------------


















# look at effect of block CV range size
ggplot(data = evals[evals$metric == "AUC" & !is.na(evals$block_cv_range), ], 
       aes(x = factor(block_cv_range), 
           y = value, 
           color = factor(model, 
                          levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf"), 
                          labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                                     "\nEnvironment +\nDOY +\nList Length\n", 
                                     "\nLat + Lon +\nDOY +\nList Length\n")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(train_data, 
                               levels = c("raw", "spat_subsamp"), 
                               labels = c("raw", "spatially\nundersampled"))) + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8)






















# plot average of predictions from all 5 folds (so 4 predictions will be to
# training data, one prediction to test data in each grid cell)

# spatial model - raw data
mill_predictions_spat_ll_rf_raw <- list(
  `Ommatoiulus sabulosus` = readRDS("mill_predictions_spat_ll_rf_noSubSamp_Ommatoiulus_sabulosus.rds"), 
  `Tachypodoiulus niger` = readRDS("mill_predictions_spat_ll_rf_noSubSamp_Tachypodoiulus_niger.rds"))
mill_predictions_spat_ll_rf_raw <- lapply(
  mill_predictions_spat_ll_rf_raw, FUN= function(x) {
    group_by(x, hectad) %>%
      summarise(mean_prediction = mean(pred, na.rm = T), 
                eastings = mean(eastings), northings = mean(northings))
  })

for(i in 1:length(mill_predictions_spat_ll_rf_raw)) {
  print(ggplot() + 
          geom_tile(data = mill_predictions_spat_ll_rf_raw[[i]], 
                    aes(x = eastings, y = northings, fill = mean_prediction)) + 
          geom_point(data = mill, aes(x = eastings, y = northings),
                     color = "light grey", size = 0.5) +
          geom_point(data = mill[mill$Genus_species ==
                                   names(mill_predictions_spat_ll_rf_raw)[i], ],
                     aes(x = eastings, y = northings), color = "orange") +
          ggtitle(paste0(names(mill_predictions_spat_ll_rf_raw)[i], 
                         " - spatial model"))
  )
}


# environmental variables model - raw data
mill_predictions_env_rf_raw <- list(
  `Ommatoiulus sabulosus` = readRDS("mill_predictions_env_rf_noSubSamp_Ommatoiulus_sabulosus.rds"), 
  `Tachypodoiulus niger` = readRDS("mill_predictions_env_rf_noSubSamp_Tachypodoiulus_niger.rds"))
mill_predictions_env_rf_raw <- lapply(
  mill_predictions_env_rf_raw, FUN= function(x) {
    group_by(x, hectad) %>%
      summarise(mean_prediction = mean(pred, na.rm = T), 
                eastings = mean(eastings), northings = mean(northings))
  })

for(i in 1:length(mill_predictions_env_rf_raw)) {
  print(ggplot() + 
          geom_tile(data = mill_predictions_env_rf_raw[[i]], 
                    aes(x = eastings, y = northings, fill = mean_prediction)) + 
          geom_point(data = mill, aes(x = eastings, y = northings),
                     color = "light grey", size = 0.5) +
          geom_point(data = mill[mill$Genus_species ==
                                   names(mill_predictions_env_rf_raw)[i], ],
                     aes(x = eastings, y = northings), color = "orange") +
          ggtitle(paste0(names(mill_predictions_env_rf_raw)[i], 
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

