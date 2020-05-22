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
## last modified: 13 May 2020
##############################
t_size <- 15

ggplot(data = evals[evals$metric == "AUC", ], 
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
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

ggplot(data = evals[evals$metric == "Kappa", ], 
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
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))





# plot average of predictions from all 5 folds (so 4 predictions will be to
# training data, one prediction to test data in each grid cell)

# spatial model - raw data
mill_predictions_spat_rf_raw <- list(
  `Ommatoiulus sabulosus` = readRDS("mill_predictions_spat_rf_noSubSamp_Ommatoiulus_sabulosus.rds"), 
  `Tachypodoiulus niger` = readRDS("mill_predictions_spat_rf_noSubSamp_Tachypodoiulus_niger.rds"))
mill_predictions_spat_rf_raw <- lapply(
  mill_predictions_spat_rf_raw, FUN= function(x) {
    group_by(x, hectad) %>%
      summarise(mean_prediction = mean(pred, na.rm = T), 
                eastings = mean(eastings), northings = mean(northings))
  })

for(i in 1:length(mill_predictions_spat_rf_raw)) {
  print(ggplot() + 
          geom_tile(data = mill_predictions_spat_rf_raw[[i]], 
                    aes(x = eastings, y = northings, fill = mean_prediction)) + 
          geom_point(data = mill, aes(x = eastings, y = northings),
                     color = "light grey", size = 0.5) +
          geom_point(data = mill[mill$Genus_species ==
                                   names(mill_predictions_spat_rf_raw)[i], ],
                     aes(x = eastings, y = northings), color = "orange") +
          ggtitle(paste0(names(mill_predictions_spat_rf_raw)[i], 
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