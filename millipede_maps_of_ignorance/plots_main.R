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
## last modified: 16 June 2020
##############################
try(rm(block_subsamp, fold_assignments, hec_names_spat, mill_fewer_vars, 
       mill_spat))
t_size <- 12

evals <- read.csv("./saved_objects/evals.csv")

# spatial evenness of training and test datasets
ggplot(data = evals[!is.na(evals$simpson_training_hectad), ], 
       aes(x = factor(block_cv_range), y = simpson_training_hectad)) + 
  geom_boxplot() + 
  facet_wrap(~train_data) + 
  ggtitle("Spatial Evenness calculated at hectad scale\nincluding hectads with zero observations")

# spatial evenness of training and test datasets
ggplot(data = evals[!is.na(evals$simpson_training_subsampBlock), ], 
       aes(x = factor(block_cv_range), y = simpson_training_subsampBlock)) + 
  geom_boxplot() + 
  facet_wrap(~train_data) + 
  ggtitle("Spatial Evenness calculated at subsample block scale\nincluding hectads with zero observations")

warning("re-do evenness for test data to include zero hectads.")
mapply(FUN = function(x, nm) hist(x, main = nm, 
                                  xlab = "Simpson Evenness - Test data"), 
       evenness_test, names(evenness_test))


# plot AUC for random CV
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "spat_ll_rf","env_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                        "\nLat + Lon +\nDOY +\nList Length\n",
                        "\nEnvironment +\nDOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("AUC\n(Cross-Validated)") + 
  ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
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
           color = factor(
             model, 
             levels = c("day_ll_rf", "spat_ll_rf", "env_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n",  
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n",
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Kappa\n(Cross-Validated)") + 
  ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
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
           color = factor(
             model, 
             levels = c("day_ll_rf", "spat_ll_rf", "env_ll_rf",  
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n",  
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Sensitivity\n(Cross-Validated)") + 
  ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
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
           color = factor(
             model, 
             levels = c("day_ll_rf", "spat_ll_rf", "env_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n",
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Specificity\n(block Cross-Validated)") + 
  ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))

# plot Brier score
ggplot(data = evals[evals$metric == "Brier" & 
                      as.character(evals$block_cv_range) == "random", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "spat_ll_rf", "env_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Brier score\n(Cross-Validated)") + 
  ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end random CV -------------------------------------------------------------



### plot variable importance only for models of interest ----------------------
# read in variable importance results
vimp <- list.files("./saved_objects/")
vimp <- vimp[grepl("var_import.*", vimp)]
vimp <- vimp[grepl(paste0(".*", analysis_resolution, ".rds"), vimp)]
vimp <- vimp[grepl(".*env_spat_ll.*", vimp)]
vimp <- lapply(vimp, function(x) readRDS(paste0("./saved_objects/", x)))
# average the variable importance from each CV fold
vimp <- lapply(vimp, FUN = function(x) {
  group_by(x, variable, cv) %>%
    summarise(MeanDecreaseGini = mean(MeanDecreaseGini), 
              species = unique(species), model = unique(model), 
              train_dat = unique(train_dat))
})

vimp <- bind_rows(vimp)

vimp_plots <- lapply(sp_to_fit, FUN = function(x, v_df) {
  dat <- v_df[v_df$cv == "random" & v_df$species == x, ]
  dat <- dat[order(dat$MeanDecreaseGini, decreasing = FALSE), ]
  ggplot(data = dat, 
         aes(x = factor(variable, levels = dat$variable, 
                        labels = dat$variable, ordered = T), 
             y = MeanDecreaseGini)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    ggtitle(paste0(dat$species[1], "\n", dat$model[1], ", trained with: ", 
                   dat$train_dat[1], "\n", "CV: ", dat$cv[1], 
                   " analysis resolution: ", analysis_resolution)) + 
    facet_wrap(~factor(train_dat), scales = "free")
}, v_df = vimp)

for(i in 1:length(vimp_plots)) print(vimp_plots[i])
### end plot variable importance ---------------------------------------------


### plot partial dependence --------------------------------------------------
# read in partial dependence files
pd <- list.files("./saved_objects/")
pd <- pd[grepl("partial_depen.*", pd)]
pd <- pd[grepl(paste0(".*", analysis_resolution, ".rds"), pd)]
pd <- pd[grepl(paste0(".*", mods_for_pd_plots, ".*", collapse = "|"), pd)]
names(pd) <- pd
pd <- lapply(pd, function(x) readRDS(paste0("./saved_objects/", x)))
# average the dependence for each variable from each CV fold
pd <- lapply(pd, FUN = function(x) {
  group_by(x, x, variable, cv) %>%
    summarise(y = mean(y), 
              species = unique(species), model = unique(model), 
              train_data = unique(train_data))
})

# make plots only using the "random" CV, which is the one I will use for results
pd_plots <- lapply(pd, FUN = function(x) {
  dat <- x[x$cv == "random", ]
  ggplot(data = dat, 
         aes(x = x, y = y)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~variable, scales = "free_x") +
    ggtitle(paste0(dat$species[1], "\n", dat$model[1], ", trained with: ", 
                   dat$train_data[1], "\n", "CV: ", dat$cv[1], 
                   "\nanalysis resolution: ", analysis_resolution))
})

for(i in 1:length(pd_plots)) print(pd_plots[i])
rm(pd)
### end plot partial dependence -----------------------------------------------


### plot predictions with standardized list length ----------------------------
# load standardized predictions
stpred <- list.files("./saved_objects/")
stpred <- stpred[grepl("standard_pre.*", stpred)]
stpred <- stpred[grepl(paste0(".*", analysis_resolution, ".rds"), stpred)]
stpred <- stpred[grepl(paste0(".*", mods_for_pd_plots, ".*", collapse = "|"), 
                       stpred)]
names(stpred) <- gsub("standard.*tions_", "", stpred)
names(stpred) <- gsub(".rds", "", names(stpred))
# stpred <- lapply(stpred, 
#                  FUN = function(x) readRDS(paste0("./saved_objects/", x)))

# plot average of predictions from all 5 folds (so 4 predictions will be to
# training data, one prediction to test data in each grid cell)
prediction_plots <- mapply(FUN = function(x, nm, mill) {
  dat <- readRDS(paste0("./saved_objects/", x)) # load object
  # map only predictions from random CV
  dat <- dat[dat$cv == "random", ] # use only random CV results
  sp <- gsub(".*Samp_", "", nm) # get species name
  sp <- gsub("1000.*", "", sp)
  sp <- gsub("_", " ", sp)
  mod <- gsub("_SubSamp.*|_noSub.*", "", nm) # get model name
  train_data <- gsub("^.*_rf_", "", nm) # get training data type name
  train_data <- gsub("_.*_.*$", "", train_data)
  
  # get average predictions for each grid cell (averaging over all folds)
  dat <- group_by(dat, en) %>%
    summarise(mean_prediction = mean(mean_pred, na.rm = T), 
              eastings = mean(eastings), northings = mean(northings))
  ggplot() + 
    geom_tile(data = dat, 
              aes(x = eastings, y = northings, fill = mean_prediction)) + 
    geom_point(data = mill, aes(x = eastings, y = northings),
               color = "light grey", size = 0.5) +
    geom_point(data = mill[mill$Genus_species == sp, ],
               aes(x = eastings, y = northings), color = "orange") +
    ggtitle(paste0(sp, "\n", mod, " trained with ", train_data)) + 
    theme_bw()
}, stpred, names(stpred), MoreArgs = list(mill = mill), SIMPLIFY = FALSE)

for(i in 1:length(prediction_plots)) print(prediction_plots[i])
### end plot standardized predictions -----------------------------------------


### save figures and tables ---------------------------------------------------
# Number of repeat visits to grid cells (though there could be multiple 
# locations visited within a grid cell, so these are not necessarily true 
# repeat visits)
repeat_visits <- data.frame(
  en = paste0(mill_wide$eastings, "_", mill_wide$northings), 
  year = mill_wide$year)
table(repeat_visits$en)[order(table(repeat_visits$en), decreasing = T)]
table(as.numeric(table(repeat_visits$en, repeat_visits$year)))
# # look at most re-visited grid cell
# data.frame(mill_wide[mill_wide$eastings == 316499 & 
#                        mill_wide$northings == 237500, ]) 

# number of checklists
nrow(mill_wide)

# proportion of lists with list lenght of 1
length(which(mill_wide$list_length == 1)) / nrow(mill_wide)
# proportion of lists with list lenght of 1
length(which(mill_wide$list_length == 2)) / nrow(mill_wide)

# Number of detections per species when using 10km resolution
n_detections_per_species_10km <- data.frame(table(mill$Genus_species))
n_detections_per_species_10km <- n_detections_per_species_10km[order(
  n_detections_per_species_10km$Freq, decreasing = FALSE), ]
colnames(n_detections_per_species_10km) <- c("species", "number_of_detections")
n_detections_per_species_10km <- n_detections_per_species_10km[
  n_detections_per_species_10km$species %in% sp_to_fit, ]
n_detections_per_species_10km

# Number of detections per species when using 1 km resolution
n_detections_per_species_1km <- data.frame(
  table(mill$Genus_species[mill$Precision <= 1000]))
n_detections_per_species_1km <- n_detections_per_species_1km[order(
  n_detections_per_species_1km$Freq, decreasing = FALSE), ]
colnames(n_detections_per_species_1km) <- c("species", "number_of_detections")
n_detections_per_species_1km <- n_detections_per_species_1km[
  n_detections_per_species_1km$species %in% sp_to_fit, ]
n_detections_per_species_1km

### end save figures and tables -----------------------------------------------
