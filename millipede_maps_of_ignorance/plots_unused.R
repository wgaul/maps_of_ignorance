#############################
## Plots for millipedes SDMs - unused plots
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 11 June 2020
## last modified: 11 June 2020
##############################
t_size <- 12

warning("re-do evenness for test data to include zero hectads.")
mapply(FUN = function(x, nm) hist(x, main = nm, 
                                  xlab = "Simpson Evenness - Test data"), 
       evenness_test, names(evenness_test))

#######################
### plot results random CV all combinations ------------------------------------
# plot AUC for random CV
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "random" & 
                      as.character(evals$test_data) == "spat_subsamp", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = species)) + 
  geom_boxplot() + 
  facet_wrap(~factor(
    model, 
    levels = c("day_ll_rf", "spat_ll_rf","env_ll_rf", 
               "env_spat_ll_rf"), 
    labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
               "\nLat + Lon +\nDOY +\nList Length\n",
               "\nEnvironment +\nDOY +\nList Length\n", 
               "\nEnvironment + \nLat + Long +\nDOY +\nList Length"))) + 
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
### end plot results random CV all combinations -------------------------------
##################

#############################################
### plot results using 100km block CV -----------------------------------------
ggplot(data = evals[evals$metric == "AUC" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year\n(DOY) +\nList Length\n", 
                        "\nEnvironment +\nDOY +\nList Length\n", 
                        "\nLat + Lon +\nDOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
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
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
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
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
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
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
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

# plot specificity
ggplot(data = evals[evals$metric == "Brier" & 
                      as.character(evals$block_cv_range) == "1e+05", ], 
       aes(x = factor(train_data, 
                      levels = c("raw", "spat_subsamp"), 
                      labels =  c("raw", "spatially\nundersampled")), 
           y = value, 
           color = factor(
             model, 
             levels = c("day_ll_rf", "env_ll_rf", "spat_ll_rf", 
                        "env_spat_ll_rf"), 
             labels = c("\nDay of Year (DOY) +\nList Length\n", 
                        "\nEnvironment + DOY +\nList Length\n", 
                        "\nLat + Lon + DOY +\nList Length\n", 
                        "\nEnvironment + \nLat + Long +\nDOY +\nList Length")))) + 
  geom_boxplot() + 
  facet_wrap(~species + factor(
    test_data, 
    levels = c("raw", "spat_subsamp"), 
    labels = c("test data - raw", "test data -\nspatially undersampled"))) + 
  xlab("Training Data") + 
  ylab("Brier score\n(block Cross-Validated)") + 
  ggtitle("100km block CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end 100km block CV ---------------------------------------------------------


#######################################################
## plot variable importance for all models ------------------------------------
# read in variable importance results
vimp <- list.files("./saved_objects/")
vimp <- vimp[grepl("var_import.*", vimp)]
vimp <- lapply(vimp, function(x) readRDS(paste0("./saved_objects/", x)))
# average the variable importance from each CV fold
vimp <- lapply(vimp, FUN = function(x) {
  group_by(x, variable, cv) %>%
    summarise(MeanDecreaseGini = mean(MeanDecreaseGini), 
              species = unique(species), model = unique(model), 
              train_dat = unique(train_dat))
})

vimp_plots <- lapply(vimp, FUN = function(x) {
  dat <- x[x$cv == "random", ]
  dat <- dat[order(dat$MeanDecreaseGini, decreasing = FALSE), ]
  ggplot(data = dat, 
         aes(x = factor(variable, levels = dat$variable, 
                        labels = dat$variable, ordered = T), 
             y = MeanDecreaseGini)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    ggtitle(paste0(dat$species[1], "\n", dat$model[1], ", trained with: ", 
                   dat$train_dat[1], "\n", "CV: ", dat$cv[1]))
})

# plot variable importance for all models
for(i in 1:length(vimp_plots)) print(vimp_plots[i])
### end variable importance ----------------------------------------------------
##################################################


### Numbers and tables for paper -------------------------------------------

# Number of detections per species when using 10km resolution
n_detections_per_species_10km <- data.frame(table(mill$Genus_species))
n_detections_per_species_10km <- n_detections_per_species_10km[order(
  n_detections_per_species_10km$Freq, decreasing = FALSE), ]
colnames(n_detections_per_species_10km) <- c("species", "number_of_detections")
n_detections_per_species_10km <- n_detections_per_species_10km[
  n_detections_per_species_10km$species %in% sp_to_fit, ]
n_detections_per_species_10km