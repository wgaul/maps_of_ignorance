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
## last modified: 12 June 2020
##############################
t_size <- 12

# spatial evenness of training and test datasets
ggplot(data = evals[!is.na(evals$simpson_training), ], 
       aes(x = factor(block_cv_range), y = simpson_training)) + 
  geom_boxplot() + 
  facet_wrap(~train_data) + 
  ggtitle("Spatial Evenness calculated at hectad scale\nincluding hectads with zero observations")

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
  ylab("AUC\n(Cross-Validated)") + 
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
  ylab("Kappa\n(Cross-Validated)") + 
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
  ylab("Sensitivity\n(Cross-Validated)") + 
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
  ggtitle("Random CV") + 
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
  ylab("Brier score\n(Cross-Validated)") + 
  ggtitle("Random CV") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size))
### end random CV -------------------------------------------------------------



### plot variable importance ---------------------------------------------------
# read in variable importance results
vimp <- list.files("./saved_objects/")
vimp <- vimp[grepl("var_import.*", vimp)]
vimp <- lapply(vimp, readRDS)
# average the variable importance from each CV fold
vimp <- lapply(vimp, FUN = function(x) {
  group_by(x, variable, cv) %>%
    summarise(MeanDecreaseGini = mean(MeanDecreaseGini), 
              species = unique(species), model = unique(model), 
              train_dat = unique(train_dat))
})

vimp_plots <- lapply(vimp, FUN = function(x) {
  dat <- x[x$cv == "random", ]
  ggplot(data = dat, 
         aes(x = variable, y = MeanDecreaseGini)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    ggtitle(paste0(dat$species[1], "\n", dat$model[1], ", trained with: ", 
                   dat$train_dat[1], "\n", "CV: ", dat$cv[1]))
})

for(i in 1:length(vimp_plots)) print(vimp_plots[i])
### end plot variable importance -----------------------------------------------


### plot partial dependence ---------------------------------------------------
# read in partial dependence files
pd <- list.files("./saved_objects/")
pd <- pd[grepl("partial_depen.*", pd)]
pd <- lapply(pd, readRDS)
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
                   dat$train_data[1], "\n", "CV: ", dat$cv[1]))
})

for(i in 1:length(pd_plots)) print(pd_plots[i])

# see effect of DOY by looking at standardized predictions
# load standardized predictions
stpred <- list.files("./saved_objects/")
stpred <- stpred[grepl("standard_pre.*", stpred)]
names(stpred) <- gsub("standard.*tions_", "", stpred)
names(stpred) <- gsub(".rds", "", names(stpred))
stpred <- lapply(stpred, 
                 FUN = function(x) readRDS(paste0("./saved_objects/", x)))
# average predictions for each day of year (taking average over all locations)
doy_plots <- lapply(stpred, FUN = function(x) {
  # keep only random CV predictions
  dat <- x[x$cv == "random", ]
  dat <- group_by(dat, day_of_year) %>%
    summarise(pred = mean(pred, na.rm = T), cos_doy = mean(cos_doy), 
                          sin_doy = mean(sin_doy))
  # # make sure 1st and last days of year are close to each other
  # ggplot(data = dat, aes(x = cos_doy, y = sin_doy, color = day_of_year, 
  #                        size = pred)) + geom_point()
  
  ggplot(data = dat, aes(x = day_of_year, y = pred)) + 
    geom_point() + geom_line() + 
    theme_bw()
})

doy_plots <- mapply(FUN = function(x, nm) {x + ggtitle(nm)}, 
                    doy_plots, names(doy_plots), SIMPLIFY = FALSE)

for(i in 1:length(doy_plots)) print(doy_plots[i])
### end plot partial dependence ------------------------------------------------


### plot predictions with standardized list length ----------------------------
# plot average of predictions from all 5 folds (so 4 predictions will be to
# training data, one prediction to test data in each grid cell)
prediction_plots <- mapply(FUN = function(x, nm, mill) {
  # map only predictions from random CV
  dat <- x[x$cv == "random", ]
  sp <- gsub(".*Samp_", "", nm)
  sp <- gsub("_", " ", sp)
  mod <- gsub("_SubSamp.*|_noSub.*", "", nm)
  train_data <- gsub("^.*_rf_", "", nm)
  train_data <- gsub("_.*_.*$", "", train_data)
  
  dat <- group_by(dat, hectad) %>%
    summarise(mean_prediction = mean(pred, na.rm = T), 
              eastings = mean(eastings), northings = mean(northings))
  ggplot() + 
    geom_tile(data = dat, 
              aes(x = eastings, y = northings, fill = mean_prediction)) + 
    geom_point(data = mill, aes(x = eastings, y = northings),
               color = "light grey", size = 0.5) +
    geom_point(data = mill[mill$Genus_species == sp, ],
               aes(x = eastings, y = northings), color = "orange") +
    ggtitle(paste0(sp, "\n", mod, " trained with ", train_data))
}, stpred, names(stpred), MoreArgs = list(mill = mill), SIMPLIFY = FALSE)

for(i in 1:length(prediction_plots)) print(prediction_plots[i])
### end plot standardized predictions -----------------------------------------


