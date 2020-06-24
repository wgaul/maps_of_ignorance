#############################
## Plots for millipedes SDMs - Supplementary materials
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 10 June 2020
## last modified: 16 June 2020
##############################
t_size <- 12

library(GGally)

## predictor variable correlations ------------------------------------------
predictors_df <- data.frame(newdata[newdata$day_of_year == 20, ])
predictors_df <- predictors_df[, colnames(predictors_df) %in% 
                                c("hectad", "mean_tn", "mean_tx", 
                                  "mean_rr", "artificial_surfaces", 
                                  "forest_seminatural_l1", 
                                  "wetlands_l1", "pasture_l2", 
                                  "arable_l2", "elev", 
                                  "eastings", "northings", 
                                  "sin_doy", "cos_doy")]
ggpairs(data = predictors_df, columns = 2:ncol(predictors_df), 
        title = "Predictor variable values\nfor all hectads")

ggpairs(data = mill, columns = which(colnames(mill) %in% 
                                       c("mean_tn", 
                                         "mean_rr", "artificial_surfaces", 
                                         "forest_seminatural_l1", 
                                         "wetlands_l1", "pasture_l2", 
                                         "arable_l2", "elev", 
                                         "eastings", "northings", 
                                         "sin_doy", "cos_doy", "day_of_year",
                                         "list_length")), 
        title = "Predictor variable values\non millipede checklists")
## end predictor variable correlatins --------------------------------------

## checklist length plots -----------------------------------------------------
hist(mill_wide$list_length, breaks = 15)
# calculated median and mean list length in each hectad
ll_df <- group_by(data.frame(mill_wide), hectad) %>%
  summarise(median_ll = median(list_length),
            mean_ll = mean(list_length), 
            n_lists = n()) %>% 
  left_join(., hec_names)

ggplot(data = ll_df, aes(x = eastings, y = northings, fill = median_ll)) + 
  geom_tile() + 
  scale_fill_continuous(name = "Median\nchecklist\nlength") + 
  theme_bw()

# ggplot(data = ll_df, aes(x = eastings, y = northings, fill = n_lists)) + 
#   geom_tile() + 
#   scale_fill_continuous(name = "Number\nof\nchecklists") + 
#   theme_bw()
## end checklist length plots -------------------------------------------------


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


# class balance - proportion of observations that are detections in training data
ggplot(data = evals[!is.na(evals$proportion_detections), ], 
       aes(x = factor(train_data), y = proportion_detections)) + 
  geom_boxplot() + 
  facet_wrap(~factor(species))

# predicton performance as a function of class balance
ggplot(data = evals[!is.na(evals$proportion_detections) & 
                      evals$metric == "AUC", ], 
       aes(x = proportion_detections, y = value, color = factor(train_data))) + 
  geom_point() + geom_smooth() + 
  ggtitle("AUC") + 
  facet_wrap(~factor(species))

ggplot(data = evals[!is.na(evals$proportion_detections) & 
                      evals$metric == "Brier", ], 
       aes(x = proportion_detections, y = value, color = factor(train_data))) + 
  geom_point() + geom_smooth() + 
  ggtitle("Brier") + 
  facet_wrap(~factor(species))

ggplot(data = evals[!is.na(evals$proportion_detections) & 
                      evals$metric == "Kappa", ], 
       aes(x = proportion_detections, y = value, color = factor(train_data))) + 
  geom_point() + geom_smooth() + 
  ggtitle("Kappa") + 
  facet_wrap(~factor(species))

ggplot(data = evals[!is.na(evals$proportion_detections) & 
                      evals$metric == "sensitivity", ], 
       aes(x = proportion_detections, y = value, color = factor(train_data))) + 
  geom_point() + geom_smooth() + 
  ggtitle("Sensitivity") + 
  facet_wrap(~factor(species))

ggplot(data = evals[!is.na(evals$proportion_detections) & 
                      evals$metric == "specificity", ], 
       aes(x = proportion_detections, y = value, color = factor(train_data))) + 
  geom_point() + geom_smooth() + 
  ggtitle("Specificity") + 
  facet_wrap(~factor(species))


## look at effect of block CV range size
ggplot(data = evals[evals$metric == "AUC" & !is.na(evals$block_cv_range), ], 
       aes(x = factor(block_cv_range), 
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
  facet_wrap(~species + factor(train_data, 
                               levels = c("raw", "spat_subsamp"), 
                               labels = c("raw", "spatially\nundersampled"))) + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8)




### plot variable importance ---------------------------------------------------
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


# graph variable importance by CV strategy
# CV strategy does not change variable importance conlcusions
vimp_plots <- lapply(vimp, FUN = function(x) {
  ggplot(data = x, 
         aes(x = variable, y = MeanDecreaseGini)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    ggtitle(paste0(x$species[1], "\n", x$model[1], ", trained with: ", 
                   x$train_dat[1], "\n")) + 
    facet_wrap(~cv)
})

for(i in 1:length(vimp_plots)) print(vimp_plots[i])
### end plot variable importance -----------------------------------------------


