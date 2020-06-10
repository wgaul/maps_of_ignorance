#############################
## Plots for millipedes SDMs - Supplementary materials
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 10 June 2020
## last modified: 10 June 2020
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
