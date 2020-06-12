#############################
## Plots for millipedes SDMs - Supplementary materials
## 
## This script is meant to be sourced from within 'main_millipede_maps.R'.
##
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 10 June 2020
## last modified: 12 June 2020
##############################
t_size <- 12

library(GGally)

## predictor variable correlations ------------------------------------------
predictors_df <- data.frame(newdata[newdata$day_of_year == 10, ])
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
                                       c("mean_tn", "mean_tx", 
                                         "mean_rr", "artificial_surfaces", 
                                         "forest_seminatural_l1", 
                                         "wetlands_l1", "pasture_l2", 
                                         "arable_l2", "elev", 
                                         "eastings", "northings", 
                                         "sin_doy", "cos_doy", 
                                         "list_length")), 
        title = "Predictor variable values\non millipede checklists")
## end predictor variable correlatins --------------------------------------


## calculate Simpson's evenness and sample size for an example training
## dataset for adding a point to the contour plot (from simulation) -----------

# make example spatial subsampling blocks 
checklists_spat <- mill_wide[ , "checklist_ID"] # make df of checklists
b_subsamp <- spatialBlock(checklists_spat, 
                          theRange = block_range_spat_undersamp,
                          k = n_folds, 
                          selection = "random", 
                          iteration = 5, 
                          # rows = 20, cols = 20, 
                          xOffset = runif(n = 1, min = 0, max = 1), 
                          yOffset = runif(n = 1, min = 0, max = 1),
                          showBlocks = FALSE, 
                          rasterLayer = pred_brick$pasture_l2, 
                          biomod2Format = FALSE)
# add spatial subsampling grid cell ID to each hectad
example_blocks <- st_join(
  st_as_sf(newdata[newdata$day_of_year == newdata$day_of_year[1], ]), 
  st_as_sf(b_subsamp$blocks[, names(
    b_subsamp$blocks) == "layer"]))
example_blocks <- data.frame(example_blocks)
colnames(example_blocks)[which(
  colnames(example_blocks) == "layer")] <- "subsamp_blocks"
# remove geometry column
example_blocks <- example_blocks[, -which(grepl(".*geomet.*", 
                                                colnames(example_blocks)))]

# subsample mill data to get example training checklists
sp_df <- data.frame(mill_wide)
sp_name <- "Ommatoiulus.sabulosus"
sp_df$spat_subsamp_cell <- block_subsamp_10k[, sample(
  2:(ncol(block_subsamp_10k)-1), size = 1)]
# spatially sub-sample absence checklists to 1 per cell
# separate presence and absence checklists.  Keep all presence checklists.
presences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 1, ]
absences <- sp_df[sp_df[colnames(sp_df) == sp_name] == 0, ]
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
sp_df <- bind_rows(absences, presences)


table_nobs <- data.frame(table(sp_df$hectad))
colnames(table_nobs) <- c("hectad", "nrec")
# add number of recs
example_blocks <- left_join(example_blocks, table_nobs, by = "hectad")
example_blocks$nrec[is.na(example_blocks$nrec)] <- 0
# sum number of recs per subsample block
example_blocks <- group_by(example_blocks, subsamp_blocks) %>%
  summarise(nrec = sum(nrec))

simpson_even(example_blocks$nrec) # print Simpson's evenness

rm(sp_df, sp_name)
## end calculate evenness and sample size from example training dataset ------



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
