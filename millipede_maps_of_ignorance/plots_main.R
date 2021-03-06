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
## last modified: 8 Dec 2020
##############################
library(patchwork)
try(rm(block_subsamp, fold_assignments, hec_names_spat, mill_fewer_vars, 
       mill_spat))
t_size <- 20

evals <- list.files("./saved_objects/")
evals <- evals[grepl("evals.*", evals)]
evals <- bind_rows(lapply(evals, 
                          FUN = function(x) read.csv(paste0("./saved_objects/", 
                                                            x))))

# order species factor by number of detections (rarest species first)
evals$species <- factor(evals$species, 
                        levels = c("Macrosternodesmus palicola", 
                                   "Boreoiulus tenuis", 
                                   "Ommatoiulus sabulosus", 
                                   "Blaniulus guttulatus", 
                                   "Glomeris marginata", 
                                   "Cylindroiulus punctatus"))

# Load Ireland coastline
ir <- readOGR(dsn='./data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))

# make annotations for N arrow and scale bar
annot <- data.frame(x1 = c(265000, 310000, 60000, 60000), 
                    x2 = c(365000, 310000, 60000, 60000), 
                    y1 = c(60000, 40000, 400000, 380000), 
                    y2 = c(60000, 40000, 455000, 380000),
                    label = c(NA, "100 km", NA, "N"), 
                    bias = "D")


### Map example spatially under-sampled data ----------------------------------
## get an example spatially under-sampled dataset
block_subsamp <- readRDS("./block_subsamp.rds") # load subsampling blocks
ex_osab <- data.frame(mill_wide)
# add a column of subsampling block assignments by chosing randomly from
# the many allocations created in "prepare_objects_for_SDM.R"
ex_osab$spat_subsamp_cell <- block_subsamp[, sample(2:(ncol(block_subsamp)-1), 
                                                  size = 1)]
# spatially sub-sample absence checklists to 1 per cell
# separate presence and absence checklists.  Keep all presence checklists.
presences <- ex_osab[ex_osab[, colnames(ex_osab) == 
                                     "Ommatoiulus.sabulosus"] == 1, ]
absences <- ex_osab[ex_osab[, colnames(ex_osab) == 
                                     "Ommatoiulus.sabulosus"] == 0, ]
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
ex_osab <- bind_rows(absences, presences)
ex_osab <- SpatialPointsDataFrame(
  coords = ex_osab[, c("eastings", "northings")], 
  data = ex_osab, proj4string = CRS("+init=epsg:29903"))

map_Osab_raw <- ggplot() + 
  geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
  geom_sf(data = st_as_sf(mill_wide[mill_wide$`Ommatoiulus sabulosus` == 0, ]), 
          color = "dark grey", size = 0.04*t_size) + 
  geom_sf(data = st_as_sf(mill_wide[mill_wide$`Ommatoiulus sabulosus` == 1, ]), 
          color = "dark orange", size = 0.04*t_size) + 
  geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label)) + 
  geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
               arrow = arrow(length = unit(0.1, "npc"))) + 
  ylab("Latitude") + xlab("Longitude") + 
  ggtitle("(a)") +
  theme_bw() + 
  theme(text = element_text(size = 0.7*t_size), 
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

map_Osab_spat_subsamp <- ggplot() + 
  geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
  geom_sf(data = st_as_sf(ex_osab[ex_osab$`Ommatoiulus.sabulosus` == 0, ]), 
          color = "dark grey", size = 0.04*t_size) + 
  geom_sf(data = st_as_sf(ex_osab[ex_osab$`Ommatoiulus.sabulosus` == 1, ]), 
          color = "dark orange", size = 0.04*t_size) + 
  xlab("Longitude") + ggtitle("(b)") + 
  theme_bw() + 
  theme(text = element_text(size = 0.7*t_size), 
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
# print maps using patchwork
map_Osab_raw + map_Osab_spat_subsamp
### end map example data -----------------------------------------------------

### plot AUC for random CV ----------------------------------------------------
auc_all_models_plot <- ggplot(
  data = evals[evals$metric == "AUC" & 
                 as.character(evals$block_cv_range) == "random" & 
                 as.character(evals$test_data) == "spat_subsamp", ], 
  aes(
    x = factor(train_data, 
               levels = c("raw", "spat_subsamp"), 
               labels =  c("raw", "spatially\nunder-sampled")), 
    y = value, 
    color = factor(
      model, 
      levels = c("day_ll_rf", "spat_ll_rf","env_ll_rf", 
                 "env_spat_ll_rf"), 
      labels = c("\nseason +\nlist length\n", 
                 "\ncoordinates +\nseason +\nlist length\n",
                 "\nenvironment +\nseason +\nlist length\n", 
                 "\nenvironment + \ncoordinates +\nseason +\nlist length")))) + 
  # geom_violin(draw_quantiles = 0.5, size = 0.9) + 
  geom_boxplot(size = 0.9) +
  facet_wrap(~factor(species, 
                     levels = c("Macrosternodesmus palicola", 
                                "Boreoiulus tenuis", 
                                "Ommatoiulus sabulosus", 
                                "Blaniulus guttulatus", 
                                "Glomeris marginata", 
                                "Cylindroiulus punctatus"),
                     labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))) + 
  xlab("Training data") + 
  ylab("AUC\n(cross-validated)") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
auc_all_models_plot
### end AUC boxplots ----------------------------------------------------------

### AUC difference plots -----------------------------------------------------
# calculate mean and se of mean for AUC
auc_summary <- evals[evals$metric == "AUC" & 
                       as.character(evals$block_cv_range) == "random" & 
                       as.character(evals$test_data) == "spat_subsamp", ]
auc_summary <- group_by(auc_summary, species, model, train_data) %>%
  summarise(mean_auc = mean(value)) %>%
  pivot_wider(names_from = train_data, values_from = mean_auc) %>%
  mutate(change_after_subsampling = spat_subsamp - raw)

auc_means_plot <- ggplot(
  data = auc_summary, 
  aes(
    x = factor(species, 
               levels = c("Macrosternodesmus palicola", 
                          "Boreoiulus tenuis", 
                          "Ommatoiulus sabulosus", 
                          "Blaniulus guttulatus", 
                          "Glomeris marginata", 
                          "Cylindroiulus punctatus"),
               labels = c("M. palicola", 
                          "Boreoiulus\ntenuis", 
                          "O. sabulosus", 
                          "Blaniulus\nguttulatus", 
                          "G. marginata", 
                          "C. punctatus")), 
    y = change_after_subsampling, 
    color = factor(
      model, 
      levels = c("day_ll_rf", "spat_ll_rf","env_ll_rf", 
                 "env_spat_ll_rf"), 
      labels = c("\nseason +\nlist length\n", 
                 "\ncoordinates +\nseason +\nlist length\n",
                 "\nenvironment +\nseason +\nlist length\n", 
                 "\nenvironment + \ncoordinates +\nseason +\nlist length")), 
    shape = factor(
      model, 
      levels = c("day_ll_rf", "spat_ll_rf","env_ll_rf", 
                 "env_spat_ll_rf"), 
      labels = c("\nseason +\nlist length\n", 
                 "\ncoordinates +\nseason +\nlist length\n",
                 "\nenvironment +\nseason +\nlist length\n", 
                 "\nenvironment + \ncoordinates +\nseason +\nlist length")))) + 
  geom_point(size = t_size/4) + 
  # facet_wrap(~factor(species, 
  #                    levels = c("Macrosternodesmus palicola", 
  #                               "Boreoiulus tenuis", 
  #                               "Ommatoiulus sabulosus", 
  #                               "Blaniulus guttulatus", 
  #                               "Glomeris marginata", 
  #                               "Cylindroiulus punctatus"),
  #                    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))) + 
  xlab("Species") + 
  ylab("Change in AUC") + 
  scale_color_viridis_d(name = "Model", #option = "magma", 
                        begin = 0.1, end = 0.8, direction = -1) + 
  scale_shape(name = "Model") + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1))
auc_means_plot
### end AUC difference plots --------------------------------------------------

## plot all performance metrics for best model --------------------------------
evals_median_best <- filter(evals, model == "env_spat_ll_rf" & 
                              block_cv_range == "random" & 
                              test_data == "spat_subsamp") %>%
  select(species, train_data, metric, value) %>%
  group_by(species, metric, train_data) %>% 
  summarise(median = median(value))

evals_median_best$species <- factor(
  evals_median_best$species, 
  labels = gsub(" ", "\n", levels(evals_median_best$species)))

performance_best_mod_plot <- ggplot(
  data = evals_median_best, 
  aes(x = factor(train_data, 
                 levels = c("raw", "spat_subsamp"), 
                 labels =  c("raw", "spatially\nundersampled")), 
      y = median, 
      group = factor(species))) + 
  geom_point(aes(color = factor(species)), size = 0.1*t_size) + 
  geom_line(aes(color = factor(species)), size = 0.05*t_size) + 
  facet_wrap(~factor(metric, 
                     levels = c("AUC", "sensitivity", "specificity", "Kappa", 
                                "Brier"), 
                     labels = c("AUC", "Sensitivity", "Specificity", "Kappa", 
                                "Brier Score"), ordered = T)) + 
  xlab("Training data") + 
  ylab("value") + 
  # ggtitle(paste0("Random CV\nmodel resolution: ", analysis_resolution)) + 
  scale_color_viridis_d(begin = 0, end = 0.8) + 
  theme_bw() + 
  theme(text = element_text(size = 0.9*t_size), 
        legend.position = c(0.84, 0.18), legend.title = element_blank(), 
        legend.key.size = unit(2*t_size, "points"), 
        legend.text = element_text(size = 0.7*t_size))
performance_best_mod_plot
## end plot performance for best model ----------------------------------------
### end random CV -------------------------------------------------------------


### plot partial dependence --------------------------------------------------
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
              train_data = unique(train_data)) %>%
    filter(variable != "cos_doy" & variable != "sin_doy")
})
pd <- bind_rows(pd)
# make eastings and northings be in km instead of m
pd$x[pd$variable == "eastings" | pd$variable == "northings"] <- 
  pd$x[pd$variable == "eastings" | pd$variable == "northings"] / 1000

## make pd plots using random CV, spatially under-sampled training data, and
## the best model
pd_plots <- lapply(sp_to_fit, FUN = function(x, dat, vimp) {
  vi <- vimp[vimp$species == x, ] # get variable importance for this sp.
  vi <- vi[order(vi$MeanDecreaseGini, decreasing = T), ]
  vi <- vi[-which(grepl(".*doy", vi$variable))[2], ]
  vi$variable <- gsub(".*doy", "day_of_year", vi$variable)
  # make better variable names
  vi$variable <- gsub("arabl.*", "arable\nland", vi$variable)
  vi$variable <- gsub("artific.*", "artificial\nsurfaces", vi$variable)
  vi$variable <- gsub("elev", "elevation", vi$variable)
  vi$variable <- gsub("fores.*", "forest and\nsemi-natural\nland", vi$variable)
  vi$variable <- gsub("list_l.*", "checklist\nlength", vi$variable)
  vi$variable <- gsub("mean_rr", "annual\nprecipitation", vi$variable)
  vi$variable <- gsub("mean_tn", "annual\nminimum\ntemperature", vi$variable)
  vi$variable <- gsub("pastu.*", "pasture", vi$variable)
  vi$variable <- gsub("wetla.*", "wetlands", vi$variable)
  vi$variable <- gsub("day_o.*", "day of\nyear", vi$variable)
  
  # make variable column a factor
  vi$variable <- factor(vi$variable, levels = vi$variable, 
                        ordered = T)
  
  pdat <- dat[dat$species == x, ] # get pd data for this sp.
  # make better variable names
  pdat$variable <- gsub("arabl.*", "arable\nland", pdat$variable)
  pdat$variable <- gsub("artific.*", "artificial\nsurfaces", pdat$variable)
  pdat$variable <- gsub("elev", "elevation", pdat$variable)
  pdat$variable <- gsub("fores.*", "forest and\nsemi-natural\nland", 
                        pdat$variable)
  pdat$variable <- gsub("list_l.*", "checklist\nlength", pdat$variable)
  pdat$variable <- gsub("mean_rr", "annual\nprecipitation", pdat$variable)
  pdat$variable <- gsub("mean_tn", "annual\nminimum\ntemperature", 
                        pdat$variable)
  pdat$variable <- gsub("pastu.*", "pasture", pdat$variable)
  pdat$variable <- gsub("wetla.*", "wetlands", pdat$variable)
  pdat$variable <- gsub("day_o.*", "day of\nyear", pdat$variable)
  # make variable column a factor
  pdat$variable <- factor(pdat$variable, levels = levels(vi$variable))
 
  ggplot(data = pdat,
         aes(x = x, y = y)) +
    geom_point() +
    geom_line() +
    facet_wrap(~variable, scales = "free_x") +
    ylab("Partial dependence") + 
    xlab("") + 
    ggtitle(pdat$species[1]) + 
    theme_bw() + 
    theme(text = element_text(size = t_size), 
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
}, dat = pd[pd$train_data == "spat_subsamp" & pd$cv == "random", ], 
vimp = vimp[vimp$train_dat == "spat_subsamp" & vimp$cv == "random", ])

for(sn in c("Macrosternodesmus palicola", "Boreoiulus tenuis",
            "Ommatoiulus sabulosus",  "Blaniulus guttulatus")) {
  print(pd_plots[sn])
}
rm(pd, vimp)
### end plot partial dependence -----------------------------------------------


### plot predictions with standardized list length ----------------------------
# load standardized predictions
osab_preds <- list(raw = readRDS("./saved_objects/standard_predictions_env_spat_ll_rf_noSubSamp_Ommatoiulus_sabulosus1000.rds"), 
                   spat_subsamp = readRDS("./saved_objects/standard_predictions_env_spat_ll_rf_SubSamp_Ommatoiulus_sabulosus1000.rds"))

# define function to calculate standard error of mean
se <- function(x) sd(x)/sqrt(length(x))

## calculate mean predictions and se of mean in each grid cell
osab_ciWidth_data <- lapply(osab_preds, function(dat) {
  # get only width of 95% CI for predictions from random CV
  dat <- dat[dat$cv == "random", ] # use only random CV results
  sp <- "Ommatoiulus sabulosus"
  mod <- "env_spat_ll_rf"

  dat <- group_by(dat, en) %>%
    summarise(mean_prediction = mean(mean_pred, na.rm = T), 
              se = se(mean_pred), 
              eastings = mean(eastings), northings = mean(northings))
  dat
})

osab_maps_ciWidth <- lapply(osab_ciWidth_data, function(dat, annot, ir_TM75) {
  # map standard error of mean
  ggplot() +
    geom_sf(data = st_as_sf(ir_TM75), fill = NA) +
    geom_tile(data = dat,
              aes(x = eastings, y = northings, fill = se)) +
    ylab("") + xlab("") +
    scale_fill_gradient(low = "white", high = "red") +  
    guides(fill = guide_colorbar(title = "",
                                 barwidth = unit(0.4 * t_size, "points"))) +
    theme_bw() +
    theme(text = element_text(size = 0.55*t_size),
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
          plot.margin = unit(c(-0.25, -0.6, -0.25, -0.5), "lines"))
}, annot = annot, ir_TM75 = ir_TM75)

osab_maps_prediction <- lapply(
  osab_ciWidth_data, function(dat, annot, ir_TM75) {
    # map mean predictions
    ggplot() + 
      geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
      geom_tile(data = dat, 
                aes(x = eastings, y = northings, fill = mean_prediction)) +
      ylab("") + xlab("") + 
      guides(fill = guide_colorbar(title = "", 
                                   barwidth = unit(0.4 * t_size, "points"))) +
      theme_bw() + 
      theme(text = element_text(size = 0.55*t_size), 
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), 
          plot.margin = unit(c(-0.25, -0.4, -0.25, -0.5), "lines"))
}, annot = annot, ir_TM75 = ir_TM75)

# map observed checklists
osab_map_observed <- ggplot() + 
  geom_sf(data = st_as_sf(ir_TM75), fill = NA) + 
  geom_sf(data = st_as_sf(mill_wide[mill_wide$`Ommatoiulus sabulosus` == 0, ]), 
          color = "dark grey", size = 0.02*t_size) + 
  geom_sf(data = st_as_sf(mill_wide[mill_wide$`Ommatoiulus sabulosus` == 1, ]), 
          color = "dark orange", size = 0.02*t_size) + 
  geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label), 
            size = 0.12*t_size) + 
  geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
               arrow = arrow(length = unit(0.1, "npc"))) + 
  ylab("") + xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 0.55*t_size), 
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), 
        plot.margin = unit(c(-0.25, -0.4, -0.5, -0.5), "lines"))

# combine maps in one ordered list
osab_maps <- list(observed = osab_map_observed, 
                  raw_pred = osab_maps_prediction$raw, 
                  spat_subsamp_pred = osab_maps_prediction$spat_subsamp, 
                  raw_ciWidth = osab_maps_ciWidth$raw, 
                  spat_subsamp_ciWidth = osab_maps_ciWidth$spat_subsamp)

# add titles
osab_maps[[1]] <- osab_maps[[1]] + ggtitle("(a)") + ylab("Latitude") 
osab_maps[[2]] <- osab_maps[[2]] + ggtitle("(b)")
osab_maps[[3]] <- osab_maps[[3]] + ggtitle("(c)")
osab_maps[[4]] <- osab_maps[[4]] + ggtitle("(d)") + 
  ylab("Latitude") + xlab("Longitude")
osab_maps[[5]] <- osab_maps[[5]] + ggtitle("(e)") + xlab("Longitude")

# print plots using patchwork
osab_maps[[1]]+ osab_maps[[2]] + osab_maps[[3]] + plot_spacer() + 
  osab_maps[[4]] + osab_maps[[5]] + plot_layout(ncol = 3)
### end plot standardized predictions -----------------------------------------


### print tables and numbers for text -----------------------------------------
# Number of repeat visits to grid cells (though there could be multiple 
# locations visited within a grid cell, so these are not necessarily true 
# repeat visits)
repeat_visits <- data.frame(
  en = paste0(mill_wide$eastings, "_", mill_wide$northings), 
  year = mill_wide$year)
# table(repeat_visits$en)[order(table(repeat_visits$en), decreasing = T)]
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

# Number of detections per species when using 1 km resolution
n_detections_per_species_1km <- data.frame(
  table(mill$Genus_species[mill$Precision <= 1000]))
n_detections_per_species_1km <- n_detections_per_species_1km[order(
  n_detections_per_species_1km$Freq, decreasing = FALSE), ]
colnames(n_detections_per_species_1km) <- c("species", "number_of_detections")
n_detections_per_species_1km <- n_detections_per_species_1km[
  n_detections_per_species_1km$species %in% sp_to_fit, ]
# add proportion of checklists that have a detection of each species
n_detections_per_species_1km$proportion_detections <- round(
  n_detections_per_species_1km$number_of_detections / nrow(mill_wide), 
  digits = 3)
n_detections_per_species_1km

# some AUC summary stats
auc_summary[auc_summary$model == "day_ll_rf", 1:4]
auc_summary[auc_summary$model == "env_spat_ll_rf", 1:4]
summary(auc_summary$spat_subsamp[auc_summary$model == "day_ll_rf"])
summary(auc_summary$spat_subsamp[auc_summary$model == "env_spat_ll_rf"])

# difference between simplest and most complex model for each sp
summary(auc_summary$spat_subsamp[auc_summary$model == "day_ll_rf"] - 
  auc_summary$spat_subsamp[auc_summary$model == "env_spat_ll_rf"])

auc_summary[auc_summary$species == "Cylindroiulus punctatus", 1:4]

data.frame(group_by(auc_summary, species) %>%
  arrange(spat_subsamp, .by_group = TRUE))
### end print tables and numbers for text  ------------------------------------

### save plots ----------------------------------------------------------------
## save as jpg
ggsave("Fig1.jpg", map_Osab_raw + map_Osab_spat_subsamp, 
       width = 20, height = 20/2, units = "cm", 
       device = "jpg")
ggsave("Fig2.jpg", auc_means_plot, width = 20, height = 20, units = "cm", 
       device = "jpg")
ggsave("Fig3.jpg", auc_all_models_plot, width = 20, height = 20, units = "cm",
       device = "jpg")
ggsave("Fig4.jpg", performance_best_mod_plot, width = 20, height = 20, 
       units = "cm", device = "jpg")
ggsave("Fig5.jpg", pd_plots[["Ommatoiulus sabulosus"]] + ggtitle(""), 
       width = 20, height = 20, units = "cm", device = "jpg")
ggsave("Fig6.jpg", osab_maps[[2]] + ylab("Latitude") + ggtitle("(a)") + 
         geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
         geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label), 
                   size = 0.2*t_size) + 
         geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
                      arrow = arrow(length = unit(0.1, "npc"))) + 
         theme(text = element_text(size = 0.8*t_size)) +
         osab_maps[[3]] + ggtitle("(b)") +  
         theme(text = element_text(size = 0.8*t_size)) +
         osab_maps[[4]] + ggtitle("(c)") + 
         theme(text = element_text(size = 0.8*t_size)) +
         osab_maps[[5]] + ggtitle("(d)") + 
         theme(text = element_text(size = 0.8*t_size)) + 
         plot_layout(ncol = 2), 
       width = 18, height = 18, units = "cm", device = "jpg")


## save as eps
ggsave("Figure 1.eps", spat_evenness_boxplot, width = 25, height = 25, 
       units = "cm", device = "eps")

### end save plots ------------------------------------------------------------
