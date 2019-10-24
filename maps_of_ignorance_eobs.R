##########################################
## Maps of ignorance
## Show distribution of samples in environmental space compared to the 
## distribution under random sampling. Also calculate mutual information.
## 
## 
## NOTES: currently uses training AND test data to get sense of overall 
## sampling distribution
## 
## TODO:  - for mutual information calculation, need to add the min and max 
##          values to the variables when discretizing. (line 300ish)
##        - re-do using number of checklists instead of number of records?  I think number of records might make biases appear where there are none in e.g. Odonata where less species rich areas at the colder end of the gradient would have fewer records even if they had appropriate sampling effort.
##        - Add Hannah's vascular plant dataset?
## 
## author: W. Gaul
## created: 31 Oct. 2017
## last modified: 13 Nov 2017
#########################################

setwd("~/Documents/Data_Analysis/UCD/maps_of_ignorance/")
rm(list = ls())
library(wgutil)
library(Hmisc)
library(ggridges)
library(parallel)
library(rgdal)
library(sp)
library(entropy)
library(tidyverse)
source("../../R_templates/biorecording-master/eastnorth2os.R")
source("../../R_templates/biorecording-master/os2eastnorth.R")

# workspace saved, so can to avoid running entire script
load("maps_of_ignorance_eobs.RData")

## ------------------------- load data ---------------------------------------
# climate data
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/annual_precip_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/elevation_hec.RData")

# bryophyte training and test data
load("~/Documents/Data_Analysis/UCD/bryophytes/data/bry_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/bryophytes/data/test_vault/bry_test_data.Rdata")

# NBDC training and test data
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/NBDC_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

hec_names <- read_csv(
  "~/Documents/Data_Analysis/UCD/mapping/data/Irish_land_hectads.csv")

# Load Ireland coastline
ir <- readOGR(dsn="../mapping/data", 
              layer="ireland_coastline")
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))

# combine test and training data for all groups
bry <- bind_rows(bry_training_data, bry_test_data)
nbdc <- test_data
for (i in 1:length(names(nbdc))) {
  taxon <- names(nbdc)[i]
  nbdc[[i]] <- bind_rows(training_data[which(names(training_data) == taxon)], 
                         test_data[which(names(test_data) == taxon)])
}

# get mean annual precipitation values for each hectad
mean_rr <- as.data.frame(krg_mean_rr_predict, xy=T)
mean_rr$x_corner <- mean_rr$s1 - 5000
mean_rr$y_corner <- mean_rr$s2 - 5000
mean_rr$en <- paste(round(mean_rr$x_corner), 
                    round(mean_rr$y_corner), sep = "_")

# get elevation values
elev_df <- as.data.frame(krg_elev_predict, xy = T)
elev_df$x_corner <- elev_df$s1 - 5000
elev_df$y_corner <- elev_df$s2 - 5000
elev_df$en <- paste(round(elev_df$x_corner), 
                    round(elev_df$y_corner), sep = "_")
colnames(elev_df)[1] <- "elevation"

# clean workspace
rm(list = c("test_data", "training_data", "bry_test_data", 
            "bry_training_data"))
## ------------------------ end load data ------------------------------------

## -------------------------- data cleaning ----------------------------------
# convert NBDC gridRefs to hectads (many currently 1km blocks)
for (i in 1:length(nbdc)) {
  nbdc[[i]]$hectad <- gridref_to_hec(nbdc[[i]]$GridReference)
}

# remove british hecs from bryophyte data
bry <- bry[which(nchar(bry$OSGR_10km) == 3), ]

# select which datasets to use
over_50 <- lapply(nbdc, FUN = nrow)
over_50 <- over_50[over_50 > 50]
nbdc <- nbdc[which(names(nbdc) %in% names(over_50))]
rm(over_50)

drop_taxa <- c("clubmoss", "quillwort", "stonewort", "hornwort", "liverwort", 
               "moss")
nbdc_bry <- nbdc[which(names(nbdc) %nin% drop_taxa)]
nbdc_bry$bry_BBS <- bry
nbdc_bry$bry_BBS$hectad <- nbdc_bry$bry_BBS$OSGR_10km
#rm(list = c("bry", "nbdc"))

add_en <- function(x){
  # add eastings and northings to each taxonomic dataset
  # ARGS: x - df of the taxanomic dataset
  # for records which have a hectad code us Jon's os2eastnorth() function
  x_h <- x[which(!is.na(x$hectad)), ]
  x_en <- os2eastnorth(x_h$hectad)
  en_df <- as_tibble(x_en$en)
  en_df$hectad <- x_en$gridref
  x_h <- left_join(x_h, en_df, by = "hectad")
  
  # for records with no hectad code us lat/long instead
  x_ll <- x[which(is.na(x$hectad)), ]
  if (nrow(x_ll) > 0) {
    x_ll <- lat_lon_to_osi(x_ll, orig_crs = "+init=epsg:4326", 
                           precision = 10000)
  }
  
  x <- full_join(x_h, x_ll)
  x$en <- paste(round(x$eastings),
                round(x$northings), sep = "_")
  x
}

nbdc_bry <- mclapply(nbdc_bry, FUN = add_en, mc.cores = 2)

add_mean_precip <- function(x, hec_precips){
  # add average temp on Jan 1 to taxonomic datasets
  # ARGS: x - df holding the taxonomic dataset
  #       hec_precips - df holding the mean annual precipitation values for 
  #                   each hectad and a column named "en" holding character 
  #                   string of e and n pasted
  x <- left_join(x, hec_precips, by = "en")
  x
}

nbdc_bry <- mclapply(nbdc_bry, FUN = add_mean_precip,
                     hec_precips = mean_rr[, c("en", "mean_annual_rr")], 
                     mc.cores = 2)

add_elev <- function(x, hec_elevs){
  # add elevation to taxonomic datasets
  # ARGS: x - df holding the taxonomic dataset
  #       hec_elevs - df holding the elevation values for each hectad and a
  #       column named "en" holding character string of e and n pasted
  browser()
  x <- left_join(x, hec_elevs, by = "en")
  colnames(x)[which(colnames(x) == "elevation")] <- "elevation"
  x
}

nbdc_bry <- mclapply(nbdc_bry, FUN = add_elev,
                     hec_elevs = elev_df[, c("en", "elevation")], 
                     mc.cores = 2)

## -------------------------- end cleaning -----------------------------------

## -------------------- randomly sample n points from IE ----------------------
# make a list of random datasets to correspond to real taxonomic datasets

sample_points <- function(n, ireland_coast, env_var) {
  # sample n random points from Ireland
  # ARGS: n - number of points to draw
  #       ireland_coast - shapefile of irish coast
  #       env_var - df holding easings/northings and env. variable values
  # OUTPUTS: a df of random points with eastings/northings and env. var. values
  
  rand_points <- sp::spsample(x = ireland_coast, n = n, type = "random")@coords
  rand_points <- as.data.frame(rand_points)
  rand_points$east_round <- rand_points$x - rand_points$x%%10000
  rand_points$north_round <- rand_points$y - rand_points$y%%10000
  rand_points$en <- paste(rand_points$east_round, rand_points$north_round, 
                          sep = "_")
  # add avg. daily temp values to random points
  rand_points <- left_join(rand_points, env_var, by = "en")
  rand_points
}

set.seed(351, kind = "L'Ecuyer-CMRG")
rand_points_mean_rr <- mclapply(lapply(nbdc_bry, FUN = nrow), 
                                  FUN = sample_points, 
                                  ireland_coast = ir_TM75, 
                                  env_var = mean_rr[, c("en", 
                                                        "mean_annual_rr")], 
                                  # mc.preschedule = T, 
                                  # mc.set.seed = T, 
                                  mc.cores = 3)

rand_points_elev <- mclapply(lapply(nbdc_bry, FUN = nrow), 
                             FUN = sample_points, 
                             ireland_coast = ir_TM75, 
                             env_var = elev_df[, c("en", 
                                                   "elevation")], 
                             # mc.preschedule = T, 
                             # mc.set.seed = T, 
                             mc.cores = 3)

## ----------------- end random sampling IE --------------------------------


## ---------------------- plot all taxonomic groups ---------------------------
make_plots <- function(x, taxon_name, random_points, env_variable) {
  # make ridge plots and overlapping histograms of sampling density in 
  # environmental space for real vs. randomly selected points
  # ARGS: x - df of real biological records with a column holding the values
  #           of the environmental variable
  #       taxon_name - character string giving taxon name for plotting
  #       random_points - dataframe of randomly selected points with env. var.
  real <- x[, which(colnames(x) == env_variable)]
  random <- random_points[, which(colnames(random_points) == env_variable)] 
  
  # make a df for ggplot 
  real_temps <- data.frame(value = real, 
                           dataset = as.character("real"))
  rand_temps <- data.frame(value = random, 
                           dataset = as.character("random"))
  real_random <- bind_rows(real_temps, rand_temps)
  rm(real_temps, rand_temps)
  
  # plot histograms of env. var. values for real and random data
  data(cbPalette)
  h <- ggplot(data = real_random, aes(x = value, 
                                      group = dataset, 
                                      color = dataset, 
                                      fill = dataset)) + 
    geom_histogram(position = "identity", alpha = 0.3, bins = 50) + 
    ggtitle(paste0("Sampling distribution - ", taxon_name, 
                   "\n", env_variable)) + 
    xlab(env_variable) + 
    scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) + 
    theme_bw()
  
  ridge <- ggplot(data = real_random, aes(x = value, 
                                          y = dataset)) + 
    geom_density_ridges() + 
    ggtitle(paste0("Sampling distribution - ", taxon_name, 
                   "\n", env_variable)) + 
    theme_ridges()
  
  plot_list <- list(h, ridge)
  plot_list
}

precip_plots <- mapply(FUN = make_plots, 
                       nbdc_bry, 
                       names(nbdc_bry), 
                       random_points = rand_points_mean_rr, 
                       MoreArgs = list(env_variable = "mean_annual_rr"))

elev_plots <- mapply(FUN = make_plots, 
                     nbdc_bry, 
                     names(nbdc_bry), 
                     random_points = rand_points_elev, 
                     MoreArgs = list(env_variable = "elevation"))

## print plots to pdf --------------------------------------------------------
pdf_out <- F
if (pdf_out == T) {
  pdf("sampling_mismatch_plots_annual_rr.pdf")
  for(i in 1:length(precip_plots)) {
    print(precip_plots[i])
  }
  dev.off()
  
  pdf("sampling_mismatch_plots_elevation.pdf")
  for(i in 1:length(elev_plots)) {
    print(elev_plots[i])
  }
  dev.off()
}


## calculate mutual information ----------------------------------------------
# try example for moths (low mutual information) and birds (high mut. inf).
moth_disc <- discretize2d(
  nbdc_bry$`insect - moth`$elevation[order(nbdc_bry$`insect - moth`$elevation)], 
  rand_points_elev$`insect - moth`$elevation[order(
    rand_points_elev$`insect - moth`$elevation)], 
  numBins1 = 30, numBins2 = 30)
mi.empirical(moth_disc)

bird_rr_disc <- discretize2d(
  nbdc_bry$bird$mean_annual_rr[order(nbdc_bry$bird$mean_annual_rr)], 
  rand_points_mean_rr$bird$mean_annual_rr[order(
    rand_points_mean_rr$bird$mean_annual_rr)], 
  numBins1 = 100, numBins2 = 100)
mi.empirical(bird_rr_disc, unit = "log2")

# this might be giving values too high because the binning goes from the min
# to max values for the variable in each independent dataset, but the real data
# has no observations from low or high ends of the spectrum.  The result is that
# the first few bins for joint distribution look like they line up well b/c they
# have similar numbers of observations, but in fact they should do poorly b/c
# the absolute values of rr for the bins are very different.  So, need to add
# the min and max values to the variables when discretizing.
spider_rr_disc <- discretize2d(
  nbdc_bry$`spider (Araneae)`$mean_annual_rr[order(
    nbdc_bry$`spider (Araneae)`$mean_annual_rr)], 
  rand_points_mean_rr$`spider (Araneae)`$mean_annual_rr[order(
    rand_points_mean_rr$`spider (Araneae)`$mean_annual_rr)], 
  numBins1 = 100, numBins2 = 100)
mi.empirical(spider_rr_disc, unit = "log2")


## ----------------------------------------------------------------------------
# ## -------------------------- Bryophytes -------------------------------------
# bry_en <- os2eastnorth(bry$OSGR_10km)
# bry$eastings <- NA
# bry$northings <- NA
# for(i in 1:nrow(bry)) {
#   bry$eastings[i] <- bry_en$en[which(bry_en$gridref == bry$OSGR_10km[i]), 1]
#   bry$northings[i] <- bry_en$en[which(bry_en$gridref == bry$OSGR_10km[i]), 2]
# }
# bry$en <- paste(round(bry$eastings), 
#                 round(bry$northings), sep = "_")
# 
# # add env. var. values to sp data
# bry$mean_jan_daily_temp <- NA
# for(i in 1:nrow(bry)) {
#   bry$mean_jan_daily_temp[i] <- jan1_mean_temp_df$var1.pred[which(
#     jan1_mean_temp_df$en == bry$en[i])]
# }
# 
# 
# # make a df for ggplot of 2n rows where each record is either a real or random
# # record, there is a col for env. var. value, and a column indicating real
# # or random
# bry_temps <- data.frame(jan1_temp = bry$mean_jan_daily_temp, 
#                         dataset = as.character("real"))
# rand_temps <- data.frame(jan1_temp = rand_points$mean_jan_daily_temp, 
#                          dataset = as.character("random"))
# bry_real_random <- bind_rows(bry_temps, rand_temps)
# rm(bry_temps, rand_temps)
# 
# ## ------------------- RESULT PLOT --------------------------------------------
# # plot histograms of env. var. values for real and random data
# ggplot(data = bry_real_random, aes(x = jan1_temp, 
#                                    group = dataset, 
#                                    color = dataset, 
#                                    fill = dataset)) + 
#   geom_histogram(position = "identity", alpha = 0.3, 
#                  bins = 30) + 
#   ggtitle("Sampling distribution - real vs. randomly selected points\nBryophytes") + 
#   theme_bw()
# 
# ggplot(data = bry_real_random, aes(x = jan1_temp, 
#                                    y = dataset)) + 
#   geom_density_ridges() + 
#   ggtitle("Sampling distribution relative to avg. temp on Jan 1\nBryophytes") + 
#   theme_ridges()
# ## ----------------------------- end Bryophytes -------------------------------
# 
# ## ------------------------------- moths -------------------------------------
# moths <- nbdc$`insect - moth`
# 
# moths_en <- os2eastnorth(moths$hectad)
# moths$eastings <- NA
# moths$northings <- NA
# for(i in 1:nrow(moths)) {
#   moths$eastings[i] <- moths_en$en[which(
#     moths_en$gridref == moths$hectad[i]), 1]
#   moths$northings[i] <- moths_en$en[which(
#     moths_en$gridref == moths$hectad[i]), 2]
# }
# moths$en <- paste(round(moths$eastings), 
#                   round(moths$northings), sep = "_")
# 
# # add env. var. values to sp data
# moths$mean_jan_daily_temp <- NA
# for(i in 1:nrow(moths)) {
#   moths$mean_jan_daily_temp[i] <- jan1_mean_temp_df$var1.pred[which(
#     jan1_mean_temp_df$en == moths$en[i])]
# }
# 
# # make a df for ggplot 
# moths_temps <- data.frame(jan1_temp = moths$mean_jan_daily_temp, 
#                           dataset = as.character("real"))
# rand_temps_moths <- data.frame(jan1_temp = rand_points$mean_jan_daily_temp, 
#                                dataset = as.character("random"))
# moths_real_random <- bind_rows(moths_temps, rand_temps_moths)
# rm(moths_temps, rand_temps_moths)
# 
# ## ------------------- RESULT PLOT --------------------------------------------
# # plot histograms of env. var. values for real and random data
# ggplot(data = moths_real_random, aes(x = jan1_temp, 
#                                      group = dataset, 
#                                      color = dataset, 
#                                      fill = dataset)) + 
#   geom_histogram(position = "identity", alpha = 0.3, bins = 50) + 
#   ggtitle("Sampling distribution - real vs. randomly selected points\nMoths") + 
#   theme_bw()
# 
# ggplot(data = moths_real_random, aes(x = jan1_temp, 
#                                      y = dataset)) + 
#   geom_density_ridges() + 
#   ggtitle("Sampling distribution relative to avg. temp on Jan 1\nMoths") + 
#   theme_ridges()
# ## ----------------------------- end Moths -------------------------------
# 


save.image("maps_of_ignorance_eobs.RData")
