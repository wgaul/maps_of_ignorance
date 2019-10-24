####################################
## Maps of Ignorance - spatial coverage
## Distance to nearest record in space
## 
## author: Willson Gaul
## created: 22 Oct 2019 
## last modified: 24 Oct 2019
#####################################

setwd("~/Documents/Data_Analysis/UCD/maps_of_ignorance/")
library(wgutil)
library(Hmisc)
library(raster)
library(spatstat)
library(GGally)
library(tidyverse)
source("~/Documents/Data_Analysis/R_templates/biorecording-master/os2eastnorth.R")

### load data --------------------------------------------------------------
# bryophyte training and test data
load("~/Documents/Data_Analysis/UCD/bryophytes/data/bry_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/bryophytes/data/test_vault/bry_test_data.Rdata")

# NBDC training and test data
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/NBDC_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

hec_names <- read_csv("~/Documents/Data_Analysis/UCD/mapping/data/Irish_land_hectads.csv")
### end load data --------------------------------------------------------------

### format data ------------------------------------------------------------
bry_test_data$split <- rep("test", nrow(bry_test_data))
bry_training_data$split <- rep("training", nrow(bry_training_data))

# combine test and training data for all groups
bry <- bind_rows(bry_training_data, bry_test_data)
nbdc <- test_data
for (i in 1:length(names(nbdc))) {
  taxon = names(nbdc)[i]
  nbdc[[i]] <- bind_rows(training_data[which(names(training_data) == taxon)], 
                         test_data[which(names(test_data) == taxon)])
}
# remove NBDC bryophyte groups, since I have the BBS bry data
nbdc <- nbdc[-c(which(names(nbdc) %in% c("liverwort", "hornwort", "moss")))]

# convert NBDC gridRefs to hectads (many currently 1km blocks)
for (i in 1:length(nbdc)) {
  nbdc[[i]]$hectad <- gridref_to_hec(nbdc[[i]]$GridReference)
}

# remove british hecs from bryophyte data
bry <- bry[which(nchar(bry$OSGR_10km) == 3), ]

# make a version of bird data without Bird Atlas 2007 - 2011
nbdc$bird_noAtlas <- nbdc$bird[nbdc$bird$DataSetHash != 
                                 "AB385E4D6E73BFA8C42CD208BD00B1E0", ]
# make a version of bird data with only Bird Atlas 2007 - 2011
nbdc$bird_Atlas <- nbdc$bird[nbdc$bird$DataSetHash == 
                               "AB385E4D6E73BFA8C42CD208BD00B1E0", ]

# make unique checklist identifier
nbdc <- lapply(nbdc, FUN = function(x) {x$checklist_ID <- paste(
  x$StartDate, x$GridReference, x$SiteName, x$RecorderHash, sep = "_"); x})
bry$checklist_ID <- paste(bry$Year, bry$Month, bry$Date, bry$OSGR,
                          bry$Locality, bry$observer_1,
                          bry$observer_2, bry$observer_3, bry$observer_4,
                          bry$observer_5, sep = "_")

rm(test_data, training_data, bry_test_data, bry_training_data) # clean workspace
### end format data ------------------------------------------------------------

stop("extend function to more dimensions?")
dist_to_nearest_record <- function(df, query_points, coords) {
  # calculate the euclidean distance to the nearest record in df from the 
  # points in query_points.
  # All dimensions should be scaled before being sent into this function.
  # ARGS: df - data frame of observations containing columns for coordinates
  #       query_points - data frame giving the points from which distances to 
  #           nearest records are desired.  This should have two columns for 
  #           coordinates with the same names as those in df and coords
  #       coords - character vector of any length giving the names of the columns
  #           containing the coordinates of each record
  # VALUE: a data frame with all the points in query_points and the distances 
  #         from those points to the nearest records in df
  df_ranges <- lapply(df[, coords], function(x) max(x) - min(x))
  if(any(df_ranges > 1)) stop("Scale dimensions of records before using this function.")
  q_ranges <- lapply(query_points[, coords], function(x) max(x) - min(x))
  if(any(q_ranges > 1.5)) stop("Scale dimensions of query points before using this function.")
  
  
  dists <- query_points # make df to hold resulting minimum distances
  dists$dist_to_nearest_rec <- NA

  for (i in 1:nrow(query_points)) {
    # format data for use in dist by adding a single query_points row to df
    d <- data.frame(bind_rows(query_points[i, coords], df[, coords]))
    # coerce all columns to numeric
    for(j in 1:ncol(d)) {
      d[, j] <- as.numeric(d[, j])
    }

    dm <- dist(d, method = "euclidean", diag = F, upper = F)
    dm <- as.matrix(dm)
    dists$dist_to_nearest_rec[i] <- min(dm[1, 2:ncol(dm)]) # find smallest distance
  }
  # make distances relative so most distant grid query point has distance = 1
  dists$dist_to_nearest_rec <- dists$dist_to_nearest_rec/max(dists$dist_to_nearest_rec)
  dists # return df with distance from each query point to nearest record
}

# test 
od <- nbdc$`insect - dragonfly (Odonata)`
query_points <- hec_names
query_points$StartDate <- as.Date("2017-10-23")

# scale dimensions to 0 to 1 range
# eastings and northings should be scaled in the same way so space isn't more
# important in one direction
spat_range <- max(c(max(od$eastings) - min(od$eastings), 
                max(od$northings) - min(od$northings)))
od$eastings_sc <- od$eastings/spat_range
od$northings_sc <- od$northings/spat_range

date_range <- max(as.numeric(od$StartDate)) - min(as.numeric(od$StartDate))
od$StartDate_sc <- as.numeric(od$StartDate)/date_range

# a hectad (10 km) is about as much "distance" as a year in this dataset if 
# dimensions are scaled to 0 to 1:
spat_range/10000 # spatial range is about 44 hectads, so 0 to 1 will span 44 hectads
year(max(od$StartDate)) - year(min(od$StartDate)) # 0 to 1 will span 47 years

# center and scale query point variables
query_points$eastings_sc <- query_points$eastings/spat_range
query_points$northings_sc <- query_points$northings/spat_range

query_points$StartDate_sc <- as.numeric(query_points$StartDate)/date_range

# measure spatial/temporal distance
test_dist_spTime <- dist_to_nearest_record(od[, ], 
                                           query_points = query_points, 
                                           coords = c("eastings_sc", 
                                                      "northings_sc", 
                                                      "StartDate_sc"))
ggplot(data = test_dist_spTime, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  # geom_point(data = od[1:1000, ], aes(x = eastings, y = northings,
  #                                     size = StartDate),
  #            color = "orange") +
  # geom_point(data = od[tp, ], aes(x = eastings, y = northings),
  #            color = "red") + 
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial / temporal distance")

# measure spatial distance only
test_dist_space <- dist_to_nearest_record(od[1:1000, ], 
                                           query_points = query_points, 
                                           coords = c("eastings_sc", 
                                                      "northings_sc"))
ggplot(data = test_dist_space, aes(x = eastings, y = northings)) + 
  geom_point(aes(color = dist_to_nearest_rec), size = 5) + 
  geom_point(data = od[1:1000, ], aes(x = eastings, y = northings),
             color = "orange") +
  scale_color_gradient(limits = c(0, 1)) + 
  ggtitle("spatial distance")
