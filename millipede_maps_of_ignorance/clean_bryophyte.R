##################
## Clean bryophyte data from NBN Atlas
## Bryophyte data for Great Britain and Ireland from the British Bryological 
## Society held by BRC: both 2014 Atlas data and data compiled post-Atlas
## from the NBN gateway - downloaded 24 Aug 2017
## data provided under a  CC-BY licence
##
## TODO: 
## 
## NOTES: This script reads in the original data downloaded from NBN, which 
## includes names of recorders.  This script then anonymises the names and saves
## the resulting dataset as .csv and .RData files.  After I run this (13 Feb
## 2018), I will delete names column from the original dataset.
##
## author: Willson Gaul
## created: 7 June 2017
## last modified: 16 Feb 2018 by wg (deleted recorder names)
###################

setwd("~/Documents/Data_Analysis/UCD/bryophytes/wg_scripts/")
rm(list = ls())
library(wgutil)
library(stringr)
library(doParallel)
library(foreach)
library(doRNG)
library(Hmisc)
library(tidyverse)
# source("../../../R_templates/wg_util.R")

## read in the datasets:
atlas <- read_csv(file = "../data/bryophyte_2014_atlas.csv", quote = "\"")
post_atlas <- read_csv(file = "../data/bryophyte_post_2014atlas.csv", 
                       quote = "\"")

bry <- bind_rows(atlas, post_atlas)
rm(list = c("atlas", "post_atlas"))


# remove white space from column names
reduce_underscores <- function(x) {
  while(any(grepl(pattern = "__", x))) {
    x <- gsub(pattern = "__", x, rep = "_")
  }
  x
}
colnames(bry) <- make.names(colnames(bry))
colnames(bry) <- gsub(colnames(bry), 
                      pattern = "\\.", 
                      rep = "_")
colnames(bry) <- reduce_underscores(colnames(bry))

bry <- select(bry, -c(Basis_Of_Record_original, Basis_of_record, Sex, 
                      Longitude_original, Geodetic_Datum, Licence))

# fix column name  
#   - 'Decimal_latitude_WGS84_1' actually gives longitude according
#     to 'headings.csv' metadata file.
colnames(bry)[which(
  colnames(bry) == "Decimal_latitude_WGS84_1")] <- "Decimal_longitude_WGS84"

## TODO: Parse collector names field to get individual collector names and to 
## identify how many collectors were on a survey. (see code in 
## 'spread_collector_names_bryophyte.R')

# parse collector names
collectors <- mclapply(bry$Collector, 
                       FUN = str_split, 
                       pattern = "; ", mc.cores = 3)
# make codes for names
names_df <- data.frame(unique_names = unique(unlist(collectors)))
names_df$code <- as.numeric(as.factor(names_df$unique_names))
names_df$code <- paste0("p_", names_df$code)
group_indicator <- c(".*group.*", ".*club.*", ".*committee.*", ".*union.*", 
                     ".*survey.*", ".*society.*", ".*unit.*", ".*meeting.*", 
                     ".*excursion.*", ".*council.*", ".*association.*", 
                     ".*BBS.*", ".*bryologists.*", ".*atlas.*", ".*team.*", 
                     ".*BSBI.*", ".*museum.*", ".*partnership.*", 
                     ".*naturalist.*", ".*centre.*", ".*university.*", 
                     ".*agency.*", ".*party.*")
names_df$group <- grepl(pattern = paste(group_indicator, collapse = "|"), 
                        names_df$unique_names, ignore.case = TRUE)

# create dataframe with individual columns to hold individual collector names 
# (from parsed "Collector" field)
obs_names <- data.frame(matrix(ncol = 0, nrow = nrow(bry)))
obs_names$observer_1 <- NA
obs_names$observer_2 <- NA
obs_names$observer_3 <- NA
obs_names$observer_4 <- NA
obs_names$observer_5 <- NA

# To run in parallel, split collectors list up into 3 parts.  
col_1 <- collectors[1:round((1/3)*length(collectors))]
col_2 <- collectors[(1 + round((1/3)*length(collectors))):round(
  (2/3)*length(collectors))]
col_3 <- collectors[(1 + round((2/3)*length(collectors))):length(collectors)]

col_list <- list(col_1, col_2, col_3) # make list to use in lapply

anonymize_names <- function(collectors, empty_df, codes_df) {
  # put the names of individual recorders into a df from the list returned by
  # mclapply above
  # ARGS: collectors - a list, with each element a list of names of recorders
  #       for a single site
  #       empty_df - a df with a column for each individual name.  To be filled.
  #       codes_df - a df with anonymous codes for each individual name.
  for (i in 1:length(collectors)) { # put names in df  
    for (j in 1:length(collectors[[i]][[1]])) { 
      if(j > 5) stop("j bigger than 5")
      nm <- collectors[[i]][[1]][[j]]
      if(is.na(nm)) {
        obs_names[i, j] <- NA
      } else {
        obs_names[i, j] <- codes_df$code[which(codes_df$unique_names == nm)]
      }
    }
  }
  obs_names[1:length(collectors), ]
}

# anonymize names
names_df_list <- mclapply(col_list, FUN = anonymize_names, empty_df = obs_names,
                          codes_df = names_df, mc.cores = 3)
obs_names <- bind_rows(names_df_list)
# join individual anonymised names to bry dataframe
bry <- cbind(bry, obs_names) 

# Add "group" column from names_df matching by 'code'.  Do this for each column
# of anonymised names
bry <- left_join(bry, names_df[, 2:3], by = c("observer_1" = "code"))
bry <- left_join(bry, names_df[, 2:3], by = c("observer_2" = "code"))
bry <- left_join(bry, names_df[, 2:3], by = c("observer_3" = "code"))
bry <- left_join(bry, names_df[, 2:3], by = c("observer_4" = "code"))
bry <- left_join(bry, names_df[, 2:3], by = c("observer_5" = "code"))

bry <- bry %>% # make a column indicating whether a group is listed as collector
  rowwise() %>%
  mutate(group_listed_as_observer = sum(group, group.y.y, group.x.x, 
                                        group.y, group.x, na.rm = TRUE)) %>%
  select(-c(group, group.y.y, group.x.x, group.y, group.x)) %>%
  select(-Collector) # remove real collector names
bry$group_listed_as_observer <- as.logical(bry$group_listed_as_observer)

# remove all objects with individual names
rm(names_df, collectors, col_1, col_2, col_3, col_list, nm)

write_csv(bry, path = "../data/bryophyte_atlasAndPostAtlas_anonymized.csv")
save(bry, file = "../data/bry.Rdata")
### ----------------------------------------------------------------------------













