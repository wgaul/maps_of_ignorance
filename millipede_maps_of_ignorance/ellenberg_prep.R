####################################
## Prepare plant records as indicator of soil acidity using Ellenberg values
## 
## author: Willson Gaul wgaul@hotmail.com 
## created: 28 Feb 2020
## last modified: 28 Feb 2020
####################################

## load BryAtt and PLANTATT data
bryatt <- read.csv("./data/BryoAtt/Bryoatt_Oct_2007.csv", 
                   stringsAsFactors = FALSE)
bryatt <- select(bryatt, "Taxon.name", "R")
plantatt <- read.csv("./data/PLANTATT/PLANTATT_19_Nov_08.csv", 
                     stringsAsFactors = FALSE)
plantatt <- select(plantatt, "Taxon.name", "R")

### load bryophyte and plant data
# NBDC training and test data
load("./data/6_Oct_download/NBDC_training_data.Rdata")
load("./data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

# combine test and training data for plants
ptest <- test_data$`flowering plant`
ptrain <- training_data$`flowering plant`
plant <- bind_rows(ptest, ptrain)
rm(ptest, ptrain, test_data, training_data)

# BBS bryophyte data
stop("TODO: look over clean_bryophyte.R for millipedes 28 Feb2020.")
source("clean_bryophyte.R")


# make 10km, 1km, and 100m square template rasters
irish_hec_raster <- raster(xmn = -60000, xmx = 450000, ymn = -70000, ymx = 550000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_hec_raster) <- 10000
irish_hec_shp <- st_as_sf(rasterToPolygons(irish_hec_raster), 
                          crs = CRS(irish_hec_raster))

irish_1km_raster <- raster(xmn = 10000, xmx = 380000, ymn = -30000, ymx = 500000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_1km_raster) <- 1000
df_1km <- data.frame(rasterToPoints(irish_1km_raster, spatial = F))
colnames(df_1km) <- c("eastings", "northings")

## Choose indicator species



## find whether there are indicator species in each grid cell