###########################
## calculate distance to nearest coast for each grid cell
## 
## Willson Gaul willson.gaul@ucdconnect.ie
## created: 12 March 2019
## last modified: 25 Oct 2019 - changed to millipedes map of ignorance project
###########################

print_plots <- F
### load data ---------------------------------------------------------------
# Load Ireland coastline
ir <- st_read(dsn='./data/', layer='ireland_coastline')
ir_TM75 <- st_transform(ir, crs = 29903)
rm(ir)

# load Irish hectad shapefile
hecs_shp <- st_read(dsn='./data/',
                    layer='IE_10km_hecs')
hecs_shp <- st_transform(hecs_shp, crs = 29903)

# make 10km square template raster
irish_hec_raster <- raster(xmn = -60000, xmx = 450000, ymn = -70000, ymx = 550000, 
                           crs = CRS("+init=epsg:29903"), vals = 1)
res(irish_hec_raster) <- 10000

# make 1km raster and shapefile templates 
ie_1km_raster <- disaggregate(irish_hec_raster, fact = 10)
ie_1km_shp <- st_as_sf(rasterToPolygons(ie_1km_raster)) %>% 
  st_transform(ie_1km_shp, crs = 29903)
### end load data -----------------------------------------------------------

### make 1km raster and shapefile templates -----------------------------------
ie_1km_raster <- disaggregate(irish_hec_raster, fact = 10)
ie_1km_shp <- st_as_sf(rasterToPolygons(ie_1km_raster)) %>% 
  st_transform(ie_1km_shp, crs = 29903)
## end make 1km raster and shapefile templates --------------------------------



### calculate distance from hectads to coast ----------------------------------
# make points at the centroid of each hectad
hec_points <- st_centroid(hecs_shp)

# make a linestring object of the coast polygons
ir_lines <- st_cast(ir_TM75, to = "MULTILINESTRING")

# calculate distance between each hectad centroid and the nearest part of each
# coast polygon.  Note that there are 18 coast polygons, some of them are small
# islands, etc.  
dist <- st_distance(hec_points, ir_lines) 
# find the minimum distance in each row, which will be the distance from the
# hectad centroid to the nearest coastline.
dist <- apply(dist, MARGIN = 1, FUN = min) 

# add distance to nearest coast to hec_points
hec_points$dist_to_coast <- dist

if(print_plots) {
  ggplot() + 
    geom_sf(data = ir_TM75) + 
    geom_sf(data = hec_points, aes(color = dist_to_coast), size = 5)
}
# rasterize hectad distances from nearest coast
coast_dist_hec_rast <- raster::rasterize(hec_points, irish_hec_raster, 
                                         field = "dist_to_coast")
### end calculate distance from hectads to coast ------------------------------


### calculate distance from 1km cells to coast --------------------------------
# make points at the centroid of each 1km grid cell
points_1km <- st_centroid(ie_1km_shp)

# calculate distance between each 1km cell centroid and the nearest part of each
# coast polygon.  Note that there are 18 coast polygons, some of them are small
# islands, etc.  
dist_1k <- st_distance(points_1km, ir_lines) 
# find the minimum distance in each row, which will be the distance from the
# 1km cell centroid to the nearest coastline.
dist_1k <- apply(dist_1k, MARGIN = 1, FUN = min) 

# add distance to nearest coast to points_1km
points_1km$dist_to_coast <- dist_1k

# rasterize hectad distances from nearest coast
coast_dist_1km_rast <- raster::rasterize(points_1km, ie_1km_raster, 
                                         field = "dist_to_coast")
### end calculate distance from 1km cells to coast ----------------------------


## save outputs ---------------------------------------------------------------
coast_dist_hec_shp <- hec_points # rename shapefile descriptively
coast_dist_1km_shp <- points_1km

save(coast_dist_hec_shp, coast_dist_hec_rast, file = "coast_dist_hec.RData")
save(coast_dist_1km_shp, coast_dist_1km_rast, file = "coast_dist_1km.RData")
## end save outputs -----------------------------------------------------------
