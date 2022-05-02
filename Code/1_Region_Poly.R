############################
# Read and organize DFO data
############################
# Load R package
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(mapdata) 
library(maptools)

## Read shapefile
gulf <- readRDS("Gulf_poly.RDS")

# Extract map of Canada
canada <- map("world", "Canada", fill = TRUE,
              col = "transparent", plot = FALSE)
IDs <- sapply(strsplit(canada$names, ":"), function(x) x[1])
canada <- map2SpatialPolygons(canada,
                                IDs = IDs,
                                proj4string = gulf@proj4string)

# Build spatialPolygon for region
gulfRegion <- SpatialPolygons(list(Polygons(list(gulf@polygons[[1]]@Polygons[[1]]),
                                            ID=1)), proj4string = gulf@proj4string)

# Build a 0.5 degree buffer around the region
Buffer <- gBuffer(gulfRegion, width = 0.25)

# Find the region of the buffer that is in water
water <- gDifference(Buffer, canada)

# Convert projection
bufferProj <- spTransform(Buffer, CRS("EPSG:3798"))
canadaProj <- spTransform(canada, CRS("EPSG:3798"))
waterProj <- spTransform(water, CRS("EPSG:3798"))

plot(waterProj)

# Save region polygons
saveRDS(bufferProj, file = "./Data/regionBuffer.RDS")
saveRDS(canadaProj, file = "./Data/Canada.RDS")
saveRDS(waterProj, file = "./Data/waterRegion.RDS")
