#################################################
# Read and organize DFO data in a standard format
# 
# explan : for samples 
#################################################
# Load R package
library(sf)
library(spacetime)
library(sp)
library(rgdal)
library(raster)
library(gslea)
library(mgcv)
library(spatialEco)
library(cleangeo)
library(rgeos)

#==========
# Read data 
#==========
#-----------------------
# Species code and names
#-----------------------
spRef <- read.csv("species_names.csv", row.names = 1)

#---------------------------------
# Raw data (site X year X Species) 
#---------------------------------
dfo <- read.csv("RVdata_adjusted_strata_ 415_439_1971_2020.csv")

#---------------
# Read shapefile
#---------------
strata <- readOGR("Gulf-Sept-RV-strata_poly")
gulf <- readRDS("Gulf_poly.RDS")

#--------------
# Load gulf DEM
#--------------
load(file = "gulf.dem.rda")
dem <- raster(gulf.dem)

# Read seal data (there is a lot of junk in this file)
# I will focus on the "seals" object
load("RVdatnoOgac.RData")
seals <- rbind(seals, data.frame(year = c(2019, 2020),
                                 sealsum = c(NA,NA)))

# Organize temporal data
dateTime <- as.POSIXlt(paste(paste(dfo$year,
                                   dfo$month,
                                   dfo$day, sep = "-"),
                             paste(dfo$start.hour,
                                   dfo$start.minute,
                                   dfo$start.second, sep = ":"),
                             sep = " "), tz = "Canada/Atlantic")

dfo <- cbind(time = dateTime, dfo)

# Organize spatial data
dfoST <- stConstruct(dfo,
                     space = c(which(colnames(dfo) == "longitude"),
                               which(colnames(dfo) == "latitude")), 
                     time = "time",
                     crs = CRS("EPSG:4326"))

## Transform projection
dfoST <- spTransform(dfoST, CRSobj = CRS("EPSG:3798"))
strata <- spTransform(strata, CRSobj = CRS("EPSG:3798"))
gulf <- spTransform(gulf, CRSobj = CRS("EPSG:3798"))
dem <- projectRaster(dem,crs = CRS("EPSG:3798"))

# Remove the non-sampled strata from the gulf polygon
regionSampled <- aggregate(strataClean[-(1:3),])

poly1 <- regionSampled@polygons[[1]]@Polygons[[1]]
poly2 <- regionSampled@polygons[[1]]@Polygons[[3]]
poly2@coords <- poly2@coords[-c(1,333),]
poly2@coords <- rbind(poly2@coords,poly2@coords[1,])

regionGulf <- SpatialPolygons(list(Polygons(list(poly1, poly2),ID=1)),
                              proj4string = regionSampled@proj4string)

# Extract and organize data to use 
#======================
# Explanatory variables
#======================
#-------------------------------------
# Organize basic explanatory variables
#-------------------------------------
explan <- dfoST
explan@data <- explan@data[,c("Cruise.code",
                              "set.number",
                              "stratum",
                              "duration",
                              "depth",
                              "distance",
                              "surface.temperature",
                              "bottom.temperature",
                              "bottom.salinity")]

# Replace 99.9 by NA
explan@data[,7:9][which(explan@data[,7:9] == 99.9, arr.ind = TRUE)] <- NA

#-------------------------------------
# Add additional explanatory variables
#-------------------------------------
#__________
# Seal data
#__________
# Impute Seal data
seals[49:50,2] <- predict.gam(gam(sealsum ~ s(year),
                                  data = seals),newdata = seals[49:50,])

# Weight of seal predation based on depth
depthWeight <- function(x, upThresh = -40) {
  weight <- exp(-0.5 * ((x - upThresh) / 25)^2)

  weight[x > upThresh & x < 0] <- 1
  weight[x > 0] <- 0
  
  return(weight)
}

# Weight of predation
demSeal <- dem
demSeal@data@values[demSeal@data@values < -200] <- NA
predWeightAll <- depthWeight(demSeal@data@values)
demSeal@data@values <- predWeightAll

demSeal@data@values[is.na(demSeal@data@values)] <- 0

predWeight <- extract(demSeal,dfoST@sp)

yearRep <- summary(as.factor(dateTime$year + 1900))
sealsRep <- numeric()
for(i in 1:nrow(seals)){
  sealsRep <- c(sealsRep, rep(seals[i,2], yearRep[i]))
}

explan@data$seal <- sealsRep * predWeight

#_________
# Ice data
#_________
# Load GSLEA map and convert projection
gsleaMap <- readOGR("gslea_shp")
gsleaMap <- spTransform(gsleaMap, CRSobj = CRS("EPSG:3798"))

gsleaMapSplit <- gDifference(gsleaMap[5,], gsleaMap[9,])

gsleaMapSplit@polygons[[1]]@ID <- "4"

gsleaSplit <- SpatialPolygonsDataFrame(gsleaMapSplit,
                                       data = gsleaMap[5,]@data)

gsleaList <- list(gsleaMap[1,],
                  gsleaMap[2,],
                  gsleaMap[3,],
                  gsleaMap[4,],
                  gsleaSplit,
                  gsleaMap[6,],
                  gsleaMap[7,],
                  gsleaMap[8,],
                  gsleaMap[9,])

gsleaSP <- SpatialPolygons(lapply(gsleaList, function(x){x@polygons[[1]]}),
                         proj4string = CRS("EPSG:3798"))

gslea <- SpatialPolygonsDataFrame(gsleaSP, data = gsleaMap@data)

#saveRDS(gslea, file = "./Data/gsleaClean.RDS")

#====================
# Species information
#====================
## Weight
weight <- dfoST
weight@data <- weight@data[,grep(".weight.", colnames(weight@data))]

### Extract code for weight
spWeightXCode <- unlist(strsplit(colnames(weight@data), split = ".weight.caught"))
spWeightCode <- as.numeric(sapply(strsplit(spWeightXCode, "X"), function(x) x[2]))

# Rename species names
colnames(weight@data) <- spRef[order(spRef[,1]),2]

## Number
number <- dfoST
number@data <- number@data[,grep(".number.", colnames(number@data))]

### Replace code by species name
spNumberXCode <- unlist(strsplit(colnames(number@data), split = ".number.caught"))
spNumberCode <- as.numeric(sapply(strsplit(spNumberXCode, "X"), function(x) x[2]))

# Rename species names
number@data <- round(number@data)
colnames(number@data) <- spRef[order(spRef[,1]),2]

## Fish richness
fishRich <- dfoST
fishRich@data <- data.frame(species.fish.number = fishRich@data[,"species.fish.number"])

## Invertebrate richness
invertRich <- dfoST
invertRich@data <- data.frame(species.invertebrate.number = invertRich@data[,"species.invertebrate.number"])

## Total weight caught
catchTotal <- dfoST
catchTotal@data <- data.frame(catch.total.weight = catchTotal@data[,"catch.total.weight"])

## Save data into .RDS object
saveRDS(explan, file = "./Data/explanRaw.RDS")
saveRDS(weight, file = "./Data/weight.RDS")
saveRDS(number, file = "./Data/number.RDS")

saveRDS(fishRich, file = "./Data/fishRich.RDS")
saveRDS(invertRich, file = "./Data/invertRich.RDS")
saveRDS(catchTotal, file = "./Data/catchTotal.RDS")

saveRDS(strata, file = "./Data/strata.RDS")
saveRDS(gulf, file = "./Data/gulf.RDS")
saveRDS(regionGulf, file = "./Data/gulfSampled.RDS")
saveRDS(dem, file = "./Data/dem.RDS")
