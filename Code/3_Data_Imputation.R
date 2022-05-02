###########################################
# Read and organize DFO data
# 
# explan : for samples 
# explanMesh : for mesh points through time
###########################################
# Load R package
library(sf)
library(spacetime)
library(sp)
library(rgdal)
library(raster)
library(insol)
library(gslea)
library(gstat)
library(mgcv)
library(spatialEco)

#==========
# Read data 
#==========
# Non-imputed explanatory variables
explan <- readRDS("./Data/explanRaw.RDS")

# Gulf DEM
dem <- readRDS("./Data/dem.RDS")

# Gulf
gulf <- readRDS("./Data/gulf.RDS")

# Spatial mesh
meshSpace <- readRDS("./Mesh/meshSpace.RDS")

# GSLEA
gslea <- readRDS("./Data/gsleaClean.RDS")

#=================
# Build explanMesh
#=================
year <- unique(as.POSIXlt(explan@endTime)$year) + 1900
nyear <- length(unique(as.POSIXlt(explan@endTime)$year))

time <- as.POSIXlt(paste(rep(year, each = nrow(meshSpace$loc)),
                         "09-15", sep = "-"), tz = "Canada/Atlantic")

yearSite <- data.frame(year = time,
                       x = rep(meshSpace$loc[,1], nyear),
                       y = rep(meshSpace$loc[,2], nyear))

explanMesh <- stConstruct(yearSite, 
                          space = c("x", "y"),
                          time = "year",
                          crs = explan@sp@proj4string)

#---------------------
# Depth for explanMesh
#---------------------
explanMesh@data$depth <- -extract(dem, explanMesh@sp)

#--------------------
# Seal for explanMesh
#--------------------
load("RVdatnoOgac.RData")
seals <- rbind(seals, data.frame(year = c(2019, 2020),
                                 sealsum = c(NA,NA)))

# Impute Seal data
seals[49:50,2] <- predict.gam(gam(sealsum ~ s(year),
                                  data = seals),newdata = seals[49:50,])

depthWeight <- function(x, upThresh = -40) {
  weight <- exp(-0.5 * ((x - upThresh) / 25)^2)
  
  weight[x > upThresh & x < 0] <- 1
  weight[x > 0] <- 0
  
  return(weight)
}
dem <- readRDS(dem, file = "./Data/dem.RDS")

# Weight of predation
demSeal <- dem
demSeal@data@values[demSeal@data@values < -200] <- NA
predWeightAll <- depthWeight(demSeal@data@values)
demSeal@data@values <- predWeightAll

demSeal@data@values[is.na(demSeal@data@values)] <- 0
predWeight <- extract(demSeal,explanMesh@sp)

yearRep <- summary(as.factor(as.POSIXlt(explanMesh@endTime)$year + 1900))
sealsRep <- numeric()
for(i in 1:nrow(seals)){
  sealsRep <- c(sealsRep, rep(seals[i,2], yearRep[i]))
}

explanMesh@data$seal <- sealsRep * predWeight

#-------------------
# Ice for explanMesh
#-------------------
# Find EAR for each point
EAR <- point.in.poly(explanMesh@sp, gslea)

land <- which(is.na(EAR$PID))
sea <- which(!is.na(EAR$PID))

# find ice data
year <- as.POSIXlt(explanMesh@endTime)$year + 1900
iceEA <- EA.query.f(variables="ice.duration",
                    years=sort(unique(year)),
                    EARs=as.numeric(unique(EAR$PID)))

# Extract ice data and match to right year and EAR
ice <- numeric()
for(i in 1:nrow(explanMesh@data)){
  if(!is.na(EAR$PID[i])){
    ice[i] <- iceEA$value[which(iceEA$year == year[i] & iceEA$EAR == EAR$PID[i])]
  }else{
    ice[i] <- 0
  }
}

# Include in explanMesh 
explanMesh@data$ice <- ice

#---------------
# Ice for explan
#---------------
# Find EAR for each point
EAR <- point.in.poly(explan@sp, gslea)

# Replace NAs by 5
EAR$PID[which(is.na(EAR$PID))] <- "5"

# find ice data
year <- as.POSIXlt(explan@endTime)$year + 1900
iceEA <- EA.query.f(variables="ice.duration",
                    years=sort(unique(year)),
                    EARs=as.numeric(unique(EAR$PID)))

# Extract ice data and match to right year and EAR
ice <- numeric()
for(i in 1:nrow(explan@data)){
  ice[i] <- iceEA$value[which(iceEA$year == year[i] & iceEA$EAR == EAR$PID[i])]
}

# Include in explan 
explan@data <- data.frame(explan@data, ice = ice)

#--------------------------------------
# North Atlantic oscillation for explan
#--------------------------------------
NAO <- find.vars.f("nao.month")

# Across the months of the year
NAOYear <- EA.query.f(NAO, EARs = -1, years= c(1971:2020))
NAOYearMean <- tapply(NAOYear$value, as.factor(NAOYear$year), mean)

# In september
NAOYearSept <- EA.query.f(NAO[12], EARs = -1, years= c(1971:2020))$value

year <- as.POSIXlt(explan@endTime)$year + 1900

NOAMean <- numeric()
NOASept <- numeric()
for(i in 1:nrow(explan@data)){
  NOAMean[i] <- NAOYearMean[(year[i] - 1970)]
  NOASept[i] <- NAOYearSept[(year[i] - 1970)]
}

# Include in explan 
explan@data <- data.frame(explan@data, NOAYear = NOAMean, NOASept = NOASept)

#------------------------------------------
# North Atlantic oscillation for explanMesh
#------------------------------------------
year <- as.POSIXlt(explanMesh@endTime)$year + 1900

NOAMean <- numeric()
NOASept <- numeric()
for(i in 1:nrow(explanMesh@data)){
  NOAMean[i] <- NAOYearMean[(year[i] - 1970)]
  NOASept[i] <- NAOYearSept[(year[i] - 1970)]
}

# Include in explanMesh
explanMesh@data <- data.frame(explanMesh@data, NOAYear = NOAMean, NOASept = NOASept)

#======================
# Impute missing values 
#======================
#-------------------
# Bottom temperature 
#-------------------
bottomTempNA <- which(is.na(explan@data$bottom.temperature))

#_______________
# For explanMesh
#_______________
bottomTemp <-  explan@data$bottom.temperature[-bottomTempNA]
yearNoNA <- (as.POSIXlt(explan@endTime)$year + 1900)[-bottomTempNA]
locNoNA <- explan@sp[-bottomTempNA,]

bottomTempMesh <- numeric()
for(i in 1:nyear){
  yearSel <- which(yearNoNA == (i + 1970))
  res <- idw(bottomTemp[yearSel] ~ 1, 
             locations = locNoNA[yearSel,],
             newdata = SpatialPoints(coords = meshSpace$loc[,1:2],
                                     proj4string = explan@sp@proj4string))
  
  bottomTempMesh <- c(bottomTempMesh,res@data$var1.pred)
}

explanMesh@data$bottom.temperature <- bottomTempMesh

#___________
# For explan
#___________
bottomTempModel <- gam(bottom.temperature ~ s(depth), data = explan@data)

explan@data$bottom.temperature[bottomTempNA] <- predict(bottomTempModel,
                                                        newdata = data.frame(depth = explan@data$depth[bottomTempNA]))

#--------------------
# Surface temperature
#--------------------
surfaceTempNA <- which(is.na(explan@data$surface.temperature))
#_________________
# For explanMesh 
#
# For all but 2008
#_________________
surfaceTemp <-  explan@data$surface.temperature[-surfaceTempNA]
yearNoNA <- (as.POSIXlt(explan@endTime)$year + 1900)[-surfaceTempNA]
yearUniqueNA <- (1971:2020)[!(1971:2020 %in% unique(yearNoNA))]
locNoNA <- explan@sp[-surfaceTempNA,]

surfaceTempMesh <- numeric()

# Year 2008 (38th year) has not data
for(i in 1:nyear){
  if(any((i + 1970) == yearUniqueNA)){
    surfaceTempMesh <- c(surfaceTempMesh, rep(NA, nrow(meshSpace$loc[,1:2])))
  }else{
    yearSel <- which(yearNoNA == (i + 1970))
    res <- idw(surfaceTemp[yearSel] ~ 1, 
               locations = locNoNA[yearSel,],
               newdata = SpatialPoints(coords = meshSpace$loc[,1:2],
                                       proj4string = explan@sp@proj4string))
    
    surfaceTempMesh <- c(surfaceTempMesh,res@data$var1.pred)
  }
}

explanMesh@data$surface.temperature <- surfaceTempMesh

# For 2008 take average across the two years before and the two years after at
# the same location
yearMesh <- as.POSIXlt(explanMesh@endTime)$year + 1900
sel <- which(yearMesh >= 2006 & yearMesh <= 2010)

yearSel <- yearMesh[sel]
xySel <- coordinates(explanMesh@sp)[sel,]
surfaceTempSel <- surfaceTempMesh[sel]

ID <- rep(which(yearSel == unique(yearSel)[1]),5)

explanMesh@data$surface.temperature[is.na(explanMesh@data$surface.temperature)] <- tapply(surfaceTempSel, as.factor(ID), mean, na.rm=TRUE)

#___________
# For explan
#___________
timelt <- as.POSIXlt(explan@endTime)
year <- timelt$year + 1900

surfaceTempNA <- which(is.na(explan@data$surface.temperature))
yearNA <- year[surfaceTempNA]
yearUnique <- unique(yearNA)
yearNAsum <- summary(as.factor(yearNA), maxsum = 100)

# * i = 25 is year 2008 where no data is available for this year
for(i in (1:length(yearUnique))[-25]){
  # Extract values for a year with missing values
  xyYearNA <- coordinates(explan@sp[which(year == yearUnique[i])])
  datYearNA <- explan@data$surface.temperature[which(year == yearUnique[i])]
  
  # Remove the missing values
  NAtoRm <- which(is.na(datYearNA))
  xyYear <- xyYearNA[-NAtoRm,]
  datYear <- data.frame(surfTemp = datYearNA[-NAtoRm])
  
  # Build SPDF
  spDF <- SpatialPointsDataFrame(coords = xyYear,
                                 data = datYear,
                                 proj4string = crs(explan@sp))
  
  sp <- SpatialPoints(coords = matrix(xyYearNA[NAtoRm,], ncol = 2),
                      proj4string = crs(explan@sp))
  
  pred <- idw(surfTemp ~ 1, location = spDF, newdata = sp)
  
  # Save prediction
  explan@data$surface.temperature[which(year == yearUnique[i])][NAtoRm] <- pred$var1.pred
}

# Plot to check everything is OK
plot(explan@data$surface.temperature, type = "l")
points((1:nrow(explan))[surfaceTempNA],
       explan@data$surface.temperature[surfaceTempNA], 
       col = "red", pch = 19, cex = 0.2)

# Imputation by regression analysis (for 2018)
set.seed(42)
index <- 1:nrow(explan@data)
reg <- lm(explan@data$surface.temperature ~ index)
predReg <- predict.lm(reg,
                      newdata = data.frame(index = surfaceTempNA),
                      interval = "prediction", level = pnorm(1))
# The pnorm(1) value is to get 1 standard deviation. It will become useful for resampling

explan@data$surface.temperature[surfaceTempNA] <- rnorm(length(surfaceTempNA),
                                                        mean = predReg[,1],
                                                        sd = (predReg[,3] - predReg[,2])/3)

# Plot to check everything is OK
plot(index, explan@data$surface.temperature, type = "l", main = "surface temperature")
points(index[surfaceTempNA], explan@data$surface.temperature[surfaceTempNA],
       pch = 19, col = "red", cex = 0.2)

# Note : there are a few imputed samples for which the bottom temperature is 
# warmer than the surface temperature. After much checking, this is the results 
# of warm observed (not imputed) bottom temperature. So, it is OK.

#----------------
# Bottom salinity 
#----------------
# Replace extreme values by NAs
explan@data$bottom.salinity[1:6400][which(explan@data$bottom.salinity[1:6400] < 26)] <- NA
explan@data$bottom.salinity[which(explan@data$bottom.salinity > 36)] <- NA

#_______________
# For explanMesh
#_______________

bottomSalNA <- which(is.na(explan@data$bottom.salinity))
#___________________________________
# For explanMesh 
#
# For all but the years without data
#___________________________________
bottomSal <-  explan@data$bottom.salinity[-bottomSalNA]
yearNoNA <- (as.POSIXlt(explan@endTime)$year + 1900)[-bottomSalNA]
yearUniqueNA <- (1971:2020)[!(1971:2020 %in% unique(yearNoNA))]
locNoNA <- explan@sp[-bottomSalNA,]

bottomSalMesh <- numeric()

# Year 2008 (38th year) has not data
for(i in 1:nyear){
  if(any((i + 1970) == yearUniqueNA)){
    bottomSalMesh <- c(bottomSalMesh, rep(NA, nrow(meshSpace$loc[,1:2])))
  }else{
    yearSel <- which(yearNoNA == (i + 1970))
    res <- idw(bottomSal[yearSel] ~ 1, 
               locations = locNoNA[yearSel,],
               newdata = SpatialPoints(coords = meshSpace$loc[,1:2],
                                       proj4string = explan@sp@proj4string))
    
    bottomSalMesh <- c(bottomSalMesh,res@data$var1.pred)
  }
}

explanMesh@data$bottom.salinity <- bottomSalMesh

#____________________________________________________________________________
# For the years without data take average across the two years before and all 
# the years before 2018
#____________________________________________________________________________
yearMesh <- as.POSIXlt(explanMesh@endTime)$year + 1900
sel <- which(yearMesh >= 1971 & yearMesh <= 2018)

yearSel <- yearMesh[sel]
xySel <- coordinates(explanMesh@sp)[sel,]
bottomSalSel <- bottomSalMesh[sel]

ID <- rep(which(yearSel == unique(yearSel)[1]),length(1971:2018))

explanMesh@data$bottom.salinity[is.na(explanMesh@data$bottom.salinity)] <- tapply(bottomSalSel, as.factor(ID), mean, na.rm=TRUE)

#___________
# For explan
#___________

# Find all the values within 27 and 35 before 2019
ref <- which(explan@data$bottom.salinity[1:6472]>27 & explan@data$bottom.salinity[1:6472]<35)

# Find all the NAs before 2019
refNA <- which(is.na(explan@data$bottom.salinity[1:6472]))

# Replace NAs by values sampled from referenced values
set.seed(42)
explan@data$bottom.salinity[refNA] <- sample(explan@data$bottom.salinity[1:6472][ref],
                                             size = length(refNA), replace = TRUE)

# Find all the values within 27 and 35 after 2019
ref2019 <- which(!is.na(explan@data$bottom.salinity[6473:nrow(explan)]))

# Find all the NAs after 2019
refNA2019 <- which(is.na(explan@data$bottom.salinity[6473:nrow(explan)]))

# Replace NAs by values sampled from referenced values
set.seed(42)
explan@data$bottom.salinity[6473:nrow(explan)][refNA2019] <- sample(explan@data$bottom.salinity[6473:nrow(explan)][ref2019],
                                                                    size = length(refNA2019), replace = TRUE)

# Visual
plot(explan@data$bottom.salinity, type = "l")
points(refNA,explan@data$bottom.salinity[refNA],
       pch = 19, cex = 0.2, col = "red")

points(6472 + refNA2019,explan@data$bottom.salinity[6472 + refNA2019],
       pch = 19, cex = 0.2, col = "red")

# Although there is clearly some autocorrelation within a sampling year, it was 
# decided to sample randomly (so without accounting for autocorrelation) because
# autocorrelation was not always in the same direction across years.

#===============================
# Generate explanatory variables
#===============================
# To do when all the variables are imputed 
#------
# Slope
#------
Slope <- raster(slope(cgrad(dem), degrees = TRUE),
                xmn = xmin(dem),
                xmx = xmax(dem),
                ymn = ymin(dem),
                ymx = ymax(dem),
                crs = CRS("EPSG:3798"))

#___________
# explanMesh
#___________
SlopeSmp <- extract(Slope, explanMesh@sp)
explanMesh@data$slope <- SlopeSmp

#_______
# explan
#_______
SlopeSmp <- extract(Slope, explan@sp)
explan@data <- data.frame(explan@data, slope = SlopeSmp)

#-------------------------------------------
# Surface temperature variation over 5 years
#-------------------------------------------
#___________
# explanMesh
#___________
timelt <- as.POSIXlt(explanMesh@endTime)
year <- timelt$year + 1900
yearUni <- unique(year)
yearSel <- yearUni[-c(1:2,length(yearUni)-1,length(yearUni))]

sTemVar <- rep(NA, length = length(yearUni))
names(sTemVar) <- yearUni

for(i in 1:length(yearSel)){
  # Sample to use
  smplToUse <- which(year %in% (yearSel[i] - 2):(yearSel[i] + 2))
  
  # Save new variable
  sTemVar[which(names(sTemVar) == yearSel[i])] <- var(explanMesh@data$surface.temperature[smplToUse], na.rm = TRUE)
}

# First year
smplToUse <- which(year %in% yearUni[1]:(yearUni[1] + 2))
sTemVar[which(names(sTemVar) == yearUni[1])] <- var(explanMesh@data$surface.temperature[smplToUse], na.rm = TRUE)

# Second year
smplToUse <- which(year %in% (yearUni[1]):(yearUni[2] + 2))
sTemVar[which(names(sTemVar) == yearUni[2])] <- var(explanMesh@data$surface.temperature[smplToUse], na.rm = TRUE)

# Last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 2)])
sTemVar[which(names(sTemVar) == yearUni[length(yearUni)])] <- var(explanMesh@data$surface.temperature[smplToUse], na.rm = TRUE)

# Second to last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 3)])
sTemVar[which(names(sTemVar) == yearUni[length(yearUni)-1])] <- var(explanMesh@data$surface.temperature[smplToUse], na.rm = TRUE)

sTemVarSmp <- numeric()
nYear <- summary(as.factor(year), maxsum = 100)
countStart <- 1
countEnd <- nYear[1]
for(i in 1:length(nYear)){
  sTemVarSmp[countStart:countEnd] <- rep(sTemVar[i], each = nYear[i])
  countStart <- countEnd + 1
  countEnd <- countEnd + nYear[i + 1]
}

explanMesh@data$surface.temperature.5.var <- sTemVarSmp

#_______
# explan
#_______
timelt <- as.POSIXlt(explan@endTime)
year <- timelt$year + 1900
yearUni <- unique(year)
yearSel <- yearUni[-c(1:2,length(yearUni)-1,length(yearUni))]

sTemVar <- rep(NA, length = length(yearUni))
names(sTemVar) <- yearUni

for(i in 1:length(yearSel)){
  # Sample to use
  smplToUse <- which(year %in% (yearSel[i] - 2):(yearSel[i] + 2))
  
  # Save new variable
  sTemVar[which(names(sTemVar) == yearSel[i])] <- var(explan@data$surface.temperature[smplToUse], na.rm = TRUE)
}

# First year
smplToUse <- which(year %in% yearUni[1]:(yearUni[1] + 2))
sTemVar[which(names(sTemVar) == yearUni[1])] <- var(explan@data$surface.temperature[smplToUse], na.rm = TRUE)

# Second year
smplToUse <- which(year %in% (yearUni[1]):(yearUni[2] + 2))
sTemVar[which(names(sTemVar) == yearUni[2])] <- var(explan@data$surface.temperature[smplToUse], na.rm = TRUE)

# Last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 2)])
sTemVar[which(names(sTemVar) == yearUni[length(yearUni)])] <- var(explan@data$surface.temperature[smplToUse], na.rm = TRUE)

# Second to last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 3)])
sTemVar[which(names(sTemVar) == yearUni[length(yearUni)-1])] <- var(explan@data$surface.temperature[smplToUse], na.rm = TRUE)

sTemVarSmp <- numeric()
nYear <- summary(as.factor(year), maxsum = 100)
countStart <- 1
countEnd <- nYear[1]
for(i in 1:length(nYear)){
  sTemVarSmp[countStart:countEnd] <- rep(sTemVar[i], each = nYear[i])
  countStart <- countEnd + 1
  countEnd <- countEnd + nYear[i + 1]
}

explan@data <- data.frame(explan@data, surface.temperature.5.var = sTemVarSmp)

#------------------------------------------
# Bottom temperature variation over 5 years
#------------------------------------------
#___________
# explanMesh
#___________
timelt <- as.POSIXlt(explanMesh@endTime)
year <- timelt$year + 1900
yearUni <- unique(year)
yearSel <- yearUni[-c(1:2,length(yearUni)-1,length(yearUni))]

bTemVar <- rep(NA, length = length(yearUni))
names(bTemVar) <- yearUni

for(i in 1:length(yearSel)){
  # Sample to use
  smplToUse <- which(year %in% (yearSel[i] - 2):(yearSel[i] + 2))
  
  # Save new variable
  bTemVar[which(names(bTemVar) == yearSel[i])] <- var(explanMesh@data$bottom.temperature[smplToUse], na.rm = TRUE)
}

# First year
smplToUse <- which(year %in% yearUni[1]:(yearUni[1] + 2))
bTemVar[which(names(bTemVar) == yearUni[1])] <- var(explanMesh@data$bottom.temperature[smplToUse], na.rm = TRUE)

# Second year
smplToUse <- which(year %in% (yearUni[1]):(yearUni[2] + 2))
bTemVar[which(names(bTemVar) == yearUni[2])] <- var(explanMesh@data$bottom.temperature[smplToUse], na.rm = TRUE)

# Last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 2)])
bTemVar[which(names(bTemVar) == yearUni[length(yearUni)])] <- var(explanMesh@data$bottom.temperature[smplToUse], na.rm = TRUE)

# Second to last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 3)])
bTemVar[which(names(bTemVar) == yearUni[length(yearUni)-1])] <- var(explanMesh@data$bottom.temperature[smplToUse], na.rm = TRUE)

bTemVarSmp <- numeric()
nYear <- summary(as.factor(year), maxsum = 100)
countStart <- 1
countEnd <- nYear[1]
for(i in 1:length(nYear)){
  bTemVarSmp[countStart:countEnd] <- rep(bTemVar[i], each = nYear[i])
  countStart <- countEnd + 1
  countEnd <- countEnd + nYear[i + 1]
}

explanMesh@data$bottom.temperature.5.var <- bTemVarSmp

#_______
# explan
#_______
timelt <- as.POSIXlt(explan@endTime)
year <- timelt$year + 1900
yearUni <- unique(year)
yearSel <- yearUni[-c(1:2,length(yearUni)-1,length(yearUni))]

bTemVar <- rep(NA, length = length(yearUni))
names(bTemVar) <- yearUni

for(i in 1:length(yearSel)){
  # Sample to use
  smplToUse <- which(year %in% (yearSel[i] - 2):(yearSel[i] + 2))
  
  # Save new variable
  bTemVar[which(names(bTemVar) == yearSel[i])] <- var(explan@data$bottom.temperature[smplToUse], na.rm = TRUE)
}

# First year
smplToUse <- which(year %in% yearUni[1]:(yearUni[1] + 2))
bTemVar[which(names(bTemVar) == yearUni[1])] <- var(explan@data$bottom.temperature[smplToUse], na.rm = TRUE)

# Second year
smplToUse <- which(year %in% (yearUni[1]):(yearUni[2] + 2))
bTemVar[which(names(bTemVar) == yearUni[2])] <- var(explan@data$bottom.temperature[smplToUse], na.rm = TRUE)

# Last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 2)])
bTemVar[which(names(bTemVar) == yearUni[length(yearUni)])] <- var(explan@data$bottom.temperature[smplToUse], na.rm = TRUE)

# Second to last year
smplToUse <- which(year %in% (yearUni[length(yearUni)]):yearUni[(length(yearUni) - 3)])
bTemVar[which(names(bTemVar) == yearUni[length(yearUni)-1])] <- var(explan@data$bottom.temperature[smplToUse], na.rm = TRUE)

bTemVarSmp <- numeric()
nYear <- summary(as.factor(year), maxsum = 100)
countStart <- 1
countEnd <- nYear[1]
for(i in 1:length(nYear)){
  bTemVarSmp[countStart:countEnd] <- rep(bTemVar[i], each = nYear[i])
  countStart <- countEnd + 1
  countEnd <- countEnd + nYear[i + 1]
}

explan@data <- data.frame(explan@data, bottom.temperature.5.var = bTemVarSmp)

#-----------
# Coarseness
#-----------
coarsnessRaster <- aggregate(dem, fact = 2, fun = var)

#___________
# explanMesh
#___________
coarsnessSmp <- extract(coarsnessRaster, explanMesh@sp)
explanMesh@data$coarsness <- coarsnessSmp

#_______
# explan
#_______
coarsnessSmp <- extract(coarsnessRaster, explan@sp)
explan@data <- data.frame(explan@data, coarsness = coarsnessSmp)

#----------------------------------------
# Bottom - surface temperature difference
#----------------------------------------
#___________
# explanMesh
#___________
explanMesh@data$surfaceBottomTempDiff <- explanMesh@data$surface.temperature - explanMesh@data$bottom.temperature

#_______
# explan
#_______
explan@data <- data.frame(explan@data,
                          surfaceBottomTempDiff = explan@data$surface.temperature - explan@data$bottom.temperature)

#-------------------------------------------
# Sea surface temperature slope - 10 km grid
#-------------------------------------------
# Build reference raster
gulfRasterBase <- raster(extent(gulf),
                         crs = gulf@proj4string,
                         resolution = 1000)

nCell <- ncell(gulfRasterBase)

# Extract adjacent cells (queen) for later calculation
adjCell <- vector("list", length = nCell)

for(i in 1:nCell){
  adjCell[[i]] <- adjacent(gulfRasterBase, i, direction = 8, pairs = FALSE)
}

# Extract coordinates for interpolation / extrapolation
coord <- coordinates(gulfRasterBase)

# Extract surface temperature
surfaceTemp <- explan@data$surface.temperature

for(i in 1:nyear){
  yearSel <- which(year == (i + 1970))
  surfaceTempIDW <- idw(surfaceTemp[yearSel] ~ 1, 
                        locations = explan@sp[yearSel,],
                        newdata = SpatialPoints(coords = coord,
                                                proj4string = explan@sp@proj4string))
  
# Use slope function of the insol package
  values(gulfRasterBase) <- surfaceTempIDW$var1.pred
  res <- raster(slope(cgrad(gulfRasterBase), degrees = TRUE),
         xmn = xmin(gulfRasterBase),
         xmx = xmax(gulfRasterBase),
         ymn = ymin(gulfRasterBase),
         ymx = ymax(gulfRasterBase),
         crs = gulfRasterBase@crs)
  
}

Slope <- raster(slope(cgrad(dem), degrees = TRUE),
                xmn = xmin(dem),
                xmx = xmax(dem),
                ymn = ymin(dem),
                ymx = ymax(dem),
                crs = CRS("EPSG:3798"))

#============
# Save object
#============
saveRDS(explanMesh, file = "./Data/explanMesh.RDS")
saveRDS(explan, file = "./Data/explan.RDS")
