# Load R package
library(INLA)
library(raster)
library(animation)
library(rgeos)

#__________
# Load data
#__________
# Data
sp <- readRDS("./Data/number.RDS")

# Transform to presence-absence data
sp@data <- as.data.frame(ifelse(sp@data>0,1,0))

# Focus only on species with 50 occurrences or more across full data
sp@data <- sp@data[,which(colSums(sp@data) >= 50)]
spName <- colnames(sp@data)
spName[64] <- "Lycodes"
spName[74] <- "Atlantic lyre crab_Greater toad crab"

#___________
# Load model 
#___________
mapGulfStackMeanList <- vector("list", length = length(spName))
mapGulfStack.025List <- vector("list", length = length(spName))
mapGulfStack.975List <- vector("list", length = length(spName))

names(mapGulfStackMeanList) <- spName
names(mapGulfStack.025List) <- spName
names(mapGulfStack.975List) <- spName

for(i in 75:length(spName)){
  type <- "occ" # "abund" or "occ"
  
  model <- readRDS(paste("occ -",spName[i],"- model.RDS"))
  ID <- readRDS(paste(spName[i],"- ID.RDS")) # Important for projection
  
  #__________________
  # Load mapping info
  #__________________
  gulf <- readRDS("./Data/gulfSampled.RDS")
  canada <- readRDS("./Data/Canada.RDS")
  
  # Load data (to draw strata)
  explan <- readRDS("./Data/explan.RDS")
  
  # Independent year
  year <- 1971:2020
  nyear <- length(year)
  
  #_____________________________
  # Load mesh and related object
  #_____________________________
  meshSpace <- readRDS("./Mesh/meshSpace.RDS")
  
  ############
  # Plot model
  ############
  # Build basis of the map
  stepsize <- 1000
  rangeX <- range(meshSpace$loc[,1])
  rangeY <- range(meshSpace$loc[,2])
  nxy <- round(c(diff(rangeX), 
                 diff(rangeY)) / stepsize)
  
  # Define basis of the map (warning is OK)
  mapBasis <- inla.mesh.projector(meshSpace,
                                  xlim = rangeX,
                                  ylim = rangeY,
                                  dim = nxy,
                                  crs = meshSpace$crs)
  
  # Result objects
  mapGulfMean <- vector("list", length = nyear)
  mapGulf.025 <- vector("list", length = nyear)
  mapGulf.975 <- vector("list", length = nyear)
  
  for(j in 1:nyear){
    if(!is.null(model[[j]])){
      # Model prediction
      fitMean <- inla.mesh.project(mapBasis, 
                                   model[[j]]$summary.fitted.values$mean[ID[,j]])
      fit.025 <- inla.mesh.project(mapBasis, 
                                   model[[j]]$summary.fitted.values$`0.025quant`[ID[,j]])
      fit.975 <- inla.mesh.project(mapBasis, 
                                   model[[j]]$summary.fitted.values$`0.975quant`[ID[,j]])
    }else{
      # Model prediction
      fitMean <- inla.mesh.project(mapBasis, 
                                   rep(0, meshSpace$n))
      fit.025 <- inla.mesh.project(mapBasis, 
                                   rep(0, meshSpace$n))
      fit.975 <- inla.mesh.project(mapBasis, 
                                   rep(0, meshSpace$n))
    }  
    # Build maps with confidence intervals
    mapMean <- raster(t(fitMean[,ncol(fitMean):1]),
                      xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                      ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                      crs = crs(gulf))
    
    map.025 <- raster(t(fit.025[,ncol(fit.025):1]),
                      xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                      ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                      crs = crs(gulf))
    
    map.975 <- raster(t(fit.975[,ncol(fit.975):1]),
                      xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                      ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                      crs = crs(gulf))
    
    # Crop and mask
    mapGulfMean[[j]] <- mask(crop(mapMean, extent(gulf)), gulf)
    mapGulf.025[[j]] <- mask(crop(map.025, extent(gulf)), gulf)
    mapGulf.975[[j]] <- mask(crop(map.975, extent(gulf)), gulf)
  }
  
  # Make stack
  mapGulfStackMean <- stack(mapGulfMean)
  mapGulfStack.025 <- stack(mapGulf.025)
  mapGulfStack.975 <- stack(mapGulf.975)
  
  # Change layers name
  names(mapGulfStackMean) <- paste0("y", year)
  names(mapGulfStack.025) <- paste0("y", year)
  names(mapGulfStack.975) <- paste0("y", year)
  
  mapGulfStackMeanList[[i]] <- mapGulfStackMean
  mapGulfStack.025List[[i]] <- mapGulfStack.025
  mapGulfStack.975List[[i]] <- mapGulfStack.975
}

saveRDS(mapGulfStackMeanList, file = "mapGulfStackMeanList.RDS")
saveRDS(mapGulfStack.025List, file = "mapGulfStack.025List.RDS")
saveRDS(mapGulfStack.975List, file = "mapGulfStack.975List.RDS")
