# Load R package
library(spacetime)
library(INLA)
library(raster)

#__________
# Load data
#__________
# Data
sp <- readRDS("./Data/number.RDS")
explan <- readRDS("./Data/explan.RDS")
explanMesh <- readRDS("./Data/explanMesh.RDS")

# Transform to presence-absence data
sp@data <- as.data.frame(ifelse(sp@data>0,1,0))

# Focus only on species with 50 occurrences or more across full data
sp@data <- sp@data[,which(colSums(sp@data) >= 50)]
spName <- colnames(sp@data)

# WAIC results from hypothesis
WAIC <- read.csv("./Results/occ_WAICDiff_Hypothese.csv", row.names = 1)
WAICprob <- round(exp(-WAIC[which(rownames(WAIC) %in% spName),1:5]/2), 5)

#__________________
# Load mapping info
#__________________
gulf <- readRDS("./Data/gulf.RDS")
canada <- readRDS("./Data/Canada.RDS")

#_____________________________
# Load mesh and related object
#_____________________________
meshSpace <- readRDS("./Mesh/meshSpace.RDS")
barrierTriangle <- readRDS("./Mesh/barrierTriangle.RDS")

#_________________________________________________
# Space only barrier model across years separately
#_________________________________________________
# Get species info to build model
year <- unique(as.POSIXlt(sp@endTime)$year)
nyear <- length(year)

# Species to model
for(i in 75:76){
  
  # result object
  model <- vector("list", length = nyear)
  
  # ID for prediction
  ID <- matrix(NA,nrow = meshSpace$n, ncol = nyear)
  
  for(j in 1:nyear){
    pointerYear <-  which(as.POSIXlt(sp@endTime)$year == year[j])
    locYear <- coordinates(sp@sp)[pointerYear,]
    spYear <- sp@data[pointerYear,i]
  
    if(sum(spYear)>0){
      
      pointerExplanMesh <- which(as.POSIXlt(explanMesh@endTime)$year == year[j])
      
      # Build stack for estimation
      AEst <- inla.spde.make.A(meshSpace, locYear)
      stackEst <- inla.stack(data = list(y = spYear),
                             A = list(AEst,1,1,1,1),
                             effect = list(s = 1:meshSpace$n,
                                           intercept = rep(1, length(pointerYear)),
                                           seal = scale(explan@data$seal)[pointerYear],
                                           depth = scale(explan@data$depth)[pointerYear],
                                           surfaceBottomTempDiff = scale(explan@data$surfaceBottomTempDiff)[pointerYear]),
                             tag = "est")
      
      # Build stack for prediction
      APred <-  inla.spde.make.A(meshSpace, meshSpace$loc[,1:2])
      stackPred <- inla.stack(data = list(y = NA),
                              A = list(APred,1,1,1,1),
                              effect = list(s = 1:meshSpace$n,
                                            intercept = rep(1, length(pointerExplanMesh)),
                                            seal = scale(explanMesh@data$seal)[pointerExplanMesh],
                                            depth = scale(explanMesh@data$depth)[pointerExplanMesh],
                                            surfaceBottomTempDiff = scale(explanMesh@data$surfaceBottomTempDiff)[pointerExplanMesh]),
                              tag = "pred")
      
      # Formula
      Formula <- y ~ 0 + intercept + seal + depth + surfaceBottomTempDiff + f(s, model = barrierModel)
      
      # Join stacks
      stack <- inla.stack(stackEst, stackPred)
      
      # ID in stack for prediction
      ID[,j] <- inla.stack.index(stack, tag="pred")$data
      
      # Build basis of barrier model
      barrierModel <- inla.barrier.pcmatern(meshSpace, 
                                            barrier.triangles = barrierTriangle,
                                            prior.range = c(100000, 0.5),
                                            prior.sigma = c(30, 0.5))
      
      # Fit model
      model[[j]] <- inla(Formula,
                         data = inla.stack.data(stack),
                         control.predictor = list(A = inla.stack.A(stack),
                                                  compute = TRUE,
                                                  link = 1),
                         family = "binomial",
                         control.inla = list(int.strategy = "eb"),
                         control.family = list(control.link = list(model = "logit")),
                         control.compute = list(waic = TRUE))
    }
    print(paste(1900 + year[j],"Done"))
  }
  
  saveRDS(model, file = paste0("occ - ",colnames(sp@data)[i]," - model.RDS"))
  saveRDS(ID, file = paste0(colnames(sp@data)[i]," - ID.RDS"))
}
