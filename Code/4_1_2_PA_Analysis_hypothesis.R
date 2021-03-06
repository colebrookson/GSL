# Load R package
library(spacetime)
library(INLA)

#__________
# Load data
#__________
sp <- readRDS("./Data/number.RDS")
explan <- readRDS("./Data/explan.RDS")

# Transform to presence-absence data
sp@data <- as.data.frame(ifelse(sp@data>0,1,0))

# Focus only on species with 50 occurrences or more across full data
sp@data <- sp@data[,which(colSums(sp@data) >= 50)]

# Identify which species is a fish and an invertebrate
fishInvert <- c(rep("fish", 68), 
                rep("invertebrate", 8))

# Number of species
nsp <- length(fishInvert)

#-------
# Priors
#-------
# PC prior for standard deviation

#------------------
# Energy hypothesis
#------------------
# Result object
energyWAIC <- matrix(NA, nrow = nsp, ncol = 2)
rownames(energyWAIC) <- colnames(sp@data)
colnames(energyWAIC) <- c("binomial","zeroinflatedbinomial0")

for(i in i:nsp){
  # Build data
  dat <- data.frame(species = sp@data[,i],
                    bottom.temperature = scale(explan@data$bottom.temperature),
                    surface.temperature = scale(explan@data$surface.temperature))
  
  # Define formula
  if(fishInvert[i] == "fish"){
    Formula <- species ~ bottom.temperature + surface.temperature
  }else{
    Formula <- species ~ bottom.temperature
  }
  
  # Estimate model
  for(j in 1:2){
    try(model <- inla(Formula,
                      data = dat,
                      family = colnames(energyWAIC)[j],
                      control.family = list(control.link = list(model = "logit")),
                      control.compute = list(waic = TRUE)))
    
    # Save results
    if(length(model) == 1){
      energyWAIC[i,j] <- NA
    }else{
      energyWAIC[i,j] <- model$waic$waic
    }
    model <- NA
  }
  print(i)
}

write.csv(energyWAIC, file = "occ_energyWAIC.csv")

#------------------------
# Productivity hypothesis
#------------------------
# Result object
productivityWAIC <- matrix(NA, nrow = nsp, ncol = 2)
rownames(productivityWAIC) <- colnames(sp@data)
colnames(productivityWAIC) <- c("binomial", "zeroinflatedbinomial0")

for(i in 1:ncol(sp@data)){
  # Build data
  dat <- data.frame(species = sp@data[,i],
                    bottom.salinity = scale(explan@data$bottom.salinity),
                    ice = scale(explan@data$ice))
  
  # Define formula
  Formula <- species ~ bottom.salinity + ice

  # Estimate model
  for(j in 1:2){
    try(model <- inla(Formula,
                  data = dat,
                  family = colnames(productivityWAIC)[j],
                  control.family = list(control.link = list(model = "logit")),
                  control.compute = list(waic = TRUE)))
    
    # Save results
    if(length(model) == 1){
      productivityWAIC[i,j] <- NA
    }else{
      productivityWAIC[i,j] <- model$waic$waic
    }
    model <- NA
  }
  print(i)
}

write.csv(productivityWAIC, file = "occ_productivityWAIC.csv")

#-----------------------------
# Climate stability hypothesis
#-----------------------------
# Result object
climateWAIC <- matrix(NA, nrow = nsp, ncol = 2)
rownames(climateWAIC) <- colnames(sp@data)
colnames(climateWAIC) <- c("binomial", "zeroinflatedbinomial0")

for(i in 1:ncol(sp@data)){
  # Build data
  dat <- data.frame(species = sp@data[,i],
                    bottom.temperature.5.var = scale(explan@data$bottom.temperature.5),
                    surface.temperature.5.var = scale(explan@data$surface.temperature.5.var))
  
  # Define formula
  if(fishInvert[i] == "fish"){
    Formula <- species ~ surface.temperature.5.var
  }else{
    Formula <- species ~ bottom.temperature.5.var
  }
  
  # Estimate model
  for(j in 1:2){
    try(model <- inla(Formula,
                      data = dat,
                      family = colnames(climateWAIC)[j],
                      control.family = list(control.link = list(model = "logit")),
                      control.compute = list(waic = TRUE)))
    
    # Save results
    if(length(model) == 1){
      climateWAIC[i,j] <- NA
    }else{
      climateWAIC[i,j] <- model$waic$waic
    }
    model <- NA
  }
  print(i)
}

write.csv(climateWAIC, file = "occ_climateWAIC.csv")

#---------------------------------
# habitat heterogeneity hypothesis
#---------------------------------
# Result object
habitatWAIC <- matrix(NA, nrow = nsp, ncol = 2)
rownames(habitatWAIC) <- colnames(sp@data)
colnames(habitatWAIC) <- c("binomial", "zeroinflatedbinomial0")

for(i in 1:ncol(sp@data)){
  # Build data
  dat <- data.frame(species = sp@data[,i],
                    slope = scale(explan@data$slope),
                    coarsness = scale(explan@data$coarsness),
                    bottom.salinity = scale(explan@data$bottom.salinity))
  
  # Define formula
  Formula <- species ~ slope + coarsness + bottom.salinity

  # Estimate model
  for(j in 1:2){
    try(model <- inla(Formula,
                      data = dat,
                      family = colnames(habitatWAIC)[j],
                      control.family = list(control.link = list(model = "logit")),
                      control.compute = list(waic = TRUE)))
    
    # Save results
    if(length(model) == 1){
      habitatWAIC[i,j] <- NA
    }else{
      habitatWAIC[i,j] <- model$waic$waic
    }
    model <- NA
  }
  print(i)
}

write.csv(habitatWAIC, file = "occ_habitatWAIC.csv")

#--------------------------------
# stress heterogeneity hypothesis
#--------------------------------
# Result object
stressWAIC <- matrix(NA, nrow = nsp, ncol = 2)
rownames(stressWAIC) <- colnames(sp@data)
colnames(stressWAIC) <- c("binomial", "zeroinflatedbinomial0")

for(i in 1:ncol(sp@data)){
  # Build data
  dat <- data.frame(species = sp@data[,i],
                    depth = scale(explan@data$depth),
                    seal = scale(explan@data$seal),
                    surfaceBottomTempDiff = scale(explan@data$surfaceBottomTempDiff))
  
  # Define formula
  Formula <- species ~ depth + seal + surfaceBottomTempDiff
  
  # Estimate model
  for(j in 1:2){
    try(model <- inla(Formula,
                      data = dat,
                      family = colnames(stressWAIC)[j],
                      control.family = list(control.link = list(model = "logit")),
                      control.compute = list(waic = TRUE)))
    
    # Save results
    if(length(model) == 1){
      stressWAIC[i,j] <- NA
    }else{
      stressWAIC[i,j] <- model$waic$waic
    }
    model <- NA
  }
  print(i)
}

write.csv(stressWAIC, file = "occ_stressWAIC.csv")

energyWAIC <- read.csv("./Results/occ_energyWAIC.csv", row.names = 1)
productivityWAIC <- read.csv("./Results/occ_productivityWAIC.csv", row.names = 1)
climateWAIC<- read.csv("./Results/occ_climateWAIC.csv", row.names = 1)
habitatWAIC <- read.csv("./Results/occ_habitatWAIC.csv", row.names = 1)
stressWAIC <- read.csv("./Results/occ_stressWAIC.csv", row.names = 1)

# Check best family 
sort((energyWAIC[,1] - energyWAIC[,2])[which((energyWAIC[,1] - energyWAIC[,2]) > 0)])
sort((productivityWAIC[,1] - productivityWAIC[,2])[which((productivityWAIC[,1] - productivityWAIC[,2]) > 0)])
sort((climateWAIC[,1] - climateWAIC[,2])[which((climateWAIC[,1] - climateWAIC[,2]) > 0)])
sort((habitatWAIC[,1] - habitatWAIC[,2])[which((habitatWAIC[,1] - habitatWAIC[,2]) > 0)])
sort((stressWAIC[,1] - stressWAIC[,2])[which((stressWAIC[,1] - stressWAIC[,2]) > 0)])

# Check non-convergence problem (no in binomial so I will focus on binomial) 
any(is.infinite(energyWAIC[,1]))
any(is.infinite(productivityWAIC[,1]))
any(is.infinite(climateWAIC[,1]))
any(is.infinite(habitatWAIC[,1]))
any(is.infinite(stressWAIC[,1]))

#----------------------
# Check best hypothesis
#----------------------
hypoBinomial <- cbind(energyWAIC[,1],
                      productivityWAIC[,1],
                      climateWAIC[,1],
                      habitatWAIC[,1],
                      stressWAIC[,1])

colnames(hypoBinomial) <- c("energy", 
                            "productivity", 
                            "climate", 
                            "habitat", 
                            "stress")

# Find best hypothesis
bestHypoSp <- apply(hypoBinomial, 1, which.min)

# Find WAIC difference between best hypothesis and other hypothesis 
hypoBinomialDiff <- matrix(NA,nrow = nrow(hypoBinomial),
                            ncol = ncol(hypoBinomial))

colnames(hypoBinomialDiff) <- colnames(hypoBinomial)
rownames(hypoBinomialDiff) <- rownames(hypoBinomial)

for(i in 1:nrow(hypoBinomialDiff)){
  hypoBinomialDiff[i,] <- hypoBinomial[i,]-hypoBinomial[i,bestHypoSp[i]]
}

# Organize and save results
abund <- colSums(sp@data)
occ <- colSums(sp@data>0)
WAICres <-  data.frame(round(hypoBinomialDiff,3), fishInver = fishInvert, abundance = abund, occurrence = occ)

write.csv(WAICres, file = "./Results/occ_WAICDiff_Hypothese.csv")
