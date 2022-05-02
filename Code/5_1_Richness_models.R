# set-up =======================================================================

# Load R package
library(spacetime)
library(INLA)
library(here)
library(raster)

#__________
# Load data
#__________
sp <- readRDS(here("./Data/number.RDS"))
sp_rich <- readRDS(here("./Data/richnessRaw.RDS"))
explan <- readRDS(here("./Data/explan.RDS"))

# Number of species
nsp_rich <- 103

# make just the year to put into the dataframes for both the sp and explan one
year = as.POSIXlt(sp_rich@endTime)$year + 1900 # year value
uniqueYears = unique(year)

# insert year 
sp_rich@data$year = year
explan@data$year = year

# model fits ===================================================================


#-------
# Priors
#-------
# PC prior for standard deviation

#------------------
# Energy hypothesis
#------------------
# Result object
energyWAIC <- matrix(NA, nrow = length(unique(uniqueYear)), ncol = 1)
rownames(energyWAIC) <- unique(uniqueYear)
colnames(energyWAIC) <- c("poisson")

waic_row = 1
for(i in uniqueYears){
  
  # Filter to just the year at hand
  sp_rich_temp = sp_rich@data[which(sp_rich@data[,"year"] == i),]
  explan_temp = explan@data[which(explan@data[,"year"] == i),]
  
  # Build data
  dat <- data.frame(richness = sp_rich_temp[,1],
                    bottom.temperature = scale(explan_temp$bottom.temperature),
                    surface.temperature = 
                      scale(explan_temp$surface.temperature))
  
  # Model formula
  Formula <- richness ~ bottom.temperature + surface.temperature

  # Estimate model
  model <- inla(Formula,
                    data = dat,
                    family = "poisson",
                    control.family = list(link='log'),
                    control.compute = list(waic = TRUE))
    
  # Save results
  energyWAIC[waic_row,1] <- model$waic$waic
  
  waic_row = waic_row + 1
}

write.csv(energyWAIC, 
          file = here("./ModelOutputs/occ_energyWAIC.csv"))

#------------------------
# Productivity hypothesis
#------------------------
# Result object
productivityWAIC <- matrix(NA, nrow = length(unique(uniqueYear)), ncol = 1)
rownames(productivityWAIC) <- unique(uniqueYear)
colnames(productivityWAIC) <- c("poisson")

waic_row = 1
for(i in uniqueYears){
  
  # Filter to just the year at hand
  sp_rich_temp = sp_rich@data[which(sp_rich@data[,"year"] == i),]
  explan_temp = explan@data[which(explan@data[,"year"] == i),]
  
  # Build data
  dat <- data.frame(richness = sp_rich_temp[,1],
                    bottom.salinity = scale(explan_temp$bottom.salinity),
                    ice = scale(explan_temp$ice))
  
  # Define formula
  Formula <- richness ~ bottom.salinity + ice

  # Estimate model
  model <- inla(Formula,
                data = dat,
                family = "poisson",
                control.family = list(link='log'),
                control.compute = list(waic = TRUE))
  
  # Save results
  productivityWAIC[waic_row,1] <- model$waic$waic
  
  waic_row = waic_row + 1
}

write.csv(productivityWAIC, 
          file = here("./ModelOutputs/occ_productivityWAIC.csv"))

#-----------------------------
# Climate stability hypothesis
#-----------------------------
# Result object
climateWAIC <- matrix(NA, nrow = length(unique(uniqueYear)), ncol = 1)
rownames(climateWAIC) <- unique(uniqueYear)
colnames(climateWAIC) <- c("poisson")

waic_row = 1
for(i in uniqueYears){
  
  # Filter to just the year at hand
  sp_rich_temp = sp_rich@data[which(sp_rich@data[,"year"] == i),]
  explan_temp = explan@data[which(explan@data[,"year"] == i),]
  
  # Build data
  dat <- data.frame(richness = sp_rich_temp[,1],
                    bottom.temperature.5.var = 
                      scale(explan_temp$bottom.temperature.5),
                    surface.temperature.5.var = 
                      scale(explan_temp$surface.temperature.5.var))
  
  # Define formula
  Formula <- richness ~ surface.temperature.5.var + bottom.temperature.5.var
  
  # Estimate model
  model <- inla(Formula,
                data = dat,
                family = "poisson",
                control.family = list(link='log'),
                control.compute = list(waic = TRUE))
    
  # Save results
  climateWAIC[waic_row,1] <- model$waic$waic
  
  waic_row = waic_row + 1
}

write.csv(climateWAIC, file = here("./ModelOutputs/occ_climateWAIC.csv"))

#---------------------------------
# habitat heterogeneity hypothesis
#---------------------------------
# Result object
habitatWAIC <- matrix(NA, nrow = length(unique(uniqueYear)), ncol = 1)
rownames(habitatWAIC) <- unique(uniqueYear)
colnames(habitatWAIC) <- c("poisson")

waic_row = 1
for(i in uniqueYears){
  
  # Filter to just the year at hand
  sp_rich_temp = sp_rich@data[which(sp_rich@data[,"year"] == i),]
  explan_temp = explan@data[which(explan@data[,"year"] == i),]
  
  # Build data
  dat <- data.frame(richness = sp_rich_temp[,1],
                    slope = scale(explan_temp$slope),
                    coarsness = scale(explan_temp$coarsness),
                    bottom.salinity = scale(explan_temp$bottom.salinity))
  
  # Define formula
  Formula <- richness ~ slope + coarsness + bottom.salinity

  # Estimate model
  model <- inla(Formula,
                data = dat,
                family = "poisson",
                control.family = list(link='log'),
                control.compute = list(waic = TRUE))
  
  # Save results
  habitatWAIC[waic_row,1] <- model$waic$waic
  
  waic_row = waic_row + 1
}

write.csv(habitatWAIC, file = here("./ModelOutputs/occ_habitatWAIC.csv"))

#--------------------------------
# stress heterogeneity hypothesis
#--------------------------------
# Result object
stressWAIC <- matrix(NA, nrow = length(unique(uniqueYear)), ncol = 1)
rownames(stressWAIC) <- unique(uniqueYear)
colnames(stressWAIC) <- c("poisson")

waic_row = 1
for(i in uniqueYears){
  
  # Filter to just the year at hand
  sp_rich_temp = sp_rich@data[which(sp_rich@data[,"year"] == i),]
  explan_temp = explan@data[which(explan@data[,"year"] == i),]
  
  # Build data
  dat <- data.frame(richness = sp_rich_temp[,1],
                    depth = scale(explan_temp$depth),
                    seal = scale(explan_temp$seal),
                    surfaceBottomTempDiff = 
                      scale(explan_temp$surfaceBottomTempDiff))
  # Define formula
  Formula <- richness ~ depth + seal + surfaceBottomTempDiff
  
  # Estimate model
  model <- inla(Formula,
                data = dat,
                family = "poisson",
                control.family = list(link='log'),
                control.compute = list(waic = TRUE))
    
  # Save results
  stressWAIC[waic_row,1] <- model$waic$waic
  
  waic_row = waic_row + 1
}

write.csv(stressWAIC, file = here("./ModelOutputs/occ_stressWAIC.csv"))

#
energyWAIC <- read.csv(here("./ModelOutputs/occ_energyWAIC.csv"), 
                            row.names = 1)
productivityWAIC <- read.csv(here("./ModelOutputs/occ_productivityWAIC.csv"), 
                                  row.names = 1)
climateWAIC<- read.csv(here("./ModelOutputs/occ_climateWAIC.csv"), 
                            row.names = 1)
habitatWAIC <- read.csv(here("./ModelOutputs/occ_habitatWAIC.csv"), 
                        row.names = 1)
stressWAIC <- read.csv(here("./ModelOutputs/occ_stressWAIC.csv"), 
                       row.names = 1)

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
