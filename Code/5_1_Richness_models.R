# set-up =======================================================================

# Load R package
library(spacetime)
library(INLA)
library(here)
library(raster)
library(tidyverse)

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

# results ======================================================================

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

# Check non-convergence problem (no in binomial so I will focus on binomial) 
ifelse(
  any(
    any(is.infinite(energyWAIC[,1])),
    any(is.infinite(productivityWAIC[,1])),
    any(is.infinite(climateWAIC[,1])),
    any(is.infinite(habitatWAIC[,1])),
    any(is.infinite(stressWAIC[,1]))),
  stop("ERROR - Non-convergence"), 
  print("No non-convergence"))


#----------------------
# Check best hypothesis
#----------------------
poisson <- cbind(energyWAIC[,1],
                  productivityWAIC[,1],
                  climateWAIC[,1],
                  habitatWAIC[,1],
                  stressWAIC[,1])

colnames(poisson) <- c("energy", 
                            "productivity", 
                            "climate", 
                            "habitat", 
                            "stress")

# Find best hypothesis
bestHypopoi <- apply(poisson, 1, which.min)

# Find WAIC difference between best hypothesis and other hypothesis 
poissonDiff <- matrix(NA,nrow = nrow(poisson),
                            ncol = ncol(poisson))

colnames(poissonDiff) <- colnames(poisson)
rownames(poissonDiff) <- rownames(poisson)

for(i in 1:nrow(poissonDiff)){
  poissonDiff[i,] <- poisson[i,]-poisson[i,bestHypopoi[i]]
}

# Organize and save results
write.csv(poissonDiff, file = here("./ModelOutputs/occ_WAICDiff_Hypothese.csv"))


waic_sig_table = matrix(NA, nrow = nrow(poisson),
                        ncol = ncol(poisson))
colnames(waic_sig_table) = colnames(poissonDiff)
rownames(waic_sig_table) = uniqueYears

for(i in 1:50) {
  for(j in 1:5) { 
    if(poissonDiff[i,j] == 0){ # best model
      waic_sig_table[i,j] = 1
    } else if((poissonDiff[i,j] > 0) & (poissonDiff[i,j] < 2)) { # close to best
      waic_sig_table[i,j] = 2
    } else if((poissonDiff[i,j] > 2) & (poissonDiff[i,j] < 10)) { # not too far
      waic_sig_table[i,j] = 3
    } else { # not at all good model
      waic_sig_table[i,j] = 4
    }
  }
}

waic_sig_table_df = data.frame(waic_sig_table)
waic_sig_table_df$year = rownames(waic_sig_table)

waic_sig_table_df_long = tidyr::pivot_longer(
  waic_sig_table_df,
  cols = c("energy", "productivity", "climate", "habitat", "stress"),
  names_to = "hypothesis",
  values_to = "waic")
waic_sig_table_df_long$waic = as.factor(waic_sig_table_df_long$waic)

waic_sig_table_df_long = waic_sig_table_df_long %>% 
  mutate(cat_waic = 
           case_when(
             waic == 1 ~ "dWAIC = 0.0",
             waic == 2 ~ "0.0 < dWAIC <= 2.0", 
             waic == 3 ~ "2.0 < dWAIC <= 10.0",
             waic == 4 ~ "dWAIC > 10.0")
  )

ggplot(waic_sig_table_df_long, 
       aes(x = hypothesis, y = year)) + 
  geom_tile(aes(fill = cat_waic), colour = "white") + 
  scale_fill_manual("wAIC levels", values=c("red", "blue", "black", "green"))









