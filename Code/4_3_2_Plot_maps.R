# Load R package
library(INLA)
library(raster)
library(animation)
library(rgeos)

type <- "occ"
#__________________
# Load mapping info
#__________________
gulf <- readRDS("./Data/gulfSampled.RDS")
canada <- readRDS("./Data/Canada.RDS")

# Independent year
year <- 1971:2020
nyear <- length(year)

#_____________________________
# Load mesh and related object
#_____________________________
meshSpace <- readRDS("./Mesh/meshSpace.RDS")

#___________________
# Load projected map
#___________________
mapGulfStackMeanList <- readRDS("mapGulfStackMeanList.RDS")
mapGulfStack.025List <- readRDS("mapGulfStack.025List.RDS")
mapGulfStack.975List <- readRDS("mapGulfStack.975List.RDS")

#--------------------------
# With confidence intervals
#--------------------------
# Select species of interest
for(i in 1:length(mapGulfStackMeanList)){
  i = 73
  spName <- names(mapGulfStackMeanList)[i]
  mapGulfStackMean <- mapGulfStackMeanList[[i]]
  mapGulfStack.025 <- mapGulfStack.025List[[i]]
  mapGulfStack.975 <- mapGulfStack.975List[[i]]
  
  ############
  # Plot model
  ############
  # Extract most extreme values 
  val <- c(values(mapGulfStack.025),
           values(mapGulfStack.975))
  
  # Define class intervals
  if(type == "occ"){
    fixBreak <- seq(0,1, length = 100)
  }
  
  color <- colorRampPalette(c("grey90", "steelblue4", "steelblue2", 
                              "steelblue1", "gold", "red1", "red4"),
                            bias = 1)
  
  
  # Plot figure for model mean
  saveGIF(
    {
      for(j in 1:nyear){
        # Base of the figure
        par(mfrow = c(1,4), mar=c(1,1,1,1))
        #=====
        # 2.5%
        #=====
        # Base of plot
        plot(mapGulfStack.025@extent,
             xlim = c(mapGulfStack.025@extent@xmin,
                      mapGulfStack.025@extent@xmax),
             ylim = c(mapGulfStack.025@extent@ymin,
                      mapGulfStack.025@extent@ymax),
             type = "n", axes=FALSE,
             xlab = "", ylab = "", asp = 1)
        
        # Background
        rect(mapGulfStack.025@extent@xmin-25000,
             mapGulfStack.025@extent@ymin-25000,
             mapGulfStack.025@extent@xmax+25000,
             mapGulfStack.025@extent@ymax+25000,
             col = rgb(193/256,236/256,250/256),
             border = rgb(193/256,236/256,250/256))
        
        # Model
        plot(mapGulfStack.025[[j]], 
             col = color(length(fixBreak)),
             add = TRUE,
             breaks = fixBreak,
             legend = FALSE)
        
        plot(canada, add = TRUE, col ="darkkhaki")
        legend("topright", "2.5%", bty = "n", cex = 3, adj = c(0,0))
        
        #=====
        # Mean
        #=====
        # Base of plot
        plot(mapGulfStackMean@extent,
             xlim = c(mapGulfStackMean@extent@xmin,
                      mapGulfStackMean@extent@xmax),
             ylim = c(mapGulfStackMean@extent@ymin,
                      mapGulfStackMean@extent@ymax),
             type = "n", axes=FALSE,
             xlab = "", ylab = "", asp = 1)
        
        # Background
        rect(mapGulfStackMean@extent@xmin-25000,
             mapGulfStackMean@extent@ymin-25000,
             mapGulfStackMean@extent@xmax+25000,
             mapGulfStackMean@extent@ymax+25000,
             col = rgb(193/256,236/256,250/256),
             border = rgb(193/256,236/256,250/256))
        
        # Model
        plot(mapGulfStackMean[[j]], 
             col = color(length(fixBreak)),
             add = TRUE,
             breaks = fixBreak,
             legend = FALSE)
  
        plot(canada, add = TRUE, col ="darkkhaki")
        legend("topright", "Moyenne", bty = "n", cex = 3, adj = c(0,0))
        
        #=======
        # Legend
        #=======
        par(mar =c(5,28,10,28))
        
        image(matrix(1:length(fixBreak), nrow = 1), col = color(length(fixBreak)),
              xaxt = "n", yaxt = "n")
        
        if(type == "abund"){
          axis(4, at = seq(0.05,0.95,by=0.1),
               labels = c(1, 3, 10, 25, 60, 150, 400, 1000, 2500, 6500), 
               las = 1, cex.axis = 2.5)
        }
        
        if(type == "occ"){
          axis(4, at = seq(0.05,0.95,by=0.1),
               labels = seq(0.05,0.95,by=0.1), 
               las = 1, cex.axis = 2.5)
        }
        
        title(main = year[j], cex.main = 5)
        
        #======
        # 97.5%
        #======
        # Base of plot
        plot(mapGulfStack.975@extent,
             xlim = c(mapGulfStack.975@extent@xmin,
                      mapGulfStack.975@extent@xmax),
             ylim = c(mapGulfStack.975@extent@ymin,
                      mapGulfStack.975@extent@ymax),
             type = "n", axes=FALSE,
             xlab = "", ylab = "", asp = 1)
        
        # Background
        rect(mapGulfStack.975@extent@xmin-25000,
             mapGulfStack.975@extent@ymin-25000,
             mapGulfStack.975@extent@xmax+25000,
             mapGulfStack.975@extent@ymax+25000,
             col = rgb(193/256,236/256,250/256),
             border = rgb(193/256,236/256,250/256))
        # Model
        plot(mapGulfStack.975[[j]], 
             col = color(length(fixBreak)),
             add = TRUE,
             breaks = fixBreak,
             legend = FALSE)
  
        plot(canada, add = TRUE, col ="darkkhaki")
        legend("topright", "97.5%", bty = "n", cex = 3, adj = c(0,0))
        
      }
    }, 
    movie.name = paste0(spName,".gif"),
    ani.height = 500,
    ani.width = 2400
  )
}


#-----------------------------
# Without confidence intervals
#-----------------------------
# Load R package
library(INLA)
library(raster)
library(animation)
library(rgeos)

#__________________
# Load mapping info
#__________________
gulf <- readRDS("./Data/gulfSampled.RDS")
canada <- readRDS("./Data/Canada.RDS")

# Independent year
year <- 1971:2020
nyear <- length(year)

#_____________________________
# Load mesh and related object
#_____________________________
meshSpace <- readRDS("./Mesh/meshSpace.RDS")

#___________________
# Load projected map
#___________________
mapGulfStackMeanList <- readRDS("mapGulfStackMeanList.RDS")
mapGulfStack.025List <- readRDS("mapGulfStack.025List.RDS")
mapGulfStack.975List <- readRDS("mapGulfStack.975List.RDS")


# Select species of interest
for(i in 1:length(mapGulfStackMeanList)){
  
  spName <- names(mapGulfStackMeanList)[i]
  mapGulfStackMean <- mapGulfStackMeanList[[i]]

  ############
  # Plot model
  ############
  # Extract most extreme values 
  val <- values(mapGulfStackMean)
  
  # Define class intervals
  fixBreak <- seq(0,1, length = 100)

  color <- colorRampPalette(c("grey90", "steelblue4", "steelblue2", 
                              "steelblue1", "gold", "red1", "red4"),
                            bias = 1)
  
  
  # Plot figure for model mean
  saveGIF(
    {
      for(j in 1:nyear){
        # Base of the figure
        par(mfrow = c(1,2), mar=c(5,5,5,5))
        #=====
        # Mean
        #=====
        # Base of plot
        plot(mapGulfStackMean@extent,
             xlim = c(mapGulfStackMean@extent@xmin,
                      mapGulfStackMean@extent@xmax),
             ylim = c(mapGulfStackMean@extent@ymin,
                      mapGulfStackMean@extent@ymax),
             type = "n", axes=FALSE,
             xlab = "", ylab = "", asp = 1)
        
        # Background
        rect(mapGulfStackMean@extent@xmin-25000,
             mapGulfStackMean@extent@ymin-25000,
             mapGulfStackMean@extent@xmax+25000,
             mapGulfStackMean@extent@ymax+25000,
             col = rgb(193/256,236/256,250/256),
             border = rgb(193/256,236/256,250/256))
        
        # Model
        plot(mapGulfStackMean[[j]], 
             col = color(length(fixBreak)),
             add = TRUE,
             breaks = fixBreak,
             legend = FALSE)
        
        plot(canada, add = TRUE, col ="darkkhaki")
        
        # Legend
        par(mar =c(5,28,10,28))
        
        image(matrix(1:length(fixBreak), nrow = 1), col = color(length(fixBreak)),
              xaxt = "n", yaxt = "n")
        
        axis(4, at = seq(0.05,0.95,by=0.1),
             labels = seq(0.05,0.95,by=0.1), 
             las = 1, cex.axis = 2.5)

        title(main = year[j], cex.main = 5)
        
      }
    }, 
    movie.name = paste0(spName," - Mean.gif"),
    ani.height = 1000,
    ani.width = 1800
  )
}
