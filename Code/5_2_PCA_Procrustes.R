# Describe / visualise environmental changes

library(dplyr)
library(lubridate)
library(vegan)
library(ggplot2)

## Prepare the environmental data ----

df <- readRDS("Data/explan.RDS")

# remove columns we don't want
df@data <- subset(df@data, select =-c(Cruise.code, set.number, stratum, duration, distance))

# standardize the variables

# make a list of environmental dataframes (one per year)
env <- vector(mode = "list", length = length(1971:2020))
years <- 1971:2020
for(i in seq_along(years)){
  
  # index to rows in the given year
  index <- which(year(as.POSIXct(df@endTime)) == years[i])
  
  # transform the data
  env[[i]] <- df@data[index,] 
  env[[i]][,-c(8,9)] <- decostand(env[[i]][,-c(8,9)], method = "standardize")
  env[[i]][,c(8,9)] <- scale(sqrt(env[[i]][,c(8,9)]), center = FALSE)
}

## groups of variables for each hypothesis ----

# energy: bottom.temperature + surface.temperature
# productivity: bottom.salinity + ice
# climate: surface.temperature.5.var + bottom.temperature.5.var
# habitat: slope + coarsness + bottom.salinity
# stress: depth + seal + surfaceBottomTempDiff

# make a palette representing these hypothesis groups
pal <- RColorBrewer::brewer.pal(5, name = "Dark2")


## PCA ----

env_pca <- lapply(env, rda)
# pdf("Outputs/pca_biplots.pdf", onefile = TRUE)
# for(i in seq_along(years)){
#   biplot(env_pca[[i]], 
#          col = pal,
#          scaling = "species", 
#          display = c("species", "sites"),
#          type = c("text"),
#          main = years[i])
# }
# dev.off()

pdf("Outputs/pca_biplots.pdf", onefile = TRUE)
for(i in seq_along(years)){
  
  ## extract scores - these are coordinates in the RDA space
  sc_si <- scores(env_pca[[i]], display="sites", choices=c(1,2), scaling=2)
  sc_sp <- scores(env_pca[[i]], display="species", choices=c(1,2), scaling=2)
  
  # get scalings per hypothesis group
  sc_energy <- sc_sp[c("bottom.temperature", "surface.temperature"),]
  sc_prod <- sc_sp[c("bottom.salinity", "ice"),]
  sc_climate <- sc_sp[c("surface.temperature.5.var", "bottom.temperature.5.var"),]
  sc_habitat <- sc_sp[c("slope", "coarsness", "bottom.salinity"),]
  sc_stress <- sc_sp[c("depth", "seal", "surfaceBottomTempDiff"),]
  sc_hyp <- list(sc_energy, sc_prod, sc_climate, sc_habitat, sc_stress)
  
  # make a biplot
  plot(env_pca[[i]],
       scaling = 2, # set scaling type 
       type = "none", # this excludes the plotting of any points from the results
       frame = FALSE,
       # # set axis limits
       xlim = c(-2,2), 
       ylim = c(-2,2),
       # label the plot (title, and axes)
       main = years[i]
  )
  # add arrows for effects of the expanatory variables
  for(n in seq_along(sc_hyp)){
    arrows(0, 0, sc_hyp[[n]][,1], sc_hyp[[n]][,2], col = pal[n])
    text(x = sc_hyp[[n]][,1] + 0.1, 
         y = sc_hyp[[n]][,2] + 0.1, 
         labels = rownames(sc_hyp[[n]]),
         col = pal[n])
    legend(x = 2.2, y = 1.5, fill = pal, bty = "n",
           legend = c("Energy","Productivity","Climate","Habitat","Stress"),
    )
  }
  
}
dev.off()

# Procrustes analysis ----

pro_res <- vector("list", length = length(years))
pdf("Outputs/pca_procrustes.pdf", onefile = TRUE)
for(i in seq_along(years)[-1]){
  
  x1 <- scores(env_pca[[i-1]], scaling = 2)$species
  x2 <- scores(env_pca[[i]], scaling = 2)$species
  
  pro <- protest(x1, x2)
  
  plot(pro, 
       main = paste(c(years[i-1], "-", years[i])))
  text(pro, display = "target")
  
  pro_res[[i]] <- pro
}
dev.off()

pro_res[[1]] <- NULL

# make a table with procrustes results

pro_df <- data.frame(
  "year1" = years[-length(years)],
  "year2" = years[-1],
  "R" = lapply(pro_res, function(x) permustats(x)$statistic) %>% unlist(),
  "p-val" = lapply(pro_res, function(x) x$signif) %>% unlist()
)

ggplot(pro_df) +
  geom_line(aes(x = year1, y = R)) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  labs(y = "Procrustes correlation", y = "Year")
