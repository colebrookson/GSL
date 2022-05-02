# Load R package
library(spacetime)
library(INLA)
library(mapdata) 
library(maptools)

#__________
# Load data
#__________
region <- readRDS("./Data/regionBuffer.RDS")
water <- readRDS("./Data/waterRegion.RDS")
canada <- readRDS("./Data/Canada.RDS")

#_______________________________________
# Build spatial mesh (for barrier model)
#_______________________________________
maxEdge <- 19500
meshSpace <- inla.mesh.2d(boundary = water, 
                          max.edge=maxEdge * c(1,5),
                          cutoff=maxEdge/4,
                          offset = c(10000,20000),
                          crs = CRS("EPSG:3798"))

# Find the triangle 
triangleWater <- inla.over_sp_mesh(x = water,
                               y = meshSpace,
                               type = "centroid")

# Number of triangles in the mesh
ntriangle <- nrow(meshSpace$graph$tv)

# Triangles within the barrier
barrierTriangle <- setdiff(1:ntriangle, triangleWater)

# 
polyBarrier <- inla.barrier.polygon(meshSpace,
                                    barrier.triangles = barrierTriangle)

# Visual check to make sure everything is OK
plot(meshSpace, asp = 1)
plot(region, add=TRUE, border = "orange", lwd = 3)
plot(canada, add=TRUE, col ="grey")
plot(polyBarrier, border = "red", add = TRUE)

# Save mesh and barrier triangle (for barrier model)
saveRDS(meshSpace,file = "./Mesh/meshSpace.RDS")
saveRDS(barrierTriangle, file = "./Mesh/barrierTriangle.RDS")
