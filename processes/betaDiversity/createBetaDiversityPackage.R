
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(terra)
library(dplyr)
library(sf)
library(csmaps)
focalTaxa <- read.csv("data/focalTaxa.csv")

betaDirectory <- "processes/betaDiversity/data/"

###---------------------###
### 1. Combine rasters ####
###---------------------###
rastersToImport <- c("vascularPlants", "birds")
rasterPackage <- rast(paste0(betaDirectory, rastersToImport, "10betaDiversity.tiff")) |>
  setNames(rastersToImport)

### Fungi import
fungiLayer <- rast(paste0(betaDirectory, "fungi10betaDiversity.tiff")) |>
  terra::project(rasterPackage, method = "mode") |>
  setNames("fungi")

# Get border to map to
norwayBorder <- st_union(nor_municip_map_b2024_default_sf)
norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), rasterPackage[[1]]))

# Crop to better border
rasterPackageCropped <- c(rasterPackage, fungiLayer) |>
  crop(norwayBorderProjected)

### Save to file 
writeRaster(rasterPackageCropped, filename = paste0(betaDirectory, "rasterPackage4.tiff"), overwrite = TRUE)
