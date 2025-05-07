
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(terra)
library(dplyr)
library(sf)
library(csmaps)
focalTaxa <- read.csv("data/focalTaxa.csv")

betaDirectory <- "processes/betaDiversity/data/"

baseRaster <- rast("processes/betaDiversity/data/fungi10betaDiversity.tiff")

###---------------------###
### 1. Combine rasters ####
###---------------------###
rastersToImport <- grep("10", list.files(betaDirectory), value = TRUE)
rasterPackage <- lapply(rastersToImport, FUN = function(x) {
  rawR <- rast(paste0(betaDirectory, x))
  if (ext(rawR) != ext(baseRaster)) {
    return(project(rawR, baseRaster, "mode"))
  } else {return(rawR)}
}) |> setNames(gsub("10betaDiversity.tiff", "", rastersToImport))
rasterPackage <- rast(rasterPackage)

# Get border to map to
norwayBorder <- st_union(nor_municip_map_b2024_default_sf)
norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), rasterPackage[[1]]))

# Crop to better border
rasterPackageCropped <- crop(rasterPackage, norwayBorderProjected)

### Save to file 
writeRaster(rasterPackageCropped, filename = paste0(betaDirectory, "rasterPackage4.tiff"), overwrite = TRUE)
