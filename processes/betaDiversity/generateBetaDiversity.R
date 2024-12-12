
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(tidyterra)

# Import ansvarsarter
ansvarsArter <- readRDS("data/ansvarsArterList.RDS")
source("functions/defineRegion.R")

# Import full ansvartsarter species data at all sites
# Direct pc to find correct files
ogDirectory <- "../BioDivMapping/data/run_2024-10-11/modelOutputs/processedOutputs"
listProbFiles <- list.files(ogDirectory, full.names = TRUE)[grep("probability" ,list.files(ogDirectory))]
listProbFiles <- listProbFiles[-grep(".xml", listProbFiles)]

###------------------------------###
### 1. Calculate beta diversity ####
###------------------------------###

# Specific species data
groupsToProcess <- c("hymenopterans", "vascularPlants", "birds")
resolutionInKm <- 10

for (focalGroup in groupsToProcess) {
  focalGroup <- "hymenopterans"
  listProbFilesGroup <- listProbFiles[grep(focalGroup, listProbFiles)]
  
  # Import probabilities for different species and cut down to ansvarsarter
  probabilites <- rast(listProbFilesGroup)
  speciesToMap <- intersect(names(probabilites), ansvarsArter$simpleScientificName)
  ansvarsarterProbs <- probabilites[[speciesToMap]]
  
  # Convert everything above 0.7 to a 1 (this is where processing time starts to blow out)
  print(paste0("Calculating local richness for ", focalGroup))
  speciesAssumed <- ifel(ansvarsarterProbs > 75, 1, 0)
  localRichness <- sum(speciesAssumed)
  
  # Get species richness in the 10 x 10km area surrounding each pixel (LONG PROCESSING TIME)
  print(paste0("Calculating beta diversity for ", focalGroup))
  focalAttempt <- terra::focal(speciesAssumed, w = resolutionInKm*4+1, fun = any, na.rm = TRUE)
  focalRichness <- sum(focalAttempt)
  focalBeta <- localRichness/focalRichness
  writeRaster(focalBeta, paste0("processes/betaDiversity/data/",focalGroup, resolutionInKm ,"betaDiversity.tiff"), 
              overwrite = TRUE)
}
