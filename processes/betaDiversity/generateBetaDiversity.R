
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(tidyterra)

args <- commandArgs(TRUE)
# THis should only run if the script is being run from the command line

if (length(args) != 0) {
  # Set arguments
  groupsToProcess <- args[1]
  # Set the working directory
  setwd("~/HotspotsDerived")
}

# Import ansvarsarter
ansvarsArter <- readRDS("data/ansvarsArterList.RDS")
source("functions/defineRegion.R")
source("functions/get_gbif_backbone.R")

# Import full ansvartsarter species data at all sites
# Direct pc to find correct files
ogDirectory <- "../BioDivMapping/data/run_2025-01-06/modelOutputs/processedOutputs"
listProbFiles <- list.files(ogDirectory, full.names = TRUE)[grep("probability" ,list.files(ogDirectory))]

norwayBorder2 <- sf::read_sf("../BioDivMapping/data/external/norge_border/Noreg_polygon.shp")

# Code for bird groups
birdSubgroups <- c("groundNesters", "waders", "woodpeckers")
if (groupsToProcess %in% birdSubgroups) {
  birdList <- read.csv("../BioDivMapping/data/external/birdTypeList.csv", sep = ",") %>%
    filter(group == groupsToProcess)
  birdGroup <- TRUE
} else {birdGroup <- FALSE}

###------------------------------###
### 1. Calculate beta diversity ####
###------------------------------###

# Specific species data
resolutionInKm <- 10

for (focalGroup in groupsToProcess) {
  
  if (birdGroup) {
    listProbFilesGroup <- listProbFiles[grep("birds", listProbFiles)]
    refinementList <- gsub(" ", "_", birdList$simpleName)
  } else {
    listProbFilesGroup <- listProbFiles[grep(focalGroup, listProbFiles)]
    refinementList <- ansvarsArter$simpleScientificName
  }
  
  # Import probabilities for different species and cut down to ansvarsarter
  probabilites <- rast(listProbFilesGroup)
  norwayBorderProjected <- terra::project(vect(norwayBorder2), probabilites)
  speciesToMap <- intersect(names(probabilites), refinementList)
  ansvarsarterProbs <- probabilites[[speciesToMap]]
  ansvarsarterProbs <- crop(ansvarsarterProbs, norwayBorderProjected, mask = T)
  
  # Convert everything above 0.7 to a 1 (this is where processing time starts to blow out)
  print(paste0("Calculating local richness for ", focalGroup))
  speciesAssumed <- ifel(ansvarsarterProbs > 75, 1, 0)
  localRichness <- sum(speciesAssumed)
  
  # Get species richness in the 10 x 10km area surrounding each pixel (LONG PROCESSING TIME)
  print(paste0("Calculating beta diversity for ", focalGroup))
  focalAttempt <- terra::focal(speciesAssumed, w = 4*resolutionInKm+1, fun = any, na.rm = TRUE)
  focalRichness <- sum(focalAttempt)
  focalBeta <- 1-localRichness/focalRichness
  writeRaster(focalBeta, paste0("processes/betaDiversity/data/",focalGroup, resolutionInKm ,"betaDiversity.tiff"), 
              overwrite = TRUE)
}
