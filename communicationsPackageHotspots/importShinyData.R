### Need to get all necessary data into a folder, and then aggregate it all

### --------------------------------------------- ###
#### 0. Set base parameters and import libraries ####
### --------------------------------------------- ###

library(dplyr)
library(terra)

relevantTaxa <- c("vascularPlants", "insects", "birds", "fungi", "lichens")
appDirectory <- "communicationsPackageHotspots/data"
importedDirectory <- paste0(appDirectory, "/importedData")
baseDirectory <- "../BioDivMapping/data/run_2025-01-06/modelOutputs/processedOutputs/"

if (!dir.exists(importedDirectory)) {
  create.dir(importedDirectory)
}

### ----------------------------------------- ###
#### 1. Import red list and Norwegian border ####
### ----------------------------------------- ###

# Import red list
source("functions/importRedList.R")
redList <- importRedList(c("VU", "EN", "CR"))
redList$simpleScientificName <- gsub(" ", "_", redList$species)
saveRDS(redList, "communicationsPackageHotspots/data/importedData/redList.RDS")

norwayBorder2 <- sf::read_sf("data/norge_border/Noreg_polygon.shp")

### ------------------------- ###
#### 2. Import richness list ####
### ------------------------- ###

richnessDirectory <- "data/modelOutputs/run_2025-01-06"

availableFiles <- list.files(richnessDirectory, full.names = TRUE)
filesToImport <- c(sapply(relevantTaxa, FUN = function(taxa) {
  availableFiles[grep(taxa, availableFiles)]
}))

for (taxaFile in filesToImport) {
  print(paste("aggregating ", taxaFile))
  starterRast <- rast(taxaFile)
  aggRast <- aggregate(starterRast, fact = 10)
  fileName <- paste0(importedDirectory, "/", gsub(paste0(richnessDirectory,"/"), "", taxaFile))
  writeRaster(aggRast, fileName, overwrite = TRUE)
}

### -------------------------- ###
#### 3. Import beta diversity ####
### -------------------------- ###

betaDirectory <- "processes/betaDiversity/data"
availableFiles <- list.files(betaDirectory, full.names = TRUE)
availableFiles <- grep("10", availableFiles, value = TRUE)
filesToImport <- c(sapply(relevantTaxa, FUN = function(taxa) {
  allFiles <- availableFiles[grep(taxa, availableFiles)]
  allFiles
}))

for (taxaFile in filesToImport) {
  print(paste("aggregating ", taxaFile))
  starterRast <- rast(taxaFile)
  starterRastCropped <- crop(starterRast, project(vect(norwayBorder2), starterRast), mask = T)
  aggRast <- aggregate(starterRastCropped, fact = 10)
  fileName <- paste0(importedDirectory, "/", gsub("processes/betaDiversity/data/", "", taxaFile))
  writeRaster(aggRast, fileName, overwrite = TRUE)
}

### -------------------- ###
#### 4. Import Hotspots ####
### -------------------- ###

# This requires you to have run processes/createSpeciesHotspots.R
hotspotFile <- rast("data/allHotspots.tiff")
aggRast <- aggregate(hotspotFile, fact = 10, fun = "modal", na.rm=TRUE)
fileName <- paste0(importedDirectory, "/allHotspots.tiff")
writeRaster(aggRast, fileName, overwrite = TRUE)

### ----------------------------------- ###
#### 5. Import individual species data ####
### ----------------------------------- ###

# Define species to plot
speciesToPlot <- read.csv("communicationsPackageHotspots/data/speciesForDisplay.csv", sep = ",")
namesFromGroups <- do.call(rbind, lapply(relevantTaxa, FUN = function(taxa) {
  taxaNames <- names(rast(paste0(baseDirectory, "speciesprobability_", 
                         taxa, ".tiff")))
  df <- data.frame(species = taxaNames, taxa = taxa)
  df <- df[df$species %in% speciesToPlot$species,]
  df$fullName <- gsub("_"," ", df$species)
  df
}))
focalSpeciesDDVector <- split(namesFromGroups$fullName, f = namesFromGroups$taxa)
focalSpeciesDDList <- lapply(focalSpeciesDDVector, FUN = function(x) {
  speciesList <- as.list(x)
  names(speciesList) <- x
  speciesList
}
)
saveRDS(focalSpeciesDDList, "communicationsPackageHotspots/data/importedData/ddList.RDS")
saveRDS(namesFromGroups, "communicationsPackageHotspots/data/importedData/namesFromGroups.RDS")


# Import all individual data
for (taxa in relevantTaxa) {
  fullList <- rast(paste0(baseDirectory, "speciesprobability_", 
                          taxa, ".tiff"))
  fullUncertaintyList <- rast(paste0(baseDirectory, "speciesuncertainty_", 
                          taxa, ".tiff"))
  narrowedList <- fullList[[namesFromGroups$species[namesFromGroups$taxa %in% taxa]]]
  narrowedListCropped <- crop(narrowedList, project(vect(norwayBorder2), narrowedList), mask = T)
  aggNarrowedList <- aggregate(narrowedListCropped, fact = 10)
  narrowedUList <- fullUncertaintyList[[namesFromGroups$species[namesFromGroups$taxa %in% taxa]]]
  narrowedUListCropped <- crop(narrowedUList, project(vect(norwayBorder2), narrowedUList), mask = T)
  aggnarrowedUList <- aggregate(narrowedUListCropped, fact = 10)
  writeRaster(aggNarrowedList, paste0(importedDirectory, "/", taxa, "_probability.tiff"), overwrite = TRUE)
  writeRaster(aggnarrowedUList, paste0(importedDirectory, "/", taxa, "_uncertainty.tiff"), overwrite = TRUE)
}

# Import images from iNaturalist

taxaFolder <- "communicationsPackageHotspots/data/importedData/photos/"
if (!dir.exists(taxaFolder)) { dir.create(taxaFolder)}
imageCredit <- data.frame(species = unlist(focalSpeciesDDList),
                          credit = NA,
                          url = NA,
                          stringsAsFactors = FALSE)
for (species in imageCredit$species) {
  simpleName <- namesFromGroups$species[namesFromGroups$fullName == species]
  # Create species folder
  speciesFolder <- paste0(taxaFolder,species)
  if (!file.exists(speciesFolder)) {
    dir.create(speciesFolder)
  }
  # Create species file (if it doesn't already exist)
  if (file.exists(paste0(speciesFolder, "/speciesImage.jpg"))) {next}
  tryCatch({inatAttempt <- rinat::get_inat_obs(taxon_name = simpleName, maxresults = 10)
  imageRow <- inatAttempt[!is.na(inatAttempt$image_url) & inatAttempt$image_url != "",]
  imageCredit$url[imageCredit$species == species] <- imageRow$url[1]
  imageCredit$credit[imageCredit$species == species] <- imageRow$user_login[1]
  download.file(imageRow$image_url[1], destfile = paste0(speciesFolder,"/speciesImage.jpg" ), mode = 'wb')
  }
  , error = function(x){print(paste0("Error for species ",species, ". ", x))})
}

imageCredit$credit[imageCredit$credit == "" | is.na(imageCredit$credit)] <- "User name not provided"
saveRDS(imageCredit, "communicationsPackageHotspots/data/imageCredit.RDS")

