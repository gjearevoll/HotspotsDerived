### Need to get all necessary data into a folder, and then aggregate it all
library(dplyr)
library(terra)

relevantTaxa <- c("vascularPlants", "hymenopterans", "birds")
appDirectory <- "communicationsPackageHotspots/data"
baseDirectory <- "../BioDivMapping/data/run_2024-10-11/modelOutputs/processedOutputs/"

# Import red list
source("functions/importRedList.R")
redList <- importRedList(c("VU", "EN", "CR"))
redList$simpleScientificName <- gsub(" ", "_", redList$species)
saveRDS(redList, "communicationsPackageHotspots/data/redList.RDS")

# Import all richness data
richnessDirectory <- "data/modelOutputs"

availableFiles <- list.files(richnessDirectory, full.names = TRUE)
filesToImport <- c(sapply(relevantTaxa, FUN = function(taxa) {
  availableFiles[grep(taxa, availableFiles)]
}))

for (taxaFile in filesToImport) {
  print(paste("aggregating ", taxaFile))
  starterRast <- rast(taxaFile)
  aggRast <- aggregate(starterRast, fact = 10)
  fileName <- paste0(appDirectory, "/", gsub("data/modelOutputs/", "", taxaFile))
  writeRaster(aggRast, fileName, overwrite = TRUE)
}

# Import all beta diversity data
betaDirectory <- "betaDiversity/data"
availableFiles <- list.files(betaDirectory, full.names = TRUE)
filesToImport <- c(sapply(relevantTaxa, FUN = function(taxa) {
  allFiles <- availableFiles[grep(taxa, availableFiles)]
  allFiles[grep("10", allFiles)]
}))

for (taxaFile in filesToImport) {
  print(paste("aggregating ", taxaFile))
  starterRast <- rast(taxaFile)
  aggRast <- aggregate(starterRast, fact = 10)
  fileName <- paste0(appDirectory, "/", gsub("betaDiversity/data/", "", taxaFile))
  writeRaster(aggRast, fileName, overwrite = TRUE)
}

# Import all hotspots data
hotspotDirectory <- "overlays/data/hotspots"
availableFiles <- list.files(hotspotDirectory, full.names = TRUE)
filesToImport <- c(sapply(relevantTaxa, FUN = function(taxa) {
  taxaGrouping <- ifelse(taxa == "birds", "fugler", ifelse(taxa == "vascularPlants", "karplanter", "insekter"))
  fileNames <- availableFiles[grep(taxaGrouping, availableFiles)]
  fileNames[!grepl("xml", fileNames)]
}))

for (taxaFile in filesToImport) {
  print(paste("aggregating ", taxaFile))
  starterRast <- rast(taxaFile)
  aggRast <- aggregate(starterRast, fact = 10, fun = "modal", na.rm=TRUE)
  fileName <- paste0(appDirectory, "/", gsub("overlays/data/hotspots/", "", taxaFile))
  writeRaster(aggRast, fileName, overwrite = TRUE)
}

# Define ansvararter to plot
ansvarsarter <- readRDS("data/ansvarsArterList.RDS")
namesFromGroups <- do.call(rbind, lapply(relevantTaxa, FUN = function(taxa) {
  taxaNames <- names(rast(paste0(baseDirectory, "speciesprobability_", 
                         taxa, ".tiff")))
  df <- data.frame(species = taxaNames, taxa = taxa)
  df <- df[df$species %in% ansvarsarter$simpleScientificName,]
  df$fullName <- ansvarsarter$GBIFName[match(df$species, ansvarsarter$simpleScientificName)]
  df
})) %>% group_by(taxa) %>%
  sample_n(4)
focalSpeciesDDVector <- split(namesFromGroups$fullName, f = namesFromGroups$taxa)
focalSpeciesDDList <- lapply(focalSpeciesDDVector, FUN = function(x) {
  speciesList <- as.list(x)
  names(speciesList) <- x
  speciesList
}
)
saveRDS(focalSpeciesDDList, "communicationsPackageHotspots/data/ddList.RDS")
saveRDS(namesFromGroups, "communicationsPackageHotspots/data/namesFromGroups.RDS")


# Import all individual data
for (taxa in relevantTaxa) {
  fullList <- rast(paste0(baseDirectory, "speciesprobability_", 
                          taxa, ".tiff"))
  fullUncertaintyList <- rast(paste0(baseDirectory, "speciesuncertainty_", 
                          taxa, ".tiff"))
  narrowedList <- fullList[[namesFromGroups$species[namesFromGroups$taxa %in% taxa]]]
  aggNarrowedList <- aggregate(narrowedList, fact = 10)
  narrowedUList <- fullUncertaintyList[[namesFromGroups$species[namesFromGroups$taxa %in% taxa]]]
  aggnarrowedUList <- aggregate(narrowedUList, fact = 10)
  writeRaster(aggNarrowedList, paste0(appDirectory, "/", taxa, "_probability.tiff"), overwrite = TRUE)
  writeRaster(aggnarrowedUList, paste0(appDirectory, "/", taxa, "_uncertainty.tiff"), overwrite = TRUE)
}

# Import images
# We need to download photos from iNaturalist, including their URL and the user name of the individual
# who took the photo to give appropriate credit

# # Check if there is already an image credit file
# if (file.exists("visualisation/hotspotMaps/data/imageCredit.RDS")) {
#   imageCredit <- readRDS("visualisation/hotspotMaps/data/imageCredit.RDS")
# } else {
#   imageCredit <- data.frame(species = redList$species[redList$valid],
#                             credit = NA,
#                             url = NA,
#                             stringsAsFactors = FALSE)
# }

taxaFolder <- "communicationsPackageHotspots/data/photos/"
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

