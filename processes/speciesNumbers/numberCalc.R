
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(dplyr)
library(terra)

sourceDirectory <- "../BioDivMapping/data/run_2024-10-11/modelOutputs/processedOutputs/"

###--------------------------------###
### 1. Set up data for processing ####
###--------------------------------###

### Get names of all taxa we currently have finished data for
taxaNames <- unique(gsub(".tiff", "",
                         gsub("final.tiff", "", 
                              gsub("allspeciesstats_", "", 
                                   list.files(sourceDirectory, pattern = "allspecies")))))

# Remove bird names
taxaNames <- taxaNames[c(1:9, 11:16)]

# Import each taxa and get names of species therein
speciesLists <- lapply(taxaNames, FUN = function(taxa) {
  focalRast <- rast(paste0(sourceDirectory, "speciesprobability_", taxa, ".tiff"))
  nameList <- names(focalRast)
  nameList
}) |> setNames(taxaNames)

# Check speciesDataProcessed file
dataMatching <- data.frame(taxa = taxaNames, grouping = c("AmphibiansReptiles", "Insects", "Bats", "Insects", "Birds", 
                                                          "Insects", "Insects", "Insects", "Fungi", "Insects", "Fungi",
                                                          "Mammals", "Mosses", "Insects", "Plants"))

###--------------------------------###
### 2. Set up data for processing ####
###--------------------------------###

# For now, do this without a loop
finalStats <- lapply(taxaNames, FUN = function(x) {
  grouping <- dataMatching$grouping[dataMatching$taxa == x]
  nameList <- speciesLists[[x]]
  
  # Download correct species data
  print(paste0("Loading ", grouping, " processed data for ", x))
  focalData <- readRDS(paste0("data/processedData/speciesDataProcessed", grouping, ".RDS"))
  
  # Dataset types
  datasetTypes <- sapply(focalData, FUN= function(dt) {
    dt$dataType[1]
  }) 
  datasetTypes <- tally(group_by(data.frame(datasetTypes), datasetTypes))
  
  # Get basic stats
  numberDatasets <- length(focalData)
  numberSpecies <- length(nameList)
  stats <- data.frame(numberSpecies, numberDatasets)
  
  # Process and return only relevant data
  print(paste0("Reducing ", grouping, " data for ", x))
  focalDataFrame <- do.call(rbind, lapply(1:length(focalData), FUN = function(d) {
    data <- focalData[[d]]
    reducedData <- data[,c("simpleScientificName", "taxa", "dataType")] %>%
      filter(simpleScientificName %in% nameList)
    if (nrow(reducedData) == 0) {return(NA)}
    reducedData$dsName <- names(focalData)[d]
    reducedData
  })) %>%
    as.data.frame()
  
  rm("focalData")
  gc()
  
  # Now produce species tallies
  focalDataSpecies <- focalDataFrame %>%
    filter(taxa %in% taxa) 
  
  focalDataCountByTypes <- focalDataSpecies %>%
    group_by(dataType) %>%
    tally()
  focalDataCountByDataset <- focalDataSpecies %>%
    group_by(dsName) %>%
    tally()
  
  
  finalList <- list(dataTypes = as.data.frame(focalDataCountByTypes), datasets = as.data.frame(focalDataCountByDataset), stats = stats, datasetTypes = datasetTypes)
  finalList
}
) |> setNames(taxaNames)

# Save everything
saveRDS(finalStats, "speciesNumbers/data/speciesDataCounts.RDS")

