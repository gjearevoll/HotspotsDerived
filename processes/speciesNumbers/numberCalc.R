
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(dplyr)
library(terra)
library(sf)

sourceDate <- "2025-01-06"
sourceDirectory <- paste0("../BioDivMapping/data/run_",sourceDate,"/modelOutputs/processedOutputs/")

###--------------------------------###
### 1. Set up data for processing ####
###--------------------------------###

### Get names of all taxa we currently have finished data for
taxaNames <- c("insects", "fungi", "lichens", "vascularPlants", "birds")
insectList <- c("aquaticInsects", "beetles", "bugs", "butterfliesMoths", "cockroaches", "earwigs", "flies", 
                "grasshopperLikeInsects", "hymenopterans", "netWingedInsects", "scorpionflies", "snakeflies",
                "spiders")

# Import each taxa and get names of species therein
speciesLists <- lapply(taxaNames, FUN = function(taxa) {
  focalRast <- rast(paste0(sourceDirectory, "speciesprobability_", taxa, ".tiff"))
  nameList <- names(focalRast)
  nameList
}) |> setNames(taxaNames)

###--------------------------------###
### 2. Set up data for processing ####
###--------------------------------###

# For now, do this without a loop
finalStats <- lapply(taxaNames, FUN = function(x) {
  nameList <- speciesLists[[x]]
  
  # Download correct species data
  print(paste0("Loading processed data for ", x))
  nameChange <- paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep="")
  focalData <- readRDS(paste0("data/processedData/speciesDataProcessed", nameChange, ".RDS"))

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
  print(paste0("Reducing data for ", x))
  focalDataFrame <- do.call(rbind, lapply(1:length(focalData), FUN = function(d) {
    data <- st_drop_geometry(focalData[[d]])
    reducedData <- data[,c("simpleScientificName", "taxa", "dataType")] %>%
      filter(simpleScientificName %in% nameList)
    if (nrow(reducedData) == 0) {return(NA)}
    reducedData$dsName <- names(focalData)[d]
    reducedData
  })) %>%
    as.data.frame()
  
  if (colnames(focalDataFrame)[1] == "V1") {return(NA)}
  
  rm("focalData")
  gc()
  
  # Now produce species tallies
  focalDataSpecies <- focalDataFrame %>%
    filter(taxa %in% x) 
  if (x == "insects") {
    focalDataSpecies <- focalDataFrame %>%
      filter(taxa %in% insectList) 
  }
  
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
saveRDS(finalStats, paste0("processes/speciesNumbers/data/speciesDataCounts_",sourceDate,".RDS"))

