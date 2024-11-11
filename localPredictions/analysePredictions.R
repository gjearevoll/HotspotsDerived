#### Prediction compilation
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(igraph)


# Import field data
fieldWorkResults <- readRDS("localPredictions/data/fieldWorkResults.RDS")

processedSpecies <- gsub(".tiff", "", list.files("localPredictions/data/predictions/resolution1000"))
species <- intersect(processedSpecies, names(fieldWorkResults))

# # # Importing predictions from model for 500m resolution (if needed)
full500Raster <- rast("../BioDivMapping/data/run_2024-10-11/modelOutputs/speciesprobability_vascularPlants.tiff")
speciesRaster <- full500Raster[[species]]
speciesRaster <- crop(speciesRaster, ext(st_transform(fieldWorkResults,crs = st_crs(speciesRaster))))
writeRaster(speciesRaster, "localPredictions/data/predictions/vascularPlantPredictions.tiff",overwrite = TRUE)

# Import relevant species
resolutions <- c("250", "500", "1000")

datasets <- lapply(resolutions, FUN = function(x) {
  if (x == "500") {
    predictions <- rast("localPredictions/data/predictions/vascularPlantPredictions.tiff") |>
      project(sf::st_crs(25833)$proj4string)
    predictions$pixel <- 1:ncell(predictions)
  } else {
    fileList <- list.files(paste0("localPredictions/data/predictions/resolution", x), full.names = TRUE)
    speciesNames <- gsub(".tiff", "" ,list.files(paste0("localPredictions/data/predictions/resolution", x)))
    
    dataList <- lapply(fileList, rast) |>
      setNames(speciesNames)
    predictions <- rast(lapply(dataList, `[`, "mean"))[[species]]
    predictions$pixel <- 1:ncell(predictions)
  }
  return(predictions)
}) |> setNames(resolutions)

### First we get two tables, one showing species richness at all points and another showing species presences
datasets[[1]]

fieldWork <- readRDS("localPredictions/data/fieldWorkResultsWide.RDS")


# Import predictive data
comparisonTables <- lapply(resolutions, FUN = function(x) {
  dataList <- datasets[[x]]
  dataList$pixel <- 1:ncell(dataList)
  dataFrame <- as.data.frame(dataList)
  speciesRichness <- c(sum(dataList[[1:(nlyr(dataList)-1)]]), dataList$pixel)
  
  # Get cumulative richness based on the raster pixel they're in
  dfToUse <- fieldWork
  dfToUse$pixel <- terra::extract(speciesRichness, dfToUse)[,"pixel"]
  dfToUseGrouped <- aggregate(. ~ pixel, data=st_drop_geometry(dfToUse[,!(colnames(dfToUse) %in% c("siteID", "predictedRichness"))]), FUN=sum)
  
  dfToUseGrouped[,-1] <- lapply(dfToUseGrouped[,-1], FUN = function(x) {
    x[x > 1] <- 1
    x
  })
  
  # Group the observations based on distance
  # Sort these observations into groups
  pixelDistance <- dfToUse[dfToUse$pixel %in% dfToUseGrouped$pixel, c("pixel")]
  adj <- st_distance(pixelDistance)
  adj <- matrix(as.numeric(as.numeric(adj)) < 25000, nrow = nrow(adj))
  
  g <- graph_from_adjacency_matrix(adj)
  pixelDistance$group <- factor(components(g)$membership)
  
  dfToUseGrouped$group <- pixelDistance$group[match(dfToUseGrouped$pixel, pixelDistance$pixel)]
  
  speciesTable <- lapply(species, FUN = function(sp) {
    if (!(sp %in% names(dfToUseGrouped)) | !(sp %in% names(dataFrame))) {return(NA)}
    
    focalSpeciesObserved <- dfToUseGrouped[,c(sp, "pixel", "group")]
    # Check whether we actually have enough presences or absences for a valid comparison
    ratioPA <- sum(focalSpeciesObserved[,1])/nrow(focalSpeciesObserved)
    if (ratioPA > 0.9 | ratioPA < 0.1) {return(NA)}
    
    focalSpeciesTable <- merge(focalSpeciesObserved, 
                               dataFrame[,c(sp, "pixel")],
                               by = "pixel", all.x = TRUE)
    names(focalSpeciesTable)[c(2,4)] <- c("observedPresence", "predictedPresence")
    focalSpeciesTable
  }) |> setNames(species)
  speciesTable <- speciesTable[!is.na(speciesTable)]
  
  # Start with richness
  dfToUseGrouped$observedRichness <- rowSums(dfToUseGrouped[,2:(ncol(dfToUseGrouped)-1)])
  
  predictedRichnessTable <- terra::extract(speciesRichness, dfToUse)
  predictedRichnessTable <- predictedRichnessTable[!duplicated(predictedRichnessTable$pixel) & !is.na(predictedRichnessTable$sum),] %>%
    rename(predictedRichness = sum) %>%
    as.data.frame()
  
  # Merged stuff
  mergedNumbers <- merge(dfToUseGrouped[,c("pixel", "observedRichness", "group")], 
                         predictedRichnessTable[,c("pixel", "predictedRichness")],
                         by = "pixel", all.x = TRUE)
  
  return(list(richness = mergedNumbers, species = speciesTable))
}) |> setNames(paste0("resolution",resolutions))



# First we can run a single species comparison
availableSpecies <- names(comparisonTables$resolution500$species)
speciesSelected <- "Anemone_nemorosa" 

focalRaster <- lapply(comparisonTables, FUN = function(x) {
  speciesTable <- x$species[[speciesSelected]]
  speciesTable <- speciesTable[complete.cases(speciesTable),]
  speciesTable
})

resolution <- "250"
summary(lm(data = focalRaster[[paste0("resolution", resolution)]], predictedPresence ~ factor(observedPresence)))
ggplot(focalRaster[[1]], aes(x=factor(observedPresence), y=predictedPresence)) + 
  geom_boxplot()

# Richness comparison
df <- as.data.frame(comparisonTables[[1]]$richness) %>%
  filter(!is.na(predictedRichness))

plot(data = df, observedRichness ~ predictedRichness)
summary(lm(data = df, predictedRichness ~ observedRichness))


