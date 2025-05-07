
#### PREDICTION VS. FIELDWORK COMPARISON ####

###-----------------###
### 1. Data import ####
###-----------------###

#### Prediction compilation
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(igraph)
library(tidyterra)
library(csmaps)


# Import field data
fieldWorkResults <- readRDS("processes/localPredictions/data/fieldWorkResults.RDS")

processedSpecies <- gsub(".tiff", "", list.files("processes/localPredictions/data/predictions/resolution500"))
species <- intersect(processedSpecies, names(fieldWorkResults))

###-----------------------###
### 2. Prediction import ####
###-----------------------###

# Import relevant species
resolutions <- c("250", "500", "1000")

# Organise different resolutiosn of predictions into list objects

datasets <- lapply(resolutions, FUN = function(x) {

    fileList <- list.files(paste0("processes/localPredictions/data/predictions/resolution", x), full.names = TRUE)
    speciesNames <- gsub(".tiff", "" ,list.files(paste0("processes/localPredictions/data/predictions/resolution", x)))
    
    dataList <- lapply(fileList, rast) |>
      setNames(speciesNames)
    validSpecies <- speciesNames[speciesNames %in% species]
    predictions <- rast(lapply(dataList, `[`, "mean"))[[validSpecies]]
    predictions$pixel <- 1:ncell(predictions)

  return(predictions)
}) |> setNames(resolutions)

modelledSpecies <- Reduce(intersect,list(names(datasets[[1]]), names(datasets[[2]]), names(datasets[[3]]), species))

# Let's map out our predictions and our datasets first, including a map of the coastline
norwayBorder <- st_union(csmaps::nor_municip_map_b2024_default_sf)
norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), datasets[[1]]))

ggplot()  +
  geom_spatvector(data = norwayBorderProjected) + 
  geom_spatraster(data = datasets$`250`, aes(fill = Acer_platanoides)) +
  scale_fill_continuous(na.value = "transparent") +
  geom_sf(data = fieldWorkResults["Acer_platanoides"]) +
  coord_sf(crs = crs(datasets$`1000`$Acer_platanoides),xlim = c(-61000, 345000), ylim = c(6584000, 7100000))


### First we get two tables, one showing species richness at all points and another showing species presences
fieldWork <- readRDS("processes/localPredictions/data/fieldWorkResultsWide.RDS")

###----------------------------------###
### 3. Produce comparative datasets ####
###----------------------------------###

# Import predictive data
comparisonTables <- lapply(resolutions, FUN = function(x) {
  dataList <- datasets[[x]][[modelledSpecies]]
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
  
  speciesTable <- lapply(modelledSpecies, FUN = function(sp) {
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
  }) |> setNames(modelledSpecies)
  speciesTable <- speciesTable[!is.na(speciesTable)]
  
  # Start with richness
  dfToUseGrouped$observedRichness <- rowSums(dfToUseGrouped[,2:(ncol(dfToUseGrouped)-1)])
  
  predictedRichnessTable <- terra::extract(speciesRichness, dfToUse)
  predictedRichnessTable <- predictedRichnessTable[!duplicated(predictedRichnessTable$pixel) & !is.na(predictedRichnessTable$sum),] %>%
    rename(predictedRichness = sum) %>%
    as.data.frame()
  
  # Scale predictedRichness
  predictedRichnessTable$scaledRichness <- (predictedRichnessTable$predictedRichness-min(predictedRichnessTable$predictedRichness))/
    (max(predictedRichnessTable$predictedRichness)-min(predictedRichnessTable$predictedRichness))
  
  # Merged stuff
  mergedNumbers <- merge(dfToUseGrouped[,c("pixel", "observedRichness", "group")], 
                         predictedRichnessTable[,c("pixel", "predictedRichness", "scaledRichness")],
                         by = "pixel", all.x = TRUE)
  # Label resolution column
  mergedNumbers$res <- x
  speciesTable$res <- x
  
  return(list(richness = mergedNumbers, species = speciesTable))
}) |> setNames(paste0("resolution",resolutions))

###-------------------------------------###
### 3. Compare fieldwork with datasets ####
###-------------------------------------###

# Check which species it actually makes a difference for
availableSpecies <- intersect(intersect(names(comparisonTables$resolution1000$species),
                              names(comparisonTables$resolution500$species)),
                              names(comparisonTables$resolution250$species))

# See if there's a species in there with a significant difference
speciesSignificanceCheck <- lapply(availableSpecies[1:34], FUN = function(sp) {
  tableForChecking <- comparisonTables[[1]][["species"]][[sp]] %>%
    filter(!is.na(predictedPresence))
  modelled <- summary(lm(data = tableForChecking, observedPresence ~ predictedPresence))
  if (nrow(modelled$coefficients) == 1) {return(NA)}
  if (modelled$coefficients[2,4] > 0.1) {return(NA)} else {return(modelled)}
}) |> setNames(availableSpecies[1:34])
viableSpecies <- names(speciesSignificanceCheck[!is.na(speciesSignificanceCheck)])

# Now do a species check
speciesSelected <- "Filipendula_ulmaria"
focalRaster <- lapply(resolutions, FUN = function(x) {
  comparisonTable <- comparisonTables[[paste0("resolution", x)]]
  speciesTable <- comparisonTable$species[[speciesSelected]]
  speciesTable <- speciesTable[!is.na(speciesTable$predictedPresence),]
  speciesTable$res <- x
  speciesTable
}) |> setNames(paste0("resolution",resolutions))
combinedStats <- do.call(rbind, focalRaster)
combinedStats$factorRes <- factor(combinedStats$res, levels=c("250", "500", "1000"))
combinedStats$presence <- as.factor(combinedStats$observedPresence)
levels(combinedStats$presence) <- c("Funn", "Ikke funn")

# Comparison
resolution <- "250"
summary(lm(data = focalRaster[[paste0("resolution", resolution)]], predictedPresence ~ factor(observedPresence)))
ggplot(combinedStats, aes(x=factorRes, y=predictedPresence, fill = presence)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("white", "grey")) +
  #ylim(c(0,1)) +
  xlab("Oppløysing (m)") + 
  ylab("Predikert førekomst") + 
  #xlab("") + 
  #ylab("") + 
  labs(fill = "") +
  ggtitle(gsub("_",  " ",speciesSelected)) 


# All species together
comparisonTablesCompiled <- do.call(rbind, lapply(comparisonTables, FUN = function(ct) {
  boundComparisons <- do.call(rbind, ct$species[names(ct$species) != "res"])
  boundComparisons <- boundComparisons[!is.na(boundComparisons$predictedPresence),]
  boundComparisons$res <- unique(ct$richness$res)
  boundComparisons
}))
comparisonTablesCompiled$presence <- as.factor(comparisonTablesCompiled$observedPresence)
comparisonTablesCompiled$factorRes <- factor(comparisonTablesCompiled$res, levels=c("250", "500", "1000"))
ggplot(comparisonTablesCompiled, aes(x=factorRes, y=predictedPresence, fill = presence)) + 
  geom_boxplot() +
  ggtitle("All species compiled") + xlab("Resolution (m)") + ylab("Predicted presence")



# Richness comparison
df <- do.call(rbind, lapply(comparisonTables, FUN = function(x) {x$richness}))
df$res1 <- factor(df$res , levels=c('250', '500', '1000'))
# df <- as.data.frame(comparisonTables[[1]]$richness) %>%
#   filter(!is.na(predictedRichness))

ggplot(df, aes(x=scaledRichness, y=observedRichness, colour=res1)) + 
  geom_point(size=2) +
  scale_colour_brewer(palette = "Set2") +
  theme_bw() + labs(colour = "Oppløysing\n(m)", x = "Predikert rikdom (skalert)", y = "Observert rikdom")
summary(lm(data = df, predictedRichness ~ observedRichness# + factor(res))
           ))


