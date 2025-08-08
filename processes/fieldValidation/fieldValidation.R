# Import relevant libraries
library(sf)
library(terra)
library(dplyr)

args <- commandArgs(TRUE)

if (length(args) != 0) {
  # Set arguments
  runDate <- args[1]
  # Set the working directory
  setwd("~/HotspotsDerived")
}

#runDate <- "2025-01-06"

# Get total species richness for an area using the NTNU Vascular Plant Field Notes

if (!file.exists("processes/fieldValidation/data/GBIFDownload/event.txt")) {
  destFile <- "processes/fieldValidation/data/GBIFDownload.zip"
  download.file("https://gbif.vm.ntnu.no/archive.do?r=vascularplantfieldnotes", destFile, mode = "wb", timeout = 400)
  unzip(destFile, exdir = "processes/fieldValidation/data/GBIFDownload")
}

events <- read.delim("processes/fieldValidation/data/GBIFDownload/event.txt")
occurrences <- read.delim("processes/fieldValidation/data/GBIFDownload/occurrence.txt")

# Get all species surveyed
speciesSurveyed <- unique(occurrences$scientificName)

cat("Downloaded occurrence and event data\n")

# Get richness per event
eventsJoined <- merge(events[,c("eventID", "locationID")], occurrences, all.x = TRUE, by = "eventID") %>%
  filter(individualCount == 1)
eventsJoined <- eventsJoined[!duplicated(eventsJoined[,c("eventID","scientificName")]),]
eventsRichness <- tally(group_by(eventsJoined, eventID))

cat("Calculated richness per event\n")

### Import shapefile for vascular plants ###
# We're still not sure whether we'll make this publicly available #
plantShapeFile <- read_sf("processes/fieldValidation/data/VascularLocalities/VascularLocalities.shp") |>
  st_transform("EPSG:25833")

# Check that all events in GBIF have a corresponding location
eventsInGBIF <- unique(eventsJoined$locationID)
eventsInSHP <- st_drop_geometry(plantShapeFile$locationID)

# summary(eventsInGBIF %in% eventsInSHP)

# There's 1 missing, so we need to remove this
events <- events[events$locationID %in% eventsInSHP,]
events$richness <- eventsRichness$n[match(events$eventID, eventsRichness$eventID)]
eventsInGBIF <- unique(events$locationID)

cat("Got all relevant events as shapefiles\n")

# Now remove any locations that are not in our event data
summary(eventsInSHP %in% eventsInGBIF)
plantShapeFile <- vect(plantShapeFile[plantShapeFile$locationID %in% eventsInGBIF,])

# Find the area of each polygon
plantShapeFile$area <- expanse(plantShapeFile)/1000000
plantShapeFile$richness <- events$richness[match(plantShapeFile$locationID, events$locationID)]

# Take sample for practice
plantShapeFileSample <- plantShapeFile

# Get species richness for each area
probabilities <- rast(paste0("../BioDivMapping/data/run_", runDate,"/modelOutputs/processedOutputs/speciesprobability_vascularPlants.tiff"))
cat("Download species probabilities\n")

# Narrow down to only species surveyed
probabilities <- probabilities[[names(probabilities) %in% gsub(" ", "_", speciesSurveyed)]]

# Get probabilities megred and as fixed values
# probabilities <- probabilities[[1:5]]
extractedProbs <- predictedPresences <- extract(probabilities, plantShapeFile, fun = "mean", na.rm = TRUE)
predictedPresences[,-1] <- ifelse(predictedPresences[,-1] > 75, 1, 0)
cat("Calculated presence/absence in location\n")
plantShapeFileSample$predictedRichness <- rowSums(predictedPresences[,-1])
cat("Writing vector\n")
writeVector(plantShapeFileSample, "processes/fieldValidation/data/predictedRichness.shp", overwrite = TRUE)

# Get species per event
cat("Producing individual species maps\n")
speciesEvent1 <- eventsJoined[,c("locationID", "scientificName")]
speciesEvent1$occurrence <- 1
speciesEvent <- reshape(speciesEvent1, idvar = "locationID", timevar = "scientificName", direction = "wide")
colnames(speciesEvent) <- gsub("occurrence." , "",colnames(speciesEvent))
speciesEvent[is.na(speciesEvent)] <- 0

# Restrict only to events that are relevant
speciesEventRelevant <- speciesEvent[speciesEvent$locationID %in% eventsInGBIF,]

# Select a species
plantShapeFileSampleSpecies <- plantShapeFileSample[c("locationID", "area", "richness")]
for (sp in names(probabilities)) {
  plantShapeFileSampleSpecies[[paste0(sp, "_predicted")]] <- extractedProbs[,sp]
  speciesEventRelevant$focalSpecies <- speciesEventRelevant[,gsub("_", " ", sp)]
  plantShapeFileSampleSpecies[[paste0(sp, "_observed")]] <- speciesEventRelevant$focalSpecies[match(plantShapeFileSampleSpecies$locationID, speciesEventRelevant$locationID)]
}
writeVector(plantShapeFileSampleSpecies, "processes/fieldValidation/data/predictedPresences.shp", overwrite = TRUE)
saveRDS(names(plantShapeFileSampleSpecies), "processes/fieldValidation/data/indSpeciesLayerNames.RDS")

### It makes no sense to automate the running the script from here onwards, hence the STOP message below
stop("Automated section of script completed")

# Richness analysis
plantRichness <- vect("processes/fieldValidation/data/predictedRichness.shp")
plantRichness <- plantRichness[plantRichness$richness >= 0,]
plot(plantRichness$richness, plantRichness$predictedR)
plantRichnessDF <- as.data.frame(plantRichness)
ggplot(plantRichnessDF, aes(x=richness, y=predictedR)) + 
  labs(x = "Observed richness", y = "Predicted richness") +
  theme_bw() +
  geom_point(alpha = 0.3)

# Individual species analysis
plantPresences <- vect("processes/fieldValidation/data/predictedPresences.shp")
names(plantPresences) <- readRDS("processes/fieldValidation/data/indSpeciesLayerNames.RDS")

# Get species names
speciesNames <- gsub("_observed", "", grep("observed", names(plantPresences), value = T))

focalSpecies <- speciesNames[5]
plantInds <- plantPresences[[grep(focalSpecies, names(plantPresences))]]
colnames(plantInds) <- c("predicted", "observed")
plantInds <- plantInds[plantInds$observed %in% c(0,1),]
library(ggplot2)
plantInds$observed <- as.factor(plantInds$observed)
ggplot(plantInds, aes(x=observed, y=predicted)) + 
  geom_boxplot()

