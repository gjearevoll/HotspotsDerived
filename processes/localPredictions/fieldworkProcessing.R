### Need to process the fieldwork data


library(sf)
library(dplyr)
library(reshape2)
library(igraph)
library(terra)

projCRS <- sf::st_crs(25833)$proj4string

# Import local functions
sapply(list.files("functions", pattern = "\\.R$", full.names = TRUE), source)

# Start with on Arvid's data


# Need to merge these two, creating locations IDs as we go
fieldworkSurveyCodeB <- read.csv("localPredictions/data/rawFieldData/fieldworkSurveyCode2.csv", sep = ";")
fieldworkSurveyCodeB$RuteID <- paste0(1:nrow(fieldworkSurveyCodeB), "B")
fieldworkSurveyCodeA <- read.csv("localPredictions/data/rawFieldData/fieldworkSurveyCode.csv", sep = ";") %>%
  dplyr::select(colnames(fieldworkSurveyCodeB))
fieldworkSurveyCodeA$RuteID <- paste0(fieldworkSurveyCodeA$RuteID, "A")

# Since this has an inconsistency, we need to remove all 5m sites
fieldworkSurveyCodeA <- fieldworkSurveyCodeA[grep("1" ,fieldworkSurveyCodeA$Størrelse),]

fieldworkSurveyCode <- rbind(fieldworkSurveyCodeA, fieldworkSurveyCodeB)

# Isolate lines to fix
badCoords <- fieldworkSurveyCode$X != ""

# Sort out messed up lines of code
fieldworkSurveyCode$X <- ifelse(badCoords, fieldworkSurveyCode$X, sub('.*\n', '', fieldworkSurveyCode$Y))
fieldworkSurveyCode$Y <- ifelse(badCoords, fieldworkSurveyCode$Y, substring(fieldworkSurveyCode$Y,1,10))

# Now convert to decimal coordinates
fieldworkSurveyCode$Xdecimal <- as.numeric(gsub(",", ".", gsub("E0", "", fieldworkSurveyCode$X)))
fieldworkSurveyCode$Ydecimal <- as.numeric(gsub(",", ".", gsub("N0", "", fieldworkSurveyCode$Y)))

fieldworkSurveyCode$size <- as.integer(substring(fieldworkSurveyCode$Størrelse, 1, 1))

# Clean up date collections
fieldworkSurveyCode$date <- gsub("\\.", "/", fieldworkSurveyCode$Dato)
fieldworkSurveyCode$date <- as.Date(fieldworkSurveyCode$date, format = c("%d/%m/%Y"))
fieldworkSurveyCode$month <- format(fieldworkSurveyCode$date, "%m")
fieldworkSurveyCode$year <- "24"
fieldworkSurveyCode$day <- format(fieldworkSurveyCode$date, "%d")
fieldworkSurveyCode$date <- as.Date(paste(fieldworkSurveyCode$day, fieldworkSurveyCode$month,
                                          fieldworkSurveyCode$year, sep="-"), "%d-%m-%Y")

# Get down to relevant rows and columns
fieldworkSurveyCode2 <- fieldworkSurveyCode[,c("RuteID", "date", "Xdecimal", "Ydecimal")]
fieldworkSurveyCode2 <- fieldworkSurveyCode2[complete.cases(fieldworkSurveyCode2),]

# Now import our results and assign the same locationIDs
fieldworkSurveyResultsA <- read.csv("localPredictions/data/rawFieldData/fieldworkSurveyResults.csv", sep = ";")
fieldworkSurveyResultsA[is.na(fieldworkSurveyResultsA)] <- 0
fieldworkSurveyResultsA$RuteID <- paste0(fieldworkSurveyResultsA$RuteID, "A")

fieldworkSurveyResultsB <- read.csv("localPredictions/data/rawFieldData/fieldworkSurveyResults2.csv", sep = ";")
fieldworkSurveyResultsB[is.na(fieldworkSurveyResultsB)] <- 0
fieldworkSurveyResultsB$RuteID <- paste0(1:52, "B")

fieldworkSurveyResults <- bind_rows(fieldworkSurveyResultsA, fieldworkSurveyResultsB[,-1])
colnames(fieldworkSurveyResults) <- gsub("\\.", "_", colnames(fieldworkSurveyResults))

# All species not sampled for at one site need NA at the other
fieldworkSurveyResults[is.na(fieldworkSurveyResults)] <- 0


fieldworkSurvey <- merge(fieldworkSurveyCode2, fieldworkSurveyResults, all.x= TRUE, by = "RuteID")

# Convert to an sf object and norrow down to one set of observations
df1 <- st_as_sf(x = fieldworkSurvey,                         
               coords = c("Xdecimal", "Ydecimal"),
               crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
df1 <- st_transform(df1, crs = projCRS)

# And now move on to Ivar's data
load("localPredictions/data/rawFieldData/trondelagFieldDataSurvey.rda")
load("localPredictions/data/rawFieldData/trondelagFieldDataCode.rda")
colnames(Y_oekosystem_troendelag)[3:ncol(Y_oekosystem_troendelag)] <- gsub(" ","_",colnames(Y_oekosystem_troendelag)[3:ncol(Y_oekosystem_troendelag)])

mergedData <- merge(StudyDesign_ET, Y_oekosystem_troendelag, all.y = TRUE, by = "pointID")

df2 <- st_transform(mergedData, crs = projCRS)

# Prepare both datasets for merging
df1ForMerging <- df1 %>%
  mutate(survey = "UiB") %>%
  dplyr::select(-date) %>%
  rename(SiteID = RuteID)

df2ForMerging <- df2 %>%
  mutate(survey = "NTNU") %>%
  dplyr::select(-pointID, -siteID)

fieldworkSurveyResults <- bind_rows(df1ForMerging, df2ForMerging)
fieldworkSurveyResults[is.na(fieldworkSurveyResults)] <- 0


# We can also group the data into clusters now
adj <- st_distance(fieldworkSurveyResults)
adj <- matrix(as.numeric(as.numeric(adj)) < 25000, nrow = nrow(adj))

g <- graph_from_adjacency_matrix(adj)
fieldworkSurveyResults$group <- factor(components(g)$membership)

saveRDS(fieldworkSurveyResults, "localPredictions/data/fieldWorkResults.RDS")

### Now we also need to aggregate species richness to different resolutions


# Get the data into a more useable format
dfLong <- melt(fieldworkSurveyResults, id.vars = c("SiteID", "geometry", "group", "survey"), variable.name = "species",
                      value.name = "occurrences") %>%
  group_by(species, geometry, group, survey) %>%
  dplyr::summarise(summedOccurrences = sum(occurrences, na.rm = TRUE))
dfLong$presence <- ifelse(dfLong$summedOccurrences > 0, 1, 0)

# Check which species are present at all and keep only them
dfSpeciesToKeep <- summarise(group_by(dfLong, species), present = sum(presence))
speciesToKeep <- dfSpeciesToKeep$species[dfSpeciesToKeep$present > 0]
dfLong <- dfLong[dfLong$species %in% speciesToKeep,]

dfWide <- st_as_sf(as.data.frame(tidyr::spread(dfLong[,c("species", "geometry", "presence")], key = species, value = presence)))
dfWide$siteID <- paste0("site",1:nrow(dfWide))

# make a species table
surveyedSpecies <- data.frame(simpleSpecies = unique(dfLong$species))
# Match to accepted names
speciesBackbones <- get_gbif_backbone(surveyedSpecies$simpleSpecies)
surveyedSpecies$GBIFName <- speciesBackbones$scientificName
completeSpecies <- surveyedSpecies[complete.cases(surveyedSpecies),]

duplicatedSpecies <- completeSpecies[duplicated(completeSpecies$GBIFName),]

saveRDS(dfWide, "localPredictions/data/fieldWorkResultsWide.RDS")

