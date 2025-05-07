### Here we add the necessary columns to turn the fielwork data into DWC acceptable tables.

library(sf)
library(dplyr)
library(reshape2)
library(igraph)
library(terra)

projCRS <- sf::st_crs(25833)$proj4string

# Import local functions
sapply(list.files("functions", pattern = "\\.R$", full.names = TRUE), source)


###-----------------------###
### Edit occurrence data ####
###-----------------------###

# Need to merge these two, creating locations IDs as we go
fieldworkSurveyCodeB <- read.csv("localPredictions/data/rawFieldData/fieldworkSurveyCode2.csv", sep = ";")
fieldworkSurveyCodeB$RuteID <- paste0(1:nrow(fieldworkSurveyCodeB), "B")
fieldworkSurveyCodeA <- read.csv("localPredictions/data/rawFieldData/fieldworkSurveyCode.csv", sep = ";") %>%
  dplyr::select(colnames(fieldworkSurveyCodeB))
fieldworkSurveyCodeA$RuteID <- paste0(fieldworkSurveyCodeA$RuteID, "A")

fieldworkSurveyCode <- rbind(fieldworkSurveyCodeA, fieldworkSurveyCodeB)

# Get rid of gaps in størrelse column
fieldworkSurveyCode$Størrelse <- gsub(" ","", fieldworkSurveyCode$Størrelse)

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

# Start adding relevant columns
fieldworkSurveyCode$fieldNumber <- fieldworkSurveyCode$RuteID
fieldworkSurveyCode$locality <- fieldworkSurveyCode$Område
fieldworkSurveyCode$language <- "en"
fieldworkSurveyCode$rightsHolder <- "UiB"
fieldworkSurveyCode$institutionID <- "https://ror.org/03zga2b32"
fieldworkSurveyCode$institutionCode <- "UiB"
fieldworkSurveyCode$datasetID <- rep(uuid::UUIDgenerate(), nrow(fieldworkSurveyCode))
fieldworkSurveyCode$ownerInstitutionCode <- "UiB"
fieldworkSurveyCode$basisOfRecord <- "humanObservation"
fieldworkSurveyCode$eventID <- uuid::UUIDgenerate(n = nrow(fieldworkSurveyCode))
fieldworkSurveyCode$eventDate <- fieldworkSurveyCode$date
fieldworkSurveyCode$samplingProtocol <- paste0(fieldworkSurveyCode$Størrelse, " quadrat")
fieldworkSurveyCode$verbatimLongitude <- fieldworkSurveyCode$Xdecimal
fieldworkSurveyCode$verbatimLatitude <- fieldworkSurveyCode$Ydecimal

# Get down to relevant columns
eventData <- fieldworkSurveyCode[,c("date", "month", "year", "day", "locality", "language",
                                                  "rightsHolder", "institutionID",
                                                 "institutionCode", "datasetID", "ownerInstitutionCode",
                                                 "basisOfRecord", "eventID", "eventDate", 
                                                 "samplingProtocol", "verbatimLongitude",
                                                 "verbatimLatitude", "fieldNumber")]

###-----------------------###
### Edit occurrence data ####
###-----------------------###

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

# Turn this into long format
fieldworkSurveyResultsLong <- melt(fieldworkSurveyResults,
             id.vars = c("RuteID"),
             variable.name = "species")

fieldworkSurveyResultsLong$individualCount <- fieldworkSurveyResultsLong$value
fieldworkSurveyResultsLong$occurrenceID <- uuid::UUIDgenerate(n = nrow(fieldworkSurveyResultsLong))
fieldworkSurveyResultsLong$scientificName <- gsub("_", " ", fieldworkSurveyResultsLong$species)
fieldworkSurveyResultsLong$fieldNumber <- fieldworkSurveyResultsLong$RuteID
fieldworkSurveyResultsLong$eventID <- eventData$eventID[match(fieldworkSurveyResultsLong$fieldNumber,
                                                              eventData$fieldNumber)]

occurrenceData <- fieldworkSurveyResultsLong[,c("fieldNumber", "scientificName", "occurrenceID",
                                                "individualCount", "eventID")]

write.table(eventData, file = "localPredictions/data/event.txt")
write.table(occurrenceData, file = "localPredictions/data/occurrence.txt")


