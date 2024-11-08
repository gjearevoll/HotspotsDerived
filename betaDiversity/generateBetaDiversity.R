### Import relevant libraries
library(terra)
library(dplyr)

# Import ansvarsarter
ansvarsArter <- readRDS("data/ansvarsArterList.RDS")

# Import full ansvartsarter species data at all sites
# Direct pc to find correct files
ogDirectory <- "../BioDivMapping/data/run_2024-10-11/modelOutputs/"
listProbFiles <- list.files(ogDirectory, full.names = TRUE)[grep("probability" ,list.files(ogDirectory))]
listProbFiles <- listProbFiles[-grep(".xml", listProbFiles)]

# Import probabilities for different species and cut down to ansvarsarter
probabilites <- rast(listProbFiles)
speciesToMap <- intersect(names(probabilites), ansvarsArter$simpleScientificName)
ansvarsarterProbs <- probabilites[[speciesToMap]]

# Convert everything above 0.7 to a 1 (this is where processing time starts to blow out)
speciesAssumed <- ifel(ansvarsarterProbs > 75, 1, 0)
localRichness <- sum(speciesAssumed)

# Get species richness in the 20 x 20km area surrounding each pixel (LONG PROCESSING TIME)
focalAttempt <- terra::focal(speciesAssumed, w = 41, fun = any, na.rm = TRUE)
focalRichness <- sum(focalAttempt)
focalBeta <- localRichness/focalRichness
writeRaster(focalBeta, "betaDiversity/data/betaDiversity.tiff")
