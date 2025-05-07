### Import libraries
library(terra)
library(dplyr)
library(sf)
focalTaxa <- read.csv("data/focalTaxa.csv")

dateUsed <- "2025-01-06"

###---------------------###
### Construct hotspots ####
###---------------------###

# Get table to merging all taxa together
finalFiles <- list.files(paste0("data/modelOutputs/run_",dateUsed), pattern = "final")
filePatterns <- unique(gsub("final.tiff", "",gsub(".*_","", finalFiles)))
allSIFiles <- list.files(paste0("data/modelOutputs/run_",dateUsed), pattern = "density")

norwayBorder2 <- vect(sf::read_sf("../BioDivMapping/data/external/norge_border/Noreg_polygon.shp"))

translations <- read.csv("data/taxaTranslations.csv")
translations$groupsNynorsk[translations$groups == "Insects"] <- "Insekt og edderkoppdyr"

firstlow <- function(x) {
  substr(x, 1, 1) <- tolower(substr(x, 1, 1))
  x
}

baseRaster <- rast(paste0("data/modelOutputs/run_", dateUsed, "/", finalFiles[1])) 
norwayBorderProjected <- terra::project(norwayBorder2, baseRaster)

# Get all available taxa
hotspotList <- list()
biasList <- list()
#for (taxa in unique(groupingDF$grouping)) {
for (taxa in filePatterns) {
  allRelevantFiles <- finalFiles[grepl(taxa, finalFiles)]
  allRelevantSIFiles <- allSIFiles[grep(taxa, allSIFiles)]
  
  types <- if (taxa %in% c("groundNesters", "waders", "woodpeckers")) "allspecies" else c("allspecies", "ansvarsarter", "threatenedspecies")
  
  for (type in types) {
    typeFiles <- allRelevantFiles[grep(type, allRelevantFiles)]
    if (length(typeFiles) == 0) {next}
    statsRaster <- rast(paste0("data/modelOutputs/run_", dateUsed, "/", typeFiles)) 
    if(ext(statsRaster) != ext(rast(baseRaster))) {
      statsRaster <- project(statsRaster, rast(baseRaster))
    }
    statsRaster <- crop(statsRaster, norwayBorderProjected, mask = T)

    # Create quantiles
    quantiles <- quantile(values(statsRaster[["skalertRikhet"]]), c(0.9, 0.95, 0.99), na.rm = TRUE)
    hotspotsRaster1 <- c(ifel(statsRaster$skalertRikhet > quantiles[1], TRUE, FALSE),
                         ifel(statsRaster$skalertRikhet > quantiles[2], TRUE, FALSE),
                         ifel(statsRaster$skalertRikhet > quantiles[3], TRUE, FALSE))
      
    names(hotspotsRaster1) <- c("HS3", "HS2", "HS1")
    print(paste0("Saving ", type, " for ", taxa))
    # writeRaster(hotspotsRaster1, paste0("processes/overlays/data/hotspots/", taxa, "_", type, "_rikhet.tiff"), 
    #             overwrite = TRUE)
    
    # Create Norwegian name and save for upload
    typeNorsk <- ifelse(type == "allspecies", "alle", ifelse(type == "ansvarsarter", "ansvar", "trua"))
    nameNorsk <- firstlow(gsub(" ", "-", translations$nynorsk[translations$engelsk == taxa]))
    rasterName <- paste0("data/forUpload/", typeNorsk, "-",nameNorsk,"-hotspots_norge_2025_32633.tiff")
    writeRaster(hotspotsRaster1, rasterName, overwrite = TRUE)
    
    # Add version to large raster for local use
    names(hotspotsRaster1) <-  paste0(taxa, "_", type, "_", c("HS3", "HS2", "HS1"))
    hotspotList <- c(hotspotList, hotspotsRaster1)
    
    # Import density file
    typeFilesSI <- allRelevantSIFiles[grep(type, allRelevantSIFiles)]
    statsRaster <- rast(paste0("data/modelOutputs/run_", dateUsed, "/", typeFilesSI)) 
    if(ext(statsRaster) != ext(rast(baseRaster))) {
      statsRaster <- project(statsRaster, rast(baseRaster))
    }
    
    # Save for upload
    print(paste0("Saving bias for ", taxa))
    SIRasterName <- paste0("data/forUpload/", typeNorsk, "-",nameNorsk,"-innsamlingsintensitet_norge_2025_32633.tiff")
    writeRaster(project(statsRaster, hotspotsRaster1), SIRasterName,
                overwrite = TRUE)
    
    # Rename and save for local use
    names(statsRaster) <-  paste0(taxa, "_", type, "_bias")
    biasList <- c(biasList, statsRaster)
  }
  
}
hotspots <- do.call(c, hotspotList)
writeRaster(hotspots, "data/allHotspots.tiff", overwrite = TRUE)

bias <- do.call(c, biasList)
writeRaster(bias, "data/allBiases.tiff", overwrite = TRUE)
