### Import libraries
library(terra)
library(dplyr)
library(sf)
focalTaxa <- read.csv("data/focalTaxa.csv")


###---------------------###
### Construct hotspots ####
###---------------------###

# Get table to merging all taxa together
finalFiles <- list.files("data/modelOutputs", pattern = "final")
filePatterns <- unique(gsub("final.tiff", "",gsub(".*_","", finalFiles)))
allSIFiles <- list.files("data/modelOutputs", pattern = "bias")
groupingDF <- data.frame(speciesGroup = filePatterns, grouping = c("insekter", NA, "insekter", "fugler",
                                                                   "insekter", "insekter", "insekter", "sopp", NA,
                                                                   "insekter", "lav", NA, NA,
                                                                   "edderkopper", "karplanter",
                                                                   NA, NA))
groupingDF <- groupingDF[!is.na(groupingDF$grouping),]



# Get all available taxa
for (taxa in unique(groupingDF$grouping)) {
  allRelavantTaxa <- groupingDF$speciesGroup[groupingDF$grouping %in% taxa]
  allRelevantFiles <- finalFiles[grepl(paste0(allRelavantTaxa, collapse = "|"), finalFiles)]
  
  for (type in c("allspecies", "ansvarsarter", "threatenedspecies")) {
    typeFiles <- allRelevantFiles[grep(type, allRelevantFiles)]
    if (length(typeFiles) == 0) {next}
    combinedRaster <- rast(paste0("data/modelOutputs/", typeFiles)) 
    statsRaster <- sum(combinedRaster[[names(combinedRaster) == "skalertRikhet"]]) |> 
      setNames("rikhet")
    
    norwayBorder <- st_union(csmaps::nor_municip_map_b2024_default_sf)
    norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), statsGrouped[[1]]$stats))
    statsRaster <- crop(statsRaster, norwayBorderProjected, mask = T)
    
    quantiles <- quantile(values(statsRaster[["rikhet"]]), c(0, 0.9, 0.95, 0.99, 1), na.rm = TRUE)
    categorisedRichness <- classify(statsRaster[["rikhet"]], 
                                    rcl=quantiles)
    levels(categorisedRichness) <- data.frame(ID=0:3, label=c("1", "2", "3", "4"))
    names(categorisedRichness) <- paste(taxa, type, sep = "_")
    print(paste0("Saving ", type, " for ", taxa))
    writeRaster(categorisedRichness, paste0("overlays/data/hotspots/", taxa, "_", type, "_rikhet.tiff"), 
                overwrite = TRUE)
  
  }
  
 
  allRelevantSIFiles <- allSIFiles[grepl(paste0(allRelavantTaxa, collapse = "|"), allSIFiles)]
  combinedSIRaster <- rast(paste0("data/modelOutputs/", allRelevantSIFiles)) 
  statsRaster <- sum(combinedSIRaster[[names(combinedSIRaster) == "skalertInnsamlingsintensitet"]]) |> 
    setNames("innsamlingsIntensitet")
  quantiles <- quantile(values(statsRaster[["innsamlingsIntensitet"]]), c(0, 0.9, 0.95, 0.99, 1), na.rm = TRUE)
  categorisedSI <- classify(statsRaster[["innsamlingsIntensitet"]], 
                                  rcl=quantiles)
  levels(categorisedSI) <- data.frame(ID=0:3, label=c("1", "2", "3", "4"))
  names(categorisedSI) <- paste(taxa, sep = "_")
  print(paste0("Saving bias for ", taxa))
  writeRaster(categorisedSI, paste0("overlays/data/hotspots/", taxa, "_innsamlingsintensitet.tiff"), 
              overwrite = TRUE)
  
}

