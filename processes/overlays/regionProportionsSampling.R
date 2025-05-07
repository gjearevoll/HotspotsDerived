
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###

library(terra)
library(dplyr)
library(sf)
library(exactextractr)

baseRaster <- rast("data/modelOutputs/run_2025-01-06/density_allspecies_insects.tiff")

sourceDir <- "data/modelOutputs/run_2025-01-06"
samplingFiles <- list.files(sourceDir, pattern = "density")
densities <- do.call(c, lapply(paste(sourceDir, samplingFiles, sep = "/"), FUN = function(g) {
  project(rast(g), baseRaster)
  
  })) |> 
  setNames(gsub(".tiff","",samplingFiles))

translationTable <- read.csv("data/taxaTranslations.csv")[c(15,18:23,25),]


###---------------###
#### WATER AREAS ####
###---------------###

waterAreas <- st_read("processes/overlays/data/waterRegions/Vannomrader_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)

#extract the area of each raster cell covered by the plot and summarize
waterAreaProportions <- lapply(seq_along(samplingFiles), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- strsplit(names(densities)[x], "_")[[1]][3]
  speciesName2 <- translationTable$nynorsk[translationTable$engelsk == speciesName]
  grouping <- ifelse(grepl("allspecies", strsplit(names(densities)[x], "_")[[1]][2]), "Alle ", 
                     ifelse(grepl("ansvar", strsplit(names(densities)[x], "_")[[1]][2]), "Ansvar ", "Trua ")) 
  
  # Get average sampling density in each region
  table1 <- terra::extract(densities[[x]], waterAreas)
  colnames(table1)[2] <- "density"
  matchingTable <- data.frame(navn = waterAreas$navn, id = 1:nrow(waterAreas))
  table1$navn <- matchingTable$navn[match(table1$ID, matchingTable$id)]
  
  table2 <- table1 %>%
    group_by(navn) %>%
    summarise(samplingIntensity = mean(density,  na.rm = TRUE))
  table2$grouping <- paste0(grouping, speciesName2)
  table2
}
)

# Collate data frame for all species groups
asDF <- as.data.frame(do.call(rbind, waterAreaProportions))
asDF$samplingIntensity <- round(asDF$samplingIntensity,4)
widVersionWaterAreas <- reshape(asDF, idvar = "navn", timevar = "grouping", direction = "wide")
colnames(widVersionWaterAreas) <- gsub("samplingIntensity.", "", colnames(widVersionWaterAreas))
widVersionWaterAreas[is.na(widVersionWaterAreas)] <- 0
colnames(widVersionWaterAreas)[1] <- "Vannområde"

write.csv(widVersionWaterAreas, "processes/overlays/data/waterAreaSamplingDensities.csv")

###-----------------###
#### WATER REGIONS ####
###-----------------###

waterRegions <- st_read("processes/overlays/data/waterRegions/Vannregioner_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)
waterRegions <- waterRegions[waterRegions$navn != "Jan Mayen",]

#extract the area of each raster cell covered by the plot and summarize
waterRegionsProportions <- lapply(seq_along(samplingFiles), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- strsplit(names(densities)[x], "_")[[1]][3]
  speciesName2 <- translationTable$nynorsk[translationTable$engelsk == speciesName]
  grouping <- ifelse(grepl("allspecies", strsplit(names(densities)[x], "_")[[1]][2]), "Alle ", 
                     ifelse(grepl("ansvar", strsplit(names(densities)[x], "_")[[1]][2]), "Ansvar ", "Trua ")) 
  
  # Get average sampling density in each region
  table1 <- terra::extract(densities[[x]], waterRegions)
  colnames(table1)[2] <- "density"
  matchingTable <- data.frame(navn = waterAreas$navn, id = 1:nrow(waterAreas))
  table1$navn <- matchingTable$navn[match(table1$ID, matchingTable$id)]
  
  table2 <- table1 %>%
    group_by(navn) %>%
    summarise(samplingIntensity = mean(density,  na.rm = TRUE))
  table2$grouping <- paste0(grouping, speciesName2)
  table2
}
)

# Collate data frame for all species groups
asDF <- as.data.frame(do.call(rbind, waterRegionsProportions))
asDF$samplingIntensity <- round(asDF$samplingIntensity,4)
widVersionWaterRegions <- reshape(asDF, idvar = "navn", timevar = "grouping", direction = "wide")
colnames(widVersionWaterRegions) <- gsub("samplingIntensity.", "", colnames(widVersionWaterRegions))
widVersionWaterRegions[is.na(widVersionWaterRegions)] <- 0
colnames(widVersionWaterRegions)[1] <- "Vannregion"

write.csv(widVersionWaterRegions, "processes/overlays/data/waterRegionSamplingDensities.csv")

###------------------###
#### MUNICIPALITIES ####
###------------------###

municipalities <- csmaps::nor_municip_map_b2024_default_sf
# municipalityRaster <- terra::rasterize(municipalitiesVect, hotspotLayers, field = "location_code")
# names(municipalityRaster) <- "municipalities"
# counties <- csmaps::nor_county_map_b2024_default_sf
# countiesVect <- terra::project(vect(counties), hotspotLayers)

#extract the area of each raster cell covered by the plot and summarize
municipalityProportions <- lapply(seq_along(samplingFiles), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- strsplit(names(densities)[x], "_")[[1]][3]
  speciesName2 <- translationTable$nynorsk[translationTable$engelsk == speciesName]
  grouping <- ifelse(grepl("allspecies", strsplit(names(densities)[x], "_")[[1]][2]), "Alle ", 
                     ifelse(grepl("ansvar", strsplit(names(densities)[x], "_")[[1]][2]), "Ansvar ", "Trua ")) 
  
  # Get average sampling density in each region
  table1 <- terra::extract(densities[[x]], municipalities)
  colnames(table1)[2] <- "density"
  matchingTable <- data.frame(navn = municipalities$location_code, id = 1:nrow(municipalities))
  table1$navn <- matchingTable$navn[match(table1$ID, matchingTable$id)]
  
  table2 <- table1 %>%
    group_by(navn) %>%
    summarise(samplingIntensity = mean(density,  na.rm = TRUE))
  table2$grouping <- paste0(grouping, speciesName2)
  table2
}
)

# Collate data frame for all species groups
asDF <- as.data.frame(do.call(rbind, municipalityProportions))
asDF$samplingIntensity <- round(asDF$samplingIntensity,4)
wideMunicipalities <- reshape(asDF, idvar = "navn", timevar = "grouping", direction = "wide")
colnames(wideMunicipalities) <- gsub("samplingIntensity.", "", colnames(wideMunicipalities))
wideMunicipalities[is.na(wideMunicipalities)] <- 0
colnames(wideMunicipalities)[1] <- "Kommune"

write.csv(wideMunicipalities, "processes/overlays/data/municipalitySamplingDensities.csv")

###------------###
#### COUNTIES ####
###------------###

counties <- csmaps::nor_county_map_b2024_default_sf

#extract the area of each raster cell covered by the plot and summarize
countyProportions <- lapply(seq_along(samplingFiles), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- strsplit(names(densities)[x], "_")[[1]][3]
  speciesName2 <- translationTable$nynorsk[translationTable$engelsk == speciesName]
  grouping <- ifelse(grepl("allspecies", strsplit(names(densities)[x], "_")[[1]][2]), "Alle ", 
                     ifelse(grepl("ansvar", strsplit(names(densities)[x], "_")[[1]][2]), "Ansvar ", "Trua ")) 
  
  # Get average sampling density in each region
  table1 <- terra::extract(densities[[x]], counties)
  colnames(table1)[2] <- "density"
  matchingTable <- data.frame(navn = counties$location_code, id = 1:nrow(counties))
  table1$navn <- matchingTable$navn[match(table1$ID, matchingTable$id)]
  
  table2 <- table1 %>%
    group_by(navn) %>%
    summarise(samplingIntensity = mean(density,  na.rm = TRUE))
  table2$grouping <- paste0(grouping, speciesName2)
  table2
}
)

# Collate data frame for all species groups
asDF <- as.data.frame(do.call(rbind, countyProportions))
asDF$samplingIntensity <- round(asDF$samplingIntensity,4)
wideCounties <- reshape(asDF, idvar = "navn", timevar = "grouping", direction = "wide")
colnames(wideCounties) <- gsub("samplingIntensity.", "", colnames(wideCounties))
wideCounties[is.na(wideCounties)] <- 0
colnames(wideCounties)[1] <- "Fylke"

write.csv(wideCounties, "processes/overlays/data/countySamplingDensities.csv")
