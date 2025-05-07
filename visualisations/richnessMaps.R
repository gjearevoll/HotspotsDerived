### Visualisation for figures
library(tidyterra)
library(ggplot2)
library(ggtext)
library(terra)
library(sf)
library(dplyr)
library(terrainr)

###------------------###
### 1. Data loading ####
###------------------###

source("functions/defineRegion.R")
sourceDate <- "2025-01-06"
sourceDirectory <- paste0("data/modelOutputs/run_",sourceDate)


### Get names of all taxa we currently have finished data for
taxaNames <- unique(gsub(".tiff", "",
                         gsub("final.tiff", "", 
                              gsub("allspeciesstats_", "", 
                                   list.files(sourceDirectory, pattern = "allspecies")))))

translations <- read.csv("data/taxaTranslations.csv")
translations$groupsNynorsk[translations$engelsk == "insects"] <- "Insekt og edderkoppdyr\n"
titleTranslations <- data.frame(norsk = c("Alle artar", "Ansvarsartar", "Trua artar"),
                                english = c("allspecies", "ansvarsarter", "threatenedspecies"))

baseRaster <- rast(list.files(sourceDirectory, full.names = TRUE)[4])

norwayBorder2 <- sf::read_sf("../BioDivMapping/data/external/norge_border/Noreg_polygon.shp")
norwayBorderProjected <- terra::project(vect(norwayBorder2), baseRaster)

# Code for cropping:

# cropX <- c(600000, 800000)
# cropY <- c(7600000, 7800000)
# croppedSpecies <- crop(richnessTIFF, ext(c(cropX, cropY)))
# croppedBorder <- crop(project(norwayBorderProjected, croppedSpecies), ext(c(cropX, cropY)))

###--------------------###
### 1. Richness plots ####
###--------------------###

# Groups to plot
groupsToPlot <- unique(translations$groups[translations$engelsk %in% taxaNames])
groupsToPlot <- c("waders", "woodpeckers")

for (group in groupsToPlot) {
  
  if (group %in% c("waders", "groundNesters", "woodpeckers")) {
    typesToPlot <- "allspecies"
    relevantSpecies <- group
    nynorskName <- translations$nynorsk[translations$engelsk == group]
  } else {
    typesToPlot <- titleTranslations$english
    relevantSpecies <- translations$engelsk[translations$groups %in% group & translations$engelsk %in% taxaNames]
    nynorskName <- unique(translations$groupsNynorsk[translations$groups %in% group & translations$engelsk %in% taxaNames])
  }
  
  if (group == "Birds") {relevantSpecies <- "birds"}
  
  for (type in typesToPlot) {
    taxaNames <- unique(gsub(".tiff", "",
                             gsub("final.tiff", "", 
                                  gsub(paste0(type, "stats_"), "", 
                                       list.files(sourceDirectory, pattern = type)))))
    groupData <- rast(paste0(sourceDirectory, "/", type,"stats_", relevantSpecies, "final.tiff"))
    
    if (ext(groupData) != ext(baseRaster)) {
      groupData <- project(groupData, baseRaster, method = "mode")
    }
    
    richnessTIFF <- sum(groupData[[names(groupData) %in% "skalertRikhet"]])
    richnessTIFF$skalertRikhet <- richnessTIFF$sum/minmax(richnessTIFF$sum)[2]
    
    cat("\nSaving raster for", type, group)
    writeRaster(richnessTIFF$skalertRikhet, paste0("visualisations/figures/rawMetrics/",group,"_", type ,"richness.tiff" ), overwrite = TRUE)

    
    cat("\nCreating plot for", type, group)
    richnessPlot <- ggplot(norwayBorderProjected) +
      geom_spatraster(data = richnessTIFF, mapping = aes(fill = skalertRikhet),
                      maxcell = 5e+07) +
      geom_sf(fill = NA, lwd = 0.17, colour = "black") +
      coord_sf(datum = st_crs(richnessTIFF)) +
      scale_fill_grass_c(palette = "forest_cover", na.value = NA, limits = c(0, 1), direction = 1) +
      guides(
        fill = guide_colourbar(position = "inside")
      ) +
      theme_void() +
      #theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(vjust = -10),
            legend.position.inside = c(0.8, 0.2)
      ) +
      labs(fill = "Skalert\nartsrikdom", title = paste0(nynorskName, " (",titleTranslations$norsk[titleTranslations$english == type] , ")")
           )
    ggsave(richnessPlot, filename = paste0("visualisations/figures/rawMetrics/",group,"_", type ,"richness.png" ), units = "px",
           width = 1600, height = 1600
    )
  }
}


###-----------------------###
### 2. Uncertainty plots ####
###-----------------------###

for (group in groupsToPlot) {
  
  if (group %in% c("waders", "groundNesters", "woodpeckers")) {
    typesToPlot <- "allspecies"
    relevantSpecies <- group
    nynorskName <- translations$nynorsk[translations$engelsk == group]
  } else {
    typesToPlot <- titleTranslations$english
    relevantSpecies <- translations$engelsk[translations$groups %in% group & translations$engelsk %in% taxaNames]
    nynorskName <- unique(translations$groupsNynorsk[translations$groups %in% group & translations$engelsk %in% taxaNames])
  }
  
  if (group == "Birds") {relevantSpecies <- "birds"}
  
  for (type in typesToPlot) {
    
    taxaNames <- unique(gsub(".tiff", "",
                             gsub("final.tiff", "", 
                                  gsub(paste0(type,"stats_"), "", 
                                       list.files(sourceDirectory, pattern = type)))))
    groupData <- rast(paste0(sourceDirectory, "/", type,"stats_", relevantSpecies, "final.tiff"))
    
    richnessTIFF <- sqrt(sum((groupData[[names(groupData) %in% "skalertUsikkerhet"]])^2))
    richnessTIFF$skalertUsikkerhet <- richnessTIFF$sum/minmax(richnessTIFF$sum)[2]
    
    # Need to truncate usikkerhet
    qg <- global(richnessTIFF$skalertUsikkerhet, quantile, probs=c(0.1, 0.9), na.rm=TRUE)
    richnessTIFF$truncatedScale <- ifel(richnessTIFF$skalertUsikkerhet > qg[1,2], qg[1,2], richnessTIFF$skalertUsikkerhet)
    richnessTIFF$truncatedScale <- ifel(richnessTIFF$truncatedScale < qg[1,1], qg[1,1], richnessTIFF$truncatedScale)
    richnessTIFF$newScale <- (richnessTIFF$truncatedScale - minmax(richnessTIFF$truncatedScale)[1])/
      (minmax(richnessTIFF$truncatedScale)[2] - minmax(richnessTIFF$truncatedScale)[1])
    
    cat("\nSaving raster for", group)
    writeRaster(richnessTIFF[[c("sum", "newScale")]], paste0("visualisations/figures/rawMetrics/",group, "_", type, "_usikkerhet.tiff" ), overwrite = TRUE)

    plotTitle <- paste0(nynorskName, " (",titleTranslations$norsk[titleTranslations$english == type] , ")")
    
    richnessPlot <- ggplot(norwayBorderProjected) + 
      geom_spatraster(data = richnessTIFF, mapping = aes(fill = newScale),
                      maxcell = 5e+07) +
      geom_sf(fill = NA, lwd = 0.17, colour = "black") + 
      coord_sf(datum = st_crs(richnessTIFF)) +
      scale_fill_viridis_c(option = "mako", na.value = NA, limits = c(0, 1), direction = -1) +
      guides(
        fill = guide_colourbar(position = "inside")) +
      theme_void() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(vjust = -10),
            legend.position.inside = c(0.8, 0.2)
      ) +
      labs(fill = "Skalert\nusikkerheit", title = plotTitle)
    ggsave(richnessPlot, filename = paste0("visualisations/figures/rawMetrics/",group, "_", type, "_usikkerhet.png" ), units = "px",
           width = 1600, height = 1600
    )
  }
}



###----------------------------###
### 3. Sampling density plots ####
###----------------------------###

for (group in groupsToPlot) {
  
  if (group %in% c("waders", "groundNesters", "woodpeckers")) {
    typesToPlot <- "allspecies"
    relevantSpecies <- group
    nynorskName <- translations$nynorsk[translations$engelsk == group]
  } else {
    typesToPlot <- titleTranslations$english
    relevantSpecies <- translations$engelsk[translations$groups %in% group & translations$engelsk %in% taxaNames]
    nynorskName <- unique(translations$groupsNynorsk[translations$groups %in% group & translations$engelsk %in% taxaNames])
  }
  
  for (type in typesToPlot) {
    
    groupData <- rast(paste0(sourceDirectory, "/density_", type, "_",relevantSpecies,".tiff"))

    densityTIFF <- sum(groupData)
    densityTIFF$scaledDensity <- densityTIFF$sum/minmax(densityTIFF$sum)[2]
    
    cat("Saving density raster for", group)
    writeRaster(densityTIFF$scaledDensity, paste0("visualisations/figures/rawMetrics/",group,"_",type,"_innsamlingsintensitet.tiff" ), overwrite = TRUE)
    
    
    densityPlot <- ggplot(norwayBorderProjected) + 
      geom_spatraster(data = densityTIFF, mapping = aes(fill = scaledDensity),
                      maxcell = 5e+07) +
      geom_sf(fill = NA, lwd = 0.1, colour = "black") + 
      coord_sf(datum = st_crs(densityTIFF)) +
      scale_fill_grass_c(palette = "haxby", na.value = NA, limits = c(0, 1), direction = -1) +
      guides(
        fill = guide_colourbar(position = "inside")
      ) +
      theme_void() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(vjust = -10),
            legend.position.inside = c(0.8, 0.2)
      ) +
      labs(fill = "Skalert\ninnsamlingstettleik", title = paste0(nynorskName, " (",titleTranslations$norsk[titleTranslations$english == type] , ")")
      )
    ggsave(densityPlot, filename = paste0("visualisations/figures/rawMetrics/",group,"_",type,"innsamlingstettleik.png" ), units = "px",
           width = 1600, height = 1600
    )
  }
}



###----------------------------------###
### 4. Sampling density in hotspots ####
###----------------------------------###

for (group in groupsToPlot) {
  
  if (group %in% c("waders", "groundNesters", "woodpeckers")) {
    typesToPlot <- "allspecies"
    relevantSpecies <- group
    nynorskName <- translations$nynorsk[translations$engelsk == group]
  } else {
    typesToPlot <- titleTranslations$english
    relevantSpecies <- translations$engelsk[translations$groups %in% group & translations$engelsk %in% taxaNames]
    nynorskName <- unique(translations$groupsNynorsk[translations$groups %in% group & translations$engelsk %in% taxaNames])
  }
  
  richnessTiffs <- list.files("visualisations/figures/rawMetrics", pattern = "richness.tiff")
  richnessTiffs <- rast(paste0("visualisations/figures/rawMetrics/", richnessTiffs[grepl(group, substring(richnessTiffs,1, 13))])) |>
    setNames(typesToPlot)
  densityTiff <- list.files("visualisations/figures/rawMetrics", pattern = "innsamlingsintensitet.tiff")
  densityTiff <- rast(paste0("visualisations/figures/rawMetrics/", densityTiff[grepl(group, densityTiff)])) |>
    setNames(typesToPlot)
  
  cat("Saving raster for", group)
  
  if (ext(densityTiff) != ext(richnessTiffs)) {
    densityTiff <- terra::project(densityTiff, richnessTiffs, method = "mode")
  }
  
  for (type in typesToPlot) {
    combRaster <- c(densityTiff[[type]], richnessTiffs[[type]]) |> setNames(c("skalertInnsamlingintensitet", "skalertRikhet"))
    hotspotThreshold <- global(combRaster$skalertRikhet, quantile, probs=0.9, na.rm=TRUE)[1,1]
    combRaster$maskedDensity <- ifel(combRaster$skalertRikhet >= hotspotThreshold, combRaster$skalertInnsamlingintensitet, NA)
    
    densityPlot <- ggplot(norwayBorderProjected) + 
      geom_sf(fill = "lightgrey", lwd = 0.07, colour = "lightgrey") + 
      geom_spatraster(data = combRaster, mapping = aes(fill = maskedDensity),
                      maxcell = 5e+07) +
      
      geom_sf(data = norwayBorderProjected, fill = NA, lwd = 0.08, colour = "lightgrey") + 
      coord_sf(datum = st_crs(combRaster)) +
      scale_fill_gradient2(low = "black", mid = "darkgreen", high = "lightgreen", midpoint = 0.5, limits = c(0, 1), na.value = NA) +
      #scale_fill_viridis_c(option = "magma", na.value = NA, limits = c(0, 1), direction = -1) +
      #theme_bw() +
      guides(
        fill = guide_colourbar(position = "inside")
      ) +
      theme_void() +
      theme(plot.title = element_text(vjust = -10),
            legend.position.inside = c(0.8, 0.2)
      ) +
      labs(fill = "Maskert\ninnsamlingstettleik", title = paste0(nynorskName, " (", 
                                                                 titleTranslations$norsk[titleTranslations$english == type], ")")
      )
    ggsave(densityPlot, filename = paste0("visualisations/figures/hotspots/hotspotsTaxa/",group,"_",type,"hotspotsinnsamlingstettleik.png" ), units = "px",
           width = 1600, height = 1600
    )
  }
  
}



###--------------------------------------------------------###
### 5. Get hotspot size in and outside inngrepsfrie areas ####
###--------------------------------------------------------###

###-----------------------###
### Inngrepsfrie omrader ####
###-----------------------###

# This one is also available online from Miljodirektoratet:
# https://kartkatalog.miljodirektoratet.no/Dataset/Details/100
inngrepsfrieOmrader <- read_sf("processes/overlays/data/inngrepsfrieOmrader/statusPolygon.shp")
ifVect <- terra::project(vect(inngrepsfrieOmrader), baseRaster)
ifRaster <- terra::rasterize(ifVect, baseRaster, field = "vsone")
translationTable <- data.frame(ints = 0:3, categories = c("Sone1", "Sone2", NA, "Vill"))
levels(ifRaster) <- translationTable
names(ifRaster) <- "intactAreas"


for (group in groupsToPlot) {
  nynorskName <- unique(translations$groupsNynorsk[translations$groups %in% group & translations$engelsk %in% taxaNames])
  
  richnessTiffs <- list.files("visualisations/figures/rawMetrics", pattern = "richness.tiff")
  richnessTiffs <- rast(paste0("visualisations/figures/rawMetrics/", richnessTiffs[grepl(group, richnessTiffs)]))[[1]]
  
  hotspotThreshold <- global(richnessTiffs$skalertRikhet, quantile, probs=0.9, na.rm=TRUE)[1,1]
  richnessTiffs$hotspot <- ifel(richnessTiffs$skalertRikhet >= hotspotThreshold, 1, NA)
  
  hotspots <- patches(richnessTiffs$hotspot, directions=8, zeroAsNA=TRUE) 
  rz <- zonal(cellSize(hotspots, unit="km"), hotspots, sum)
  rz$type <- "Alle hotspots"
  
  # Inngrepsfrie hotspots
  richnessTiffs$hotspotsVill <- ifel(ifRaster$intactAreas == "Vill", richnessTiffs$hotspot, NA)
  hotspotsVill <- patches(richnessTiffs$hotspotsVill, directions=8, zeroAsNA=TRUE) 
  rzVill <- zonal(cellSize(hotspotsVill, unit="km"), hotspotsVill, sum) 
  rzVill$type <- "Villmarksprega\nhotspots"
  
  richnessTiffs$hotspotsInngrepsfrie <- ifel(ifRaster$intactAreas %in% "Sone1", richnessTiffs$hotspot, NA)
  hotspotsInngrepsfrie <- patches(richnessTiffs$hotspotsInngrepsfrie, directions=8, zeroAsNA=TRUE) 
  rzInngrespfrie <- zonal(cellSize(hotspotsInngrepsfrie, unit="km"), hotspotsInngrepsfrie, sum) 
  rzInngrespfrie$type <- "Villmarksprega\nhotspots"
  
  combinedSizes <- rbind(rz, rzVill)
  combinedSizes2 <- rbind(rz, rzInngrespfrie)
  
  histPlot <- ggplot(combinedSizes, aes(area, fill = type)) + 
    geom_histogram(position = 'identity', bins = 20, colour = "#333333", width = 0.7) + xlab("Areal (kvadratkilometer)") +
    scale_fill_manual(values=c("#99FF99", "#006600")) +
    ylab("") + labs(title = nynorskName, fill = "") +
    theme_minimal() +
    scale_x_continuous(trans="log",breaks = scales::trans_breaks("log", function(x) round(exp(x))))
  
  ggsave(histPlot, filename = paste0("visualisations/figures/rawMetrics/",group,"_villprgeaHotspotAreal.png" ), units = "px",
         width = 2200, height = 1400
  )
  
  histPlot2 <- ggplot(combinedSizes2, aes(area, fill = type)) + 
    geom_histogram(position = 'identity', bins = 20, colour = "#333333", width = 0.7) + xlab("Areal (kvadratkilometer)") +
    scale_fill_manual(values=c("#CCCC33", "#999933")) +
    ylab("") + labs(title = nynorskName, fill = "") +
    theme_minimal() +
    scale_x_continuous(trans="log",breaks = scales::trans_breaks("log", function(x) round(exp(x))))
  
  ggsave(histPlot2, filename = paste0("visualisations/figures/rawMetrics/",group,"_villprgeaHotspotAreal2.png" ), units = "px",
         width = 2200, height = 1400
  )
  
  
}





###------------------------###
### Supplementary figures ####
###------------------------###

# Here are maps that are useful to have easy access to but are not necessary for the report


###-------------------------###
### S1. Individual species ####
###-------------------------###

# Show individual species
speciesGroup <- "vascularPlants"
indSpeciesTIFF <- rast(paste0("../BioDivMapping/data/run_2025-01-06/modelOutputs/processedOutputs/speciesprobability_", 
                              speciesGroup, ".tiff"))

availableSpecies <- names(indSpeciesTIFF)

speciesToPlot <- "Antennaria_porsildii"

speciesPredictions <- indSpeciesTIFF[[speciesToPlot]]
speciesPredictions$scaledSpecies <- speciesPredictions/100
speciesPredictions <- crop(speciesPredictions, norwayBorderProjected, mask = T)

# 628 by 492
ggplot(norwayBorderProjected) + 
  geom_spatraster(data = speciesPredictions, mapping = aes(fill = scaledSpecies),
                  maxcell = 5e+07) +
  geom_sf(fill = NA, lwd = 0.1, colour = "black") + 
  #geom_sf(data = speciesObservations, fill = NA, colour = "black") + 
  coord_sf(datum = st_crs(speciesPredictions)) +
  scale_fill_grass_c(palette = "grass", na.value = NA, limits = c(0, 1), direction = -1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.subtitle = ggtext::element_markdown(),
        plot.title = element_text()
  ) +
  labs(fill = "Forekomst\nsannsyn", title = gsub("_", " ", speciesToPlot))
