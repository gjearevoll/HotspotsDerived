#### Analyse overlap between diffferent layers ####

library(terra)
library(tidyterra)
library(ggtext)
library(ggplot2)
library(sf)


# Take species richness layers from rawMetrics
translations <- read.csv("data/taxaTranslations.csv")
translations$groupsNynorsk[translations$groups == "Insects"] <- "Insekt og edderkoppdyr\n"
titleTranslations <- data.frame(norsk = c("Alle artar", "Ansvarsartar", "Trua artar"),
                                english = c("allspecies", "ansvarsarter", "threatenedspecies"))

groupsToPlot <- c("Birds", "Insects", "Fungi", "Lichens", "Plants", "groundNesters", "woodpeckers", "waders")
groupsToPlot <- c("woodpeckers", "waders")

# Bring all relevant taxa together
finalFiles <- list.files("visualisations/figures/rawMetrics", pattern = "richness.tiff", full.names = TRUE)

norwayBorder2 <- sf::read_sf("data/norge_border/Noreg_polygon.shp")
norwayBorderProjected <- terra::project(vect(norwayBorder2), rast(finalFiles))

###-------------------------###
### Hotspot Visualisations ####
###-------------------------###

for (group in groupsToPlot) {
  
  if (group %in% c("waders", "groundNesters", "woodpeckers")) {
    nynorskName <- translations$nynorsk[translations$engelsk == group]
    typesToPlot <- "Alle artar"
  } else {
    nynorskName <- unique(translations$groupsNynorsk[translations$groups %in% group])
    typesToPlot <- c("Alle artar", "Ansvarsartar", "Trua artar")
  }
  
  availableLayers <- rast(finalFiles[grepl(group, finalFiles)]) |> setNames(typesToPlot)
  hotspots <- lapply(availableLayers, FUN = function(lyr) {
    quantiles <- quantile(values(lyr), c(0, 0.9, 0.95, 0.99, 1), na.rm = TRUE)
    categorisedRichness <- classify(lyr, 
                                    rcl=quantiles)
    levels(categorisedRichness) <- data.frame(ID=0:3, label=c(NA, "10%", "5%", "1%"))
    categorisedRichness
  }) |> setNames(names(availableLayers))
  hotspotsAgg <- rast(hotspots)
  
  cat("Saving raster for", group)
  writeRaster(hotspotsAgg, paste0("visualisations/figures/rawMetrics/",group,"_hotspots.tiff" ), overwrite = TRUE)
  
  
  # cropX <- c(600, 800)
  # cropY <- c(7600, 7800)
  # hotspotMap <- hotspots[[layerName]]
  # hotspotMap <- crop(hotspotMap, ext(c(cropX, cropY)))
  
  for (hp in seq_along(names(hotspots))) {
    hotspotType <- names(hotspots)[[hp]]
    hotspotMap <- hotspots[[hp]]
    hotspotsMap <- ggplot(norwayBorderProjected) + 
      geom_sf(fill = "beige", lwd = 0.1, colour = "black") + 
      #geom_sf(fill = "beige", lwd = 0, colour = "black") + 
      geom_spatraster(data = hotspotMap, mapping = aes(fill = label),
                      maxcell = 5e+07) +
      scale_fill_manual(values = c("orange", "red", "black"), na.value = NA, na.translate = F)+
      guides(
        fill = guide_legend(position = "inside")
      )  +
      #theme_bw() +
      theme_void() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(vjust = -10),
            legend.position.inside = c(0.8, 0.2)
      ) +
      labs(fill = "Hotspot\nstatus", title =paste0(nynorskName, " (",hotspotType , ")")
      ) +
      coord_sf(datum = st_crs(hotspots$`Alle artar`)
       #         , xlim = cropX, ylim =cropY
      )

      filenameToUse <- paste0("visualisations/figures/hotspots/hotspotsTaxa/hotspots_", hotspotType, "_",  nynorskName,".png")

    ggsave(filename =  filenameToUse,
           plot = hotspotsMap, units = "px",
           width = 1600, height = 1600
    )
  }
}

###----------------------------------------###
### Probability of Hotspot Visualisations ####
###----------------------------------------###

# Get masks ready
mask100 <- rast("data/mask100.tiff")  
mask100 <- crop(mask100, project(norwayBorderProjected, crs(mask100)), mask = T)

groupsToPlot <- c("vascularPlants", "birds", "fungi", "lichens", "insects", "woodpeckers", "groundNesters", "waders")
### Workflow for producing hotspots
for (group in groupsToPlot) {
  
  norskNavn <- translations$nynorsk[translations$engelsk == group]
  relevantFiles <- list.files("../BioDivMapping/data/run_2025-01-06/modelOutputs/processedOutputs", 
                              full.names = TRUE, pattern = "stats")
  relevantFiles <- relevantFiles[!grepl("final", relevantFiles)]
  relevantFiles <- relevantFiles[sapply(relevantFiles, FUN = function(x) {any(sapply(group, grepl, x))})]
  
  if (group %in% c("waders", "groundNesters", "woodpeckers")) {
    typesToPlot <- "allspecies"
  } else {
    typesToPlot <- c("allspecies", "threatenedspecies", "ansvarsarter")
  }
  
  
  for (type in typesToPlot) {
    rightStats <- rast(relevantFiles[grepl(type, relevantFiles)])
    nonScaledVersion <- rightStats
    # mask things
    nonScaledVersionCropped <- crop(nonScaledVersion, project(norwayBorderProjected, crs(mask100)), mask = T)
    if (ext(mask100) != ext(nonScaledVersionCropped)) {mask100 <- terra::project(mask100, nonScaledVersionCropped)}
    nonScaledVersionCropped$probability <- nonScaledVersionCropped$probability * mask100
    nonScaledVersionCropped$SD <- nonScaledVersionCropped$uncertainty * 10
    
    
    nonScaledVersionCropped$threshold <- quantile(values(nonScaledVersionCropped$probability), 0.9, na.rm = TRUE)
    normFu <- function(x,y,z) {1-pnorm(q=z, mean=x, sd=y)}
    probs <- lapp(nonScaledVersionCropped[[c(1,3,4)]], normFu)
    
    # # If cropped version necessary
    # cropX <- c(600, 800)
    # cropY <- c(7600, 7800)
    # cropX <- c(200, 400)
    # cropY <- c(7000, 7200)
    # probs <- crop(probs, ext(c(cropX, cropY)))
    
    cat("\nCreating plot for", type, group)
    norskTypeName <- unique(titleTranslations$norsk[titleTranslations$english %in% type])
    
    hotspotsProbMap <- ggplot(norwayBorderProjected) + 
      geom_spatraster(data = probs, mapping = aes(fill = lyr1),
                      maxcell = 5e+07) +
      geom_sf(fill = NA, lwd = 0.1, colour = "black") + 
      #geom_sf(fill = NA, lwd = 0, colour = "black") + 
      coord_sf(datum = st_crs(probs)
               #, xlim = cropX, ylim = cropY
               ) +
      scale_fill_grass_c(palette = "viridis", na.value = NA, limits = c(0, 1), direction = -1) +
      guides(
        fill = guide_colourbar(position = "inside")
      )  +
      theme_void() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(vjust = -10),
            legend.position.inside = c(0.8, 0.2)
      ) +
      labs(fill = "Hotspot\nsannsyn", title = paste0(norskNavn, " (",norskTypeName , ")")
      )
    ggsave(filename =  paste0("visualisations/figures/hotspots/hotspotprobs_", group, "_",  type,".png"),
           plot = hotspotsProbMap, units = "px",
           width = 1600, height = 1600)
  }
  
}
