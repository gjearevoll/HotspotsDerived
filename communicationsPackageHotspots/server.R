#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(shinydashboard)
library(shinyjs)
library(sf)
library(ggplot2)
library(dplyr)
library(terra)
library(tidyterra)

terraOptions(memmax = 0.9)

# Also need to mask out areas outside of Norway
norwayBorder <- st_union(csmaps::nor_municip_map_b2024_default_sf)

# # Import dropdown list for taxa
ddList <- readRDS("data/importedData/ddList.RDS")
namesFromGroups <- readRDS("data/importedData/namesFromGroups.RDS")

creditList <- readRDS("data/importedData/imageCredit.RDS")
ansvarsarter <- readRDS("data/ansvarsArterList.RDS")

hotspotListAll <- mappedData <- rast("data/importedData/allHotspots.tiff")

# projCRS <- sf::st_crs(25833)$proj4string
# 
choiceList <- c("Norway", "50", "54", "46", "34") |> 
  setNames(c("Norge", "Trøndelag", "Finnmark", "Vestland", "Innlandet"))
taxaChoiceList <- c("insects", "birds", "vascularPlants", "lichens", "fungi") |> 
  setNames(c("Insekt", "Fugler", "Karplanter", "Lav", "Sopp"))

source("data/textFunctions.R")
source("data/defineRegion.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  # Input of species for intensity figure
  observe({
    updateSelectInput(inputId = "species", choices = ddList[[input$taxaIndivid]])
  })
  
  # Species intensity map
  output$speciesMap <- renderPlot({
    regionToUse <- ifelse(nchar(input$region) == 2, "county", "country")
    regionName <- names(choiceList)[choiceList == input$region]
    regionGeometry <- defineRegion(regionToUse, input$region)
    
    dataToUse <- input$speciesGroup
    taxaToUse <- input$taxa
    metricToUse <- input$statistic
    
    taxaName <- names(taxaChoiceList)[taxaChoiceList == taxaToUse]
    
    focalData <- rast(paste0("data/importedData/", dataToUse, "stats_", taxaToUse, "final.tiff"))
    focalData$layer <- focalData[[metricToUse]]
    
    if (metricToUse == "skalertUsikkerhet") {
      qg <- global(focalData$layer, quantile, probs=c(0.1, 0.9), na.rm=TRUE)
      truncatedScale <- ifel(focalData$layer > qg[1,2], qg[1,2], focalData$layer)
      truncatedScale <- ifel(truncatedScale < qg[1,1], qg[1,1], truncatedScale)
      focalData$layer <- (truncatedScale - minmax(truncatedScale)[1])/
        (minmax(truncatedScale)[2] - minmax(truncatedScale)[1])
    }
    
    if (regionToUse == "country") {
      norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), focalData))
      lwdF = 0.1
    } else {
      norwayBorderProjected <- terra::project(vect(regionGeometry), crs(focalData))
      focalData <- crop(focalData, norwayBorderProjected, mask = T)
      lwdF = 0
    }
    
    if (input$statistic == "skalertRikhet") {
      scaledPar <- scale_fill_grass_c(palette = "forest_cover", na.value = NA, limits = c(0, 1), direction = 1)
    } else {
      scaledPar <- scale_fill_viridis_c(option = "mako", na.value = NA, limits = c(0, 1), direction = -1)
    }
    
    returnPlot <- ggplot(norwayBorderProjected) + 
      geom_spatraster(data = focalData, mapping = aes(fill = layer)) +
      geom_sf(fill = NA, lwd = lwdF, colour = "black") + 
      coord_sf(datum = st_crs(focalData)) +
      scaledPar +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.subtitle = ggtext::element_markdown(),
            plot.title = ggtext::element_textbox_simple(
              size = 16, lineheight = 1, linetype = 1, padding = margin(5, 5, 5, 5), # background fill color
              r = grid::unit(3, "pt"), # radius for rounded corners,
              margin = margin(5,5,15,2)),
      ) +
      labs(
        title = paste0("Artsrikhet<br><span style = 'font-size:10pt;'>Artsrikhet per 10 x 10 km piksel for ", 
                       taxaName,"</span>"),
        fill = ifelse(metricToUse == "skalertRikhet", "Skalert\nartsrikhet", "Skalert\nusikkerhet")
      )
    
    return(returnPlot)
  }, height = 500)
  
  # Text box for beta diversity tab
  output$betaMap <- renderPlot({
    regionToUse <- ifelse(nchar(input$regionBeta) == 2, "county", "country")
    regionName <- names(choiceList)[choiceList == input$regionBeta]
    regionGeometry <- defineRegion(regionToUse, input$regionBeta)
    taxaToUse <- input$taxaBeta
    taxaName <- names(taxaChoiceList)[taxaChoiceList == taxaToUse]
    
    mappedData <- rast(paste0("data/importedData/",taxaToUse, "10betaDiversity.tiff"))
    
    if (regionToUse == "country") {
      norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), mappedData))
      mappedData <- crop(mappedData, norwayBorderProjected, mask = T)
      lwdF = 0.1
    } else {
      norwayBorderProjected <- terra::project(vect(regionGeometry), crs(mappedData))
      mappedData <- crop(mappedData, norwayBorderProjected, mask = T)
      lwdF = 0
    }
    
    
    plotToReturn <- ggplot(norwayBorderProjected) +  
      geom_spatraster(data = mappedData,
                      aes(fill = sum)) +
      scale_fill_grass_c(palette = "forest_cover", na.value = NA, limits = c(0, 1), direction = 1) +
      geom_sf(fill = NA, lwd = lwdF, colour = "black") +
      theme_bw() +
      coord_sf(datum = st_crs(mappedData)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.subtitle = ggtext::element_markdown(),
            plot.title = ggtext::element_textbox_simple(
              size = 16, lineheight = 1, linetype = 1, padding = margin(5, 5, 5, 5), # background fill color
              r = grid::unit(3, "pt"), # radius for rounded corners,
              margin = margin(5,5,15,2)),
      ) +
      labs(
        title = paste0("Betadiversitet<br><span style = 'font-size:10pt;'>Andel av alle arter funnet i regional område som er funnet i piksel for ", 
                       taxaName," i ", regionName,".</span>"),
        fill = "Betadiversitet"
      )
    return(plotToReturn)
  }, height = 500)
  
  
  # Text box for species intensity tab
  output$individuals <- renderPlot({
    regionToUse <- ifelse(nchar(input$regionIndivid) == 2, "county", "country")
    regionName <- names(choiceList)[choiceList == input$regionIndivid]
    regionGeometry <- defineRegion(regionToUse, input$regionIndivid)
    taxaName <- input$taxaIndivid
    speciesSimpleName <- namesFromGroups$species[namesFromGroups$fullName == input$species]
    
    mappedData <- rast(paste0("data/importedData/",input$taxaIndivid, "_", input$statisticIndivid,".tiff"))[[speciesSimpleName]]
    mappedData$layer <- mappedData[[1]]
    
    if (input$statisticIndivid == "uncertainty") {
      mappedData$layer <- 100 * (mappedData$layer - minmax(mappedData$layer)[1])/(minmax(mappedData$layer)[2] - minmax(mappedData$layer)[1])
      fillName <- "Usikkerhet"
      palDir <- 1
    } else {fillName <- "Forekomst\nsannsynlighet"
    palDir <- -1}
    
    if (regionToUse == "country") {
      norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), mappedData))
      mappedData <- crop(mappedData, norwayBorderProjected, mask = T)
      lwdF = 0.1
    } else {
      norwayBorderProjected <- terra::project(vect(regionGeometry), crs(mappedData))
      mappedData <- crop(mappedData, norwayBorderProjected, mask = T)
      lwdF = 0
    }

    paletteUsedIndivid <- ifelse(input$statisticIndivid == "probability", "grass", "haxby")
    
    basePlot <- ggplot(norwayBorderProjected) +  
      geom_spatraster(data = mappedData, aes(fill = layer)) +
      geom_sf(fill = NA, lwd = lwdF, colour = "black") +
      coord_sf(datum = st_crs(mappedData)) +
      scale_fill_grass_c(palette = paletteUsedIndivid, na.value = NA, limits = c(0, 100), direction = palDir) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.subtitle = ggtext::element_markdown(),
            plot.title = ggtext::element_textbox_simple(
              size = 16, lineheight = 1, linetype = 1, padding = margin(5, 5, 5, 5), # background fill color
              r = grid::unit(3, "pt"), # radius for rounded corners,
              margin = margin(5,5,15,2)),
      ) +
      labs(
        title = paste0("Artsforekomst<br><span style = 'font-size:10pt;'>Forekomstsannsynlighet for ", 
                       gsub("_", " ", speciesSimpleName)," i ", regionName,".</span>"),
        fill = fillName
      )
    
    return(basePlot)
  }, height = 500)
  
  
  # Text box for hotspots tab
  output$hotspots <- renderPlot({
    
    regionToUse <- ifelse(nchar(input$regionHotspots) == 2, "county", "country")
    regionGeometry <- defineRegion(regionToUse, input$regionHotspots)
    regionName <- names(choiceList)[choiceList == input$regionHotspots]
    taxaToUse <- input$taxaHotspots
    taxaName <- names(taxaChoiceList)[taxaChoiceList == taxaToUse]
    speciesGroupToUse <- input$speciesGroupHotspots
    
    hotspotLayers <- grep(speciesGroupToUse, grep(taxaToUse, names(hotspotListAll), value = TRUE), value = TRUE)
    mappedData <- hotspotListAll[[hotspotLayers]]
    names(mappedData) <- gsub("_", "", gsub(taxaToUse, "", gsub(speciesGroupToUse, "", names(mappedData))))
    hotspotSimple <- ifel(mappedData$HS1 == 1, 1, 
                          ifel(mappedData$HS2 == 1, 2, 
                               ifel(mappedData$HS3 == 1, 3, NA)))
    
    levels(hotspotSimple) <- data.frame(ID=3:1, label=c("HS3", "HS2", "HS1"))
    
    if (regionToUse == "country") {
      norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), hotspotSimple))
      hotspotSimple <- crop(hotspotSimple, norwayBorderProjected, mask = T)
      lwdF = 0.1
    } else {
      norwayBorderProjected <- terra::project(vect(regionGeometry), crs(hotspotSimple))
      hotspotSimple <- crop(hotspotSimple, norwayBorderProjected, mask = T)
      lwdF = 0
    }
    
    basePlot <- ggplot(norwayBorderProjected) +  
      geom_spatraster(data = hotspotSimple,
                      aes(fill = label)) +
      geom_sf(fill = NA, lwd = lwdF, colour = "black")
    
    if (input$overlayHotspots != "none") {
      overlays1 <- vect(readRDS("data/overlays.RDS")$intactAreas)
      overlays <- overlays1[overlays1$arealkm2 > input$overlayLimitHotspots] |> 
        project(hotspotSimple)
      mappedDataMasked <- terra::crop(hotspotSimple, overlays, mask = T)
    } else {mappedDataMasked <- hotspotSimple}
      
    hotspotsMap <- ggplot(norwayBorderProjected) + 
      geom_sf(fill = "beige", lwd = lwdF, colour = "black") + 
      geom_spatraster(data = mappedDataMasked, mapping = aes(fill = label)) +
      scale_fill_manual(values = c("orange", "red", "black"), na.value = NA, na.translate = F) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.subtitle = ggtext::element_markdown(),
            plot.title = ggtext::element_textbox_simple(
              size = 16, lineheight = 1, linetype = 1, padding = margin(5, 5, 5, 5), # background fill color
              r = grid::unit(3, "pt"), # radius for rounded corners,
              margin = margin(5,5,15,2)),
      ) +
      labs(
        title = paste0("Hotspots<br><span style = 'font-size:10pt;'>Hotspots for  ", 
                       taxaName," i ", regionName,".</span>"),
        fill = "Hotspot\nstatus"
      ) +
      coord_sf(datum = st_crs(mappedDataMasked))
    
    return(hotspotsMap)
  }, height = 500)
  
  # Image box for species intensity tab
  output$imageBox1 <- renderImage({
    ImgTxt <- paste0("data/importedData/photos/", input$species,"/speciesImage.jpg")
    if (!file.exists(ImgTxt)) {ImgTxt <- "data/photos/imageNotAvailable.png"}
    list(src = ImgTxt,
         contentType = "image/jpg",
         width = "80%"
    )
  }, deleteFile = FALSE)

  # Text box for species intensity tab
  output$textBox1 <- renderUI({
    simpleName <- namesFromGroups$species[namesFromGroups$fullName == input$species]
    redListStatus <- ansvarsarter$Kategori2021[ansvarsarter$simpleScientificName %in% simpleName]
    redListStatusFull <- ifelse(redListStatus == "EN", "Endangered", ifelse(redListStatus == "VU", "Sarbar", "Kritisk truet"))
    imageUser <- creditList$credit[creditList$species == input$species]
    imageURL <- creditList$url[creditList$species == input$species]
    HTML(paste0("<strong>Vitenskapelig navn:</strong> ", gsub("_"," ",simpleName)  ,
                "<br/><strong>Rødliststatus:</strong> ", redListStatusFull,
                "<br/><strong>Foto:</strong> <a href = ", imageURL, ">", imageUser, "<a/>"))
  })
  
  # Display of metadata report
  getPage<-function() {
    return(includeHTML(paste0("data/speciesMetadata_",input$metaTaxa,".html")))
  }
  output$inc<-renderUI({getPage()})
  
})



