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
ddList <- readRDS("data/ddList.RDS")
namesFromGroups <- readRDS("data/namesFromGroups.RDS")

creditList <- readRDS("data/imageCredit.RDS")
ansvarsarter <- readRDS("data/ansvarsArterList.RDS")

# projCRS <- sf::st_crs(25833)$proj4string
# 
choiceList <- c("Norway", "50", "54") |> 
  setNames(c("Norway", "Trondelag", "Finnmark"))
taxaChoiceList <- c("hymenopterans", "birds", "vascularPlants") |> 
  setNames(c("Hymenopterans", "Birds", "Vascular Plants"))

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
    dataToUse <- input$speciesGroup
    taxaToUse <- input$taxa
    metricToUse <- input$statistic
    
    taxaName <- names(taxaChoiceList)[taxaChoiceList == taxaToUse]
    
    focalData <- rast(paste0("data/", dataToUse, "stats_", taxaToUse, "final.tiff"))
    focalData$layer <- focalData[[metricToUse]]
    
    norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), focalData))
    
    paletteUsed <- ifelse(input$statistic == "skalertRikhet", "forest_cover", "haxby")
    
    returnPlot <- ggplot(norwayBorderProjected) + 
      geom_spatraster(data = focalData, mapping = aes(fill = layer)) +
      geom_sf(fill = NA, lwd = 0.25, colour = "black") + 
      coord_sf(datum = st_crs(focalData)) +
      scale_fill_grass_c(palette = paletteUsed, na.value = NA, limits = c(0, 1), direction = 1) +
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
        title = paste0("Species richness<br><span style = 'font-size:10pt;'>Species richness per 10 x 10 km pixel for ", 
                       taxaName," across Norway.</span>"),
        fill = ifelse(metricToUse == "skalertRikhet", "Scaled\nspecies\nrichness", "Scaled\nuncertainty")
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
    
    mappedData <- rast(paste0("data/",taxaToUse, "10betaDiversity.tiff"))
    
    regionGeometryVec <- terra::project(vect(regionGeometry), crs(mappedData))
    mappedDataCropped <- mappedData |> 
      crop(regionGeometryVec, mask = T)
    
    lineW <- ifelse(regionToUse == "country", 0.25,0)
    
    plotToReturn <- ggplot(regionGeometryVec) +  
      geom_spatraster(data = mappedDataCropped,
                      aes(fill = sum)) +
      scale_fill_grass_c(palette = "forest_cover", na.value = NA, limits = c(0, 1), direction = 1) +
      geom_sf(fill = NA, lwd = lineW, colour = "black") +
      theme_bw() +
      coord_sf(datum = st_crs(mappedDataCropped)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.subtitle = ggtext::element_markdown(),
            plot.title = ggtext::element_textbox_simple(
              size = 16, lineheight = 1, linetype = 1, padding = margin(5, 5, 5, 5), # background fill color
              r = grid::unit(3, "pt"), # radius for rounded corners,
              margin = margin(5,5,15,2)),
      ) +
      labs(
        title = paste0("Beta diversity<br><span style = 'font-size:10pt;'>Percentage of total species found in larger area that are present in given 10 x 10 km pixel for  ", 
                       taxaName," across ", regionName,".</span>"),
        fill = "Species beta\ndiversity"
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
    
    mappedData <- rast(paste0("data/",input$taxaIndivid, "_", input$statisticIndivid,".tiff"))[[speciesSimpleName]]
    mappedData$layer <- mappedData[[1]]
    
    regionGeometryVec <- terra::project(vect(regionGeometry), crs(mappedData))
    mappedData <- mappedData |> 
      crop(regionGeometryVec, mask = TRUE)
    
    
    basePlot <- ggplot(regionGeometryVec) +  
      geom_spatraster(data = mappedData, aes(fill = layer)) +
      geom_sf(fill = NA, lwd = 0.25, colour = "black") +
      coord_sf(datum = st_crs(mappedData)) +
      scale_fill_grass_c(palette = "grass", na.value = NA, limits = c(0, 100), direction = -1) +
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
        title = paste0("Species occurrence<br><span style = 'font-size:10pt;'>Probability of occurrence for ", 
                       gsub("_", " ", speciesSimpleName)," across ", regionName,".</span>"),
        fill = "Occurrence\nprobability"
      )
    
    return(basePlot)
  }, height = 500)
  
  
  # Text box for species intensity tab
  output$hotspots <- renderPlot({
    
    regionToUse <- ifelse(nchar(input$regionHotspots) == 2, "county", "country")
    regionGeometry <- defineRegion(regionToUse, input$regionHotspots)
    regionName <- names(choiceList)[choiceList == input$regionBeta]
    taxaToUse <- input$taxaHotspots
    taxaName <- names(taxaChoiceList)[taxaChoiceList == taxaToUse]
    
    taxaGrouping <- ifelse(input$taxaHotspots == "birds", "fugler", ifelse(input$taxaHotspots == "vascularPlants", "karplanter", "insekter"))
    mappedData <- rast(paste0("data/", taxaGrouping, "_", input$speciesGroupHotspots, "_rikhet.tiff"))
      
    regionGeometryVec <- terra::project(vect(regionGeometry), crs(mappedData))
    mappedData <- mappedData |> 
      crop(regionGeometryVec, mask = TRUE)
    mappedData$layer <- mappedData[[1]]
    
    basePlot <- ggplot(regionGeometryVec) +  
      geom_spatraster(data = mappedData,
                      aes(fill = layer)) +
      geom_sf(fill = NA, lwd = 0.5, colour = "black")
    
    if (input$overlayHotspots != "none") {
      overlays <- vect(readRDS("data/overlays.RDS")$intactAreas)
      overlays <- overlays[overlays$arealkm2 > input$overlayLimitHotspots] |> 
        mask(regionGeometryVec)
      mappedDataMasked <- terra::crop(mappedData, overlays, mask = T)
    } else {mappedDataMasked <- mappedData}
      
    hotspotsMap <- ggplot(regionGeometryVec) + 
      geom_sf(fill = "white", lwd = 0.1, colour = "black") + 
      geom_spatraster(data = mappedDataMasked, mapping = aes(fill = layer)) +
      scale_fill_manual(values = c("beige", "orange", "red", "black"), na.value = NA, na.translate = F) +
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
        title = paste0("Hotspots<br><span style = 'font-size:10pt;'>Locations of hotspots for  ", 
                       taxaName," across ", regionName,".</span>"),
        fill = "Hotspot\nstatus"
      ) +
      coord_sf(datum = st_crs(mappedDataMasked))
    
    return(hotspotsMap)
  }, height = 500)
  
  # Image box for species intensity tab
  output$imageBox1 <- renderImage({
    ImgTxt <- paste0("data/photos/", input$species,"/speciesImage.jpg")
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
    redListStatusFull <- ifelse(redListStatus == "EN", "Endangered", ifelse(redListStatus == "VU", "Vulnerable", "Critical"))
    imageUser <- creditList$credit[creditList$species == input$species]
    imageURL <- creditList$url[creditList$species == input$species]
    HTML(paste0("<strong>Scientific name:</strong> ", gsub("_"," ",simpleName)  ,
                "<br/><strong>Red list status:</strong> ", redListStatusFull,
                "<br/><strong>Image Credit:</strong> <a href = ", imageURL, ">", imageUser, "<a/>"))
  })
  
  # Display of metadata report
  getPage<-function() {
    return(includeHTML(paste0("data/speciesMetadata_",input$metaTaxa,".html")))
  }
  output$inc<-renderUI({getPage()})
  
})



