
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###

library(terra)
library(dplyr)
library(sf)
library(exactextractr)

hotspotFiles <- list.files("data/modelOutputs/run_2025-01-06", pattern = "final", full.names = TRUE)

translationTable <- read.csv("data/taxaTranslations.csv")
#translationTable$nynorsk[25] <- "Insekt"

# produce 10% hotspots
hotspots <- lapply(hotspotFiles, FUN  = function(x) {
  rasterToProcess <- rast(x)
  quantiles <- quantile(values(rasterToProcess[["skalertRikhet"]]), c(0, 0.9, 1), na.rm = TRUE)
  categorisedRichness <- classify(rasterToProcess[["skalertRikhet"]], 
                                  rcl=quantiles)
  levels(categorisedRichness) <- data.frame(ID=0:1, label=c("Not hotspot", "Hotspot"))
  return(categorisedRichness)
}) |> setNames(gsub("data/modelOutputs/run_2025-01-06", "",hotspotFiles))

### Summarizing function
sum_cover <- function(x){
  list(x %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)) %>%
         mutate(proportion = total_area/sum(total_area)))
  
}

###---------------###
#### WATER AREAS ####
###---------------###

waterAreas <- st_read("processes/overlays/data/waterRegions/Vannomrader_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)

#extract the area of each raster cell covered by the plot and summarize
waterAreaProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], waterAreas, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- waterAreas$navn
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Vannområde") %>%
    filter(!is.na(value)) %>%
    group_by(Vannområde) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea) %>%
    filter(value == 1) %>%
    dplyr::select(Vannområde, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  if (nrow(test) > 0) {test$group <- paste0(grouping, speciesName)}
  
  test
  
}
)

# Collate data frame for all species groups
asDF <- do.call(bind_rows, waterAreaProportions)
asDF$proportion <- round(asDF$proportion,4)
widVersionWaterAreas <- reshape(asDF, idvar = "Vannområde", timevar = "group", direction = "wide")
colnames(widVersionWaterAreas) <- gsub("proportion.", "", colnames(widVersionWaterAreas))
widVersionWaterAreas[is.na(widVersionWaterAreas)] <- 0

write.csv(widVersionWaterAreas, "processes/overlays/data/waterAreaProportions.csv")

###-----------------###
#### WATER REGIONS ####
###-----------------###

waterRegions <- st_read("processes/overlays/data/waterRegions/Vannregioner_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)
waterRegions <- waterRegions[waterRegions$navn != "Jan Mayen",]

#extract the area of each raster cell covered by the plot and summarize
waterRegionsProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], waterRegions, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- waterRegions$navn
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Vannregion") %>%
    filter(!is.na(value)) %>%
    group_by(Vannregion) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea) %>%
    filter(value == 1) %>%
    dplyr::select(Vannregion, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  if (nrow(test) > 0) {test$group <- paste0(grouping, speciesName)}
  
  test
}
)

# Collate data frame for all species groups
asDF <- do.call(bind_rows, waterRegionsProportions)
asDF$proportion <- round(asDF$proportion,4)
widVersionWaterRegions <- reshape(asDF, idvar = "Vannregion", timevar = "group", direction = "wide")
colnames(widVersionWaterRegions) <- gsub("proportion.", "", colnames(widVersionWaterRegions))
widVersionWaterRegions[is.na(widVersionWaterRegions)] <- 0


write.csv(widVersionWaterRegions, "processes/overlays/data/waterRegionProportions.csv")

###------------------###
#### MUNICIPALITIES ####
###------------------###

municipalities <- csmaps::nor_municip_map_b2024_default_sf
# municipalityRaster <- terra::rasterize(municipalitiesVect, hotspotLayers, field = "location_code")
# names(municipalityRaster) <- "municipalities"
# counties <- csmaps::nor_county_map_b2024_default_sf
# countiesVect <- terra::project(vect(counties), hotspotLayers)

#extract the area of each raster cell covered by the plot and summarize
municipalityProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], municipalities, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- gsub("municip_nor", "", municipalities$location_code)
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Kommune") %>%
    filter(!is.na(value)) %>%
    group_by(Kommune) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea) %>%
    filter(value == 1) %>%
    dplyr::select(Kommune, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  if (nrow(test) > 0) {test$group <- paste0(grouping, speciesName)}
  
  test
}
)

# Collate data frame for all species groups
asDF <- do.call(bind_rows, municipalityProportions)
asDF$proportion <- round(asDF$proportion,4)
wideMunicipalities <- reshape(asDF, idvar = "Kommune", timevar = "group", direction = "wide")
colnames(wideMunicipalities) <- gsub("proportion.", "", colnames(wideMunicipalities))
wideMunicipalities[is.na(wideMunicipalities)] <- 0


write.csv(wideMunicipalities, "processes/overlays/data/municipalityProportions.csv")

###------------###
#### COUNTIES ####
###------------###

counties <- csmaps::nor_county_map_b2024_default_sf

#extract the area of each raster cell covered by the plot and summarize
countyProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], counties, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- gsub("municip_nor", "", counties$location_code)
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Fylke") %>%
    filter(!is.na(value)) %>%
    group_by(Fylke) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea) %>%
    filter(value == 1) %>%
    dplyr::select(Fylke, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  if (nrow(test) > 0) {test$group <- paste0(grouping, speciesName)}
  
  test
}
)

# Collate data frame for all species groups
asDF <- do.call(bind_rows, countyProportions)
asDF$proportion <- round(asDF$proportion,4)
wideCounties <- reshape(asDF, idvar = "Fylke", timevar = "group", direction = "wide")
colnames(wideCounties) <- gsub("proportion.", "", colnames(wideCounties))
wideCounties[is.na(wideCounties)] <- 0

write.csv(wideCounties, "processes/overlays/data/countyProportions.csv")

###-----------------------------------------###
#### MAPPED NATURE TYPES - MAPPED/UNMAPPED ####
###-----------------------------------------###

# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
natureTypes<- st_read("processes/overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml") %>%
  dplyr::select(prosjektområdenavn, dekningskartverdi)
talliedGroups <- tally(group_by(natureTypes, dekningskartverdi)) %>%
  filter(dekningskartverdi != "ikkeSystematiskKartlagt")

#extract the area of each raster cell covered by the plot and summarize
natureTypeProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], talliedGroups, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- talliedGroups$dekningskartverdi
  
  # Make table to fill
  gridToFill <- data.frame(names = names(table1))
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Plot_buffer") %>%
    filter(!is.na(value)) %>%
    group_by(Plot_buffer) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea)%>%
    filter(value == 1) %>%
    dplyr::select(Plot_buffer, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  
  # Here we have a slightly different end product, since we just want to compare the two types
  gridToFill$proportion <- test$proportion[match(gridToFill$names, test$Plot_buffer)]
  gridToFill$proportion[is.na(gridToFill$proportion)] <- 0
  
  finalDF <- data.frame(grundigKartlagtMedFunn = gridToFill$proportion[gridToFill$names == "grundigKartlagtMedFunn"], 
                        kartlagtUtenFunn = gridToFill$proportion[gridToFill$names == "kartlagtUtenFunn"], 
                        group = paste0(grouping, speciesName))
  finalDF
}
)

# Collate data frame for all species groups
asDF <- do.call(bind_rows, natureTypeProportions)
asDF <- asDF %>% 
  mutate(across(1:2, function(x) {round(x,4)})) %>%
  select(3,1,2)
colnames(asDF)[1] <- "Artsgrupper"

write.csv(asDF, "processes/overlays/data/mappedNatureTypesProportions.csv")

# We also want to produce a plot for this one statistic
# Need to prepare the dataset for plotting
x2 <- reshape2::melt(asDF, id = "Artsgrupper", variable_name = "type")
x2$proportion <- x2$value * 100
x2$grouping <- trimws(x2$Artsgrupper , whitespace = "\\ .*")

x2$variableSimple <- plyr::revalue(x2$variable, c("grundigKartlagtMedFunn" = "Kartlagt med funn", 
                                                            "kartlagtUtenFunn" = "Kartlagt uten funn"))

for (t in unique(x2$grouping)) {
  allSpecies <- x2[x2$grouping == t,]
  allSpecies$groupShort <- gsub(t, "", allSpecies$Artsgrupper)
  allSpecies$groupShort <- ifelse(grepl("Bakkehekkande", allSpecies$groupShort), "Bakkehekkande\nfuglar", allSpecies$groupShort)
  
  natureTypePlot <- ggplot(allSpecies, aes(x=groupShort, y=proportion, fill=variableSimple)) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")+
    scale_fill_brewer(palette="Paired")+
    theme_minimal() +
    labs(fill = "", y = "Andel av arealet som er hotspot", x = "Artsgruppe") +
    ylim(c(0,50))
  
  filenameToUse <- paste0("visualisations/figures/overlayAnalysis/kartlagteNatur_", t, "Arter.jpg")
  ggsave(filename =  filenameToUse,
         plot = natureTypePlot, units = "px",
         width = 2000, height = 1400
  )
  
}


###------------------------------------###
#### MAPPED NATURE TYPES - CATEGORIES ####
###------------------------------------###

# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
natureTypesCats<- st_read("processes/overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml", layer = "Naturtype_nin") %>%
  dplyr::select(tilstand)
talliedGroups <- tally(group_by(natureTypesCats, tilstand)) %>%
  filter(!is.na(tilstand))

#extract the area of each raster cell covered by the plot and summarize
natureTypeProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], talliedGroups, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- st_drop_geometry(talliedGroups)$tilstand
  
  # Make table to fill
  gridToFill <- data.frame(names = names(table1))
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Plot_buffer") %>%
    filter(!is.na(value)) %>%
    group_by(Plot_buffer) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea)%>%
    filter(value == 1) %>%
    dplyr::select(Plot_buffer, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  if (nrow(test) > 0) {test$group <- paste0(grouping, speciesName)}
  
  test
})

asDF <- do.call(bind_rows, natureTypeProportions)
asDF$proportion <- round(asDF$proportion,4)
wideNatureCategories <- reshape(asDF, idvar = "Plot_buffer", timevar = "group", direction = "wide")
colnames(wideNatureCategories) <- gsub("proportion.", "", colnames(wideNatureCategories))
wideNatureCategories[is.na(wideNatureCategories)] <- 0
colnames(wideNatureCategories)[1] <- "Økologisk tilstand"

write.csv(wideNatureCategories, "processes/overlays/data/ninNatureStandardsProportions.csv")


###--------------------------###
#### SELECTED CULTURE TYPES ####
###--------------------------###

# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
selectedCultureTypes<- st_read("processes/overlays/data/utvalgteKulturtyper/KulturlandskapUtvalgte_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)

#extract the area of each raster cell covered by the plot and summarize
selectedCultureTypesProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  # Define species group to work with
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], selectedCultureTypes, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- selectedCultureTypes$navn
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Plot_buffer") %>%
    filter(!is.na(value)) %>%
    group_by(Plot_buffer) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea) %>%
    filter(value == 1) %>%
    dplyr::select(Plot_buffer, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  if (nrow(test) > 0) {test$group <- paste0(grouping, speciesName)}
  
  test
}
)

asDF <- do.call(bind_rows, selectedCultureTypesProportions)
asDF$proportion <- round(asDF$proportion,4)
wideSlectedCultureTypes <- reshape(asDF, idvar = "Plot_buffer", timevar = "group", direction = "wide")
colnames(wideSlectedCultureTypes) <- gsub("proportion.", "", colnames(wideSlectedCultureTypes))
wideSlectedCultureTypes[is.na(wideSlectedCultureTypes)] <- 0
colnames(wideSlectedCultureTypes)[1] <- "Kulturlandskap"


write.csv(wideSlectedCultureTypes, "processes/overlays/data/selectedCultureLandscapesProportions.csv")



###-------------------------###
#### SELECTED NATURE TYPES ####
###-------------------------###

# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
selectedNatureTypes<- st_read("processes/overlays/data/utvalgteNaturtyper/NaturtyperUtvalgte_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(områdenavn, utvalgtNaturtypeTekst)
selectedNatureTypes <- selectedNatureTypes %>%
  filter(utvalgtNaturtypeTekst != "Hule eiker") %>%
  group_by(utvalgtNaturtypeTekst) %>%
  tally()

#extract the area of each raster cell covered by the plot and summarize
selectedNatureTypesProportions <- lapply(seq_along(hotspots), FUN = function(x) {
  
  speciesName <- translationTable$nynorsk[translationTable$engelsk == gsub("final.tiff","",strsplit(names(hotspots[x]), "_")[[1]][2])]
  grouping <- strsplit(names(hotspots[x]), "_")[[1]][1]
  grouping <- ifelse(grepl("allspecies", grouping), "Alle ", ifelse(grepl("ansvar", grouping), "Ansvar ", "Trua "))
  
  # Get area of each area which is inside hotspots
  table1 <- exact_extract(hotspots[[x]], selectedNatureTypes, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
  names(table1) <- selectedNatureTypes$utvalgtNaturtypeTekst
  
  idTable <- data.frame(names = selectedNatureTypes$utvalgtNaturtypeTekst) 
  
  #merge the list elements into a df
  test <- bind_rows(table1, .id = "Naturtype") %>%
    filter(!is.na(value)) %>%
    group_by(Naturtype) %>%
    mutate(totalArea = sum(total_area)) %>%
    mutate(proportion2 = total_area/totalArea) %>%
    filter(value == 1) %>%
    dplyr::select(Naturtype, proportion2) %>%
    as.data.frame()
  colnames(test)[2] <- "proportion"
  
  idTable$proportion <- test$proportion[match(idTable$names, test$Naturtype)]
  idTable$proportion[is.na(idTable$proportion)] <- 0
  idTable$group <- paste0(grouping, speciesName)


  idTable
}
)

asDF <- do.call(bind_rows, selectedNatureTypesProportions)
asDF$proportion <- round(asDF$proportion,4)
wideSlectedNatureTypes <- reshape(asDF, idvar = "names", timevar = "group", direction = "wide")
colnames(wideSlectedNatureTypes) <- gsub("proportion.", "", colnames(wideSlectedNatureTypes))
wideSlectedNatureTypes[is.na(wideSlectedNatureTypes)] <- 0
colnames(wideSlectedNatureTypes)[1] <- "Naturtyp"

write.csv(wideSlectedNatureTypes, "processes/overlays/data/selectedNatureTypesProportions.csv")


