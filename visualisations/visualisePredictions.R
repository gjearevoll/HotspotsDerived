

boxedLocations <-list(c(182000, 6852000, 190000, 6869000) |>
                        setNames(c("west", "south", "east", "north")),
                      c(218000, 6585000, 245000, 6600000) |>
                        setNames(c("west", "south", "east", "north")),
                      c(-60000, 6640000, -40000, 6670000) |>
                        setNames(c("west", "south", "east", "north")),
                      c(233000, 7070000, 260000, 7100000) |>
                        setNames(c("west", "south", "east", "north")),
                      c(334000, 6951000, 345000, 6990000) |>
                        setNames(c("west", "south", "east", "north")),
                      c(260000, 7024000, 295000, 7040000) |>
                        setNames(c("west", "south", "east", "north")))

# Import field data
fieldWorkResults <- readRDS("processes/localPredictions/data/fieldWorkResults.RDS")

processedSpecies <- gsub(".tiff", "", list.files("processes/localPredictions/data/predictions/resolution500"))
species <- intersect(processedSpecies, names(fieldWorkResults))

# Import relevant species
resolutions <- c("250", "500", "1000")

datasets <- lapply(resolutions, FUN = function(x) {
  
  fileList <- list.files(paste0("processes/localPredictions/data/predictions/resolution", x), full.names = TRUE)
  speciesNames <- gsub(".tiff", "" ,list.files(paste0("processes/localPredictions/data/predictions/resolution", x)))
  
  dataList <- lapply(fileList, rast) |>
    setNames(speciesNames)
  validSpecies <- speciesNames[speciesNames %in% species]
  predictions <- rast(lapply(dataList, `[`, "mean"))[[validSpecies]]
  return(predictions[['Acer_platanoides']])
  # summedPreds <- sum(predictions)
  # summedPreds$scaledRichness <- (summedPreds$sum-minmax(summedPreds$sum)[1])/
  #   (minmax(summedPreds$sum)[2]-minmax(summedPreds$sum)[1])
  # summedPreds
}) |> setNames(resolutions)


# Let's map out our predictions and our datasets first, including a map of the coastline
norwayBorder <- st_union(csmaps::nor_municip_map_b2024_default_sf)
norwayBorderProjected <- fillHoles(terra::project(vect(norwayBorder), datasets[[1]]))

# First show whole area
ggplot()  +
  geom_spatvector(data = norwayBorderProjected) + 
  geom_spatraster(data = datasets$'250', aes(fill = Acer_platanoides)) +
  scale_fill_grass_c(palette = "grass", na.value = NA, limits = c(0, 1), direction = -1) +
  coord_sf(datum = crs(datasets$`1000`),
           xlim = c(-100000, 400000), 
           ylim = c(6400000, 7200000)) +
  labs(fill = "Occurrence\nprobability"
  )


# Level 1
croppingArea <- c(182000, 6852000, 190000, 6869000)
croppingArea2 <- c(184000, 6856000, 188000, 6865000)
croppingArea3 <- c(185000, 6856500, 187500, 6861500)
resolution <- resolutions[3]
croppedData <- datasets[[resolution]] |>
  crop(ext(c(croppingArea[1], croppingArea[3], croppingArea[2], croppingArea[4])))

ve <- vect(ext(c(croppingArea2[1], croppingArea2[3], croppingArea2[2], croppingArea2[4])), crs=crs(datasets[[resolution]]))
ve2 <- vect(ext(c(croppingArea3[1], croppingArea3[3], croppingArea3[2], croppingArea3[4])), crs=crs(datasets[[resolution]]))

# Create map
ggplot()  +
  geom_spatvector(data = norwayBorderProjected) + 
  geom_spatraster(data = croppedData, aes(fill = Acer_platanoides)) +
  scale_fill_grass_c(palette = "grass", na.value = NA, limits = c(0, 1), direction = -1) +
  geom_spatvector(data = ve, colour = "red", fill = NA, lwd = 1) + 
  geom_sf(data = fieldWorkResults["Acer_platanoides"]) +
  coord_sf(datum = crs(datasets$`1000`),
           xlim = c(croppingArea[1], croppingArea[3]), 
           ylim = c(croppingArea[2], croppingArea[4]))+
  labs(fill = "Skalert\nrikhet"
  )

# Level 2
resolution <- resolutions[2]
croppedData <- datasets[[resolution]] |>
  crop(ext(c(croppingArea2[1], croppingArea2[3], croppingArea2[2], croppingArea2[4])))

# Create map
ggplot()  +
  geom_spatvector(data = norwayBorderProjected) + 
  geom_spatraster(data = croppedData, aes(fill = Acer_platanoides)) +
  scale_fill_grass_c(palette = "grass", na.value = NA, limits = c(0, 1), direction = -1) +
  geom_spatvector(data = ve2, colour = "red", fill = NA, lwd = 1) + 
  geom_sf(data = fieldWorkResults["Acer_platanoides"]) +
  coord_sf(datum = crs(datasets$`1000`),
           xlim = c(croppingArea2[1], croppingArea2[3]), 
           ylim = c(croppingArea2[2], croppingArea2[4]))+
  labs(fill = "Occurrence\nprobability"
  )




# Level 3
resolution <- resolutions[1]
croppedData <- datasets[[resolution]] |>
  crop(ext(c(croppingArea3[1], croppingArea3[3], croppingArea3[2], croppingArea3[4])))

# Create map
ggplot()  +
  geom_spatvector(data = norwayBorderProjected) + 
  geom_spatraster(data = croppedData, aes(fill = Acer_platanoides)) +
  scale_fill_grass_c(palette = "grass", na.value = NA, limits = c(0, 1), direction = -1) +
  geom_sf(data = fieldWorkResults["Acer_platanoides"]) +
  coord_sf(datum = crs(datasets$`1000`),
           xlim = c(croppingArea3[1], croppingArea3[3]), 
           ylim = c(croppingArea3[2], croppingArea3[4]))+
  labs(fill = "Occurrence\nprobability"
  )
