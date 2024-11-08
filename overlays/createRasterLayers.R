### Import libraries
library(terra)
library(dplyr)
library(sf)
focalTaxa <- read.csv("data/focalTaxa.csv")

###------------------------------------###
### Get folders ready to receive data ####
###------------------------------------###

# List of variables which require folders for storage
overlays <- c("carbonStorage", "hovedokosystemer", "inngrepfrieOmrader", "kulturminner", "naturtyper",
              "protectedAreas", "waterRegions")
for (folder in overlays) {
  if (!dir.exists(paste0("overlays/data/", folder))) {dir.create(paste0("overlays/data/", folder))}
}

###---------------------###
### Construct hotspots ####
###---------------------###

filePatterns <- gsub("_.*","",list.files("data/modelOutputs"))
filePatterns <- filePatterns[grep("stats", filePatterns)]
hotspotsCompiled <- lapply(filePatterns, FUN = function(x) {
  allFilesToImport <- list.files("data/modelOutputs")[grep(x, list.files("data/modelOutputs"))]
  compiledRichness <- rast(paste0("data/modelOutputs/",allFilesToImport))
  compiledRichness <- sum(compiledRichness[[names(compiledRichness) %in% "probability"]])
  
  quantiles <- quantile(values(compiledRichness$sum), c(0, 0.9, 0.95, 0.99, 1), na.rm = TRUE)
  categorisedRichness <- classify(compiledRichness$sum, 
                     rcl=quantiles)
  levels(categorisedRichness) <- data.frame(ID=0:3, label=c("1", "2", "3", "4"))
  names(categorisedRichness) <- gsub("stats", "richness", x)
  categorisedRichness
})
hotspots <- rast(hotspotsCompiled)

###-----------------###
### Construct bias ####
###-----------------###

# Now let's get high richness and low bias
biasFiles <- list.files("data/modelOutputs", full.names = TRUE)[grep("bias", list.files("data/modelOutputs"))]
biasCompiled <- sum(rast(biasFiles))
names(biasCompiled) <- "bias"


###------------------###
### Protected areas ####
###------------------###

protected_areas_url <- paste0("https://nedlasting.geonorge.no/geonorge/Natur/ProtectedSites/GML/")
keyDirectory <- "overlays/data/protectedAreas"

# scrape the data URL to obtain a list of zip files
zip_links <- protected_areas_url %>%
  rvest::read_html() %>%
  rvest::html_nodes(xpath = "//a[contains(@href, '.zip')]") %>% # Extract all links that end with '.zip'
  rvest::html_attr("href") %>%  # Get the href attribute of these links
  paste0(url, .) # Concatenate the base URL with each zip link to form the full URL

# zip_links[1] is what we want. For all of Norway at 25833

# Else (quietly) download each zip file, unzip it, and then delete the zip file
invisible(lapply(zip_links[1], function(link) {
  # specify the temporary zip file name
  zip_name <- file.path("overlays", "temp.zip")
  
  # download the zip file to the specified location
  download.file(link, destfile = zip_name, mode = "wb")
  
  # unzip the downloaded file to the target directory
  unzip(zip_name, exdir = keyDirectory)
  
  # delete the zip file
  file.remove(zip_name)
}))

protectedAreasFile <- list.files(keyDirectory, full.names = TRUE)[grep(".gml", list.files(keyDirectory, full.names = TRUE))]
a<- st_read(protectedAreasFile) %>%
  dplyr::select(designation)
protectedVector <- terra::project(vect(a), hotspots)

protectedAreas <- terra::rasterize(protectedVector, hotspots, field = "designation")
names(protectedAreas) <- "protectedAreas"

###--------------------###
### Municipal regions ####
###--------------------###

if (!("csmaps" %in% installed.packages())) {install.packages("csmaps")}
library(csmaps)

municipalities <- csmaps::nor_municip_map_b2024_default_sf
municipalitiesVect <- terra::project(vect(municipalities), hotspots)
municipalityRaster <- terra::rasterize(municipalitiesVect, hotspots, field = "location_code")
names(municipalityRaster) <- "municipalities"
counties <- csmaps::nor_county_map_b2024_default_sf
countiesVect <- terra::project(vect(counties), hotspots)
countyRaster <- terra::rasterize(countiesVect, hotspots, field = "location_code")
names(countyRaster) <- "counties"

###-----------------###
### Cultural areas ####
###-----------------###

# construct the URL containing the list of zip files for the selected format
keyDirectory <- "overlays/data/kulturminner"
url <- "https://nedlasting.geonorge.no/geonorge/Kulturminner/Kulturmiljoer/GML/"

zip_links <- url %>%
  rvest::read_html() %>%
  rvest::html_nodes(xpath = "//a[contains(@href, '.zip')]") %>% # Extract all links that end with '.zip'
  rvest::html_attr("href") %>%  # Get the href attribute of these links
  paste0(url, .) # Concatenate the base URL with each zip link to form the full URL

zip_links[1]

# Else (quietly) download each zip file, unzip it, and then delete the zip file
invisible(lapply(zip_links[1], function(link) {
  # specify the temporary zip file name
  zip_name <- file.path("overlays", "temp.zip")
  
  # download the zip file to the specified location
  download.file(link, destfile = zip_name, mode = "wb")
  
  # unzip the downloaded file to the target directory
  unzip(zip_name, exdir = keyDirectory)
  
  # delete the zip file
  file.remove(zip_name)
}))
culturalAreasFile <- list.files(keyDirectory, full.names = TRUE)[grep( ".gml",list.files(keyDirectory, full.names = TRUE))]
a<- st_read(culturalAreasFile) %>%
  dplyr::select(navn)
culturalVector <- terra::project(vect(a), hotspots)

culturalAreas <- terra::rasterize(culturalVector, hotspots, field = "navn")
names(culturalAreas) <- "culturalAreas"

###---------------###
### Nature types ####
###---------------###

# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
a<- st_read("overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml") %>%
  dplyr::select(prosjektområdenavn, dekningskartverdi)
natureTypesVector <- terra::project(vect(a), hotspots)

natureTypes <- terra::rasterize(natureTypesVector, hotspots, field = "dekningskartverdi")
names(natureTypes) <- "natureTypes"

###----------------------###
### Water regions/areas ####
###----------------------###

a<- st_read("overlays/data/waterRegions/Vannomrader_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)
waterAreaVector <- terra::project(vect(a), hotspots)

waterAreas <- terra::rasterize(waterAreaVector, hotspots, field = "navn")
names(waterAreas) <- "waterAreas"

a<- st_read("overlays/data/waterRegions/Vannregioner_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)
castedA <- st_cast(a, "MULTIPOLYGON") %>% st_collection_extract("POLYGON")
waterRegionVector <- terra::project(vect(castedA), hotspots)

waterRegions <- terra::rasterize(waterRegionVector, hotspots, field = "navn")
names(waterRegions) <- "waterRegions"


###-----------------###
### Carbon storage ####
###-----------------###

belowGroundCarbon <- rast("overlays/data/carbonStorage/belowground_biomass_carbon_2010.tif")
belowGroundCarbonNorway <- terra::project(belowGroundCarbon, hotspots, method = "average")
belowGroundCarbonNorway <- crop(belowGroundCarbonNorway, hotspots[[1]], mask = T)
names(belowGroundCarbonNorway) <- "belowgroundcarbon"
aboveGroundCarbon <- rast("overlays/data/carbonStorage/aboveground_biomass_carbon_2010.tif")
aboveGroundCarbonNorway <- terra::project(aboveGroundCarbon, hotspots, method = "average")
aboveGroundCarbonNorway <- crop(aboveGroundCarbonNorway, hotspots[[1]], mask = T)
names(aboveGroundCarbonNorway) <- "abovegroundcarbon"

###-------------------###
### Hovedokosystemer ####
###-------------------###

# This is a massive file so once you've downloaded and rasterised it, for god's sake save it

if (length(list.files("overlays/data/hovedokosystemer")) == 0) {
  a <- st_read("overlays/data/hovedokosystemer/Hovedokosystem_nedlasting/Hovedokosystem.gdb") %>%
    dplyr::select("Hovedøkosystem")
  
  
  hovedokosystemVector <- terra::project(vect(a), hotspots)
  hovedokosystems <- terra::rasterize(hovedokosystemVector, hotspots, field = "ecotype")
  translationTable <- data.frame(ints = 1:12, categories = c("Bebyggelse-samferdsel", "Dyrket mark", "Grasmark", "Skog",
                                                             "Hei og apen vegetasjon", "Lite vegetert mark", "Vatmark",  "Elver-bekker",
                                                             "","Innsjoer-tern", "Svaberg kyststrender og dyner", "Apent hav"))
  levels(hovedokosystems) <- translationTable
  names(hovedokosystems) <- "Hovedokosystems"
  writeRaster(hovedokosystems, "overlays/data/hovedokosystemer/hovedokosystemer.tiff", overwrite = TRUE)
} else {
  hovedokosystems <- rast("overlays/data/hovedokosystemer/hovedokosystemer.tiff")
}

###-----------------------###
### Inngrepsfrie omrader ####
###-----------------------###

inngrepsfrieOmrader <- read_sf("overlays/data/inngrepsfrieOmrader/statusPolygon.shp")
ifVect <- terra::project(vect(inngrepsfrieOmrader), hotspots)
ifRaster <- terra::rasterize(ifVect, hotspots, field = "vsone")
translationTable <- data.frame(ints = 0:3, categories = c("Sone1", "Sone2", NA, "Vill"))
levels(ifRaster) <- translationTable
names(ifRaster) <- "intactAreas"

###------------------------###
### Import beta diversity ####
###------------------------###

betaDiversity <- rast("betaDiversity/data/betaDiversity.tiff")
betaDiversity <- crop(betaDiversity, hotspots)


###----------------------------------###
### Collate what we have up to here ####
###----------------------------------###

rasterForDelivery <- c(hotspots, biasCompiled, protectedAreas, culturalAreas, natureTypes, municipalityRaster, 
                       countyRaster, waterRegions, waterAreas, belowGroundCarbonNorway, aboveGroundCarbonNorway,
                       hovedokosystems, ifRaster, betaDiversity)
writeRaster(rasterForDelivery, file = "overlays/data/rasterPackage3.tiff", overwrite = TRUE)
