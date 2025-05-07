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
  if (!dir.exists(paste0("processes/overlays/data/", folder))) {dir.create(paste0("processes/overlays/data/", folder))}
}

translations <- read.csv("data/taxaTranslations.csv")

###---------------------###
### Download hotspots ####
###---------------------###

hotspotLayers <- rast("data/allHotspots.tiff")


###-----------------###
### Construct bias ####
###-----------------###

bias <- rast("data/allBiases.tiff")

###------------------###
### Protected areas ####
###------------------###

# Download verneomrader from geonorge
protected_areas_url <- paste0("https://nedlasting.geonorge.no/geonorge/Natur/ProtectedSites/GML/")
keyDirectory <- "processes/overlays/data/protectedAreas"

# scrape the data URL to obtain a list of zip files
zip_links <- protected_areas_url %>%
  rvest::read_html() %>%
  rvest::html_nodes(xpath = "//a[contains(@href, '.zip')]") %>% # Extract all links that end with '.zip'
  rvest::html_attr("href") %>%  # Get the href attribute of these links
  paste0(protected_areas_url, .) # Concatenate the base URL with each zip link to form the full URL

# zip_links[1] is what we want. For all of Norway at 25833

# Else (quietly) download each zip file, unzip it, and then delete the zip file
invisible(lapply(zip_links[1], function(link) {
  # specify the temporary zip file name
  zip_name <- file.path("processes/overlays", "temp.zip")
  
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
protectedVector <- terra::project(vect(a), hotspotLayers)

protectedAreas <- terra::rasterize(protectedVector, hotspotLayers, field = "designation")
names(protectedAreas) <- "protectedAreas"

###--------------------###
### Municipal regions ####
###--------------------###

# These regions are available from the csmaps package
if (!("csmaps" %in% installed.packages())) {install.packages("csmaps")}
library(csmaps)

municipalities <- csmaps::nor_municip_map_b2024_default_sf
municipalitiesVect <- terra::project(vect(municipalities), hotspotLayers)
municipalityRaster <- terra::rasterize(municipalitiesVect, hotspotLayers, field = "location_code")
names(municipalityRaster) <- "municipalities"
counties <- csmaps::nor_county_map_b2024_default_sf
countiesVect <- terra::project(vect(counties), hotspotLayers)
countyRaster <- terra::rasterize(countiesVect, hotspotLayers, field = "location_code")
names(countyRaster) <- "counties"

###-----------------###
### Cultural areas ####
###-----------------###

# construct the URL containing the list of zip files for the selected format
keyDirectory <- "processes/overlays/data/kulturminner"
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
  zip_name <- file.path("processes/overlays", "temp.zip")
  
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
culturalVector <- terra::project(vect(a), hotspotLayers)

culturalAreas <- terra::rasterize(culturalVector, hotspotLayers, field = "navn")
names(culturalAreas) <- "culturalAreas"

###---------------###
### Nature types ####
###---------------###

# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
a<- st_read("processes/overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml") %>%
  dplyr::select(prosjektområdenavn, dekningskartverdi)
natureTypesVector <- terra::project(vect(a), hotspotLayers)

natureTypes <- terra::rasterize(natureTypesVector, hotspotLayers, field = "dekningskartverdi")
names(natureTypes) <- "natureTypes"

###----------------------###
### Water regions/areas ####
###----------------------###

a<- st_read("processes/overlays/data/waterRegions/Vannomrader_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)
waterAreaVector <- terra::project(vect(a), hotspotLayers)

waterAreas <- terra::rasterize(waterAreaVector, hotspotLayers, field = "navn")
names(waterAreas) <- "waterAreas"

a<- st_read("processes/overlays/data/waterRegions/Vannregioner_0000_norge_25833_FILEGDB.gdb") %>%
  dplyr::select(navn)
castedA <- st_cast(a, "MULTIPOLYGON") %>% st_collection_extract("POLYGON")
waterRegionVector <- terra::project(vect(castedA), hotspotLayers)

waterRegions <- terra::rasterize(waterRegionVector, hotspotLayers, field = "navn")
names(waterRegions) <- "waterRegions"


###-----------------###
### Carbon storage ####
###-----------------###

# These are large files and need to be downloaded locally.
# Here is the link to the data: https://cmr.earthdata.nasa.gov/search/concepts/C2764708636-ORNL_CLOUD.html  
belowGroundCarbon <- rast("processes/overlays/data/carbonStorage/belowground_biomass_carbon_2010.tif")
belowGroundCarbonNorway <- terra::project(belowGroundCarbon, hotspotLayers, method = "average")
belowGroundCarbonNorway <- crop(belowGroundCarbonNorway, hotspotLayers[[1]], mask = T)
names(belowGroundCarbonNorway) <- "belowgroundcarbon"

aboveGroundCarbon <- rast("processes/overlays/data/carbonStorage/aboveground_biomass_carbon_2010.tif")
aboveGroundCarbonNorway <- terra::project(aboveGroundCarbon, hotspotLayers, method = "average")
aboveGroundCarbonNorway <- crop(aboveGroundCarbonNorway, hotspotLayers[[1]], mask = T)
names(aboveGroundCarbonNorway) <- "abovegroundcarbon"

###-------------------###
### Hovedokosystemer ####
###-------------------###

# This is a massive file so once you've downloaded and rasterised it, for god's sake save it.
# The file is available through Miljodirektoratet at
# https://kartkatalog.miljodirektoratet.no/Dataset/Details/3069

if (length(list.files("processes/overlays/data/hovedokosystemer")) == 0) {
  a <- st_read("processes/overlays/data/hovedokosystemer/Hovedokosystem_nedlasting/Hovedokosystem.gdb") %>%
    dplyr::select("ecotype")
  
  
  hovedokosystemVector <- terra::project(vect(a), hotspotLayers)
  hovedokosystems <- terra::rasterize(hovedokosystemVector, hotspotLayers, field = "ecotype")
  translationTable <- data.frame(ints = 1:12, categories = c("Bebyggelse-samferdsel", "Dyrket mark", "Grasmark", "Skog",
                                                             "Hei og apen vegetasjon", "Lite vegetert mark", "Vatmark",  "Elver-bekker",
                                                             "","Innsjoer-tern", "Svaberg kyststrender og dyner", "Apent hav"))
  levels(hovedokosystems) <- translationTable
  names(hovedokosystems) <- "Hovedokosystems"
  writeRaster(hovedokosystems, "processes/overlays/data/hovedokosystemer/hovedokosystemer.tiff", overwrite = TRUE)
} else {
  hovedokosystems <- rast("processes/overlays/data/hovedokosystemer/hovedokosystemer.tiff")
}

###-----------------------###
### Inngrepsfrie omrader ####
###-----------------------###

# This one is also available online from Miljodirektoratet:
# https://kartkatalog.miljodirektoratet.no/Dataset/Details/100
inngrepsfrieOmrader <- read_sf("processes/overlays/data/inngrepsfrieOmrader/statusPolygon.shp")
ifVect <- terra::project(vect(inngrepsfrieOmrader), hotspotLayers)
ifRaster <- terra::rasterize(ifVect, hotspotLayers, field = "vsone")
translationTable <- data.frame(ints = 0:3, categories = c("Sone1", "Sone2", NA, "Vill"))
levels(ifRaster) <- translationTable
names(ifRaster) <- "intactAreas"

###-----------------------------------###
### Import CORINE land cover changes ####
###-----------------------------------###

# For this one you'll need to sign up to Copernicus' data service. one you've done that you'll find the file at
# the link below, which you'll need to unzip locally.
# https://land.copernicus.eu/en/products/corine-land-cover/lcc-2012-2018

#unzip("overlays/data/landChangeCorine/Results/u2018_cha1218_v2020_20u1_raster100m.zip", exdir = "overlays/data/landChangeCorine")
corineLandChange <- rast("processes/overlays/data/landChangeCorine/u2018_cha1218_v2020_20u1_raster100m/DATA/U2018_CHA1218_12_V2020_20u1.tif")
corineLandChangeProjected <- terra::project(corineLandChange, hotspotLayers)
corineLandChangeProjected <- crop(corineLandChangeProjected, hotspotLayers[[1]], mask = T)
names(corineLandChangeProjected) <- "corineLandChange"

###----------------------------------###
### Collate what we have up to here ####
###----------------------------------###

rasterPackage3 <- c(hotspotLayers, bias, protectedAreas, culturalAreas, natureTypes, municipalityRaster, 
                       countyRaster, waterRegions, waterAreas, belowGroundCarbonNorway, aboveGroundCarbonNorway,
                       hovedokosystems, ifRaster, corineLandChangeProjected)
namesToSet <- names(rasterPackage3)
writeRaster(rasterPackage3, file = "processes/overlays/data/rasterPackage3.tiff", overwrite = TRUE, names = namesToSet)


importCheck <- rast("overlays/data/rasterPackage3.tiff")
