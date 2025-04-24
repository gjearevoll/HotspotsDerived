###-----------------###
### 1. Import data ####
###-----------------###

cat("Preparing data for model prediction.\n")

library(PointedSDMs)
library(dplyr)
library(purrr)
library(terra)
# library(tidyterra)
library(stringr)
library(intSDM)


args <- commandArgs(TRUE)

if (length(args) != 0) {
  # Set arguments
  predRes <- as.integer(args[1])
  # Set the working directory
  setwd("~/HotspotsDerived")
}

# predRes <- 250

# Import local functions
sapply(list.files("../BioDivMapping/functions", pattern = "\\.R$", full.names = TRUE), source)

# model output folder
modelFolderName <- paste0("processes/localPredictions/data/modelObjects")

# import project control parameters into the environment
focalTaxa <- read.csv(file.path("data/focalTaxa.csv"), header = T) %>%
  filter(taxa %in% "vascularPlants")
taxaFolder <- paste0("processes/localPredictions/data/predictions/resolution", predRes, "/")
for (folder in taxaFolder) {
  if (!file.exists(folder)) {
    dir.create(folder)
  }}

cat("Importing and editing environmental data\n")

environmentalDataListAll <- rast("processes/localPredictions/data/environmentalDataImported.tiff")

levels(environmentalDataListAll$land_cover_corine)[[1]][,2][is.na(levels(environmentalDataListAll$land_cover_corine)[[1]][,2])] <- "Water bodies"
levels(environmentalDataListAll$land_cover_corine)[[1]][,2][28] <- "Moors and heathland"
landCover <- environmentalDataListAll$land_cover_corine

values(environmentalDataListAll$land_cover_corine)[,1][is.nan(values(environmentalDataListAll$land_cover_corine)[,1])] <- 48
levels(environmentalDataListAll$land_cover_corine) <- levels(landCover)

environmentalDataListAll$summer_precipitation_squared <- environmentalDataListAll$summer_precipitation^2
environmentalDataListAll$summer_temperature_squared <- environmentalDataListAll$summer_temperature^2

# Import fitted models
models <- list.files(modelFolderName, pattern = "reducedModel", recursive = TRUE, full.names = TRUE)

cat("Creating mesh for model\n")
mesh <- list(cutoff = 176, max.edge=c(26850, 175903) * 2, offset= c(1760, 1200) * 10)

# Import species which we need data for
cat("Defining species to run\n")
fieldWorkResults <- readRDS("processes/localPredictions/data/fieldWorkResults.RDS") %>%
  select(-survey)
speciesMeans <- colMeans(st_drop_geometry(fieldWorkResults[,2:435]))
speciesToRun <- names(speciesMeans)[speciesMeans > 0.01]

# Check extents of various groups
groupExtents <- lapply(unique(fieldWorkResults$group), FUN = function(g) {
  ext(fieldWorkResults[fieldWorkResults$group == g,])
})

# Everything starts here
# Give list of locals
boxedLocations <-list(c(183000, 6852000, 190000, 6870000) |> #7 by 18
                        setNames(c("west", "south", "east", "north")),
                      c(218000, 6585000, 245000, 6600000) |> #27 by 15
                        setNames(c("west", "south", "east", "north")),
                      c(-60000, 6640000, -40000, 6670000) |> #20 by 30
                        setNames(c("west", "south", "east", "north")),
                      c(233000, 7070000, 260000, 7100000) |> #27 by 30
                        setNames(c("west", "south", "east", "north")),
                      c(334000, 6951000, 345000, 6990000) |> #10 by 40
                        setNames(c("west", "south", "east", "north")),
                      c(260000, 7024000, 295000, 7040000) |> #35 by 16
                        setNames(c("west", "south", "east", "north")))

# Get env data in
cat("Importing and cropping environmental data\n")
environmentalDataList <- lapply(boxedLocations, FUN = function(box) {
  regionGeometry <- defineRegion("box", extentCoords =  box)
  lapply(environmentalDataListAll, FUN = function(x) {
    crop(x, st_transform(regionGeometry, st_crs(environmentalDataListAll))) 
  })|>  
    rast()
}
)

###-----------------###
### 2. Prep objects ###
###-----------------###
# res = 2
# original covariate names
origCovs <- names(environmentalDataList[[1]])
# which factor
types <- sapply(seq(nlyr(environmentalDataList[[1]])), function(x){
  environmentalDataList[[1]][[x]][1] %>% unlist %>% class
})
# define template prediction raster
cat("Defining prediction grid\n")
predGrid <- lapply(environmentalDataList, FUN = function(envData) {
  
  predRast <-  terra::rast(ext(envData), res = c(predRes, predRes), crs = crs(st_crs(25833)$wkt))
  # Define prediction raster grid at target resolution
  if(any(types == "factor")){
    predGridfactor <- envData[[types == "factor"]] 
    # Define prediction raster grid for continuous covs (interpolate when predRes <= res, else average) 
    predGrid <- terra::project(envData[[types != "factor"]], predRast, 
                               method = if(predRes <= 1000) "bilinear" else "average") 
    # Define prediction raster grid catagorical covs
    factorRasters <- terra::project(predGridfactor, predRast, method = "mode") 
    levels(factorRasters) <- levels(predGridfactor) # reassign levels 
    # Combine binary rasters into a SpatRaster object
    predGrid <- c(predGrid, factorRasters) 
    #names(predGrid) <- origCovs
  } else {
    # Define prediction raster grid (interpolate when predRes <= res, else average) 
    predGrid <- terra::project(envData, predRast, method = if(predRes <= res) "bilinear" else "average") 
  }
})

# define geometries to combine with prediction 
geometries <- lapply(predGrid, FUN = function(pGrid) {
  xyFromCell(pGrid, seq(ncell(pGrid))) %>% 
    as.data.frame() %>% 
    st_as_sf(coords = c("x", "y"), crs = 25833) 
})
# Define model outputs based on modelRun
modelOutputs <- 'Richness'

###-------------------------###
### 3. Generate predictions ###
###-------------------------###

# timeStart <- Sys.time()

# Have done up to 138

cat("Starting model run")
for(i in seq_along(models)){
  
  predictionInterceptDataset <- "ANOData"
  focalSampleSize <- 0.25
  # import model
  model <- readRDS(models[i])
  if (!(predictionInterceptDataset %in% names(model$dataType))) next
  # identify species in model
  speciesIn <- model$species$speciesIn %>% unlist %>% unique
  speciesUsed <- speciesIn[speciesIn %in% speciesToRun]
  if (length(speciesUsed) == 0) {
    cat("no species at run", i) 
    next} 
  # indentify if bias field
  biasField <- !is.null(model$summary.random$sharedBias_biasField) | !is.null(model$spatCovs$biasFormula)
  # identify covariates used in model 
  covs <- model$spatCovs$name
  # identify categorical covariate factors 
  catCovCats <- model$summary.random[model$summary.random %>% names %>% 
                                       stringr::str_subset(paste0("^(", paste(speciesIn, collapse = "|"), ")"))]  %>% 
    sapply(function(cov){
      cov[,1]
    }) %>% unlist %>% #names %>% 
    stringr::str_remove(paste0("^(", paste(speciesIn, collapse = "|"), ")_")) %>% unique %>% 
    str_subset(paste0("^(", str_c(names(environmentalDataList[[1]][[types == "factor"]]), collapse = "|"), ")"))
  catCovs <- origCovs[sapply(origCovs, function(name) {
    any(str_detect(catCovCats, paste0("^", name)))
  })]
  
  # get complete list of covariate columns from which to predict
  covs <- unique(c(covs, catCovs))
  
  # prep prediction data 
  predData <- lapply(seq_along(predGrid), FUN = function(p) {
    pGrid <- predGrid[[p]]
    pData <- pGrid %>% 
      dplyr::select(all_of(covs)) %>% 
      as.data.frame(na.rm = FALSE) %>% 
      replicate(length(speciesUsed) + 1, ., simplify = FALSE) %>% 
      reduce(cbind) %>% 
      bind_cols(geometries[[p]], .name_repair = "unique_quiet") %>% 
      st_sf() %>% 
      na.omit()
    
    # update names
    names(pData) <- c(covs,
                      paste(rep(speciesUsed, each = length(covs)), 
                            covs, sep = "_"), "geometry")
    return(pData)
  })
  # identify bias covs
  if(!is.null(model$spatCovs$biasFormula)){
    biasCovs <- covs[covs %in% attributes(terms(model$spatCovs$biasFormula))$term.labels]
    # remove bias covariates from list
    covs <- covs[!covs %in% attributes(terms(model$spatCovs$biasFormula))$term.labels]
  }
  
  # Generate & convert & Save model/predicts (currently for all species grouped)
  
  ret <- lapply(seq_along(predData), FUN = function(preds) {
    split(1:nrow(predData[[preds]]), seq(1, ceiling(nrow(predData[[preds]]) / 5000)))
  })
  predictionsJoined <- lapply(seq_along(predData), FUN = function(pd) {
    pred <- lapply(ret[[pd]], function(x){
      richnessEst <- intSDM:::obtainRichness(model, predictionData = predData[[pd]][x, ], 
                                             predictionIntercept = predictionInterceptDataset, sampleSize = focalSampleSize)
      richnessEst$Probabilities
    })
  })
  
  for(sp in speciesUsed){
    print(sp)
    predictionsMade <- lapply(seq_along(predictionsJoined) , FUN = function(pd) {
      spPred <- lapply(predictionsJoined[[pd]], function(x){
        res <-   x[[sp]]
      })%>%
        do.call("rbind", .)%>%
        select("mean", "sd")
      spPred <- rasterize(spPred, predGrid[[pd]], names(spPred)[!names(spPred) %in% names(predData[[pd]])])
    })
    predictionsMerged <- do.call(merge, predictionsMade)
    writeRaster(predictionsMerged, paste0("processes/localPredictions/data/predictions/resolution", predRes, "/", sp, ".tiff"), overwrite = TRUE)
  }
  
}
