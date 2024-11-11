###-----------------###
### 1. Import data ####
###-----------------###

print("Preparing data for model prediction.")

library(PointedSDMs)
library(dplyr)
library(purrr)
library(terra)
# library(tidyterra)
library(stringr)
library(intSDM)

# Import local functions
sapply(list.files("../BioDivMapping/functions", pattern = "\\.R$", full.names = TRUE), source)

# model output folder
modelFolderName <- paste0("localPredictions/data/modelObjects")

# import project control parameters into the environment
focalTaxa <- read.csv(file.path("data/focalTaxa.csv"), header = T) %>%
  filter(taxa %in% "vascularPlants")
predRes <- 250
taxaFolder <- paste0("localPredictions/data/predictions/resolution", predRes, "/")
for (folder in taxaFolder) {
  if (!file.exists(folder)) {
    dir.create(folder)
  }}

environmentalDataListAll <- rast("localPredictions/data/environmentalDataImported.tiff")

levels(environmentalDataListAll$land_cover_corine)[[1]][,2][is.na(levels(environmentalDataListAll$land_cover_corine)[[1]][,2])] <- "Water bodies"
levels(environmentalDataListAll$land_cover_corine)[[1]][,2][28] <- "Moors and heathland"
landCover <- environmentalDataListAll$land_cover_corine

values(environmentalDataListAll$land_cover_corine)[,1][is.nan(values(environmentalDataListAll$land_cover_corine)[,1])] <- 48
levels(environmentalDataListAll$land_cover_corine) <- levels(landCover)

environmentalDataListAll$summer_precipitation_squared <- environmentalDataListAll$summer_precipitation^2
environmentalDataListAll$summer_temperature_squared <- environmentalDataListAll$summer_temperature^2

# Import fitted models
models <- list.files(modelFolderName, pattern = "reducedModel", recursive = TRUE, full.names = TRUE)


mesh <- list(cutoff = 176, max.edge=c(26850, 175903) * 2, offset= c(1760, 1200) * 10)

# Import species which we need data for
fieldWorkResults <- readRDS("localPredictions/data/fieldWorkResults.RDS") %>%
  select(-survey)
speciesToRun <- colnames(fieldWorkResults)[2:435]

# Check extents of various groups
lapply(unique(fieldWorkResults$group), FUN = function(g) {
  ext(fieldWorkResults[fieldWorkResults$group == g,])
})

# Everything starts here
# Give list of locals
boxedLocations <-list(c(183000, 6852000, 190000, 6870000) |>
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

# Get env data in
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
  # covs <- rownames(model$summary.fixed) %>%
  #     stringr::str_subset(paste0("^(", paste(speciesIn, collapse = "|"), ")")) %>%
  #     stringr::str_remove(paste0("^(", paste(speciesIn, collapse = "|"), ")_")) %>% 
  #     unique
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
      replicate(length(speciesIn) + 1, ., simplify = FALSE) %>% 
      reduce(cbind) %>% 
      bind_cols(geometries[[p]]) %>% 
      st_sf() %>% 
      na.omit()
    
    # update names
    names(pData) <- c(covs,
                      paste(rep(speciesIn, each = length(covs)), 
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
  for(type in modelOutputs){
    ret <- lapply(seq_along(predData), FUN = function(preds) {
      split(1:nrow(predData[[preds]]), seq(1, ceiling(nrow(predData[[preds]]) / 5000)))
    })
    
    if(type == "Bias" && biasField) {
      
      biasJoined <- lapply(seq_along(predData), FUN = function(pd) {
        pred <- lapply( ret[[pd]], function(x) {
          biasEst <- predict(model, data = predData[[pd]][x,], bias = TRUE, mesh = mesh)
          return(biasEst)
        })
        spPred <- lapply(pred, function(x){
          res <-   x[[1]][[1]]
        })%>%
          do.call("rbind", .)%>%
          select("mean")
        spPred <- rasterize(spPred, predGrid[[pd]], names(spPred)[!names(spPred) %in% names(predData[[pd]])])
      })
      
      biasMerged <- do.call(merge, biasJoined)
      writeRaster(biasMerged, paste0("data/run_", dateSaved, "/", oneTaxa, "/",type,i,"_",predRes, ".tiff"), overwrite = TRUE)
      
    } else if(type == "Richness") {
      predictionsJoined <- lapply(seq_along(predData), FUN = function(pd) {
        pred <- lapply(ret[[pd]], function(x){
          richnessEst <- intSDM:::obtainRichness(model, predictionData = predData[[pd]][x, ], 
                                                 predictionIntercept = predictionInterceptDataset, sampleSize = focalSampleSize)
          richnessEst$Probabilities
        })
      })
      
      species <- names(predictionsJoined[[1]][[1]])
      for(sp in species){
        print(sp)
        predictionsMade <- lapply(seq_along(predictionsJoined) , FUN = function(pd) {
          spPred <- lapply(predictionsJoined, function(x){
            res <-   x[[1]][[sp]]
          })%>%
            do.call("rbind", .)%>%
            select("mean", "sd")
          spPred <- rasterize(spPred, predGrid[[pd]], names(spPred)[!names(spPred) %in% names(predData[[pd]])])
        })
        predictionsMerged <- do.call(merge, predictionsMade)
        writeRaster(predictionsMerged, paste0("localPredictions/data/predictions/resolution", predRes, "/", sp, ".tiff"), overwrite = TRUE)
        
      }
      
    }
    
  }
}

