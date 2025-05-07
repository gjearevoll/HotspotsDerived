
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###

library(ggplot2)
library(dplyr)
library(terra)
library(tidyr)
library(ggridges)

source("functions/importRedList.R")

args <- commandArgs(TRUE)

if (length(args) != 0) {
  # Set arguments
  taxaList <- args[1]
  # Set the working directory
  setwd("~/HotspotsDerived")
}

sourceDirectory <- "../BioDivMapping/data/run_2025-01-06/modelOutputs"

# Get necessary species lists
redList <- importRedList(c("VU", "EN", "CR"))
redList$simpleScientificName <- gsub(" ", "_", redList$species)
modelNameList <- list.files(sourceDirectory, full.names = TRUE, recursive = TRUE, pattern = "richnessModel")
#taxaList <- unique(gsub('[[:digit:]]+', "", list.dirs(sourceDirectory, full.names = FALSE, recursive = FALSE)))
ansvarsArter <- readRDS("data/ansvarsArterList.RDS")
threatenedSpecies <- readRDS("data/redList.RDS")

insectList <- c("aquaticInsects", "beetles", "bugs", "butterfliesMoths", "cockroaches", "earwigs", "flies", 
                "grasshopperLikeInsects", "hymenopterans", "netWingedInsects", "scorpionflies", "snakeflies",
                "spiders")

insectDF <- do.call(cbind, lapply(insectList, FUN = function(x) {
  dfToUse <- readRDS(paste0("processes/covariateAnalysis/kalkinnholdEffects/", x, "CovariateList.RDS"))
  dfToUse
}))
saveRDS(insectDF, "processes/covariateAnalysis/kalkinnholdEffects/insectCovariateList.RDS")

###------------------------------------------###
### 1. Produce continuous covariate effects ####
###------------------------------------------###

# Here we need to join all our covariate tables together in order to create density plots
continuousCovariateList <- lapply(taxaList, FUN = function(taxa) {
  if (taxa == "insects") {
    modelNameListTaxa <-  grep(paste(insectList, collapse="|"), modelNameList, value = TRUE) 
  } else {
    modelNameListTaxa <- grep(paste0("/", taxa), modelNameList, value = TRUE)
  }
  
  segmentListed <- do.call(cbind, lapply(modelNameListTaxa, FUN = function(segment) {
    modelName <- segment
    cat("\nRunning model ",gsub("../BioDivMapping/data/run_2025-01-06/modelOutputs/", "", segment))
    oneModel <- readRDS(modelName)
    oneModelSummary <- oneModel$summary.fixed
    speciesUsed <- oneModel$species$speciesTable$species
    rm("oneModel")
    gc()
    
    tableSampleFull <- lapply(speciesUsed, FUN = function(sp) {
      if (stringr::str_count(sp, "_") > 1) {return(NA)}
      tableSample <- oneModelSummary[sub("^(([^_]*_){1}[^_]*).*", "\\1", row.names(oneModelSummary)) == sp,] 
      tableSample$var <- gsub(paste0(sp, "_"), "", row.names(tableSample))
      tableSample <- dplyr::select(arrange(tableSample, var), mean)
      tableSample
    }) 
    tableSampleFull2 <- do.call(cbind, tableSampleFull[!is.na(tableSampleFull)])
    #row.names(tableSampleFull)
    colnames(tableSampleFull2) <- speciesUsed[!is.na(tableSampleFull)]
    tableSampleFull2
  }))
  fileName <- paste0("processes/covariateAnalysis/continuousEffects/", taxa, "CovariateList.RDS")
  saveRDS(segmentListed, fileName)
  print(taxa)
  taxa
}) |> setNames(taxaList)

###------------------------------###
### 2. Produce group covariates ####
###------------------------------###

# Define categories
catNamesCorine <- c("Bare rocks", "Broad-leaved forest", "Built up area", "Coniferous forest","Constructed green space", 
                    "Glaciers and perpetual snow", "Marsh/bog/fen", "Mixed forest", "Moors and heathland", "Natural grasslands",
                    "Transitional woodland-shrub", "Water bodies") 
catNamesKalk <- c("svæk intermediær", "svært kalkfattig", "svært kalkrik", "temmelig kalkfattig") 

# Collate effects
corineCovariateList <- lapply(taxaList, FUN = function(taxa) {
  modelNameListTaxa <- modelNameList[grep(paste0("/", taxa), modelNameList)]
  
  # Produce a list with all effects of covariate factor levels
  segmentListed <- lapply(modelNameListTaxa, FUN = function(segment) {
    
    cat("\nRunning model ",segment)
    modelName <- segment
    oneModel <- readRDS(modelName)
    oneModelSummaryCorine <- oneModel$summary.random[grep("land_cover_corine", names(oneModel$summary.random))]
    oneModelSummaryKalk <- oneModel$summary.random[grep("kalkinnhold", names(oneModel$summary.random))]
    speciesUsed <- oneModel$species$speciesTable$species
    rm("oneModel")
    gc()
    
    # Put all effects in one table for this model
    if (length(oneModelSummaryCorine) == 0) { 
      tableSampleCorine <- NULL} else {
        tableSampleCorine <- do.call(cbind, lapply(speciesUsed, FUN = function(sp) {
          tableSample <- oneModelSummaryCorine[grep(sp, names(oneModelSummaryCorine))][[1]]
          row.names(tableSample) <- tableSample$ID
          tableSample <- merge(data.frame(catNamesCorine), tableSample, all.x=TRUE, by.x = "catNamesCorine", by.y = "ID")
          tableSample["mean"]
        })) %>% as.data.frame()
        rownames(tableSampleCorine) <- catNamesCorine
        colnames(tableSampleCorine) <- speciesUsed    
      }
    
    # Same for kalkinnhold
    if (length(oneModelSummaryKalk) == 0) {
      tableSampleKalk <- NULL} else { 
        tableSampleKalk <- do.call(cbind, lapply(speciesUsed, FUN = function(sp) {
          tableSample <- oneModelSummaryKalk[grep(sp, names(oneModelSummaryKalk))][[1]]
          row.names(tableSample) <- tableSample$ID
          tableSample <- merge(data.frame(catNamesKalk), tableSample, all.x=TRUE, by.x = "catNamesKalk", by.y = "ID")
          tableSample["mean"]
        })) %>% as.data.frame()
        
        rownames(tableSampleKalk) <- catNamesKalk
        colnames(tableSampleKalk) <- speciesUsed
      }
    
    tableSampleFull <- list(corine = tableSampleCorine, kalkinnhold = tableSampleKalk)
    tableSampleFull
  })
  
  # Clean up and finalise the lists 
  segmentListedKalk <- do.call(cbind, lapply(segmentListed, FUN = function(kalk) {
    kalk$kalkinnhold
  }))
  segmentListedKalk <- segmentListedKalk[,!duplicated(colnames(segmentListedKalk))]
  segmentListedCorine <- do.call(cbind, lapply(segmentListed, FUN = function(corine) {
    corine$corine
  }))
  segmentListedCorine <- segmentListedCorine[,!duplicated(colnames(segmentListedCorine))]
  
  # Produce file names and save data
  fileNameCorine <- paste0("processes/covariateAnalysis/corineEffects/", taxa, "CovariateList.RDS")
  fileNameKalkinnhold <- paste0("processes/covariateAnalysis/kalkinnholdEffects/", taxa, "CovariateList.RDS")
  saveRDS(segmentListedCorine, fileNameCorine)
  saveRDS(segmentListedKalk, fileNameKalkinnhold)
  print(taxa)
  taxa
})

