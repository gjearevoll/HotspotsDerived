### Covariate Analysis ####
# We need to get an overview of how each environmental covariate affects species groups as a whole


library(ggplot2)
library(dplyr)
library(terra)
source("functions/importRedList.R")

redList <- importRedList(c("VU", "EN", "CR"))
redList$simpleScientificName <- gsub(" ", "_", redList$species)
modelNameList <- list.files("localPredictions/data/modelObjects", full.names = TRUE)
ansvarsArter <- readRDS("data/ansvarsArterList.RDS")
threatenedSpecies <- readRDS("data/redList.RDS")

# Here we need to join all our covariate tables together in order to create density plots
modelList <- lapply(modelNameList, FUN = function(ml) {
  oneModel <- readRDS(ml)
  oneModelSummary <- oneModel$summary.fixed
  speciesUsed <- oneModel$species$speciesTable$species
  rm("oneModel")
  gc()
  
  tableSampleFull <- do.call(cbind, lapply(speciesUsed, FUN = function(sp) {
    tableSample <- oneModelSummary[grep(sp, row.names(oneModelSummary)),] 
    tableSample$var <- gsub(paste0(sp, "_"), "", row.names(tableSample))
    tableSample <- dplyr::select(arrange(tableSample, var), mean)
    tableSample
  }) 
  )
  colnames(tableSampleFull) <- speciesUsed
  tableSampleFull
})

# Now we have a list of all our covariate effects for each species, bind them together
modelsBound <- do.call(cbind,modelList)
row.names(modelsBound) <- gsub(paste0(colnames(modelsBound)[1], "_"), "" ,row.names(modelsBound))

# Now calculate densities across the different groupings
densities <- as.data.frame(do.call(cbind, lapply(1:nrow(modelsBound), FUN = function(r) {
  as.numeric(modelsBound[r,])
}) |> setNames(row.names(modelsBound)))) %>%
  mutate(set = "All species")

ansvarDensities <- as.data.frame(do.call(cbind, lapply(1:nrow(modelsBound), FUN = function(r) {
  threatenedTable <- modelsBound[,colnames(modelsBound) %in% ansvarsArter$simpleScientificName]
  as.numeric(threatenedTable[r,])
}) |> setNames(row.names(modelsBound)))) %>%
  mutate(set = "Ansvarsarter")

threatenedDensities <- as.data.frame(do.call(cbind, lapply(1:nrow(modelsBound), FUN = function(r) {
  threatenedTable <- modelsBound[,colnames(modelsBound) %in% redList$simpleScientificName]
  as.numeric(threatenedTable[r,])
}) |> setNames(row.names(modelsBound)))) %>%
  mutate(set = "Threatened species")

# Bind together
densityTest <- rbind(densities, ansvarDensities, threatenedDensities)
densityTestLong <- gather(densityTest, variable, effect, aspect:summer_temperature_squared, factor_key=TRUE)


#Plot
covToPlot <- "slope"
ggplot(densityTestLong[densityTestLong$variable == covToPlot,], aes(x = effect, fill = set)) + 
  geom_density(alpha = 0.5) +
  #xlim(-20,5) +
  theme_bw()

