# Quick function carrying over all finished data used to create raster packages 3 or 4

dateUsed <- "2025-01-06"
sourceDirectory <- paste0("../BioDivMapping/data/run_", dateUsed,"/modelOutputs/processedOutputs")
densitySourceDirectory <- paste0("../BioDivMapping/data/run_", dateUsed,"/samplingDensities")

filesToImport <- c(list.files(sourceDirectory, pattern = "final", recursive = FALSE, full.names = TRUE),
                   list.files(densitySourceDirectory, pattern = "density", recursive = FALSE, full.names = TRUE))
dirName <- paste0("data/modelOutputs/run_", dateUsed)
if (!dir.exists(dirName)) {
  dir.create(dirName)
}
file.copy(filesToImport, dirName)

# library(terra)
# lapply(list.files(paste0("data/modelOutputs/run_", dateUsed), full.names = T), FUN = function(x) {
#   ext(rast(x))
# }) |> setNames(list.files(paste0("data/modelOutputs/run_", dateUsed)))
