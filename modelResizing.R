library(lobstr)
# We need to transfer every vascular plants model from BioDivMapping over to here
hotspotsDir <- "../BioDivMapping/data/run_2024-10-11/modelOutputs"

source("functions/reset_environments.R")

# Get every file in there with vascular plants involved
dirsWithModel <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds", full.names = TRUE)
dirsWithModelPlants <- dirsWithModel[grep("vascularPlants", dirsWithModel)]


for (i in seq_along(dirsWithModelPlants)[92]) {
  modelName <- dirsWithModelPlants[[i]]
  richnessModel <- readRDS(modelName)
  #obj_size(richnessModel)
  reducedModel <- reset_environments(richnessModel)
  rm("richnessModel")
  fileName <- paste0("localPredictions/data/modelObjects/reducedModel", i, ".rds")
  saveRDS(reducedModel, fileName)
  if (i %% 10 == 0) cat(i, " models complete")
}

