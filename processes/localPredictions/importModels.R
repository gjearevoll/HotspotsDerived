
# We need to transfer every vascular plants model from BioDivMapping over to here
hotspotsDir <- "../BioDivMapping/data/run_2024-10-11/modelOutputs"

# Get every file in there with vascular plants involved
dirsWithModel <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds", full.names = TRUE)
dirsWithModelPlants <- dirsWithModel[grep("vascularPlants", dirsWithModel)]
newNames <- gsub("vascularPlants", "", gsub("/", "", gsub(hotspotsDir, "", dirsWithModelPlants)))

file.copy(dirsWithModelPlants, paste0("localPredictions/data/modelObjects/", newNames))
