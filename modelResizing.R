library(lobstr)
# We need to transfer every vascular plants model from BioDivMapping over to here
hotspotsDir <- "../BioDivMapping/data/run_2025-01-06/modelOutputs"

source("functions/reset_environments.R")

# Get every file in there with vascular plants involved
dirsWithModel <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds", full.names = TRUE)
dirsWithModelShort <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds")
taxaWithModels <- gsub(paste0(hotspotsDir, "/"), "", unique(gsub('[[:digit:]]+', '', sub("\\/.*", "", dirsWithModelShort))))
#taxaWithModels <- "birds"

for (taxa in taxaWithModels) {
  dirsWithModelTaxa <- dirsWithModel[grep(paste0("/",taxa), dirsWithModel)]
  for (i in seq_along(dirsWithModelTaxa)) {
    modelName <- dirsWithModelTaxa[[i]]
    richnessModel <- readRDS(modelName)
    #obj_size(richnessModel)
    reducedModel <- reset_environments(richnessModel)
    saveRDS(reducedModel, modelName)
    if (i %% 10 == 0) cat(i, " models complete for ", taxa, "\n")
  }
}
