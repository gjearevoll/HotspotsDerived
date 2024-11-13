library(lobstr)
# We need to transfer every vascular plants model from BioDivMapping over to here
hotspotsDir <- "../BioDivMapping/data/run_2024-10-11/modelOutputs"

source("functions/reset_environments.R")

# Get every file in there with vascular plants involved
dirsWithModel <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds", full.names = TRUE)
taxaWithModels <- gsub(paste0(hotspotsDir, "/"), "", unique(gsub('[[:digit:]]+', '', sub("\\/.*", "", taxaWithModels))))


for (taxa in taxaWithModels) {
  dirsWithModelTaxa <- dirsWithModel[grep(taxa, dirsWithModel)]
  for (i in seq_along(dirsWithModelTaxa)) {
    modelName <- dirsWithModelTaxa[[i]]
    richnessModel <- readRDS(modelName)
    #obj_size(richnessModel)
    reducedModel <- reset_environments(richnessModel)
    saveRDS(reducedModel, modelName)
    if (i %% 10 == 0) cat(i, " models complete")
  }
}
  
