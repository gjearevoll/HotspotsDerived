

#### MODEL RESIZING ####

# Ron has created a function which drastically reduces the size of the models without losing any of the key ingredients.
# We apply this to every model that comes out of the pipeline.

library(lobstr)
# We need to tidentify the correct directory for finding our models first.
hotspotsDir <- "../BioDivMapping/data/run_2025-01-06/modelOutputs"

# Import Ron's function
source("functions/reset_environments.R")

# Get every richness model listed
dirsWithModel <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds", full.names = TRUE)
dirsWithModelShort <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds")
taxaWithModels <- gsub(paste0(hotspotsDir, "/"), "", unique(gsub('[[:digit:]]+', '', sub("\\/.*", "", dirsWithModelShort))))
#taxaWithModels <- "spider"

# Now loop through each group and reduce the model size
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
