
###----------------------###
### 0. Bash preparation ####
###----------------------###
args <- commandArgs(trailingOnly = TRUE)

start <- Sys.time()

dateUsed <- as.character(args[1])
taxaWithModels <- as.character(args[2])
setwd("~/HotspotsDerived")

#### MODEL RESIZING ####

# Ron has created a function which drastically reduces the size of the models without losing any of the key ingredients.
# We apply this to every model that comes out of the pipeline.

library(lobstr)
# We need to tidentify the correct directory for finding our models first.
hotspotsDir <- paste0("../BioDivMapping/data/run_",dateUsed,"/modelOutputs")

# Import Ron's function
source("../HotspotsDerived/functions/reset_environments.R")

# Get every richness model listed
dirsWithModel <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds", full.names = TRUE)
dirsWithModelShort <- list.files(hotspotsDir, recursive = TRUE, pattern = "*richnessModel.rds")
#taxaWithModels <- gsub(paste0(hotspotsDir, "/"), "", unique(gsub('[[:digit:]]+', '', sub("\\/.*", "", dirsWithModelShort))))
#taxaWithModels <- "vascularPlantsB"

# Now loop through each group and reduce the model size
for (taxa in taxaWithModels) {
  dirsWithModelTaxa <- dirsWithModel[grep(paste0("/",taxa), dirsWithModel)]
  for (i in seq_along(dirsWithModelTaxa)) {
    modelName <- dirsWithModelTaxa[[i]]
    if (file.info(modelName)[1,"size"] > 300000000) {
      richnessModel <- readRDS(modelName)
      #obj_size(richnessModel)
      reducedModel <- reset_environments(richnessModel)
      saveRDS(reducedModel, modelName)
    }
    
    if (i %% 3 == 0) cat(i, " models complete for ", taxa, "\n")
  }
}
