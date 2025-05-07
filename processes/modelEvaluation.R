

###------------------###
### 1. Data loading ####
###------------------###

sourceDirectory <- "../BioDivMapping/data/run_2025-01-06/modelOutputs"

modelNameList <- list.files(sourceDirectory, full.names = TRUE, recursive = TRUE, pattern = "richnessModel")

taxaToRun <- c("fungi", "insects", "vascularPlants", "lichens", "birds")
insectGroups <- unique(read.csv("../BioDivMapping/data/run_2025-01-02/focalTaxa.csv")$taxa)

###----------------------###
### 2. Get data metrics ####
###----------------------###

fullMetrics <- lapply(taxaToRun, FUN = function(taxa) {
  
  # Collapse insects into one category if need bed
  if (taxa == "insects") {
    relevantModels <- grep(paste0(insectGroups, collapse = "|"), modelNameList, value = T)
  } else {
    relevantModels <- grep(taxa, modelNameList, value = T)
  }
  
  speciesTable <- data.frame(varSpat = as.character(), 
                             varFixedCont = as.character(),
                             varFixedCorine = as.character(),
                             varFixedKalk = as.character(),
                             varBias = as.character(),
                             species = as.character())
  
  # Process models one by one (this takes forever)
  speciesRanges <- c()
  for (models in relevantModels) {
    cat("\nDownloading", gsub(sourceDirectory, "", models))
    model1 <- readRDS(models)
    
    # Get species used
    speciesUsed <- model1$species$speciesTable$species
    
    # Get variable effects for each species
    speciesMetrics <- do.call(rbind, lapply(speciesUsed, FUN = function(x) {
      
      # Variation explained by continuous variables
      varFixedCont <- model1$summary.fixed[grepl(x, row.names(model1$summary.fixed)), 'sd']^2 
      varFixedCont <- sum(varFixedCont)
      
      # Variation  explained by bias
      varBias <- sum((model1$summary.random$Bias__Effects__Comps$sd)^2)
      
      # Variation explained by Kalkinnhold (not relevant for birds)
      if (taxa != "birds") {
        tableKalk <- model1$summary.random[grepl(x, names(model1$summary.random)) & 
                                             grepl("kalkinnhold", names(model1$summary.random))][[1]]
        tableKalkFiltered <- tableKalk[tableKalk$mean < 30,]
        varFixedKalk <- sum(tableKalkFiltered$sd^2)
      } else {
        varFixedKalk <- 0
      }
      
      # Variation explained by Corine
      tableCorine <- model1$summary.random[grepl(x, names(model1$summary.random)) & 
                                             grepl("corine", names(model1$summary.random))][[1]]
      tableCorineFiltered <- tableCorine[tableCorine$mean < 30,]
      varFixedCorine <- sum(tableCorineFiltered$sd^2)
      
      # Variation explained by spatial fields
      varSpat <- model1$summary.hyper[2, 'mean']^2
      
      c(varSpat, varFixedCont, varFixedCorine, varFixedKalk, varBias, x)
    }))
    speciesTable <- rbind(speciesTable, speciesMetrics)
  }
  
  # Collate and save all data
  speciesTable$V5 <- as.numeric(speciesTable$V5)
  speciesTable$V4 <- as.numeric(speciesTable$V4)
  speciesTable$V3 <- as.numeric(speciesTable$V3)
  speciesTable$V2 <- as.numeric(speciesTable$V2)
  speciesTable$V1 <- as.numeric(speciesTable$V1)
  colnames(speciesTable) <- c("varSpat", "varFixedCont", "varFixedCorine", "varFixedKalk", "varBias", "species")
  saveRDS(speciesTable, paste0("data/modelOutputs/modelStats/vars_", taxa, ".RDS"))
  saveRDS(speciesRanges, paste0("data/modelOutputs/modelStats/ranges_", taxa, ".RDS"))
  cat("\nRun stats for", taxa)
  speciesTable
  
}) |> setNames(taxaToRun)


stop("Script done")


###-----------------------###
### 3. Collate this data ####
###-----------------------###

results <- as.data.frame(do.call(rbind, lapply(taxaToRun, FUN = function(x) {
  
  results <- readRDS(paste0("data/modelOutputs/modelStats/vars_", x, ".RDS"))
  results$totalVar <- apply(results[,1:5], 1, sum)
  results$varFixed <- apply(results[,2:4], 1, sum)
  results$varExplained <- results$varFixed/results$totalVar
  meanResults <- mean(results$varExplained)
  sdResults <- sd(results$varExplained)
  seResults <- sd(results$varExplained)/sqrt(length(results$varExplained))
  results2 <- c(meanResults, sdResults, seResults)
  names(results2) <- c("mean", "sd", "se")
  results2
})))
results$species <- taxaToRun
saveRDS(results, "modelEvalutaion.RDS")

# 
# ### Measure spatial autocorrelation
# 
# sourceDirectory <- "../BioDivMapping/data/run_2025-01-06/modelOutputs"
# 
# modelNameList <- list.files(sourceDirectory, full.names = TRUE, recursive = TRUE, pattern = "richnessModel")
# 
# taxaToRun <- c("aquaticInsects", "beetles", "hymenopterans", "spiders", "fungi", "lichens", "birds", "flies", "bugs", "butterfliesMoths",
#                "cockroaches", "earwigs", "grasshopperLikeInsects", "netWingedInsects", "scorpionflies", "snakeflies")
# #taxaToRun <- "vascularPlants"
# 
# relevantModels <- grep("beetles", modelNameList, value = T)
# model1 <- readRDS(relevantModels[1])
# model2 <- readRDS(relevantModels[2])
# dataCheck <- inlabru::spde.posterior(result = model1, name = "speciesShared", what = 'matern.correlation')
# 
# dataCheck2 <- inlabru::spde.posterior(result = model2, name = "speciesShared", what = 'matern.correlation')
# 
# p<-ggplot(data=dataCheck, aes(x=x, y=q0.5)) + geom_line()
# p<-p+geom_ribbon(aes(ymin=dataCheck$q0.025, ymax=dataCheck$q0.975), linetype=2, alpha=0.1)

