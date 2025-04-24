
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###

library(dplyr)
library(ggplot2)
library(ggridges)
library(stringr)

source("functions/importRedList.R")
ansvarsArter <- readRDS("data/ansvarsArterList.RDS")

redList <- importRedList(c("VU", "EN", "CR"))
redList$simpleScientificName <- gsub(" ", "_", redList$species)
# Now we have a list of all our covariate effects for each species, bind them together

# Import species richness for a group
namesLayers <- c("Alle artar", "Trua artar", "Ansvarsartar")

if (!dir.exists("visualisations/figures/covariateAnalysis")) {
  dir.create("visualisations/figures/covariateAnalysis")
}

# Bring all relevant taxa together
finalFiles <- list.files("processes/covariateAnalysis", pattern = "CovariateList", recursive = TRUE)
availableTaxa <- c("vascularPlants", "birds", "fungi", "insects", "lichens")
covariateTypes <- c("continuous", "corine", "kalkinnhold")

variableTranslations <- read.csv("data/variableTranslations.csv")
taxaTranslations <- read.csv("data/taxaTranslations.csv") 
taxaTranslations$groupsNynorsk[taxaTranslations$engelsk == "insects"] <- "Insekt og edderkoppdyr"

# For bird groups
birdNames <- read.csv("../BioDivMapping/data/external/birdTypeList.csv")
birdNames$simpleName2 <- gsub(" ", "_", birdNames$simpleName)

###---------------------------------------###
### 1. Collate and set data for plotting ####
###---------------------------------------###

for (taxa in availableTaxa) {
  norskNavn <- taxaTranslations$groupsNynorsk[taxaTranslations$engelsk %in% taxa]
  
  for (cov in covariateTypes) {
    if (cov == "kalkinnhold" & taxa == "birds") {next}
    
    indModels <- readRDS(paste0("processes/covariateAnalysis/", cov, "Effects/", 
                                taxa, "CovariateList.RDS"))
    row.names(indModels) <- gsub("^[^_]*_", "", gsub("^[^_]*_", "", row.names(indModels)))
    
    # Now calculate densities across the different groupings
    densities <- as.data.frame(do.call(cbind, lapply(1:nrow(indModels), FUN = function(r) {
      as.numeric(indModels[r,])
    }) |> setNames(row.names(indModels)))) %>%
      mutate(set = "Alle artar")
    
    ansvarDensities <- as.data.frame(do.call(cbind, lapply(1:nrow(indModels), FUN = function(r) {
      threatenedTable <- indModels[,colnames(indModels) %in% ansvarsArter$simpleScientificName]
      as.numeric(threatenedTable[r,])
    }) |> setNames(row.names(indModels)))) %>%
      mutate(set = "Ansvarsartar")
    
    threatenedDensities <- as.data.frame(do.call(cbind, lapply(1:nrow(indModels), FUN = function(r) {
      threatenedTable <- indModels[colnames(indModels) %in% redList$simpleScientificName]
      as.numeric(threatenedTable[r,])
    }) |> setNames(row.names(indModels)))) %>%
      mutate(set = "Trua artar")
    
    # Bind together
    densityTest <- rbind(densities, ansvarDensities, threatenedDensities)
    densityTestLong <- tidyr::gather(densityTest, variable, effect, 1:(ncol(densityTest)-1), factor_key=TRUE) %>%
      filter(!is.na(effect))
    pointShift <- median(densityTestLong$effect)
    densityTestLong$effectShift <- densityTestLong$effect - pointShift
    getRange <- quantile(densityTestLong$effectShift, c(0.01,0.99))
    getRange <- c(-40, 20)
    
    # Need to make changes specific to each data type
    if (cov == "kalkinnhold") {
      # Set up groups in logical order
      densityTestLong$variable <- plyr::revalue(densityTestLong$variable, c("svæk intermediær" = "Svæk intermediær", 
                                                                            "temmelig kalkfattig" = "Temmelig kalkfattig",
                                                                            "svært kalkfattig" = "Svært kalkfattig",
                                                                            "svært kalkrik" = "Svært kalkrik"
      ))
      orderDF <- data.frame(variable = c("Svæk intermediær", "Svært kalkfattig", "Svært kalkrik", "Temmelig kalkfattig"),
                            order = c(3,1,4,2))
      densityTestLong$ordered <- orderDF$order[match(densityTestLong$variable, orderDF$variable)]
    } else if (cov == "continuous") {
      # Polish covariate names
      variableName <- gsub("_", " ", densityTestLong$variable)
      densityTestLong$variable <- variableTranslations$covariateNorsk[match(densityTestLong$variable, 
                                                                            gsub(" ","_",variableTranslations$covariateEnglish))]
    } else if (cov == "corine") {
      # Polish names and get rid of the one empty factor level
      densityTestLong$variable <- variableTranslations$covariateNorsk[match(densityTestLong$variable, variableTranslations$covariateEnglish)]
      densityTestLong <- densityTestLong[densityTestLong$variable != "Naturleg grasmark",]
    }
    
    ###--------------------------###
    ### 2. Create species plots ####
    ###--------------------------###
    
    # If continuous covariates
    if (cov == "continuous") {
      densityPlot <- ggplot(densityTestLong, aes(x = effectShift, y = variable, linetype = set)) + 
        geom_density_ridges(alpha = 0.1) +
        scale_linetype_manual(values=c("solid", "longdash", "dotted"))+
        xlim(floor(getRange[1]),ceiling(getRange[2])) +
        theme_bw() +
        labs(linetype = "", x = "Effekt av variabel", y = "",
             title = paste0(toupper(substr(norskNavn, 1, 1)), substr(norskNavn, 2, nchar(norskNavn))))
      ggsave(filename =  paste0("visualisations/figures/covariateAnalysis/covariateAnalysis_",  taxa,"_", cov, ".jpg"),
             plot = densityPlot, height = 1600, width = 2300, units = "px")
    } else if (cov == "kalkinnhold") {
      # If factorial covariates
      densityPlot <- ggplot(densityTestLong, aes(y = reorder(variable,ordered), x = effectShift, linetype = set)) + 
        geom_density_ridges(alpha = 0.1) +
        scale_linetype_manual(values=c("solid", "longdash", "dotted"))+
        xlim(floor(getRange[1]),ceiling(getRange[2])) +
        theme_bw() +
        theme(text = element_text(size = 12),
              axis.text.y = element_text(hjust = 1, size = 8),
              legend.text = element_text(size = 8))  +
        labs(linetype ="", x = "Effekt av variabel", y = "",
             title = paste0(toupper(substr(norskNavn, 1, 1)), substr(norskNavn, 2, nchar(norskNavn))))
      ggsave(filename =  paste0("visualisations/figures/covariateAnalysis/covariateAnalysis_",  taxa,"_", cov, ".jpg"),
             plot = densityPlot, height = 1600, width = 2300, units = "px")
    } else {
      densityPlot <- ggplot(densityTestLong, aes(y = variable, x = effectShift, linetype = set)) + 
        geom_density_ridges(alpha = 0.1) +
        scale_linetype_manual(values=c("solid", "longdash", "dotted"))+
        xlim(floor(getRange[1]),ceiling(getRange[2])) +
        theme_bw() +
        theme(text = element_text(size = 12),
              axis.text.y = element_text(hjust = 1, size = 8),
              legend.text = element_text(size = 8))  +
        labs(linetype ="", x = "Effekt av variabel", y = "",
             title = paste0(toupper(substr(norskNavn, 1, 1)), substr(norskNavn, 2, nchar(norskNavn))))
      ggsave(filename =  paste0("visualisations/figures/covariateAnalysis/covariateAnalysis_",  taxa,"_", cov, ".jpg"),
             plot = densityPlot, height = 1600, width = 2300, units = "px")
    }
  }}





###---------------------------------------###
### 1. Collate and set data for plotting ####
###---------------------------------------###


for (cov in c("continuous", "corine")) {
  
  indModels <- readRDS(paste0("processes/covariateAnalysis/", cov, "Effects/birdsCovariateList.RDS"))
  mgmtGroups <- c("woodpeckers", "waders", "groundNesters")
  row.names(indModels) <- gsub("^[^_]*_", "", gsub("^[^_]*_", "", row.names(indModels)))
  
  
  # Now calculate densities across the different groupings
  waderModels <- indModels[,colnames(indModels) %in% birdNames$simpleName2[birdNames$group == "waders"]]
  waders <- as.data.frame(do.call(cbind, lapply(1:nrow(waderModels), FUN = function(r) {
    as.numeric(waderModels[r,])
  }) |> setNames(row.names(waderModels)))) %>%
    mutate(set = "Vadefuglar")
  
  groundNesterModels <- indModels[,colnames(indModels) %in% birdNames$simpleName2[birdNames$group == "groundNesters"]]
  groundNesters <- as.data.frame(do.call(cbind, lapply(1:nrow(groundNesterModels), FUN = function(r) {
    as.numeric(groundNesterModels[r,])
  }) |> setNames(row.names(groundNesterModels)))) %>%
    mutate(set = "Bakkehekkende fuglar")
  
  hakkespettModels <- indModels[,colnames(indModels) %in% birdNames$simpleName2[birdNames$group == "woodpeckers"]]
  hakkespett <- as.data.frame(do.call(cbind, lapply(1:nrow(hakkespettModels), FUN = function(r) {
    as.numeric(hakkespettModels[r,])
  }) |> setNames(row.names(hakkespettModels)))) %>%
    mutate(set = "Hakkespettar")
  
  # Bind together
  densityTest <- rbind(waders, groundNesters, hakkespett)
  densityTestLong <- tidyr::gather(densityTest, variable, effect, 1:(ncol(densityTest)-1), factor_key=TRUE) %>%
    filter(!is.na(effect))
  pointShift <- median(densityTestLong$effect)
  densityTestLong$effectShift <- densityTestLong$effect - pointShift
  getRange <- quantile(densityTestLong$effectShift, c(0.01,0.99))
  getRange <- c(-40, 20)
  
  # Need to make changes specific to each data type
  if (cov == "continuous") {
    # Polish covariate names
    variableName <- gsub("_", " ", densityTestLong$variable)
    densityTestLong$variable <- variableTranslations$covariateNorsk[match(densityTestLong$variable, 
                                                                          gsub(" ","_",variableTranslations$covariateEnglish))]
  } else if (cov == "corine") {
    # Polish names and get rid of the one empty factor level
    densityTestLong$variable <- variableTranslations$covariateNorsk[match(densityTestLong$variable, variableTranslations$covariateEnglish)]
    densityTestLong <- densityTestLong[densityTestLong$variable != "Naturleg grasmark",]
  }
  
  ###--------------------------###
  ### 2. Create species plots ####
  ###--------------------------###
  
  # If continuous covariates
  if (cov == "continuous") {
    densityPlot <- ggplot(densityTestLong, aes(x = effectShift, y = variable, linetype = set)) + 
      geom_density_ridges(alpha = 0.1) +
      scale_linetype_manual(values=c("solid", "longdash", "dotted"))+
      xlim(floor(getRange[1]),ceiling(getRange[2])) +
      theme_bw() +
      labs(linetype = "", x = "Effekt av variabel", y = "",
           title = "Fuglarguppar")
    ggsave(filename =  paste0("visualisations/figures/covariateAnalysis/covariateAnalysis_birdGroups_", cov, ".jpg"),
           plot = densityPlot, height = 1600, width = 2300, units = "px")
  } else {
    densityPlot <- ggplot(densityTestLong, aes(y = variable, x = effectShift, linetype = set)) + 
      geom_density_ridges(alpha = 0.1) +
      scale_linetype_manual(values=c("solid", "longdash", "dotted"))+
      xlim(floor(getRange[1]),ceiling(getRange[2])) +
      theme_bw() +
      theme(text = element_text(size = 12),
            axis.text.y = element_text(hjust = 1, size = 8),
            legend.text = element_text(size = 8))  +
      labs(linetype ="", x = "Effekt av variabel", y = "",
           title = "Fuglargruppar")
    ggsave(filename =  paste0("visualisations/figures/covariateAnalysis/covariateAnalysis_birdGroups_", cov, ".jpg"),
           plot = densityPlot, height = 1600, width = 2300, units = "px")
  }
}
