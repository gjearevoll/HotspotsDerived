#### Analyse overlap between different layers ####

library(terra)
library(dplyr)
library(ggplot2)
library(sf)



namesLayers <- c("allspecies", "ansvarsarter", "threatenedspecies")
norskNames <- c("Alle artar", "Ansvarsartar", "Trua artar")

sourceDate <- "2025-01-06"
sourceDirectory <- paste0("data/modelOutputs/run_",sourceDate)


### Get names of all taxa we currently have finished data for
taxaNames <- unique(gsub(".tiff", "",
                         gsub("final.tiff", "", 
                              gsub("allspeciesstats_", "", 
                                   list.files(sourceDirectory, pattern = "allspecies")))))

translations <- read.csv("data/taxaTranslations.csv")
translations$groupsNynorsk[translations$groups == "Insects"] <- "Insekt og edderkoppdyr"
# Groups to plot

groupsToPlot <- c(unique(translations$groups[translations$engelsk %in% taxaNames]), "birdGroups")

# Import species richness for a group
statsGrouped <- lapply(groupsToPlot, FUN = function(group) {
  if (group == "birdGroups") {
    groupNames <- c("woodpeckers", "groundNesters", "waders")
    rastImport <- rast(paste0("visualisations/figures/rawMetrics/",groupNames,"_allspeciesrichness.tiff" )) |>
      setNames(translations$nynorsk[match(groupNames, translations$engelsk)])
  } else {
    rastImport <- rast(paste0("visualisations/figures/rawMetrics/",group,"_",namesLayers,"richness.tiff" )) |>
      setNames(norskNames)
  }
}) |> setNames(groupsToPlot)


# Import layers to analyse
layers <- rast("processes/overlays/data/rasterPackage3.tiff") %>%
  project(statsGrouped[[1]])

###------------------###
### Protected Areas ####
###------------------###

for (group in groupsToPlot) {
  nynorskNavn <- if (group == "birdGroups") {"Fuglargruppar"} else {unique(translations$groupsNynorsk[translations$groups %in% group])}
  protectedComparison <- c(statsGrouped[[group]], layers$protectedAreas)
  protectedComparison$protectedAreas <- ifel(protectedComparison$protectedAreas %in% c("habitatSpeciesManagementArea",
                                                                                       "strictNatureReserve",
                                                                                       "protectedLandscapeOrSeascape",
                                                                                       "wildernessArea"), TRUE, FALSE)
  protectedComparisonDF <- as.data.frame(protectedComparison)
  protectedComparisonDF <- protectedComparisonDF[!is.na(protectedComparisonDF[,1]),]
  
  protectedComparisonList <- do.call(rbind, lapply(names(statsGrouped[[group]]), FUN = function(lyr) {
    relevantDF <- protectedComparisonDF[,c(lyr, "protectedAreas")]
    relevantDF$scaledRichness <- (protectedComparisonDF[,lyr]-min(protectedComparisonDF[,lyr]))/
      (max(protectedComparisonDF[,lyr])-min(protectedComparisonDF[,lyr]))*100
    relevantDF$set <- lyr
    relevantDF[,c("protectedAreas", "scaledRichness", "set")]
  }))
  protectedComparisonList$protectedStatus <- factor(ifelse(protectedComparisonList$protectedAreas == 1, 
                                                           "Naturvernområde", "Ikkje verna område"))
  groupedList <- protectedComparisonList %>%
    group_by(set, protectedStatus) %>%
    summarise(meanR = mean(scaledRichness, na.rm = TRUE),
              seR = sd(scaledRichness, na.rm = TRUE))
  groupedList$seR <- ifelse(groupedList$meanR - groupedList$seR >= 0, groupedList$seR, groupedList$meanR)
  
  #Produce plot
  barplotProtectedAreas <- ggplot(groupedList, aes(x=set, y=meanR/100, fill=protectedStatus)) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")  +
    geom_errorbar(aes(ymin=meanR/100-seR/100, ymax=meanR/100+seR/100), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_brewer(palette="Oranges")+
    theme_minimal() +
    labs(fill = "", y = "Skalert rikdom", x = "Artsgruppe", title = nynorskNavn) +
    ylim(0,1)
  ggsave(filename =  paste0("visualisations/figures/overlayAnalysis/barplots/protectedAreaComparison_",  group,".png"),
         plot = barplotProtectedAreas, units = "px",
         width = 2600, height = 1400)
}


###--------------------------###
### Protected Areas By Type ####
###--------------------------###

for (group in groupsToPlot) {
  nynorskNavn <- if (group == "birdGroups") {"Fuglargruppar"} else {unique(translations$groupsNynorsk[translations$groups %in% group])}
  protectedComparison <- c(statsGrouped[[group]], layers$protectedAreas)  
  protectedComparisonDF <- as.data.frame(protectedComparison)
  protectedComparisonDF <- protectedComparisonDF[!is.na(protectedComparisonDF[,1]),] %>%
    mutate(protectedAreas = as.character(protectedAreas))
  
  
  # Turn all NAs to Annet vern
  protectedComparisonDF$protectedAreas[is.na(protectedComparisonDF$protectedAreas)] <- "nonProtected"
  
  protectedComparisonList <- do.call(rbind, lapply(names(statsGrouped[[group]]), FUN = function(lyr) {
    relevantDF <- protectedComparisonDF[,c(lyr, "protectedAreas")]
    relevantDF$rikhet <- (protectedComparisonDF[,lyr]-min(protectedComparisonDF[,lyr]))/
      (max(protectedComparisonDF[,lyr])-min(protectedComparisonDF[,lyr]))*100
    relevantDF$set <- lyr
    relevantDF[,c("protectedAreas", "rikhet", "set")]
  }))
  
  # Only keep hovedokosystemer we care about
  verneomraderToUse <- c("Naturreservat","Landskapsvernområde", "Nasjonalpark", "Anna vern", "Anna vern", "Ikkje verna")
  nynorskVersion <- data.frame(nynorsk = verneomraderToUse, english = c("strictNatureReserve", "protectedLandscapeOrSeascape",
                                                                        "nationalPark", "naturalMonument", 
                                                                        "habitatSpeciesManagementArea", "nonProtected"))
  protectedComparisonList <- protectedComparisonList[protectedComparisonList$protectedAreas %in% nynorskVersion$english,]
  #hovedokosystemComparisonList$englishNames <- englishVersion$english[match(hovedokosystemComparisonList$Hovedokosystems, englishVersion$norsk )]
  protectedComparisonList$nynorskNames <- nynorskVersion$nynorsk[match(protectedComparisonList$protectedAreas, nynorskVersion$english )]
  orderDF <- data.frame(variable =c("Naturreservat","Nasjonalpark", "Landskapsvernområde", "Anna vern", "Ikkje verna"),
                        order = c(1,2,3,4,5))
  
  groupedList <- protectedComparisonList %>%
    group_by(set, nynorskNames) %>%
    summarise(meanR = mean(rikhet, na.rm = TRUE),
              seR = sd(rikhet, na.rm = TRUE))
  groupedList$seR <- ifelse(groupedList$meanR - groupedList$seR >= 0, groupedList$seR, groupedList$meanR)
  groupedList$seR <- ifelse(groupedList$meanR + groupedList$seR <= 100, groupedList$seR, 100-groupedList$meanR)
  groupedList$ordered <- orderDF$order[match(groupedList$nynorskNames, orderDF$variable)]
  
  #Produce plot
  barplotProtectedAreas <- ggplot(groupedList, aes(x=set, y=meanR/100, fill=reorder(nynorskNames,ordered))) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")  +
    geom_errorbar(aes(ymin=meanR/100-seR/100, ymax=meanR/100+seR/100), width=.2,
                  position=position_dodge(.7)) +
    scale_fill_brewer(palette="Greens", direction = -1)+
    theme_minimal()+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(fill = "", y = "Skalert rikdom", x = "", title = nynorskNavn) +
    ylim(0,1)
  
  
  ggsave(filename =  paste0("visualisations/figures/overlayAnalysis/barplots/protectedAreaComparisonGrouped_",  group,".png"),
         plot = barplotProtectedAreas, units = "px",
         width = 2600*0.8, height = 1400*0.7)
}

###---------------------###
### inngrepsfrie areas ####
###---------------------###

for (group in groupsToPlot) {
  nynorskNavn <- if (group == "birdGroups") {"Fuglargruppar"} else {unique(translations$groupsNynorsk[translations$groups %in% group])}
  
  
  wildernessComparison <- c(statsGrouped[[group]], layers$intactAreas)
  wildernessComparison$intactAreas <- ifel(wildernessComparison$intactAreas == "Vill", TRUE, FALSE)
  wildernessComparisonDF <- as.data.frame(wildernessComparison)
  wildernessComparisonDF <- wildernessComparisonDF[!is.na(wildernessComparisonDF[,1]),]
  
  wildernessComparisonList <- do.call(rbind, lapply(names(statsGrouped[[group]]), FUN = function(lyr) {
    relevantDF <- wildernessComparisonDF[,c(lyr, "intactAreas")]
    relevantDF$scaledRichness <- (wildernessComparisonDF[,lyr]-min(wildernessComparisonDF[,lyr]))/
      (max(wildernessComparisonDF[,lyr])-min(wildernessComparisonDF[,lyr]))*100
    relevantDF$set <- lyr
    relevantDF[,c("intactAreas", "scaledRichness", "set")]
  }))
  wildernessComparisonList$intactStatus <- factor(ifelse(wildernessComparisonList$intactAreas == 1, "Innanfor", "Utanfor"))
  wildernessComparisonList$intactStatus[is.na(wildernessComparisonList$intactStatus)] <- "Utanfor"
  
  groupedList <- wildernessComparisonList %>%
    group_by(set, intactStatus) %>%
    summarise(meanR = mean(scaledRichness, na.rm = TRUE),
              seR = sd(scaledRichness, na.rm = TRUE))
  groupedList$seR1 <- ifelse(groupedList$meanR - groupedList$seR >= 0, groupedList$seR, groupedList$meanR)
  groupedList$seR2 <- ifelse(groupedList$meanR + groupedList$seR <= 100, groupedList$seR, 100-groupedList$meanR)
  
  
  # barplotWildernessAreas <- ggplot(wildernessComparisonList, aes(x=set, y=scaledRichness/100, fill = intactStatus)) + 
  #   geom_boxplot(outlier.alpha = 0.01) + 
  #   scale_fill_brewer(palette="Blues") +
  #   theme_minimal() +
  #   labs(fill="Villmarksprega\nområde", y = "Skalert artsrikdom", x = "", title = nynorskNavn)
  barplotWildernessAreas <- ggplot(groupedList, aes(x=set, y=meanR/100, fill=intactStatus)) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")  +
    geom_errorbar(aes(ymin=meanR/100-seR1/100, ymax=meanR/100+seR2/100), width=.2,
                  position=position_dodge(.7)) +
    scale_fill_brewer(palette="Blues")+
    theme_minimal()+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(fill = "", y = "Skalert rikdom", x = "", title = nynorskNavn) +
    ylim(0,1)
  ggsave(filename =  paste0("visualisations/figures/overlayAnalysis/barplots/wildernessAreaComparison_",  group,".png"),
         plot = barplotWildernessAreas, units = "px",
         width = 2600*0.7, height = 1400*0.7)
}


###---------------------###
### Hovedokosystems areas ####
###---------------------###

for (group in groupsToPlot) {
  nynorskNavn <- if (group == "birdGroups") {"Fuglargruppar"} else {unique(translations$groupsNynorsk[translations$groups %in% group])}
  hovedokosystemComparison <- c(statsGrouped[[group]],layers$Hovedokosystems)
  hovedokosystemComparisonDF <- as.data.frame(hovedokosystemComparison)
  hovedokosystemComparisonDF <- hovedokosystemComparisonDF[!is.na(hovedokosystemComparisonDF[,1]),] %>%
    filter(!is.na(`Hovedokosystems`))
  
  
  hovedokosystemComparisonList <- do.call(rbind, lapply(names(statsGrouped[[group]]), FUN = function(lyr) {
    relevantDF <- hovedokosystemComparisonDF[,c(lyr, "Hovedokosystems")]
    relevantDF$rikhet <- hovedokosystemComparisonDF[,lyr]
    relevantDF$set <- lyr
    relevantDF[,c("Hovedokosystems", "rikhet", "set")]
  }))
  
  # Only keep hovedokosystemer we care about
  hovedokosystemsToUse <- c("Hei og apen vegetasjon","Lite vegetert mark", "Vatmark", "Skog", "Dyrket mark", "Grasmark")
  # englishVersion <- data.frame(norsk = hovedokosystemsToUse, english = c("Sparsely vegetated ecosystems", "Heathland and shrub",
  #                                                                        "Wetlands", "Forest and woodland", "Cropland", "Grassland"))
  nynorskVersion <- data.frame(norsk = hovedokosystemsToUse, nynorsk = c("Hei og open vegetasjon", "Lite vegetert mark",
                                                                         "Myr/våtmark", "Skog", "Dyrka mark", "Grasmark"))
  hovedokosystemComparisonList <- hovedokosystemComparisonList[hovedokosystemComparisonList$Hovedokosystems %in% hovedokosystemsToUse,]
  #hovedokosystemComparisonList$englishNames <- englishVersion$english[match(hovedokosystemComparisonList$Hovedokosystems, englishVersion$norsk )]
  hovedokosystemComparisonList$nynorskNames <- nynorskVersion$nynorsk[match(hovedokosystemComparisonList$Hovedokosystems, nynorskVersion$norsk )]
  
  groupedList <- hovedokosystemComparisonList %>%
    group_by(set, nynorskNames) %>%
    summarise(meanR = mean(rikhet, na.rm = TRUE),
              seR = sd(rikhet, na.rm = TRUE))
  groupedList$seR1 <- ifelse(groupedList$meanR - groupedList$seR >= 0, groupedList$seR, groupedList$meanR)
  groupedList$seR2 <- ifelse(groupedList$meanR + groupedList$seR <= 1, groupedList$seR, 1-groupedList$meanR)
  
  # barplotHovedokosystem <- ggplot(hovedokosystemComparisonList, aes(x=set, y=rikhet, fill = nynorskNames)) + 
  #   geom_boxplot(outlier.alpha = 0.01) + 
  #   scale_fill_brewer(palette="Set3") +
  #   theme_minimal() +
  #   labs(fill="Hovudøkosystem", y = "Skalert artsrikdom", x = "", title = norskNavn)
  barplotHovedokosystem <- ggplot(groupedList, aes(x=set, y=meanR, fill=nynorskNames)) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")  +
    geom_errorbar(aes(ymin=meanR-seR1, ymax=meanR+seR2), width=.2,
                  position=position_dodge(.7)) +
    scale_fill_brewer(palette="Set3")+
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(fill = "", y = "Skalert rikdom", x = "", title = nynorskNavn) +
    ylim(0,1)
  ggsave(filename =  paste0("visualisations/figures/overlayAnalysis/barplots/hovedokosystemComparison_",  nynorskNavn,".png"),
         plot = barplotHovedokosystem, units = "px",
         width = 2600*0.8, height = 1400*0.7)
}




###-------------------###
### NIN Nature areas ####
###-------------------###

# Have to import the areas for this one
# Locally downloaded file, can be found at http://kartkatalog.miljodirektoratet.no
natureTypesCats<- st_read("processes/overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml", layer = "Naturtype_nin") %>%
  dplyr::select(tilstand)
talliedGroups <- tally(group_by(natureTypesCats, tilstand)) %>%
  filter(!is.na(tilstand))


for (group in groupsToPlot) {
  nynorskNavn <- if (group == "birdGroups") {"Fuglargruppar"} else {unique(translations$groupsNynorsk[translations$groups %in% group])}
  relevantStats <- statsGrouped[[group]]
  
  
  tilstandComparisonList <- do.call(rbind, lapply(names(statsGrouped[[group]]), FUN = function(lyr) {
    relevantDF <- relevantStats[lyr]
    extractAttempt <- extract(relevantDF, talliedGroups, fun = "mean")
    extractAttemptSD <- extract(relevantDF, talliedGroups, fun = "sd")
    valueDF <- merge(extractAttempt, extractAttemptSD, by = "ID")
    valueDF$set <- lyr
    valueDF$tilstand <- st_drop_geometry(talliedGroups)$tilstand
    colnames(valueDF)[2:3] <- c("mean", "sd") 
    valueDF[,c("mean", "sd", "set", "tilstand")]
  }))
  
  nynorskVersion <- data.frame(norsk = st_drop_geometry(talliedGroups)$tilstand, nynorsk = c("Dårleg", "God", "Moderat", "Svært redusert"))
  tilstandComparisonList$nynorskNames <- nynorskVersion$nynorsk[match(tilstandComparisonList$tilstand, nynorskVersion$norsk )]
  orderDF <- data.frame(variable =nynorskVersion$nynorsk,
                        order = c(3,1,2,4))
  tilstandComparisonList$ordered <- orderDF$order[match(tilstandComparisonList$nynorskNames, orderDF$variable)]
  
  tilstandComparisonList$seR1 <- ifelse(tilstandComparisonList$mean - tilstandComparisonList$sd >= 0, tilstandComparisonList$sd, tilstandComparisonList$mean)
  tilstandComparisonList$seR2 <- ifelse(tilstandComparisonList$mean + tilstandComparisonList$sd <= 1, tilstandComparisonList$sd, 1-tilstandComparisonList$mean)
  
  # barplotHovedokosystem <- ggplot(hovedokosystemComparisonList, aes(x=set, y=rikhet, fill = nynorskNames)) + 
  #   geom_boxplot(outlier.alpha = 0.01) + 
  #   scale_fill_brewer(palette="Set3") +
  #   theme_minimal() +
  #   labs(fill="Hovudøkosystem", y = "Skalert artsrikdom", x = "", title = norskNavn)
  barplotNinTilstand <- ggplot(tilstandComparisonList, aes(x=set, y=mean, fill=reorder(nynorskNames,ordered))) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")  +
    geom_errorbar(aes(ymin=mean-seR1, ymax=mean+seR2), width=.2,
                  position=position_dodge(.7)) +
    scale_fill_brewer(palette="YlOrRd")+
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(fill = "Tilstand", y = "Skalert rikdom", x = "", title = nynorskNavn) +
    ylim(0,1)
  ggsave(filename =  paste0("visualisations/figures/overlayAnalysis/barplots/ninTilstand_",  nynorskNavn,".png"),
         plot = barplotNinTilstand, units = "px",
         width = 2600*0.7, height = 1400*0.7)
}


###----------------------###
### NIN Nature kvalitet ####
###----------------------###

# More complicated system here
# First need to import all areas and get their qualities
natureTypesKvalitet<- st_read("processes/overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml", layer = "Naturtype_nin") %>%
  dplyr::select(lokalitetskvalitet, prosjektnavn)

# Now need to find the projects NOT in this list, which should be the ones that area "kartlagt uten funn"
natureTypesDekningStatus <- st_read("processes/overlays/data/naturtyper/Naturtyper_nin_0000_norge_25833_GML.gml") %>%
  dplyr::select(dekningskartverdi, prosjektnavn) %>%
  filter(!(prosjektnavn %in% natureTypesKvalitet$prosjektnavn)) %>%
  filter(dekningskartverdi == "kartlagtUtenFunn")
natureTypesDekningStatus$lokalitetskvalitet <- "kartlagtUtenFunn"

natureTypes <- rbind(natureTypesKvalitet, natureTypesDekningStatus[,c("lokalitetskvalitet", "prosjektnavn")])

talliedGroupsKvalitet <- tally(group_by(natureTypes, lokalitetskvalitet)) %>%
  filter(!is.na(lokalitetskvalitet)) %>%
  filter(lokalitetskvalitet != "ikkeKvalitetsvurdert")


for (group in groupsToPlot) {
  nynorskNavn <- if (group == "birdGroups") {"Fuglargruppar"} else {unique(translations$groupsNynorsk[translations$groups %in% group])}
  relevantStats <- statsGrouped[[group]]
  
  
  kvalitetComparisonList <- do.call(rbind, lapply(names(statsGrouped[[group]]), FUN = function(lyr) {
    relevantDF <- relevantStats[lyr]
    extractAttempt <- extract(relevantDF, talliedGroupsKvalitet)
    colnames(extractAttempt)[2] <- "species"
    extractAttempt2 <- extractAttempt %>%
      group_by(ID) %>% summarise(mean = mean(species, na.rm = TRUE),
                                 sd = sd(species, na.rm = TRUE))
    extractAttempt2$set <- lyr
    extractAttempt2$kvalitet <- st_drop_geometry(talliedGroupsKvalitet)$lokalitetskvalitet
    extractAttempt2[,c("mean", "sd", "set", "kvalitet")]
  }))
  
  nynorskVersion <- data.frame(norsk = st_drop_geometry(talliedGroupsKvalitet)$lokalitetskvalitet, nynorsk = c("Høg", "Kartlagt uten funn", "Låg", "Moderat", "Svært høg", "Svært låg"))
  kvalitetComparisonList$nynorskNames <- nynorskVersion$nynorsk[match(kvalitetComparisonList$kvalitet, nynorskVersion$norsk )]
  orderDF <- data.frame(variable =nynorskVersion$nynorsk,
                        order = c(2,6,4,3,1,5))
  kvalitetComparisonList$ordered <- orderDF$order[match(kvalitetComparisonList$nynorskNames, orderDF$variable)]
  
  kvalitetComparisonList$seR1 <- ifelse(kvalitetComparisonList$mean - kvalitetComparisonList$sd >= 0, kvalitetComparisonList$sd, kvalitetComparisonList$mean)
  kvalitetComparisonList$seR2 <- ifelse(kvalitetComparisonList$mean + kvalitetComparisonList$sd <= 1, kvalitetComparisonList$sd, 1-kvalitetComparisonList$mean)
  
  # barplotHovedokosystem <- ggplot(hovedokosystemComparisonList, aes(x=set, y=rikhet, fill = nynorskNames)) + 
  #   geom_boxplot(outlier.alpha = 0.01) + 
  #   scale_fill_brewer(palette="Set3") +
  #   theme_minimal() +
  #   labs(fill="Hovudøkosystem", y = "Skalert artsrikdom", x = "", title = norskNavn)
  barplotNinKvalitet <- ggplot(kvalitetComparisonList, aes(x=set, y=mean, fill=reorder(nynorskNames,ordered))) +
    geom_bar(stat = "identity", width = 0.7, position=position_dodge(), colour = "#333333")  +
    geom_errorbar(aes(ymin=mean-seR1, ymax=mean+seR2), width=.2,
                  position=position_dodge(.7)) +
    scale_fill_manual(values=c(rev(brewer.pal(n = 5, name = "YlOrRd")), "gray"))+
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(fill = "Kvalitet", y = "Skalert rikdom", x = "", title = nynorskNavn) +
    ylim(0,1)
  ggsave(filename =  paste0("visualisations/figures/overlayAnalysis/barplots/ninKvalitet_",  nynorskNavn,".png"),
         plot = barplotNinKvalitet, units = "px",
         width = 2600*0.7, height = 1400*0.7)
}


