
###----------------------------------------###
### 0. Import relevant libraries and data ####
###----------------------------------------###
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(tidyterra)

###--------------------------------###
### 1. Define relevant parameters ####
###--------------------------------###

# Visualisation
focalGroup <- "waders"
chosenResolution <- "10"
focalBeta <- rast(paste0("processes/betaDiversity/data/",focalGroup, chosenResolution, "betaDiversity.tiff"))


translations <- read.csv("data/taxaTranslations.csv")
translations$groupsNynorsk[translations$engelsk == "insects"] <- "Insekt og edderkoppdyr"
focalGroupNorsk <-  translations$nynorsk[translations$engelsk == focalGroup]

norwayBorder2 <- sf::read_sf("../BioDivMapping/data/external/norge_border/Noreg_polygon.shp")
norwayBorderProjected <- terra::project(vect(norwayBorder2), focalBeta)
focalBetaCropped <- crop(focalBeta, norwayBorderProjected, mask = TRUE)

# # If cropped version necessary
# cropX <- c(600, 800)
# cropY <- c(7600, 7800)
# focalBetaCropped <- crop(focalBetaCropped, ext(c(cropX, cropY)))

###-----------------###
### 2. Create plot ####
###-----------------###

betaPlot <- ggplot(norwayBorderProjected) +
  
  geom_spatraster(focalBetaCropped, mapping = aes(fill = sum), maxcell = 5e+07
  ) +
  #geom_sf(fill = NA, lwd = 0, colour = "black") + 
  geom_sf(fill = NA, lwd = 0.1, colour = "black") + 
  coord_sf(datum = st_crs(focalBetaCropped)#,
           #xlim = cropX, ylim = cropY
  ) +
  scale_fill_grass_c(palette = "forest_cover", na.value = NA, limits = c(0, 1), direction = 1)+
  guides(
    fill = guide_colourbar(position = "inside")
  )  +
  theme_void() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(vjust = -10),
        legend.position.inside = c(0.8, 0.2)
  ) +
  labs(fill = "Beta-\ndiversitet", title = focalGroupNorsk
  )
ggsave(betaPlot, filename = paste0("visualisations/figures/",focalGroup,"_betaDiversity.png" ), units = "px",
       width = 1600, height = 1600)
