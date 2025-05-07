landingPageText1 <- function() {
  mainPanel(
    h1("What is the Open Data Biodiversity Mapper?"),
    p("The Hotspots project was commissioned in early 2023 to map Norwegian biodiversity, in an effort to enable data-driven 
      decision making regarding Norwegian conservation. The main output of the project is a series of datasets depicting 
      various aspects of Norwegian biodiversity on a national scale, at a resolution of 500 metres by 500 metres."),
    p("This is the Open Biodiversity Mapper, a Shiny app constructed by the Hotspots team at the Gjaervoll Centre. This tool 
      has two primary aims:"),
    p("     1. Communication of the possibilities of use for the data produced from the Hotspots project."),
    p("     2. A summarised overview of the data production, and all inherent caveats involved in its use."),
    p("If you have any questions regarding the data, please do not hesitate to get in touch with Sam Wenaas Perrin at the Gjaervoll
      Centre.")
  )
}

landingPageText2 <- function() {
  mainPanel(width = 12,
    p("All code used in data production can be found on GitHub at the ", strong(a("repository BioDivMapping",
                                                                                  href = "https://github.com/gjearevoll/BioDivMapping",
                                                                                  target = "_blank")),
                                                                                  ". There is a description of the 
      pipeline both on the repositories homepage and in the file masterscript.R."),
    p("The models that are used to produce our estimates in the Hotspot project are Integrated Species Distribution Models (iSDMs), 
      so-called because they integrate different types of data, accounting for different collection methods in doing so. Every 
      dataset has its own unique biases and errors, and using one overarching model for every dataset simultaneously fails to 
      account for these. Our ISDM assigns a unique sampling model to each significantly different dataset, before then taking 
      these results and putting them through a more standard species distribution model."),
    p("These models are of course more complex than I can describe here, but you can get a comprehensive overview of how ISDMs work 
      in ", strong(a("Philip Mostert's paper on the iSDM package",
                     href = "https://www.biorxiv.org/content/10.1101/2022.09.15.507996v3",
                     target = "_blank")),
      ".")
  )
}

speciesRichnessText1<- function() {
  mainPanel(width = 12,
    p("The final species richness estimate you see in the figure above you is a scaled version of the sum of the estimated likelihoods of occurrence 
      for all species evaluated at a given pixel. For instance, if species 1 has a likelihood of occurrence of 0.7, species 2 a 
      likelihood of 0.5 and species 3 of 0.05, the resulting species richness metric for this pixel will be 1.25. This number is then
      scaled to between 0 and 1, as summing likelihoods does not give a viable total richness estimate, instead an indiction of richness relative
      to other pixels, which we suggest is best represented by a scaled metric."),
    p("As with almost all tabs on this app, you can choose to switch to a definition of species richness which only includes threatened
    species or species of national responsibility ('Ansvaraster' in Norwegian).")
  )
}

speciesRichnessText2 <- function() {
  mainPanel(width = 12,
    p("Uncertainty here is represented by the standard error.")
  )
}

speciesRichnessText3 <- function() {
  mainPanel(width = 12,
    p("The level of precision here is lower than as calculated by our data. While we have produced estimates on a 500 x 500 metre level, 
      here they have been aggregated to a 5 by 5 kilometre level in order to not to slow the Shiny app down to a snail’s pace.")
  )
}

betaDiversityText1 <- function() {
  mainPanel(width = 12,
    p("Beta diversity here is represented by looking at how much of the regional species pool we find in any given pixel. 
      The regional species pool is here represented with a 20 by 20 kilometre box, centred on the pixel. We determined whether a 
      species is present based on whether it has more than a 70% likelihood of presence or not, and add these presences to
      determine regional species richness. We then divide the species richness  of our focal pixel by the regional species richness.
      This gives us a beta diversity value between 0 and 1.",
      strong("Note"), ": This is an arbitrary value.")
  )
}

betaDiversityText2 <- function() {
  mainPanel(width = 12,
    p("While this metric makes sense on a small scale, obviously given the general 
      geometric shape of the Norway national boundary the relevance of this metric will vary from region to region. The regional 
      species richness for a 150 by 150km area in the Dovre mountains will naturally draw from a much larger land area than a pixel 
      in northern Nordland. Thus, beta diversity needs to be considered on relevant regional scales. Here we have used a 20 by 20 
      kilometre scale so that the relevance of the metric displayed does not diminished too much for any one region, though
      areas close to Norway's borders will still have slightly lower accuracies than those further inland."))
} 

individualSpeciesText1 <- function() {
  mainPanel(width = 12,
    p("Using the multi-species richness model, we are able to estimate species-level likelihood of occurrence. This requires taking 
      the species-specific linear predictor (consisting of species-specific effects) and transforming them using the inverse of the 
      link function considered for the structured data. We considered a Poisson point process framework to estimate the model; therefore 
      the link function required to connect the structured data to this model is a cloglog link function (Miller et al. 2019)."),
    p("Estimating the true probability of occurrence in a cell is difficult given that the intensity function of the point process model 
      is adjusted by a thinning process, which is confounded in the intercept term of the mode (Fithian and Hastie, 2013). To account for 
      this, we scale the intensity function of our model by the sampling area of one of the structured datasets where we assume biases 
      are minimal.")
  )
} 

hotspotsText1 <- function() {
  mainPanel(width = 12,
    p("Defining a hotspot (an area of high biodiversity) is a completely arbitrary business and completely dependent on the goals of a
      natural resource manager. Which types of species you’re looking to conserve, their red-listed status, how likely their presence 
      is at a given site and many more factors will determine the definition of a hotspot for someone looking to use this data for 
      conservation purposes."),
    p("Here we've shown hotspots as belonging to 3 different categories, with another (level 1, beige) showing pixels not belonging to hotspots.
      Level 2 (orange) pixels represents species whose scaked species richness lands them in the 90th to 94th percentile of all 
      pixels. Level 3 represents species in the 95th to 98th percentile of all pixels, and level 4 (black) shows pixels 
      only in the 99th percentile."),
    p("This map also gives the user the option to see the intersection of wilderness areas in Norway (those more than four kilometres
      from human infrastructure) and hotspots. The 'Select Overlay' option allows you to turn this on, while the 'Minimum size' slider
      gives you the option to remove wilderness areas below a certain total area from the map.")
  )
} 

environmentalCovariateText <- function() {
  mainPanel(width = 12,
    p("The environmental data used to run our ISDMs was pulled from a range of external sources, which can be found in more detail  ",
      strong(a("at this link", href = "https://github.com/gjearevoll/BioDivMapping/tree/main/data/external/environmentalCovariates", target = "_blank")),
    "The data seen here has (as with the species richness data) been aggregated to a 5 by 5 kilometre level.")
)
}