#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


relevantPackages <- c("shiny", "shinydashboard", "shinyjs", "sf", "ggplot2", "raster", "ggtext", "terra")
installedPackages <- rownames(installed.packages())
uninstalledPackages <- relevantPackages[!(relevantPackages %in% installedPackages)]

library(shiny)
library(shinydashboard)
library(shinyjs)
library(sf)
library(ggplot2)
library(raster)
library(ggtext)
library(stringr)
library(terra)

# # Import dropdown list for taxa
ddList <- readRDS("data/importedData/ddList.RDS")
namesFromGroups <- readRDS("data/importedData/namesFromGroups.RDS")

choiceList <- c("Norway", "50", "54", "46", "34") |> 
  setNames(c("Norge", "Trøndelag", "Finnmark", "Vestland", "Innlandet"))
taxaChoiceList <- c("insects", "birds", "vascularPlants", "lichens", "fungi") |> 
  setNames(c("Insekt", "Fugler", "Karplanter", "Lav", "Sopp"))

terraOptions(memmax = 0.9)

source("data/textFunctions.R")


# Define UI for application that draws a histogram
shinyUI(
  dashboardPage(
    
    # HEADER
    dashboardHeader(title = "Hotspots", titleWidth = 400),
    
    # SIDEBAR
    dashboardSidebar(
      sidebarMenu(id = "sidebar",
                  menuItem("Hjemmeside", tabName = "home", icon = icon("house")),
                  menuItem("Artsrikhet", tabName = "speciesRichness", icon = icon("binoculars")),
                  menuItem("Betadiversitet", tabName = "betaDiversity", icon = icon("worm")),
                  menuItem("Individuelle artsforekomst", tabName = "individuals", icon = icon("cloud-sun")),
                  menuItem("Hotspots", tabName = "hotspots", icon = icon("cloud-sun")),
                  menuItem("Metadata", tabName = "meta", icon = icon("cloud-sun"))
      )
    )
    ,
    
    dashboardBody(
      tabItems(
        tabItem(tabName = "home",
                fluidRow(box(landingPageText1(), width = 12)),
                fluidRow(box(landingPageText2(), title = "Teknisk ressurser", width = 12, 
                             collapsible = TRUE, collapsed = TRUE))),
        
        tabItem(tabName = "speciesRichness",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Velg artsgruppe",
                                selectInput(inputId = "taxa", label = "Artsgruppe:",
                                            selected = "birds",
                                            choices = taxaChoiceList),
                                selectInput(inputId = "speciesGroup", label = "Artsutvalg:",
                                            selected = "allspecies",
                                            choices = c("Alle arter" = "allspecies",
                                                        "Trua arter" = "threatenedspecies",
                                                        "Ansvarsarter" = "ansvarsarter")),
                                selectInput(inputId = "statistic", label = "Statistikk:",
                                            selected = "mean",
                                            choices = c("Artsrikdom" = "skalertRikhet",
                                                        "Usikkerhet" = "skalertUsikkerhet")),
                                selectInput(inputId = "region", label = "Region:",
                                            selected = "Norge",
                                            choices = choiceList)))
                  ),
                  column(
                    width = 5, offset = 0,
                    fluidRow(
                      box(width = 12,title = NULL,
                          status = "primary",
                          plotOutput("speciesMap", height = '100%'),
                          p(),
                          textOutput("figureCaption1"))),
                    fluidRow(
                      box(width = 12,
                          title = "Hva betyr ‘artsrikhet’?",
                          speciesRichnessText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    ),
                    fluidRow(
                      box(width = 12,
                          title = "Hvordan måles usikkerhet?",
                          speciesRichnessText2(), collapsible = TRUE,
                          collapsed = TRUE)
                    ),
                    fluidRow(
                      box(width = 12,
                          title = "Hvor nøyaktig er disse produktene?",
                          speciesRichnessText3(), collapsible = TRUE,
                          collapsed = TRUE)
                    ))
                )
        ),
        
        tabItem(tabName = "betaDiversity",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Velg artsgruppe",
                                selectInput(inputId = "taxaBeta", label = "Artsgruppe:",
                                            selected = "birds",
                                            choices = taxaChoiceList[-1]),
                                selectInput(inputId = "regionBeta", label = "Region:",
                                            selected = "Norge",
                                            choices = choiceList)))
                  ),
                  column(
                    width = 5, offset = 0,
                    fluidRow(
                      box(width = 12,title = NULL,
                          status = "primary",
                          plotOutput("betaMap", height = '100%'),
                          p(),
                          textOutput("figureCaption1"))),
                    fluidRow(
                      box(width = 12,
                          title = "Hva betyr betadiversitet?",
                          betaDiversityText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    ),
                    fluidRow(
                      box(width = 12,
                          title = "Om betadiversitet i Norge",
                          betaDiversityText2(), collapsible = TRUE,
                          collapsed = TRUE)
                    ))
                )
        ),
        
        tabItem(tabName = "individuals",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Velg artsgruppe",
                                selectInput(inputId = "taxaIndivid", label = "Artsgruppe:",
                                            selected = "birds",
                                            choices = taxaChoiceList),
                                selectInput(inputId = "species", label = "Art:",
                                            selected = "Calcarius lapponicus (Linnaeus, 1758)",
                                            choices = ddList[[1]]),
                                selectInput(inputId = "regionIndivid", label = "Region:",
                                            selected = "Norway",
                                            choices = choiceList),
                                selectInput(inputId = "statisticIndivid", label = "Statistikk:",
                                            selected = "speciesRichness",
                                            choices = c("Forekomst sannsynlighet" = "probability",
                                                        "Usikkerhet" = "uncertainty"))))
                          ,
                          fluidRow(
                            box(width = 12, title = "Artsinformasjon",
                                collapsible = TRUE,
                                htmlOutput("textBox1"),
                                p(),
                                imageOutput("imageBox1"))
                          )
                  ),
                  column(
                    width = 5, offset = 0,
                    fluidRow(
                      box(width = 12,title = NULL,
                          status = "primary",
                          plotOutput("individuals", height = '100%'),
                          p(),
                          textOutput("figureCaption1"))),
                    fluidRow(
                      box(width = 12,
                          title = "Hvordan beregnes sannsynlighet for artsforekomst?",
                          individualSpeciesText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    )
                  )
                )
        ),
        
        tabItem(tabName = "hotspots",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Select species",
                                selectInput(inputId = "taxaHotspots", label = "Artsgruppe:",
                                            selected = "birds",
                                            choices = taxaChoiceList),
                                selectInput(inputId = "speciesGroupHotspots", label = "Artsutvalg:",
                                            selected = "allspecies",
                                            choices = c("Alle arter" = "allspecies",
                                                        "Trua arter" = "threatenedspecies",
                                                        "Ansvarsarter" = "ansvarsarter")),
                                selectInput(inputId = "regionHotspots", label = "Region:",
                                            selected = "Norge",
                                            choices = choiceList))),
                          fluidRow(
                            box(collapsible = TRUE, collapsed = TRUE, width = 12, title = "Velg ekstra lag",
                                selectInput(inputId = "overlayHotspots", label = "Lag:",
                                            selected = "none",
                                            choices = c("None" = "none", 
                                                        "Inngrepsfrie områder" = "intactAreas")),
                                sliderInput(inputId = "overlayLimitHotspots", label = "Minste areal", 
                                            min = 0, max = 4500, value = 0))
                          )
                  ),
                  column(
                    width = 5, offset = 0,
                    fluidRow(
                      box(width = 12,title = NULL,
                          status = "primary",
                          plotOutput("hotspots", height = '100%'),
                          p(),
                          textOutput("figureCaption1"))),
                    fluidRow(
                      box(width = 12,
                          title = "Definisjon av en Hotspot",
                          hotspotsText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    ))
                )
        ),
        tabItem(tabName = "meta",
                fluidRow(column(width =3,
                                box(width = 12,title = "Select species",
                                    selectInput(inputId = "metaTaxa", label = "Taxa:",
                                                selected = "birds",
                                                choices = taxaChoiceList))),
                         
                         column(9, htmlOutput("inc")
                         ))
        )
        
        
      )
    )
  )
)


