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
ddList <- readRDS("data/ddList.RDS")
namesFromGroups <- readRDS("data/namesFromGroups.RDS")

choiceList <- c("Norway", "50", "54") |> 
  setNames(c("Norway", "Trondelag", "Finnmark"))
taxaChoiceList <- c("hymenopterans", "birds", "vascularPlants") |> 
  setNames(c("Hymenopterans", "Birds", "Vascular Plants"))

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
                  menuItem("Home", tabName = "home", icon = icon("house")),
                  menuItem("Species Richness", tabName = "speciesRichness", icon = icon("binoculars")),
                  menuItem("Species Beta Diversity", tabName = "betaDiversity", icon = icon("worm")),
                  menuItem("Individual Species Data", tabName = "individuals", icon = icon("cloud-sun")),
                  menuItem("Hotspots", tabName = "hotspots", icon = icon("cloud-sun")),
                  menuItem("Metadata", tabName = "meta", icon = icon("cloud-sun"))
      )
    )
    ,
    
    dashboardBody(
      tabItems(
        tabItem(tabName = "home",
                fluidRow(box(landingPageText1(), width = 12)),
                fluidRow(box(landingPageText2(), title = "Technical Resources", width = 12, 
                             collapsible = TRUE, collapsed = TRUE))),
        
        tabItem(tabName = "speciesRichness",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Select species",
                                selectInput(inputId = "taxa", label = "Taxa:",
                                            selected = "hymenopterans",
                                            choices = taxaChoiceList),
                                selectInput(inputId = "speciesGroup", label = "Species to use:",
                                            selected = "allspecies",
                                            choices = c("All species" = "allspecies",
                                                        "Threatened species" = "threatenedspecies",
                                                        "Species of national responsibility" = "ansvarsarter")),
                                selectInput(inputId = "statistic", label = "Statistic:",
                                            selected = "mean",
                                            choices = c("Species richness" = "skalertRikhet",
                                                        "Uncertainty" = "skalertUsikkerhet"))))
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
                          title = "What does species richness mean?",
                          speciesRichnessText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    ),
                    fluidRow(
                      box(width = 12,
                          title = "How do you measure uncertainty?",
                          speciesRichnessText2(), collapsible = TRUE,
                          collapsed = TRUE)
                    ),
                    fluidRow(
                      box(width = 12,
                          title = "How precise is this data?",
                          speciesRichnessText3(), collapsible = TRUE,
                          collapsed = TRUE)
                    ))
                )
        ),
        
        tabItem(tabName = "betaDiversity",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Select species group",
                                selectInput(inputId = "taxaBeta", label = "Taxa:",
                                            selected = "birds",
                                            choices = taxaChoiceList[-1]),
                                selectInput(inputId = "regionBeta", label = "Region",
                                            selected = "Norway",
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
                          title = "What does beta diversity mean?",
                          betaDiversityText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    ),
                    fluidRow(
                      box(width = 12,
                          title = "A note regarding beta diversity in Norway",
                          betaDiversityText2(), collapsible = TRUE,
                          collapsed = TRUE)
                    ))
                )
        ),
        
        tabItem(tabName = "individuals",
                fluidRow(
                  column( width = 3, offset = 0,
                          fluidRow(
                            box(width = 12,title = "Select species group",
                                selectInput(inputId = "taxaIndivid", label = "Taxa:",
                                            selected = "hymenopterans",
                                            choices = taxaChoiceList),
                                selectInput(inputId = "species", label = "Species:",
                                            selected = "Bombus alpinus (Linnaeus, 1758)",
                                            choices = ddList[[1]]),
                                selectInput(inputId = "regionIndivid", label = "Region:",
                                            selected = "Norway",
                                            choices = choiceList),
                                selectInput(inputId = "statisticIndivid", label = "Statistic:",
                                            selected = "speciesRichness",
                                            choices = c("Species richness" = "probability",
                                                        "Uncertainty" = "uncertainty"))))
                          ,
                          fluidRow(
                            box(width = 12, title = "Species Info",
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
                          title = "How do we calculate individual species probabilities?",
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
                                selectInput(inputId = "taxaHotspots", label = "Taxa:",
                                            selected = "beetles",
                                            choices = taxaChoiceList),
                                selectInput(inputId = "speciesGroupHotspots", label = "Species to use:",
                                            selected = "allspecies",
                                            choices = c("All species" = "allspecies",
                                                        "Threatened species" = "threatenedspecies",
                                                        "Species of national responsibility" = "ansvarsarter")),
                                selectInput(inputId = "regionHotspots", label = "Region:",
                                            selected = "Norway",
                                            choices = choiceList))),
                          fluidRow(
                            box(collapsible = TRUE, collapsed = TRUE, width = 12, title = "Select overlay",
                                selectInput(inputId = "overlayHotspots", label = "Overlay:",
                                            selected = "none",
                                            choices = c("none", "intactAreas")),
                                sliderInput(inputId = "overlayLimitHotspots", label = "Minimum size of area to include", 
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
                          title = "Defining a hotspot",
                          hotspotsText1(), collapsible = TRUE,
                          collapsed = TRUE)
                    ))
                )
        ),
        tabItem(tabName = "meta",
                fluidRow(column(width =3,
                                box(width = 12,title = "Select species",
                                    selectInput(inputId = "metaTaxa", label = "Taxa:",
                                                selected = "hymenopterans",
                                                choices = taxaChoiceList))),
                                
                                column(9, htmlOutput("inc")
                                ))
                )
                
                
        )
      )
    )
  )
  
  
  