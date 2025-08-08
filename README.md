---
title: "README"
author: "Sam Perrin"
date: "2024-11-08"
output: html_document
---


## R Markdown

The HotspotsDerived repository is designed to process all the data produced by the Hotspots project (specifically the BioDivMapping
repository) and produce a series of derived products to be used in the reports and data delivered to Miljodirektoratet.

Unlike the BioDivMapping repository, this is not designed as a singular pipeline. Each folder contains a script or series of script 
which need to be used manually to produce the desired outputs. They all correspond to different products. There are many functions
from the BioDivMapping repository which have been imported here, as well as some data which is required for these 
scripts.

One thing you will need to make effective use of this script is to have a link of some sort to the hotspots repository and its data
outputs. I've written this repo from the perspective of having the two repos in one accessible environment, with several scripts having
the pathway from one to the other defined early on.

The folders include:

### Processes

This folder contains sub-folders, each of which is a series of scripts and external data which corresponds to a derived product.

#### Beta diversity

Here we generate beta diversity values for all species. generateBetaDiversity.R creates the raster layers for each species group,
and createBetaDiversityPackage.R brings them together in one tiff file. It's worth noting that this is a very computationally heavy 
function, and currently we only do this for ansvarsarter.

#### Field validation (localPredictions)

This folder includes four scripts. One which processes the field results initially (fieldWorkProcessing.R), one which uses
our models (imported via importModels.R) to produce predictions at different resolutions in the regions where the fieldwork took place 
(predictionGenerations.R), and another which comparies the observed species data from the fieldwork and the predicted species data 
(analysePredictions.R).

#### Overlays

This directory involves one script which imports our species richness estimates and then turns them into hotspot layers, to be exported
and then combined into a raster with other environmental layers for visual and statistical comparison. The data folder here includes 
different subfolders, some with external data that is added as layers of our final raster. Some of the raw data has been added, some is 
simply too big to fit on GitHub.

#### Covariate analysis

Here we download all our reduced models for a species group in order to find the overall effect of each environmental covariate
on the entire species group. The script simply involves compiling the effects across models and applying a density function to each
environmental effect.

#### Species numbers

This is a quick calculation which gives us statistics on the number of species used in modelling and the number of observations for these
species.

### Functions

This folder contains a series of functions (almsot all carried over from the BioDivMapping repo) that are necessary for some of the
process scripts.

### Data

Here we store all external data necessary for creating the derived products, from vernacular names, Norwegian translations and
alien species lists to Norwegian borders shapefiles and lake/city maps. We also stor data that needs to be uploaded in this folder,
and we import data produced by the BioDivMapping repository here.

### communicationsPackageHotspots

All code and data necessary for the Shiny app is stored here. All data for the app is initially imported using the importShinyData.R 
app.

### Visualisations

These R scripts produce visualisations for the Hotspot reports. A full list is given below.

- analyseOverlays.R
  - All bar plots comparing species richnesses in areas ofinterest (protected areas, hovedøkosystemer, etc.)
- betaDiversity.R
  - Beta diversity maps
- covariateEffectPlots.R
  - Density plots showing effects of both continuous and covariate variables on different species groups
- hotspotsFigures.R
  - Hotspot maps for each species groups
  - Probability of hotspot maps for each species group
- richnessMaps.R
  - Richness maps
  - Uncertainty maps
  - Sampling density maps
  - Sampling density maps inside hotspots


