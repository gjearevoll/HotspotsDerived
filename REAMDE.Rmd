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

### Field validation

This folder includes three scripts. One which processes the filed results initially (fieldWorkProcessing.R), one which uses
our models to produce predictions at different resolutions in the regions where the fieldwork took place (localPredictions.R),
and another which comparies the observed species data from the fieldwork and the predicted species data.

### Overlays

This directory involves one script which imports our species richness estimates and then turns them into layers of a raster. It
then combines them with a large range of other raster data (see the individual README file for more details) so that these
can easily be visualised in conjunction with our Hotspots data. The data folder here includes different subfolders, each with external
data that is added as layers of our final raster.

### Model analysis

Here we download all our reduced models for a species group in order to find the overall effect of each environmental covariate
on the entire species group. The script simply involves compiling the effects across models and applying a density function to each
environmental effect.

### Beta diversity

This is one script which produces an estimate for beta diversity across all ansvarsarter. The results are then fed into our overlays 
script above as a layer of this output.

### Model resizing

The model outputs which we get from the INLA-reliant models are large, often around half a gigabyte. When one of these model outputs
is required for every 10 species  segment across 1500 species, they start making any model-wide operations excessively slow. As such,
here we apply a function that imports each model, cuts out a huge amount of unnecessary data, and 
