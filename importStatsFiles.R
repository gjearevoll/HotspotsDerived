# Quick function carrying over all finished data used to create raster packages 3 or 4

sourceDirectory <- "../BioDivMapping/data/run_2024-10-11/modelOutputs/processedOutputs"
filesToImport <- c(list.files(sourceDirectory, pattern = "final", recursive = TRUE, full.names = TRUE),
                   list.files(sourceDirectory, pattern = "bias", recursive = TRUE, full.names = TRUE))

file.copy(filesToImport, "data/modelOutputs")
