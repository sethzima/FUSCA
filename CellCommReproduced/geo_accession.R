if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

#error solved by not updating packages when prompted

gse <- getGEO("GSE160526", GSEMatrix = TRUE)

scData <- getGEOSuppFiles("GSE160526")

