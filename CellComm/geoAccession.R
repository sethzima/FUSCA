if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

#error solved by not updating packages when prompted

gse <- getGEO("GSE144240", GSEMatrix = TRUE)
show(gse)

filePaths = getGEOSuppFiles("GSE144240")
filePaths

expressionSet <- gse[[1]]
pData(expressionSet)
exprs(expressionSet)

