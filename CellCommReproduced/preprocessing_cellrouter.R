list.of.packages <- c('cccd', 'grid', 'tsne', 'Rtsne', 'igraph', 'mclust', 'ggplot2', 'pheatmap', 'reshape', 'reshape2')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages, repos = c("http://cran.rstudio.com/", "https://bioconductor.org/biocLite.R"))

# In case of a missing package after running the steps above, please install them by typing:
#source('http://bioconductor.org/biocLite.R')
#biocLite('package_name')

# The package Vennerable also needs to be installed:
install.packages("Vennerable", repos = "http://R-Forge.R-project.org")

#required to fix error in multiprocessing
library(doParallel)

# Install FUSCA
library(devtools)
devtools::install_github("edroaldo/fusca")

list.of.packages<-c("cccd", "grid", "Rcpp", "tsne", "uwot", "dplyr", "fusca", "Rtsne", "igraph", "Matrix", "mclust", "tibble", 
                    "cowplot", "ggplot2", "reshape", "reshape2", "pheatmap")
lapply(list.of.packages, require, character.only = TRUE)

set.seed(42)
setwd

data_dir = "C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526"

list.files(data_dir)

data = readMM("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526/GSE160526_matrix.mtx")
genes <- make.unique(read.delim("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526/GSE160526_genes.tsv",
                                header = FALSE, stringsAsFactors = FALSE)$V2, sep =".")
rownames(data) = genes
colnames(data) <- read.delim("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526/GSE160526_barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
data = as(data, "dgCMatrix")
data[1:5, 1:5]


cellrouter.het <- CreateCellRouter(data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)


mito.genes <- grep(pattern = "^MT-", x = rownames(x = cellrouter.het@assays$RNA@ndata), value = TRUE)
percent.mito <- Matrix::colSums(cellrouter.het@assays$RNA@ndata[mito.genes, ]) / Matrix::colSums(cellrouter.het@assays$RNA@ndata)
ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(x = cellrouter.het@assays$RNA@ndata), value = TRUE)
percent.ribo <- Matrix::colSums(cellrouter.het@assays$RNA@ndata[ribo.genes, ]) / Matrix::colSums(cellrouter.het@assays$RNA@ndata)

cellrouter.het@assays$RNA@sampTab$percent.mito = percent.mito
cellrouter.het@assays$RNA@sampTab$percent.ribo = percent.ribo

p1 = ggplot(cellrouter.het@assays$RNA@sampTab, aes(x = "percent.mito", y = percent.mito)) + 
  geom_violin(fill = "grey80", colour = "#FF0000") + theme(legend.position = "none") + xlab("") 
p2 = ggplot(cellrouter.het@assays$RNA@sampTab, aes(x = "percent.ribo", y = percent.ribo)) + 
  geom_violin(fill = "grey80", colour = "#FF9900") + theme(legend.position = "none") + xlab("") 
p3 = ggplot(cellrouter.het@assays$RNA@sampTab, aes(x = "nGene", y = nGene)) + 
  geom_violin(fill = "grey80", colour = "#3366FF") + theme(legend.position = "none") + xlab("") 
plot_grid(p1, p2, p3, nrow = 1) # labels = "")

cellrouter.het <- filterCells(cellrouter.het, assay.type = "RNA", 
                               variables = c("nGene", "percent.mito"), 
                               thresholds.low = c(200, -Inf, 0.10), thresholds.high = c(2500, 0.05, Inf))


#normalize the data
cellrouter.het <- Normalize(cellrouter.het)

#identify highly variable genes and rank, take top 2000
var.genes <- FindVariableGenes(cellrouter.het, assay.type = "RNA", method = "vst", loess.span = 0.3, pvalue = 0.05)
cellrouter.het@var.genes <- rownames(var.genes[1:2000, ])

#scale the data
cellrouter.het <- scaleData(cellrouter.het, genes.use = cellrouter.het@var.genes)


#dimensional reduction (t-SNE)
cellrouter.het <- computePCA(cellrouter.het, assay.type = "RNA", seed = 42, num.pcs = 50, genes.use = cellrouter.het@var.genes)
plot(cellrouter.het@pca$sdev)

cellrouter.het <- computeTSNE(cellrouter.het, num.pcs = 15, seed = 42, max_iter = 1000)

cellrouter.het <- customSpace(cellrouter.het, cellrouter.het@tsne$cell.embeddings)

plotReducedDimension(cellrouter.het, reduction.type = 'tsne', dims.use = c(1,2), annotation = "celltype", annotation.color = 'celltype_color', showlabels = TRUE)


#dimensional reduction (UMAP)
cellrouter.het <- computePCA(cellrouter.het, assay.type = "RNA", seed = 42, num.pcs = 50, genes.use = cellrouter.het@var.genes) 
umap.done <- uwot::umap(cellrouter.het@pca$cell.embeddings[, 1:15], spread = 1, min_dist = 0.3, n_neighbors = 30, metric = "cosine")
rownames(umap.done) <- rownames(cellrouter.het@pca$cell.embeddings)
colnames(umap.done) <- c("UMAP1", "UMAP2")
cellrouter.het <- customSpace(object = cellrouter.het, matrix = umap.done)


#cluster identification
cellrouter.het <- findClusters(cellrouter.het, k = 15, num.pcs = 15, nn.type = "snn")
plotReducedDimension(cellrouter.het, reduction.type = "custom", annotation = "population", annotation.color = "colors",
                     showlabels = T, dotsize = 0.01, labelsize = 5, convex = FALSE)


#identification of cluster-specific genes
markers <- findSignatures(cellrouter.het, assay.type = "RNA", column = "population", 
                          test.use = "wilcox", min.pct = 0.25, fc.threshold = 0.25, pos.only = T) 
markers.top5 <- as.data.frame(markers %>% 
                                group_by(population) %>% 
                                top_n(5, fc))
plot <- plotSignaturesHeatmap(cellrouter.het, assay.type = "RNA", markers.top5, genes.show = as.vector(markers.top5$gene), threshold =  3, 
                              column.ann = "population", column.color = "colors")
grid::grid.draw(plot$gtable)


# Plot gene expression
genelist = c("Cldn6", "Tshz2", "Wnt16", "Mpped2", "Col2a1", "Hes5", "Tfap2b", "Fcer1g", "Ramp2",
             "Lhx3", "Onecut2", "Cdkn1a", "Rhox9")
plotDRExpression(cellrouter.het, assay.type = "RNA", genelist = genelist, reduction.type = 'custom', 
                 dims.use = c(1, 2), threshold = 3, columns = 4, title = "")
