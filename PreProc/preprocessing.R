list.of.packages <- c('cccd', 'grid', 'tsne', 'Rtsne', 'igraph', 'mclust', 'ggplot2', 'pheatmap', 'reshape', 'reshape2')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages, repos = c("http://cran.rstudio.com/", "https://bioconductor.org/biocLite.R"))

# In case of a missing package after running the steps above, please install them by typing:
#source('http://bioconductor.org/biocLite.R')
#biocLite('package_name')

# The package Vennerable also needs to be installed:
install.packages("Vennerable", repos = "http://R-Forge.R-project.org")

# Install FUSCA
library(devtools)
devtools::install_github("edroaldo/fusca")

list.of.packages<-c("cccd", "grid", "Rcpp", "tsne", "uwot", "dplyr", "fusca", "Rtsne", "igraph", "Matrix", "mclust", "tibble", 
                    "cowplot", "ggplot2", "reshape", "reshape2", "pheatmap")
lapply(list.of.packages, require, character.only = TRUE)

set.seed(42)
setwd

data = readMM("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/filtered_gene_bc_matrices/hg19/matrix.mtx")
genes <- make.unique(read.delim("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/filtered_gene_bc_matrices/hg19/genes.tsv", header = FALSE, stringsAsFactors = FALSE)$V2, sep = ".");
rownames(data) = genes;
colnames(data) <- read.delim("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/filtered_gene_bc_matrices/hg19/barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
data = as(data, "dgCMatrix")
data[1:5, 1:5]

cellrouter.pbmc <- CreateCellRouter(data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)


mito.genes <- grep(pattern = "^MT-", x = rownames(x = cellrouter.pbmc@assays$RNA@ndata), value = TRUE)
percent.mito <- Matrix::colSums(cellrouter.pbmc@assays$RNA@ndata[mito.genes, ]) / Matrix::colSums(cellrouter.pbmc@assays$RNA@ndata)
ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(x = cellrouter.pbmc@assays$RNA@ndata), value = TRUE)
percent.ribo <- Matrix::colSums(cellrouter.pbmc@assays$RNA@ndata[ribo.genes, ]) / Matrix::colSums(cellrouter.pbmc@assays$RNA@ndata)

cellrouter.pbmc@assays$RNA@sampTab$percent.mito = percent.mito
cellrouter.pbmc@assays$RNA@sampTab$percent.ribo = percent.ribo

p1 = ggplot(cellrouter.pbmc@assays$RNA@sampTab, aes(x = "percent.mito", y = percent.mito)) + 
  geom_violin(fill = "grey80", colour = "#FF0000") + theme(legend.position = "none") + xlab("") 
p2 = ggplot(cellrouter.pbmc@assays$RNA@sampTab, aes(x = "percent.ribo", y = percent.ribo)) + 
  geom_violin(fill = "grey80", colour = "#FF9900") + theme(legend.position = "none") + xlab("") 
p3 = ggplot(cellrouter.pbmc@assays$RNA@sampTab, aes(x = "nGene", y = nGene)) + 
  geom_violin(fill = "grey80", colour = "#3366FF") + theme(legend.position = "none") + xlab("") 
plot_grid(p1, p2, p3, nrow = 1) # labels = "")

cellrouter.pbmc <- filterCells(cellrouter.pbmc, assay.type = "RNA", 
                               variables = c("nGene", "percent.mito"), 
                               thresholds.low = c(200, -Inf, 0.10), thresholds.high = c(2500, 0.05, Inf))

#normalizing data
cellrouter.pbmc <- Normalize(cellrouter.pbmc)

#identify high variable genes and rank them using top 2000 genes
#scaling the data
#dimensional reduction
#cluster identification

var.genes <- FindVariableGenes(cellrouter.pbmc, assay.type = "RNA", method = "vst", loess.span = 0.3, pvalue = 0.05)
cellrouter.pbmc@var.genes <- rownames(var.genes[1:2000, ])

cellrouter.pbmc <- scaleData(cellrouter.pbmc, genes.use = cellrouter.pbmc@var.genes)

cellrouter.pbmc <- computePCA(cellrouter.pbmc, assay.type = "RNA", seed = 42, num.pcs = 50, genes.use = cellrouter.pbmc@var.genes) 
umap.done <- uwot::umap(cellrouter.pbmc@pca$cell.embeddings[, 1:15], spread = 1, min_dist = 0.3, n_neighbors = 30, metric = "cosine")
rownames(umap.done) <- rownames(cellrouter.pbmc@pca$cell.embeddings)
colnames(umap.done) <- c("UMAP1", "UMAP2")
cellrouter.pbmc <- customSpace(object = cellrouter.pbmc, matrix = umap.done)

#cluster identification
cellrouter.pbmc <- findClusters(cellrouter.pbmc, k = 15, num.pcs = 15, nn.type = "snn")
plotReducedDimension(cellrouter.pbmc, reduction.type = "custom", annotation = "population", annotation.color = "colors",
                     showlabels = T, dotsize = 0.01, labelsize = 5, convex = FALSE)

#to fix dopar error
# install.packages("doSNOW")
# 
# install.packages("doParallel") 
# 
# install.packages("doMPI")

library(doParallel)

library(dplyr)

#identification of cluster-specific gene signatures and marker genes
markers <- findSignatures(cellrouter.pbmc, assay.type = "RNA", column = "population", 
                          test.use = "wilcox", min.pct = 0.25, fc.threshold = 0.25, pos.only = T) 


markers.top5 <- as.data.frame(markers %>% 
                                group_by(population) %>% 
                                top_n(5, fc))


plot <- plotSignaturesHeatmap(cellrouter.pbmc, assay.type = "RNA", markers.top5, genes.show = as.vector(markers.top5$gene), threshold =  3, 
                              column.ann = "population", column.color = "colors")
grid::grid.draw(plot$gtable)

# Plot gene expression
genelist = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
plotDRExpression(cellrouter.pbmc, assay.type = "RNA", genelist = genelist, reduction.type = 'custom', 
                 dims.use = c(1, 2), threshold = 3, columns = 3, title = "")

#Gene markers
markers.top10 <- as.data.frame(markers %>% 
                                 group_by(population) %>% 
                                 top_n(10, fc))

tmp <- recode(as.character(cellrouter.pbmc@assays$RNA@sampTab$population), 
              "1" = "CD8", 
              "2" = "B", 
              "3" = "Memory CD4", 
              "4" = "CD14 Monocytes", 
              "5" = "NK", 
              "6" = "Naive CD4",
              "7" = "FCGR3A Monocytes", 
              "8" = "DC",
              "9" = "Megakaryocytes")
names(tmp) <- rownames(cellrouter.pbmc@assays$RNA@sampTab)

df <- data.frame(cluster = cellrouter.pbmc@assays$RNA@sampTab$population, celltype = tmp)
rownames(df) <- rownames(cellrouter.pbmc@assays$RNA@sampTab)
cellrouter.pbmc <- addInfo(cellrouter.pbmc, assay.type = "RNA", metadata = df, colname = "celltype", metadata.column = "celltype")
plotReducedDimension(cellrouter.pbmc, reduction.type = "custom", annotation = "celltype", annotation.color = "colors", 
                     showlabels = T, dotsize = 0.01, labelsize = 5, convex = FALSE)
