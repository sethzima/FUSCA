#install.packages("devtools")
library("devtools")
#devtools::install_github("edroaldo/fusca")
#required to fix error in multiprocessing
library(doParallel)

list.of.packages<-c("cccd", "scales", "grid", "Rcpp", "tsne", "uwot", "dplyr", "fusca", "Rtsne", "igraph", "Matrix", "mclust", "tibble", 
                    "cowplot", "ggplot2", "reshape", "reshape2", "pheatmap")
lapply(list.of.packages, require, character.only = TRUE)

#attempting to run the tutorial without the spatial transcriptomics data
set.seed(1)
data_dir = "C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526"

#create the cellrouter object
data = readMM("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526/GSE160526_matrix.mtx.gz");
genes <- make.unique(read.delim("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526/GSE160526_genes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)$V2, sep = ".");
rownames(data) = toupper(genes);
colnames(data) <- read.delim("C:/Users/Patron/Documents/MyDocuments/R/FUSCA/CellCommReproduced/GSE160526/GSE160526_barcodes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)$V1
data = as(data, "dgCMatrix")
cellrouter.eht <- CreateCellRouter(data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)

#find variable genes, normailize and scale data
var.genes <- FindVariableGenes(cellrouter.eht, assay.type = "RNA", method = "coefficient_variation", pvalue = 0.05)
var.genes <- var.genes[order(var.genes$adj.pvalue, decreasing = FALSE), ]
cellrouter.eht@var.genes <- rownames(var.genes[1:3000, ])

cellrouter.eht <- Normalize(cellrouter.eht, assay.type = "RNA")

cellrouter.eht <- scaleData(cellrouter.eht, assay.type = "RNA", 
                        genes.use = cellrouter.eht@var.genes, blocksize = nrow(cellrouter.eht@assays$RNA@ndata))


#Dimensionality reduction and clustering
cellrouter.eht <- computePCA(cellrouter.eht, assay.type = "RNA", 
                         seed = 42, num.pcs = 50, genes.use = cellrouter.eht@var.genes) 
plot(cellrouter.eht@pca$sdev)

cellrouter.eht <- computeTSNE(cellrouter.eht, 
                          seed = 1, num.pcs = 8, max_iter = 1000)

cellrouter.eht <- computeUMAP(cellrouter.eht, 
                          seed = 1, num.pcs = 8, metric = "euclidean", n_neighbors = 15, spread = 1, min_dist = 0.1)

cellrouter.eht <- findClusters(cellrouter.eht, assay.type = "RNA", 
                           num.pcs = 8, nn.type = "knn", k = 20)

plotReducedDimension(cellrouter.eht, assay.type = "RNA", reduction.type = "tsne", annotation = "population", annotation.color = "colors",
                     dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)

plotReducedDimension(cellrouter.eht, assay.type = "RNA", reduction.type = "umap", annotation = "population", annotation.color = 'colors',
                     dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)

#identification of marker genes
# markers <- findSignatures(cellrouter.eht, assay.type = "RNA", 
#                           column = "celltype", pos.only = TRUE, fc.threshold = 0.2, nCores = 10); gc();

markers <- findSignatures(cellrouter.eht, assay.type = "RNA", column = "population", 
                          test.use = "wilcox", min.pct = 0.25, fc.threshold = 0.25, pos.only = T) 

markers.top <- markers %>% group_by(population) %>% top_n(10, fc)
markers.top <- as.data.frame(markers.top)


#error when run
# Error in split.default(x = seq_len(nrow(x)), f = f, drop = drop, ...) : 
#   group length is 0 but data length > 0
plot <- plotSignaturesHeatmap(cellrouter.eht, assay.type = 'RNA', 
                              markers.top, genes.show = as.vector(markers.top$gene), 
                              column.ann = 'celltype', column.color  = 'celltype_color', 
                              num.cells = 700, threshold =  2)

#this block works
plot <- plotSignaturesHeatmap(cellrouter.eht, assay.type = "RNA", markers.top, genes.show = as.vector(markers.top$gene), threshold =  3, 
                              column.ann = "population", column.color = "colors")


grid::grid.draw(plot$gtable)


#mean expression of ligands and receptors per cluster
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))# from NicheNet
head(lr_network)

pairs <- lr_network
pairs$Pair.Name <- paste(pairs$from, pairs$to, sep = "_")

ligands <- unique(lr_network$from)
ligands <- intersect(ligands, rownames(cellrouter.eht@assays$RNA@ndata))

receptors <- unique(lr_network$to)
receptors <- intersect(receptors, rownames(cellrouter.eht@assays$RNA@ndata))

ligands.receptors <- unique(c(ligands, receptors))



#assay type error, changed from celltype to population
mean.expr <- computeValue(cellrouter.eht, assay.type = "RNA", 
                          genelist = ligands.receptors, column = "population", fun = "mean"); gc();

interactions <- population.pairing(mean.expr = mean.expr, pairs = pairs, ligands = ligands, receptors = receptors, threshold = 0.25)

interactions <- calculateObservedMean(mean.expr = mean.expr, interactions = interactions)
head(interactions)

markers <- findSignatures(cellrouter.eht, assay.type = "RNA", 
                          column = "population", pos.only = TRUE, fc.threshold = 0.2, nCores = 10); gc();


#calculate null distribution of intracluster gene expression means
genelist <- unique(c(interactions$ligand, interactions$receptor))

p <- clusterPermutation(cellrouter.eht, assay.type = "RNA", 
                        genelist = genelist, interactions = interactions, cluster.label = "population", nPerm = 1000, nCores = 10)

interactions.p <- calculatePvalue(p, nPerm = 1000, interactions2 = interactions)

#saveRDS(interactions.p, file = paste(bd, "/interactions_with_Pvalue_1000_cluster.rds", sep = ""))
#interactions.p <- readRDS(paste(basedir, "/interactions_with_Pvalue_1000_cluster.rds", sep = ""))
tmp <- interactions.p[which(interactions.p$pvalue < 0.01),]
head(tmp)



#network calculation
my_matrix <- interactionmatrix(tmp)
head(matrix)

my_graph <- cellnetwork3(tmp, threshold = 5)

cellrouter.eht <- calculateCentroids(cellrouter.eht, assay.type = "RNA", sample.name = "Sample1", 
                                 cluster.column = "population", cluster.type = "Cluster")

cellrouter.eht <- calculateDistanceMatrix(cellrouter.eht, assay.type = "RNA", sample.name = "Sample1", 
                                      cluster.type = "Cluster", spot_distance = 100, normalize = FALSE)
function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL) 
   {                                                                        
         if (is.object(data) || !is.atomic(data))                             
               data <- as.vector(data)                                          
           .Internal(matrix(data, nrow, ncol, byrow, dimnames, missing(nrow),   
                                      missing(ncol)))
}
      
plots <- predictCellInteractions(cellrouter.eht, assay.type = "RNA", sample.name = "Sample1", 
                                       cluster.type = "Cluster", graph = my_graph, distance.threshold = 0.75)
options(repr.plot.width = 22, repr.plot.height = 7)
gridExtra::grid.arrange(grobs = plots, ncol = 1)



