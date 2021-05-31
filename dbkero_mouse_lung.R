# Load packages -------------
library(Seurat)
library(dplyr)

getwd()
setwd('C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/DBKERO scRNA scATAC')
getwd()

# scRNA Seq -------------
rna <- Read10X(data.dir = './scRNA')

dim(rna)
head(rna)

rna <- CreateSeuratObject(rna)

dim(rna)
head(rna)

# QC:
VlnPlot(rna, features = c('nCount_RNA', 'nFeature_RNA'))
rna[['percent.mt']] <- PercentageFeatureSet(rna, pattern = '^mt-')
VlnPlot(rna, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

rna <- subset(rna, subset = nFeature_RNA > 200 & percent.mt < 25)
dim(rna)

# Normalization
rna <- SCTransform(rna, vars.to.regress = 'percent.mt')

# Clustering
DefaultAssay(rna) <- 'SCT'
rna <- RunPCA(rna)
DimPlot(rna, reduction = 'pca')
ElbowPlot(rna, ndims = 50)

rna <- RunUMAP(rna, dims = 1:40)

DimPlot(rna)

rna <- FindNeighbors(rna, dims = 1:40, verbose = TRUE)
rna <- FindClusters(rna)

DimPlot(rna)

# Find Genes for Each Cluster
markers <- FindAllMarkers(rna, assay = 'SCT')
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(top10, file = "japan.markers.top10.csv", row.names = FALSE)

