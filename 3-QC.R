# Run this script after 2-DSB.R

# Load Packages and Parallelize ---------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(Matrix)
library(BiocManager) # For a more efficient WRST (optional):
library(future)

plan()

# Load Data  ---------------------------
setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
# setwd("C:/Users/UPDATE/Desktop/COVID Lung CITE-Seq")

# Read filtered RNA data
path_data <- "Y:/Cole Keenum/CITE-seq files" # will be slow

mp4 <- Read10X(paste0(path_data,"/mp4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]
mp24 <- Read10X(paste0(path_data,"/mp24hr10xdata/mp24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]
p4 <- Read10X(paste0(path_data,"/p4hr10xdata/p4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]
p24 <- Read10X(paste0(path_data,"/p24hr10xdata/p24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]
naive <- Read10X(paste0(path_data,"/naive10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "RNA")
mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "RNA")
p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "RNA")
p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "RNA")
naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "RNA")

seuratObj <- merge(naive, c(p4, p24, mp4, mp24),
                   add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))
dim(seuratObj) # [1] 32285 29389
rm(mp24, mp4, naive, p4, p24)

saveRDS(seuratObj, 'seuratObj_rna.rds')

# Bring in Normalized DSB values with names for QC-filtered cells
dsb_data <- as.sparse(readRDS('dsb_singlets.rds'))

dim(dsb_data) # [1]   120 26935
length(intersect(colnames(dsb_data), colnames(seuratObj))) # [1] 26935

seuratObj <- subset(seuratObj, cells = colnames(dsb_data))
dim(seuratObj) # [1] 32285 26935
seuratObj[['ADT']] <- CreateAssayObject(data = dsb_data)

# Rearranging metadata for downstream plotting:
seuratObj@meta.data$orig.ident <-
  factor(x = seuratObj@meta.data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))

# Bring in scrublet_score metadata and apply to seurat object:
scrub_md <- read.csv('scrublet_metadata.csv') %>% filter(X %in% colnames(seuratObj))
all(scrub_md$X==colnames(seuratObj)) # [1] TRUE
any(scrub_md$predicted_doublet=='True') # [1] FALSE

seuratObj <- AddMetaData(object = seuratObj, 
                         metadata = scrub_md$doublet_score,
                         col.name = 'doublet_score')

seurat_list <- SplitObject(seuratObj, split.by = 'orig.ident')
rm(dsb_data, scrub_md, seuratObj)

# SCTransform ---------------------------
for (i in 1:length(x = seurat_list)){
  sample_id <- names(seurat_list)[[i]]
  print(sample_id)
  
  DefaultAssay(seurat_list[[i]]) = "RNA"
  # Calculating QC features
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^mt-")
  seurat_list[[i]][["percent.globin"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^Hb[abq]")
  seurat_list[[i]][["percent.Rpsl"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^Rp[sl]")
  
  # Plotting and saving QC figures
  qc_plot <-  VlnPlot(seurat_list[[i]], features = c("nFeature_RNA",
                                                     "percent.globin",
                                                     "percent.mt",
                                                     "percent.Rpsl"), ncol = 4)
  ggsave(paste(sample_id, "_QC_v3", ".png", sep = ""), plot = qc_plot)
  # Subsetting based on QC metrics:
  seurat_list[[i]] <- subset(seurat_list[[i]], subset = nFeature_RNA > 200 &
                               percent.globin < 10)
  seurat_list[[i]] <- SCTransform(seurat_list[[i]],
                                  method = "glmGamPoi",
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = 4000,
                                  verbose = TRUE)
}

# Cell Cycle Scoring ---------------------------
# From: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# BiocManager::install("biomaRt")
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
# Gene List from: https://science.sciencemag.org/content/352/6282/189
# Tirosh et al. Science 2019
s_genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
g2m_genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

# Save for later:
write.csv(s_genes, 's_genes.csv'); write.csv(g2m_genes, 'g2m_genes.csv')

for (i in 1:length(x = seurat_list)){ # I can move this into the above for loop...
  seurat_list[[i]] <- CellCycleScoring(seurat_list[[i]], g2m.features=g2m_genes, s.features=s_genes)
}

# RNA Integration ---------------------------
# Finding integration anchors for SCT
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")
# Integrate 6000 genes, but only use the anchors from the 3000 as anchorset
to_integrate <- SelectIntegrationFeatures(object.list = seurat_list, 
                                          nfeatures = 6000 , verbose = TRUE)

# Workspace saved here and moved to workstation computer
# save.image("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/pre_integration_workspace_07022021.RData")

combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                          features.to.integrate = to_integrate, verbose = T,
                          new.assay.name = "integratedSCT_")

# Regress out S and G2M scores ---------------------------
# This will take some time. 
combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"))

saveRDS(combined, 'combined_integrated_SCT_07022021.rds')

# ADT Integration ---------------------------
seurat_list <- anchors@object.list

# Identify the rows with isotype controls: exclude these from VariableFeatures
rownames(seurat_list[[1]][["ADT"]])
rownames(seurat_list[[1]][["ADT"]])[15:17] # [1] "IgG1-TotalA"  "IgG2-TotalA"  "IgG2b-TotalA"

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- c(rownames(seurat_list[[i]][["ADT"]])[1:14],
                                          rownames(seurat_list[[i]][["ADT"]])[18:120])
  # Data already normalized with DSB
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
combinedADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I am not sure if this is a good idea:
combined[["integratedADT_"]] <- combinedADT[["integratedADT_"]] 

rm(combinedADT, anchorsADT, anchors, seurat_list)
saveRDS(combined, "combined_integrated_SCT_DSB_07022021.rds")

# Separate RNA and ADT Clustering ---------------------------
# PCA on SCT
DefaultAssay(combined) <- "integratedSCT_"
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)
combined <- RunPCA(combined, verbose = T)
ElbowPlot(combined, ndims = 50, reduction = 'pca') # will use 40 because SCT is more accurate 
ggsave("elbow_plot_sct.png")

# aPCA on ADT
DefaultAssay(combined) <- 'integratedADT_'
combined <- ScaleData(combined)
combined <- RunPCA(combined, reduction.name = 'apca')
ElbowPlot(combined, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("elbow_plot_dsb.png")

# Assess clustering qualities on SCT and ADT umap separately:
combined <- RunUMAP(combined, reduction = 'pca', dims = 1:40, assay = 'SCT', 
                    reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
combined <- RunUMAP(combined, reduction = 'apca', dims = 1:11, assay = 'ADT', 
                    reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

combined <- FindNeighbors(combined, reduction = 'pca', dims = 1:40, 
                          graph.name = 'sct.snn')
combined <- FindClusters(combined, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

combined <- FindNeighbors(combined, reduction = 'apca', dims = 1:15,
                          graph.name = 'adt.snn')
combined <- FindClusters(combined, graph.name = 'adt.snn', resolution = 1.0, verbose = T)

Idents(combined) <- 'sct.snn_res.1'
p1 <- DimPlot(combined, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
Idents(combined) <- 'adt.snn_res.1'
p2 <- DimPlot(combined, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering.png", width = 10, height = 5)

# Plot metrics to see if celltypes cluster together:
DefaultAssay(combined) <- 'ADT'
FeaturePlot(combined, features = c('CD4-TotalA', 'CD8a-TotalA', 'CD8b-TotalA', 'TCRB-TotalA'), 
            reduction = 'adt.umap', cols = c("lightgrey","darkgreen"),
            min.cutoff = 0, max.cutoff = 'q99')
ggsave('tcell_markers_adt.png', width = 10, height = 10)

FeaturePlot(combined, features = c('CD4-TotalA', 'CD8a-TotalA', 'CD8b-TotalA', 'TCRB-TotalA'), 
            reduction = 'sct.umap', cols = c("lightgrey","darkgreen"),
            min.cutoff = 0, max.cutoff = 'q99')
ggsave('tcell_markers_sct_proj120.png', width = 10, height = 10)
DefaultAssay(combined) <- 'integratedSCT_'

# Plot doublet_score: NOTE-Each sample has a different threshold, so use caution
FeaturePlot(combined, features = c('doublet_score'), reduction = 'sct.umap')
ggsave('doublet_score_after_scrublet_sct_umap.png', width = 5, height = 5)

FeaturePlot(combined, features = c('doublet_score'), reduction = 'adt.umap')
ggsave('doublet_score_after_scrublet_adt_umap.png', width = 5, height = 5)

DimPlot(combined, reduction = 'sct.umap', label = T) + NoLegend()

saveRDS(combined, 'combined_PCA.rds')
combined@assays$RNA <- NULL
saveRDS(combined, 'combined_PCA_noRNA.rds')

# WNN clustering ---------------------------
combined <- FindMultiModalNeighbors(
  combined, reduction.list = list("pca", "apca"), 
  dims.list = list(1:40, 1:15), modality.weight.name = "SCT.weight")

combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Find clusters using a range of resolutions
resolution.range <- seq(from = 0, to = 3.0, by = 0.1)
combined <- FindClusters(object = combined, graph.name = "wsnn",
                         reduction.name = "wnn.umap", algorithm = 3, 
                         resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("wsnn_res.", res, sep = "")
  Idents(combined) <- md
  p <- DimPlot(combined, label = T, repel = T, label.size = 3, reduction = 'wnn.umap') + NoLegend()
  ggsave(paste("dimplot_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(combined, prefix = "wsnn_res.")
ggsave("clustree_output_wnn.pdf", width = 9.5, height = 11*3)

saveRDS(combined, 'combined_wnn.rds')
# combined <- readRDS('combined_wnn.rds')

res <- 1.5
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)

# Data visualization:
DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
ggsave(paste("dimplot_", "res_", res, "_used.png", sep = ""), width = 5, height = 5)

FeaturePlot(combined, features = c('S.Score', 'G2M.Score'), reduction = 'wnn.umap', 
            min.cutoff = 0)
ggsave('cell_cycle_scores.wnn.png', width = 10, height = 5)

# Run: 4.1-cell_markers.R

# Marker Genes for Unlabeled Clsuters ---------------------------
## Exporting top 10 marker genes:
DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(combined.markers, file = "gene_biomarkers_unlab.csv", row.names = FALSE)
write.csv(top10, file = "gene_biomarkers_unlab_top10.csv", row.names = FALSE)

DefaultAssay(combined) <- "ADT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(combined.markers, file = "adt_biomarkers_unlab.csv", row.names = FALSE)
write.csv(top10, file = "adt_biomarkers_unlab_top10.csv", row.names = FALSE)

# Doublets: Cluster 39
dim(combined)
combined <- subset(combined, idents = c('39'), invert = TRUE)
dim(combined)

saveRDS(combined, 'combined_res1.5.rds')

# Tables:
unlab_freq <- table(combined@active.ident, col.names = combined$orig.ident)
col_order <- c('naive', 'p4', 'mp4', 'p24', 'mp24')
unlab_freq <- unlab_freq[ ,col_order]
write.csv(unlab_freq, file = "unlab_freq.csv", row.names = T)

# Subsetting of T/NK cells ---------------------------
tcells <- subset(combined, idents = c(7,8,13,15,17,27,30))
DimPlot(tcells, reduction = 'wnn.umap', label = T) + NoLegend()
# Go deeper into clustering resolution to separate CD4 and CD8 cells based on ADT
Idents(tcells) <- "wsnn_res.3"
DimPlot(tcells, reduction = 'wnn.umap', label = T) + NoLegend()
ggsave('tcells_res3_wnnUMAP.png', width = 5, height = 5)

DefaultAssay(tcells) <- 'SCT'
plot <- FeaturePlot(tcells, features = "Tcrg-C1", reduction = 'wnn.umap')
HoverLocator(plot = plot, information = FetchData(tcells, vars = c("ident")))
DefaultAssay(tcells) <- 'integratedSCT_'

## Exporting top 10 marker genes:
DefaultAssay(tcells) <- "SCT"
tcells.markers <- FindAllMarkers(tcells, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
write.csv(tcells.markers, file = "tcells_gene_biomarkers_res3.csv", row.names = FALSE)

tcells$celltype <- Idents(tcells)
tcells <- RenameIdents(tcells, 
                       '7' = 'CD4 T',
                       '10' = 'NK',
                       '14' = 'CD4 T',
                       '17' = 'CD8 T',
                       '29' = 'Nuocyte',
                       '31' = 'CD4 NKT',
                       '37' = 'Treg',
                       '38' = 'CD8 NKT',
                       '39' = 'CD4 T',
                       '50' = 'NK',
                       '52' = 'gd T',
                       '55' = 'CD4 T') # proliferating cells, are CD4 and not CD8 and not NK/NKT
DimPlot(tcells, reduction = 'wnn.umap', label = T) + NoLegend()
ggsave('tcells_wnnUMAP_labeled.png', width = 5, height = 5)

saveRDS(tcells, 'tcells_sub.rds')
rm(tcells)

# Determine identity of cluster 31 ----
# Enriched for many TFs: Let's see how it maps onto adt umap alone
# Some doublets present: Will remove with CellSelector
plot <- DimPlot(subset(combined, ident = '31'), reduction = 'adt.umap', label = T) + NoLegend()
select.cells <- CellSelector(plot = plot)
DimPlot(subset(combined, ident = '31'), reduction = 'wnn.umap', label = T) + NoLegend()
DimPlot(subset(combined, cells = select.cells), reduction = 'wnn.umap', label = T) + NoLegend()

plot <- DimPlot(subset(combined, ident = '31'), reduction = 'wnn.umap', label = T) + NoLegend()
select.cells.2 <- CellSelector(plot = plot)

# Subset out Macrophage/Monocyte/DC Cluster
mamodc <- subset(combined, idents = c(10,18,23,24,31))
DimPlot(mamodc, reduction = 'wnn.umap', label = T) + NoLegend()

Idents(mamodc) <- "wsnn_res.3"
DimPlot(mamodc, reduction = 'wnn.umap', label = T) + NoLegend()
ggsave('mamodc_res3_wnnUMAP.png', width = 5, height = 5)

DefaultAssay(mamodc) <- 'SCT'
FeaturePlot(mamodc, features = c('Ms4a6d', 'Cybb'), reduction = 'wnn.umap')

DefaultAssay(mamodc) <- 'ADT'
FeaturePlot(subset(mamodc, cells = c(select.cells, select.cells.2), invert = T), features = c("CD19-TotalA", "CD2-TotalA", "Ly6g.Ly6c-TotalA"), reduction = 'wnn.umap')
plot <- FeaturePlot(mamodc, features = "Ly-6G-TotalA", reduction = 'wnn.umap')
select.cells.3 <- CellSelector(plot = plot)
select.cells.4 <- CellSelector(plot = plot)
DefaultAssay(mamodc) <- 'integratedSCT_'

# Remove more doublets:
dim(combined) # [1]   120 26795
combined <- subset(combined, cells = c(select.cells, 
                                       select.cells.2,
                                       select.cells.3, 
                                       select.cells.4), invert = T)
dim(combined) # [1]   120 26752
DimPlot(combined, reduction = 'wnn.umap', label = T) + NoLegend()

# Compute marker genes without doublets now:
## Exporting top 10 marker genes:
DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
write.csv(combined.markers, file = "combined_gene_biomarkers_no_dub.csv", row.names = FALSE)

# Compute marker genes without doublets now:
## Exporting top 10 marker genes:
DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
write.csv(combined.markers, file = "combined_gene_biomarkers_no_dub.csv", row.names = FALSE)

# Remove droplet clusters with migh mt% ---------------------------
# Remove clusters with very high mt gene markers, sometimes these have very high 
# transcription factor expression as well. 
FeaturePlot(combined, features = 'percent.mt', reduction = 'wnn.umap')
ggsave('percent.mt_wnn.png', width = 5, height = 5)

# Overcluster to isolate might mt% regions:
Idents(combined) <- 'wsnn_res.3'
VlnPlot(combined, features = 'percent.mt', group.by = 'orig.ident')
VlnPlot(combined, features = 'percent.mt') + NoLegend()

DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
write.csv(combined.markers, file = "gene_biomarkers_unlab_res3.csv", row.names = FALSE)

dim(combined) # [1] 18037 26752
combined <- subset(combined, idents = c(13, 39, 41), invert = TRUE)
# combined <- subset(combined, idents = 41, invert = TRUE)
dim(combined) # [1] 18037 25746

Idents(combined) <- 'wsnn_res.1.5'
DimPlot(combined, reduction = 'wnn.umap', label = T) + NoLegend()
ggsave('high_mt_removed_res1.5.png', width = 5, height = 5)

saveRDS(combined, 'combined_qc.rds')
