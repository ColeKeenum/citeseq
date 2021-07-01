# NOTES ---------------------------
# CITE-Seq Analysis of Lung Cells w/ PUUC and MPLA+PUUC at 4 and 24 hour
# 1. Used SCTransform to get better sequencing results
# 2. Scrublet to remove putative doublets
# 3. Applied outlier analysis to ADT reads as well (removed ~250 cells)
# 4. Used DSB method of normalization for ADT data (V2 as of 04152021)
# 5. Regressed out cell cycle gene effect for SCT scaling

# NEED TO  ---------------------------
# 1. Consider regressing out dissociation-induced genes

# Load Packages  ---------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(Matrix)
## For a more efficient WRST (optional):
library(BiocManager)
## For plotting
library("dittoSeq")
# paralelizing
library(future)
plan("multicore")

# Load Data  ---------------------------

# new PC-based QC
setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
# setwd("C:/Users/UPDATE/Desktop/COVID Lung CITE-Seq")

# seuratObj <- readRDS("seuratObj_RNA_ADT.rds")
seuratObj <- readRDS("seuratObj_scrub_mtQC_adtQC.rds")

# Add dsb-Normalized ADTs ---------------------------
dsb <- readRDS('dsb_seurat_obj.rds')
dsb <- dsb@assays[["CITE"]]@data
dsb <- as(dsb, 'sparseMatrix')
overlap_cells <- intersect(colnames(dsb), colnames(seuratObj))
dsb <- dsb[ , overlap_cells]
seuratObj <- SetAssayData(seuratObj, slot = 'data', new.data = dsb, assay = 'ADT')

# Remove bad ADTs (does not remove high-responders) ---------------------------
adt_counts <- seuratObj[['ADT']]@counts
adt_data <- seuratObj[['ADT']]@data

all(rownames(adt_counts) == rownames(adt_data))
all(colnames(adt_counts) == colnames(adt_data))

library(readxl)
bad_adts <- read_excel("Y:/Cole Keenum/CITE-seq files/2021-04-17 Bad ADT List V4.xlsx")
bad_adts <- bad_adts$`DSB-Final`
for (i in 1:length(bad_adts)){bad_adts[[i]] <- paste(bad_adts[[i]], "-TotalA", sep = "")}

setdiff(bad_adts, rownames(adt_counts))

idx <- match(rownames(adt_counts), bad_adts)
idx <- which(is.na(idx))

adt_counts <- adt_counts[idx, ]
adt_data <- adt_data[idx, ]

all(rownames(adt_counts) == rownames(adt_data))
all(colnames(adt_counts) == colnames(adt_data))

# re-creating ADT assay
seuratObj[['ADT']] <- NULL
seuratObj[['ADT']] <- CreateAssayObject(counts = adt_counts)
seuratObj <- SetAssayData(seuratObj, slot = 'data', 
                          new.data = adt_data, assay = 'ADT')
dim(seuratObj[['ADT']])

rm(adt_counts, adt_data, bad_adts, i, idents, idx, test)

# Make a list of Seurat objects ---------------------------
seurat_list <- SplitObject(seuratObj, split.by = "orig.ident")
rm(seuratObj, dsb, overlap_cells)

dim(seurat_list[[1]][['ADT']])

# Run Quick ADT Integration and Visualization ---------------------------
# Determine if the removel of ADTs were successful, see if there are any very
# high responsing

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- rownames(seurat_list[[i]][["ADT"]]) #select all ADT as variable features
  # DSB ALREADY NORMALIZED
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
combinedADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I am not sure if this is a good idea:
# combined_clean[["integratedADT_"]] <- combinedADT[["integratedADT_"]] 

rm(anchorsADT, seurat_list)

# TEST Visualize ADT Clustering ---------------------------
# aPCA on ADT
DefaultAssay(combinedADT) <- 'integratedADT_'
combinedADT <- ScaleData(combinedADT)
combinedADT <- RunPCA(combinedADT, reduction.name = 'apca')
ElbowPlot(combinedADT, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("TEST120_elbow_plot_dsb_adt_scrub.png")

combinedADT <- RunUMAP(combinedADT, reduction = 'apca', dims = 1:15, assay = 'ADT', 
                          reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

combinedADT <- FindNeighbors(combinedADT, reduction = 'apca', dims = 1:15,
                                graph.name = 'adt.snn')
combinedADT <- FindClusters(combinedADT, graph.name = 'adt.snn', resolution = 1.0, verbose = T)

Idents(combinedADT) <- 'adt.snn_res.1'
DimPlot(combinedADT, reduction = 'adt.umap', label = TRUE, repel = TRUE, 
        label.size = 2.5) + NoLegend()
ggsave("TEST120_SCT_or_ADT_only_clustering_DSB_SCT_scrub.png", width = 5, height = 5)

DefaultAssay(combinedADT) <- "ADT"
FeaturePlot(combinedADT, features = c("CD4-TotalA", "CD8a-TotalA","CD8b-TotalA"),
            reduction = 'adt.umap', 
            cols = c("lightgrey","darkgreen"), ncol = 3, min.cutoff=0,
            max.cutoff = "q99")
ggsave("TEST120_T_cell_markers_DSB_120_0q99.png", width = 15, height = 5)

VlnPlot(combinedADT, features = 'nFeature_ADT') + NoLegend()
ggsave('vln_plot_some_high_responders.png', width = 10, height = 5)

DefaultAssay(combinedADT) <- "ADT"
FeaturePlot(combinedADT, features = c("ESAM-TotalA", "CD309-TotalA","CD326-TotalA"),
            reduction = 'adt.umap', 
            cols = c("lightgrey","darkgreen"), ncol = 3, min.cutoff=0,
            max.cutoff = "q99")
ggsave("TEST120_EpiEndo_cell_markers_DSB_120_0q99.png", width = 15, height = 5)

DefaultAssay(combinedADT) <- "ADT"
FeaturePlot(combinedADT, features = c("IgD-TotalA", "CD45R-TotalA","CD21/35-TotalA"),
            reduction = 'adt.umap', 
            cols = c("lightgrey","darkgreen"), ncol = 3, min.cutoff=0,
            max.cutoff = "q99")
ggsave("TEST120_B_cell_markers_DSB_120_0q99.png", width = 15, height = 5)

DefaultAssay(combinedADT) <- "ADT"
FeaturePlot(combinedADT, features = c("CD24-TotalA", "CD11b-TotalA", "Ly6g.Ly6c-TotalA"),
            reduction = 'adt.umap', 
            cols = c("lightgrey","darkgreen"), ncol = 3, min.cutoff=0,
            max.cutoff = "q99")
ggsave("TEST120_neutro_markers_DSB_120_0q99.png", width = 15, height = 5)

dim(combinedADT)
dim(subset(combinedADT, idents = c('16', '19'))) # putative high-responders

# TEST Visualize ADT Clustering ---------------------------
DefaultAssay(combinedADT) <- "ADT"
combinedADT.markers <- FindAllMarkers(combinedADT, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
top10 <- combinedADT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(top10, file = "combinedADT_cluster_biomarkers_v2.csv", row.names = FALSE)

## Plotting top 5 marker genes on heatmap: 

# sketchy?
combinedADT <- ScaleData(combinedADT, assay = 'ADT')

Idents(combinedADT) <- 'adt.snn_res.1'
top5 <- combinedADT.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#uninfected_lung <- subset(combinedADT, subset = trt == "Uninfected")
p <- DoHeatmap(subset(combinedADT, downsample = 100),
               features = top5$gene, size = 3, angle = 30) + NoLegend()
ggsave("combinedADT_top5_markers_heatmap_pre_high_remove.png", plot = p, width = 10, height = 6)

top2 <- combinedADT.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#uninfected_lung <- subset(combinedADT, subset = trt == "Uninfected")
p <- DoHeatmap(subset(combinedADT, downsample = 100),
               features = top2$gene, size = 3, angle = 30) + NoLegend()
ggsave("combinedADT_top2_markers_heatmap_pre_high_remove.png", plot = p, width = 10, height = 6)

# Remove High Responders and prep for final integration ---------------------------
idents <- as.character(0:(length(unique(Idents(combinedADT)))-1))
idents <- idents[-17] # removes ident 16
idents <- idents[-19] # removes ident 19
clean <- subset(combinedADT, idents = idents)

dim(combinedADT)
dim(clean)

# Undo integratedADT assay
clean[['integratedADT_']] <- NULL

dim(seurat_list[[1]])
seurat_list <- SplitObject(clean, split.by = "orig.ident")
dim(seurat_list[[1]])

# clean up memory
rm(clean, combinedADT, combinedADT.markers, p, top10, top2, top5, idents, anchorsADT)

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
  ggsave(paste(sample_id, "_QC", ".png", sep = ""), plot = qc_plot)
  # Subsetting based on QC metrics:
  seurat_list[[i]] <- subset(seurat_list[[i]], subset = nFeature_RNA > 200 &
                               percent.globin < 10 & predicted_doublet == 'False')
  seurat_list[[i]] <- SCTransform(seurat_list[[i]],
                                  method = "glmGamPoi",
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = 4000,
                                  verbose = TRUE)
}

# Cell Cycle Scoring ---------------------------
g2m_genes <- readRDS('g2m_genes.rds')
s_genes <- readRDS('s_genes.rds')
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

# WORKSPACE SAVED HERE 04182021
rm(seurat_list, qc_plot) # to conserve memory

combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                          features.to.integrate = to_integrate, verbose = T,
                          new.assay.name = "integratedSCT_")

# Regress out S and G2M scores ---------------------------
# This will take some time. 
combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"))

saveRDS(combined, 'combined_SCT_DSB120clean_pre_adt_integration.rds')

# ADT Integration ---------------------------
seurat_list <- anchors@object.list

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- rownames(seurat_list[[i]][["ADT"]]) #select all ADT as variable features
  # DSB ALREADY NORMALIZED
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
combinedADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I am not sure if this is a good idea:
combined[["integratedADT_"]] <- combinedADT[["integratedADT_"]] 

rm(combinedADT, anchorsADT, anchors, seurat_list)
saveRDS(combined, "combined_integrated_DSB120clean.rds")
# combined <- readRDS("combined_integrated_DSB_SCT_scrub.rds")

# Separate RNA and ADT Clustering ---------------------------
# PCA on SCT
DefaultAssay(combined) <- "integratedSCT_"
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)
combined <- RunPCA(combined, verbose = T)
ElbowPlot(combined, ndims = 50, reduction = 'pca') # will use 40 because SCT is more accurate 
ggsave("elbow_plot_sct_dsb120clean.png")

# aPCA on ADT
DefaultAssay(combined) <- 'integratedADT_'
combined <- ScaleData(combined)
combined <- RunPCA(combined, reduction.name = 'apca')
ElbowPlot(combined, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("elbow_plot_dsb120clean.png")

#Idents(combined) <- "old.ident"
combined <- RunUMAP(combined, reduction = 'pca', dims = 1:40, assay = 'SCT', 
                    reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
combined <- RunUMAP(combined, reduction = 'apca', dims = 1:15, assay = 'ADT', 
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
ggsave("SCT_or_ADT_only_clustering_DSB120.png", width = 10, height = 5)

DefaultAssay(combined) <- 'ADT'
FeaturePlot(combined, features = c('CD4-TotalA', 'CD8a-TotalA', 'CD8b-TotalA', 'TCRB-TotalA'), 
            reduction = 'adt.umap', cols = c("lightgrey","darkgreen"),
            min.cutoff = 0, max.cutoff = 'q99')
ggsave('tcell_markers_adt_proj120.png', width = 10, height = 10)

FeaturePlot(combined, features = c('CD4-TotalA', 'CD8a-TotalA', 'CD8b-TotalA', 'TCRB-TotalA'), 
            reduction = 'sct.umap', cols = c("lightgrey","darkgreen"),
            min.cutoff = 0, max.cutoff = 'q99')
ggsave('tcell_markers_sct_proj120.png', width = 10, height = 10)
DefaultAssay(combined) <- 'integratedSCT_'

DimPlot(subset(combined, idents = '17'), reduction = 'sct.umap')
ggsave('ident17_on_sct_umap.png', width = 5, height = 5)

FeaturePlot(combined, features = c('doublet_score'), reduction = 'sct.umap')
ggsave('doublet_score_after_scrublet_sct_umap.png', width = 5, height = 5)

FeaturePlot(combined, features = c('doublet_score'), reduction = 'adt.umap')
ggsave('doublet_score_after_scrublet_adt_umap.png', width = 5, height = 5)

DimPlot(combined, reduction = 'sct.umap', label = T) + NoLegend()

saveRDS(combined, 'combined_dsb120.rds')
#combined <- readRDS('combined_dsb120.rds')
combined@assays$RNA <- NULL
saveRDS(combined, 'combined_dsb120_NORNA.rds')
# combined <- readRDS('combined_dsb120_NORNA.rds')


# # NEED TO FIGURE OUT WHICH ADTs ARE GOOD
# # TEST Subset out cluster 17 (17 based on adt only) ---------------------------
# idents <- as.character(0:(length(unique(Idents(combined)))-1))
# idents <- idents[-18] # removes ident 17
# combined_clean <- subset(combined, idents = idents)

## RUN: TEST_combined_clean.R script

# WNN clustering ---------------------------
# On the integrated RNA and integrated ADT data
# Identify multimodal neighbors:

combined <- FindMultiModalNeighbors(
  combined, reduction.list = list("pca", "apca"), 
  dims.list = list(1:40, 1:15), modality.weight.name = "SCT.weight")

resolution.range <- seq(from = 0, to = 2.0, by = 0.1)

# MATTER TO RUN UMAP FIRST? this is what they do in vignette
# combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Find clusters using a range of resolutions
combined <- FindClusters(object = combined, graph.name = "wsnn",
                         reduction.name = "wnn.umap", algorithm = 3, 
                         resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("wsnn_res.", res, sep = "")
  Idents(combined) <- md
  p <- DimPlot(combined, label = T, repel = T, label.size = 3) + NoLegend()
  ggsave(paste("dimplot_", "res_", res, "_dsb120_v1.png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(combined, prefix = "wsnn_res.")
ggsave("clustree_output_wnn_dsb120.pdf", width = 9.5, height = 11*2)

combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
res <- 1.5
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)

# Data visualization:
DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
ggsave(paste("dimplot_", "res_", res, "_dsb120_v1.png", sep = ""), width = 5, height = 5)

# Saving file to work off of:
# saveRDS(combined, file = "combined.WNN.SCT.CLR_78_withRNA.rds")
# combined@assays$RNA <- NULL
# saveRDS(combined, file = "combined.WNN.SCT.CLR_78.rds")
# combined <- readRDS("02-05-2020 Clustering.V2.rds")
# combined <- readRDS("combined.WNN.V2.rds") # From April 2021
# combined <- readRDS("combined.WNN.V2.1.rds") # With dims=1:10 for adt PCs
# combined <- readRDS("combined.WNN.SCT.CLR_78.rds") 

### ADD 0 CUTOFF TO ALL ADT VISUALIZATIONS ON FEATURE PLOT

DefaultAssay(combined) <- 'ADT'
FeaturePlot(combined, features = c('CD4-TotalA', 'CD8a-TotalA', 'CD8b-TotalA', 'TCRB-TotalA'), 
            reduction = 'wnn.umap', cols = c("lightgrey","darkgreen"),
            min.cutoff = 0, max.cutoff = 'q99')
ggsave('tcell_markers_adt_wnn120.png', width = 10, height = 10)
DefaultAssay(combined) <- 'SCT'
FeaturePlot(combined, features = c('Cd3e', 'Cd4', 'Cd8a', 'Trdv4'), 
            reduction = 'wnn.umap', 
            min.cutoff = 0, max.cutoff = 'q99')
ggsave('tcell_markers_sct_wnn120.png', width = 10, height = 10)

VlnPlot(combined, features = "integratedSCT_.weight", sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(paste("SCT_weight_", "res", res, "_dsb120.png", sep = ""), width = 10, height = 5)

VlnPlot(combined, features = "nCount_RNA", sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(paste("nCount_RNA_", "res", res, "_dsb120.png", sep = ""), width = 10, height = 5)

VlnPlot(combined, features = "nFeature_RNA", sort = TRUE, pt.size = 0) +
  NoLegend()
ggsave(paste("nFeature_RNA", "res", res, "_dsb120.png", sep = ""), width = 10, height = 5)

VlnPlot(combined, features = "nCount_ADT", sort = TRUE, pt.size = 0.1) +
  NoLegend()
ggsave(paste("nCount_ADT", "res", res, "_dsb120.png", sep = ""), width = 10, height = 5)

VlnPlot(combined, features = "nFeature_ADT", sort = TRUE, pt.size = 0.1) +
  NoLegend()
ggsave(paste("nFeature_ADT_V2", "res", res, "_dsb120.png", sep = ""), width = 10, height = 5)

FeaturePlot(combined, features = c('doublet_score'))

FeaturePlot(combined, features = c('S.Score', 'G2M.Score'), reduction = 'wnn.umap', 
            min.cutoff = 0)
ggsave('cell_cycle_scores.wnn.png', width = 10, height = 5)

DefaultAssay(combined) <- 'SCT'
FeaturePlot(combined, features = c('Cd209a', 'Itgae'), reduction = 'wnn.umap')

# Celltype Identification ---------------------------
# Run that script

# Marker Genes ---------------------------
## Exporting top 10 marker genes:
DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(combined.markers, file = "gene_combined_cluster_biomarkers_unlab.csv", row.names = FALSE)
write.csv(top10, file = "gene_combined_cluster_biomarkers_unlab_top10.csv", row.names = FALSE)

DefaultAssay(combined) <- "ADT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(combined.markers, file = "adt_combined_cluster_biomarkers_unlab.csv", row.names = FALSE)
write.csv(top10, file = "adt_combined_cluster_biomarkers_unlab_top10.csv", row.names = FALSE)

## Plotting top 5 marker genes on heatmap: 
# top5 <- combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# p <- DoHeatmap(subset(combined, downsample = 100),
#                features = top5$gene, size = 3, angle = 30) + NoLegend()
# ggsave("top5_markers_heatmap_unlab.png", plot = p, width = 10, height = 6)

saveRDS(combined, 'combined_04182021.rds')
# combined <- readRDS('combined_04182021.rds')

# Rename Idents ---------------------------
# Saving old labels:
combined[["old.ident"]] <- Idents(object = combined)

# N = Neutrophil, T = T cell, B = B cell, gCap = gCap Endothelial, 
# AM = Alveolar Macrophage, IM = Intersitial Macro, Vein = Vein Endo,
# Macro = UNKNOWN CURRENTLY MACROPHAGE

# Renaming:
combined <- RenameIdents(object = combined, 
                         '0' = "N 1", 
                         '1' = "N 2", 
                         '2' = "gCap 1",
                         '3' = "gCap 2",
                         '4' = "B 1",
                         '5' = "AM 1",
                         '6' = "AT2 1",
                         '7' = "gCap 3",
                         '8' = "CD4 T 1",
                         '9' = "NK",
                         '10' = "N 3",
                         '11' = "AM 2",
                         '12' = "C Mono",
                         '13' = "CD4 T 2",
                         '14' = "CD8 T",
                         '15' = "Myofib",
                         '16' = "NC Mono ",
                         '17' = "Lipofib 1",
                         '18' = "AT1",
                         '19' = "N 4",
                         '20' = "N 5",
                         '21' = "IM",
                         '22' = "CD4 T 3",
                         '23' = "aCap",
                         '24' = "cDC 1",
                         '25' = "Nuocyte",
                         '26' = "B 2",
                         '27' = "Efb1 Fib",
                         '28' = "Treg",
                         '29' = "CD4 T 4",
                         '30' = "AM 3",
                         '31' = "Vein",
                         '32' = "Baso",
                         '33' = "Macro",
                         '34' = "pDC",
                         '35' = "AM 4",
                         '36' = "Lipofib 2",
                         '37' = "N 6",
                         '38' = "AT2 2",
                         '39' = "Lymph Fib",
                         '40' = "B 3",
                         '41' = "AM 5",
                         '42' = "cDC 2",
                         '43' = "Clara",
                         '44' = "Mtx Fib")
combined[["celltype"]] <- Idents(combined)
# p4 <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
# p4

p <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
ggsave("umap_res1.5_lab6.png", height = 6, width = 6)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_res1.5_lab_clean.png", height = 5, width = 5)

# define a function that gets the colors of the clean UMAP
get_colors <- function(p, obj){
  # Formatting a vector of colors for DittoSeq plot coloring
  colors <- ggplot_build(p)$data[[1]]$colour
  unique_colors <- unique(colors)
  color_freq <- 1:length(unique_colors)
  for (i in 1:length(unique_colors)){
    color <- unique_colors[[i]]
    color_freq[[i]] <- length(grep(color, colors))
  }
  color_freq_desc <- sort(color_freq, decreasing = T)
  color_names_by_freq <- as.character(1:length(color_freq))
  for (i in 1:length(color_freq_desc)){
    color_count <- color_freq_desc[[i]]
    idx <- which(color_count == color_freq)
    color_names_by_freq[[i]] <- unique_colors[[idx]]
  }
  
  t3 <- table(obj@active.ident)
  #write.csv(t2, file = "t3.csv")
  
  cell_df <- as.data.frame(t3)
  color_df <- cbind(data.frame(unique_colors), data.frame(color_freq))
  color_df <- merge(cell_df, color_df, by.x = "Freq", by.y = "color_freq")
  colnames(color_df) <- c("Freq", "Celltype", "Color")
  color_df <- color_df[order(tolower(color_df$Celltype)),]
  
  return(color_df)
}

color_df <- get_colors(p, combined)

# plotting 100% stacked bar chart for cell type distribution by group:
p <- dittoBarPlot(object = combined, var = combined@active.ident, group.by = "orig.ident",
                  scale = c('percent'), color.panel = color_df$Color, 
                  x.labels.rotate = F, ylab = 'Fraction of Cells')
p$data$grouping <- factor(x = p$data$grouping, levels = c("naive", "p4", "mp4", "p24", "mp24"))
p + theme(axis.title.x=element_blank())
ggsave("celltype_distribution_by_trt_leg.png", width = 6, height = 4.2)
p + NoLegend()
ggsave("celltype_distribution_by_trt_noLeg.png", width = 4.2, height = 4.2)

# Clusters with Legend
p <- DimPlot(combined, reduction = 'wnn.umap', label = FALSE, repel = TRUE, pt.size = 0.1)+ 
  theme(text = element_text(size=8, family = "sans"),  
        axis.title=element_text(size=8,family = "sans", face="bold"))
ggsave("rawClusters.png", plot = p, width = 5, height = 4)

# Dimplot Split by Treatment
p <- DimPlot(combined, reduction = 'wnn.umap', split.by = "orig.ident") + NoLegend()
# change order in plot: 
p$data$orig.ident <- factor(x = p$data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))
p
ggsave("dimplot_splitby_trt.png", plot = p, width = 15, height = 5)

# Cluster by Group
p5 <- DimPlot(combined, reduction = 'wnn.umap', group.by = "orig.ident") + 
  ggtitle(NULL) + NoLegend() + 
  theme(text = element_text(size=8, family = "sans"),  
        axis.title=element_text(size=8,family = "sans", face="bold"))
ggsave("clusterByGroup.png", plot = p5, width = 5, height = 5)

# rearrange orig.ident
combined$orig.ident <- factor(combined$orig.ident, 
                              levels = c('naive', 'p4', 'mp4', 'p24', 'mp24'))
DimPlot(combined, split.by = "orig.ident") + NoLegend() # works

# Tables:
# number of cells per treatment
t1 <- table(combined$orig.ident)
# number of cells per type per treatment
t2 <- table(combined$celltype, col.names = combined$orig.ident)
col_order <- c('naive', 'p4', 'mp4', 'p24', 'mp24')
t2 <- t2[ ,col_order]
# number of cells per cluster
t3 <- table(combined@active.ident)

write.csv(t1, file = "t1.csv", row.names = F)
write.csv(t2, file = "t2.csv", row.names = T)
write.csv(t3, file = "t3.csv", row.names = F)

rm(p, p4, p5, p6)

saveRDS(combined, 'combined_04192021.rds')
# combined <- readRDS('combined_04192021.rds')

# Condensed Renaming  ---------------------------
Idents(combined) <- 'old.ident'
combined <- RenameIdents(object = combined, 
                         '0' = "N B", 
                         '1' = "N A", 
                         '2' = "gCap 1",
                         '3' = "gCap 2",
                         '4' = "B 1",
                         '5' = "AM 1",
                         '6' = "AT2 1",
                         '7' = "gCap 3",
                         '8' = "CD4 T 1",
                         '9' = "NK",
                         '10' = "N C",
                         '11' = "AM 2",
                         '12' = "C Mono",
                         '13' = "CD4 T 2",
                         '14' = "CD8 T",
                         '15' = "Myofib",
                         '16' = "NC Mono ",
                         '17' = "Lipofib 1",
                         '18' = "AT1",
                         '19' = "N C",
                         '20' = "N A",
                         '21' = "IM",
                         '22' = "CD4 T 3",
                         '23' = "aCap",
                         '24' = "cDC 1",
                         '25' = "Nuocyte",
                         '26' = "B 2",
                         '27' = "Efb1 Fib",
                         '28' = "Treg",
                         '29' = "CD4 T 4",
                         '30' = "AM 3",
                         '31' = "Vein",
                         '32' = "Baso",
                         '33' = "Macro",
                         '34' = "pDC",
                         '35' = "AM 4",
                         '36' = "Lipofib 2",
                         '37' = "N B",
                         '38' = "AT2 2",
                         '39' = "Lymph Fib",
                         '40' = "B 3",
                         '41' = "AM 5",
                         '42' = "cDC 2",
                         '43' = "Clara",
                         '44' = "Mtx Fib")
combined[["condensed"]] <- Idents(combined)

p <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
ggsave("umap_condensed.png", height = 6, width = 6)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_condensed_clean.png", height = 5, width = 5)

# Neutrophil (N) Subclustering  ---------------------------
Idents(combined) <- 'celltype'
neutro <- subset(combined, 
                 idents = c('N 1', 'N 2', 'N 3', 'N 4', 'N 5', 'N 6'))

DefaultAssay(neutro) <- 'RNA'
neutro <- DietSeurat(neutro, assays = c("RNA", "ADT"))
dim(neutro)

# Confirm subsetted metadata
table(neutro@meta.data$seurat_clusters)

# split data to rerun workflow
seurat_list <- SplitObject(neutro, split.by = "orig.ident")

# Rerun analysis on subsetted data:
for (i in 1:length(x = seurat_list)){
  sample_id <- names(seurat_list)[[i]]
  print(sample_id)
  
  seurat_list[[i]] <- SCTransform(seurat_list[[i]],
                                  method = "glmGamPoi",
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = 4000,
                                  verbose = TRUE)
}

# > N Cycle Scoring ---------------------------
# g2m_genes <- readRDS('g2m_genes.rds')
# s_genes <- readRDS('s_genes.rds')
# for (i in 1:length(x = seurat_list)){ # I can move this into the above for loop...
#   seurat_list[[i]] <- CellCycleScoring(seurat_list[[i]], g2m.features=g2m_genes, s.features=s_genes)
# }

# > N SCT Integration ---------------------------
# NOTE TO SELF_NEED TO MAKE SURE Cd3e, Cd3d, Cd4, Cd8a, Cd8b1 are included!
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")
# Integrate 6000 genes, but only use the anchors from the 3000 as anchorset

save.image(file='neutro_pre_int04192021.RData')

to_integrate <- SelectIntegrationFeatures(object.list = seurat_list, 
                                          nfeatures = 6000 , verbose = TRUE,
                                          new.assay.name = "integratedSCT_")
rm(seurat_list)

neutro <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                        features.to.integrate = to_integrate, verbose = T)
#saveRDS(neutro, file = "neutro_sct_integrated.rds")

# > N Regress out S and G2M scores ---------------------------
# This will take some time. 
neutro <- ScaleData(neutro, vars.to.regress = c("S.Score", "G2M.Score"))

# > N ADT Integration ---------------------------
seurat_list <- anchors@object.list

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- rownames(seurat_list[[i]][["ADT"]]) #select all ADT as variable features
  # DSB ALREADY NORMALIZED
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
neutroADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I am not sure if this is a good idea:
neutro[["integratedADT_"]] <- neutroADT[["integratedADT_"]] 

rm(neutroADT, anchorsADT, anchors, seurat_list)
saveRDS(neutro, "neutro_integrated.rds")

neutro[['integratedSCT_']] <- neutro[['integrated']]

# > N Select PCs ---------------------------
DefaultAssay(neutro) <- "integratedSCT_"
all.genes <- rownames(neutro)
neutro <- ScaleData(neutro, features = all.genes)
neutro <- RunPCA(neutro, verbose = T)
ElbowPlot(neutro, ndims = 50, reduction = 'pca') # will use 40 because SCT is more accurate 
ggsave("neutro_elbow_plot_sct_dsb120clean.png")

# aPCA on ADT
DefaultAssay(neutro) <- 'integratedADT_'
neutro <- ScaleData(neutro)
neutro <- RunPCA(neutro, reduction.name = 'apca')
ElbowPlot(neutro, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("neutro_elbow_plot_dsb120clean.png")

#Idents(neutro) <- "old.ident"
neutro <- RunUMAP(neutro, reduction = 'pca', dims = 1:30, 
                  assay = 'integratedSCT_', 
                  reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
neutro <- RunUMAP(neutro, reduction = 'apca', dims = 1:10, 
                  assay = 'integratedADT_', 
                  reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

neutro <- FindNeighbors(neutro, reduction = 'pca', dims = 1:30, 
                          graph.name = 'sct.snn')
neutro <- FindClusters(neutro, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

neutro <- FindNeighbors(neutro, reduction = 'apca', dims = 1:10,
                          graph.name = 'adt.snn')
neutro <- FindClusters(neutro, graph.name = 'adt.snn', resolution = 1.0, verbose = T)

Idents(neutro) <- 'sct.snn_res.1'
p1 <- DimPlot(neutro, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
Idents(neutro) <- 'adt.snn_res.1'
p2 <- DimPlot(neutro, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_neutro.png", width = 10, height = 5)

# Using original labels from combined.rds
Idents(neutro) <- 'old.ident'
p1 <- DimPlot(neutro, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(neutro, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_neutro_origlabel.png", width = 10, height = 5)

# > N WNN clustering ---------------------------
# On the integrated SCT and integrated ADT data
# Identify multimodal neighbors:

neutro <- FindMultiModalNeighbors(
  neutro, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:10), modality.weight.name = "SCT.weight")

resolution.range <- seq(from = 0, to = 1.0, by = 0.1)

# WNN - Find clusters using a range of resolutions
neutro <- FindClusters(object = neutro, graph.name = "wsnn",
                         reduction.name = "wnn.umap", algorithm = 3, 
                         resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("wsnn_res.", res, sep = "")
  Idents(neutro) <- md
  p <- DimPlot(neutro, label = T, repel = T, label.size = 3) + NoLegend()
  ggsave(paste("neutro_dimplot_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(neutro, prefix = "wsnn_res.")
ggsave("clustree_output_neutro_WNN.pdf", width = 9.5, height = 11)

# neutro <- RunUMAP(neutro, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# res <- 1.5
# neutro <- FindClusters(neutro, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)
# 
# # Data visualization:
# DimPlot(neutro, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
# ggsave(paste("dimplot_", "res_", res, "neutro.png", sep = ""), width = 5, height = 5)

# > N SCT clustering  ---------------------------
# WNN CLUSTERING IS UNSTABLE - USING SCT CLUSTERING ONLY FOR NEUTROPHILS

resolution.range <- seq(from = 0, to = 1.5, by = 0.1)

# SCT - Find clusters using a range of resolutions
neutro <- FindClusters(object = neutro, graph.name = "sct.snn",
                       resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("sct.snn_res.", res, sep = "")
  Idents(neutro) <- md
  p <- DimPlot(neutro, label = T, repel = T, label.size = 3, 
               reduction = 'sct.umap') + NoLegend()
  ggsave(paste("neutro_dimplot_sct_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(neutro, prefix = "sct.snn_res.")
ggsave("clustree_output_neutro_WNN.pdf", width = 9.5, height = 11)

neutro <- FindClusters(neutro, graph.name = "sct.snn",
                       resolution = 0.5, verbose = T)

# try a 3D plot: 
# run 3D_UMAP_plots.R
# run clustering part of citeseq_cell_markers.R script

# > N alternate projections ---------------------------
load('pre_webGL_workspace.RData')

# JUST KILL THE ASSAY BRO: SO RUN THIS
# neutro[['pca']] <- NULL
# DefaultAssay(neutro) <- "integratedSCT_"
# all.genes <- rownames(neutro)
# neutro <- ScaleData(neutro, features = all.genes)
# neutro <- RunPCA(neutro, verbose = T)

Idents(neutro) <- 'old.ident'
neutro <- RenameIdents(object = neutro, 
                         '0' = "N B", 
                         '1' = "N A", 
                         '2' = "gCap 1",
                         '3' = "gCap 2",
                         '4' = "B 1",
                         '5' = "AM 1",
                         '6' = "AT2 1",
                         '7' = "gCap 3",
                         '8' = "CD4 T 1",
                         '9' = "NK",
                         '10' = "N C",
                         '11' = "AM 2",
                         '12' = "C Mono",
                         '13' = "CD4 T 2",
                         '14' = "CD8 T",
                         '15' = "Myofib",
                         '16' = "NC Mono ",
                         '17' = "Lipofib 1",
                         '18' = "AT1",
                         '19' = "N C",
                         '20' = "N A",
                         '21' = "IM",
                         '22' = "CD4 T 3",
                         '23' = "aCap",
                         '24' = "cDC 1",
                         '25' = "Nuocyte",
                         '26' = "B 2",
                         '27' = "Efb1 Fib",
                         '28' = "Treg",
                         '29' = "CD4 T 4",
                         '30' = "AM 3",
                         '31' = "Vein",
                         '32' = "Baso",
                         '33' = "Macro",
                         '34' = "pDC",
                         '35' = "AM 4",
                         '36' = "Lipofib 2",
                         '37' = "N B",
                         '38' = "AT2 2",
                         '39' = "Lymph Fib",
                         '40' = "B 3",
                         '41' = "AM 5",
                         '42' = "cDC 2",
                         '43' = "Clara",
                         '44' = "Mtx Fib")
neutro[["condensed"]] <- Idents(neutro)

# pca_ key is the one to use for Neutrophils, PC_ was already calcualted

DimPlot(neutro, reduction = 'pca')

DimPlot(neutro, reduction = 'pca', group.by = 'celltype')
ggsave('neutro_pca_OG_celltype.png', width = 5, height = 5)

DimPlot(neutro, reduction = 'pca', group.by = 'condensed')
ggsave('neutro_pca_condensed_celltype.png', width = 5, height = 5)

p <- DimPlot(neutro, reduction = 'pca', group.by = 'celltype', split.by = 'orig.ident')
p$data$orig.ident <- factor(x = p$data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))
p
ggsave('neutro_pca_OG_celltype_orig.ident.png', width = 15, height = 3.5)

p <- DimPlot(neutro, reduction = 'pca', group.by = 'condensed', split.by = 'orig.ident')
p$data$orig.ident <- factor(x = p$data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))
p
ggsave('neutro_pca_condensed_orig.ident.png', width = 15, height = 3.5)

FeaturePlot(neutro, reduction = 'pca', features = 'Mmp8')
FeaturePlot(neutro, reduction = 'pca', features = 'Ccl3')

DefaultAssay(neutro) <- 'SCT'
FeaturePlot(neutro, reduction = 'pca', features = c('Mmp8', 'Ccl3'), blend = T)
ggsave('neutro_pca_co_expression.png', width = 16, height = 4)

VizDimLoadings(neutro, dims = 1:2, reduction = "pca")
ggsave('neutro_pca_dim_loadings.png', width = 10, height = 5)

# pca_1 pos
FeaturePlot(neutro, reduction = 'pca', ncol = 3,
            features = c('B2m', 'Ctsb', 'Ctsz', 'Bri3', 'Atp6v0c'))
ggsave('pca_1_pos_neutro.png', width = 5*3, height = 5*2)

# pca_1 neg
FeaturePlot(neutro, reduction = 'pca', ncol = 3,
            features = c('S100a8', 'Alox5ap', 'S100a11', 'S100a9', 'Retnlg', 'Anxa1'))
ggsave('pca_1_neg_neutro.png', width = 5*3, height = 5*2)

FeaturePlot(neutro, reduction = 'pca', features = 'Camp') # Cathelicidin

FeaturePlot(neutro, features = 'Cxcr2', reduction = 'pca')


# Structural Cell Subclustering  ---------------------------
# AKA Endothelial / Epithelial / Fibroblast
Idents(combined) <- 'celltype'
struct <- subset(combined, idents = c('Lymph Fib',
                                      'Myofib',
                                      'Lipofib 2',
                                      'gCap 1',
                                      'Vein',
                                      'gCap 3',
                                      'gCap 2',
                                      'aCap',
                                      'Lipofib 1',
                                      'Efb1 Fib',
                                      'Mtx Fib',
                                      'AT1',
                                      'AT2 2',
                                      'AT2 1'))

DefaultAssay(struct) <- 'RNA'
struct <- DietSeurat(struct, assays = c("RNA", "ADT"))
dim(struct)

# Confirm subsetted metadata
table(struct@meta.data$seurat_clusters)

# split data to rerun workflow
seurat_list <- SplitObject(struct, split.by = "orig.ident")

# Rerun analysis on subsetted data:
for (i in 1:length(x = seurat_list)){
  sample_id <- names(seurat_list)[[i]]
  print(sample_id)
  
  seurat_list[[i]] <- SCTransform(seurat_list[[i]],
                                  method = "glmGamPoi",
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = 4000,
                                  verbose = TRUE)
}

# > Struct SCT Integration ---------------------------
# NOTE TO SELF_NEED TO MAKE SURE Cd3e, Cd3d, Cd4, Cd8a, Cd8b1 are included!
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")

# save.image(file='struct_pre_int04192021.RData')

# to_integrate <- SelectIntegrationFeatures(object.list = seurat_list, 
#                                           nfeatures = 6000 , verbose = TRUE,
#                                           new.assay.name = "integratedSCT_")
rm(seurat_list, combined)

struct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                        new.assay.name = 'integratedSCT_', verbose = T)
saveRDS(struct, file = "struct_sct_integrated.rds")

# > Struct Regress out S and G2M scores ---------------------------
# This will take some time. 
struct <- ScaleData(struct, vars.to.regress = c("S.Score", "G2M.Score"))

# > Struct ADT Integration ---------------------------
seurat_list <- anchors@object.list

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- rownames(seurat_list[[i]][["ADT"]]) #select all ADT as variable features
  # DSB ALREADY NORMALIZED
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
structADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I am not sure if this is a good idea:
struct[["integratedADT_"]] <- structADT[["integratedADT_"]] 

rm(structADT, anchorsADT, anchors, seurat_list)
saveRDS(struct, "struct_integrated.rds")
# struct <- readRDS("struct_integrated.rds")

# > Struct Select PCs ---------------------------
DefaultAssay(struct) <- "integratedSCT_"
all.genes <- rownames(struct)
struct <- ScaleData(struct, features = all.genes)
struct <- RunPCA(struct, verbose = T)
ElbowPlot(struct, ndims = 50, reduction = 'pca') # will use 30 because SCT is more accurate 
ggsave("struct_elbow_plot_sct_dsb120clean.png")

# aPCA on ADT
DefaultAssay(struct) <- 'integratedADT_'
struct <- ScaleData(struct)
struct <- RunPCA(struct, reduction.name = 'apca')
ElbowPlot(struct, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("struct_elbow_plot_dsb120clean.png")

#Idents(struct) <- "old.ident"
struct <- RunUMAP(struct, reduction = 'pca', dims = 1:30, 
                  assay = 'integratedSCT_', 
                  reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
struct <- RunUMAP(struct, reduction = 'apca', dims = 1:10, 
                  assay = 'integratedADT_', 
                  reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

struct <- FindNeighbors(struct, reduction = 'pca', dims = 1:30, 
                        graph.name = 'sct.snn')
struct <- FindClusters(struct, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

struct <- FindNeighbors(struct, reduction = 'apca', dims = 1:10,
                        graph.name = 'adt.snn')
struct <- FindClusters(struct, graph.name = 'adt.snn', resolution = 1.0, verbose = T)

Idents(struct) <- 'sct.snn_res.1'
p1 <- DimPlot(struct, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
Idents(struct) <- 'adt.snn_res.1'
p2 <- DimPlot(struct, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_struct.png", width = 10, height = 5)

# Using original labels from combined.rds
Idents(struct) <- 'celltype'
p1 <- DimPlot(struct, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(struct, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_struct_origlabel.png", width = 10, height = 5)

# > Struct WNN clustering ---------------------------
# On the integrated SCT and integrated ADT data
# Identify multimodal neighbors:

struct <- FindMultiModalNeighbors(
  struct, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:10), modality.weight.name = "SCT.weight")

struct <- RunUMAP(struct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

resolution.range <- seq(from = 0, to = 2.0, by = 0.1)

# WNN - Find clusters using a range of resolutions
struct <- FindClusters(object = struct, graph.name = "wsnn",
                       reduction.name = "wnn.umap", algorithm = 3, 
                       resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("wsnn_res.", res, sep = "")
  Idents(struct) <- md
  p <- DimPlot(struct, label = T, repel = T, label.size = 3, reduction = 'wnn.umap') + NoLegend()
  ggsave(paste("struct_dimplot_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(struct, prefix = "wsnn_res.")
ggsave("clustree_output_struct_WNN.pdf", width = 9.5, height = 11)

# struct <- RunUMAP(struct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# res <- 1.5
# struct <- FindClusters(struct, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)
# 
# # Data visualization:
# DimPlot(struct, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
# ggsave(paste("dimplot_", "res_", res, "struct.png", sep = ""), width = 5, height = 5)

# > Struct SCT clustering  ---------------------------
# WNN CLUSTERING IS UNSTABLE - USING SCT CLUSTERING ONLY FOR structPHILS

resolution.range <- seq(from = 0, to = 1.0, by = 0.1)

# SCT - Find clusters using a range of resolutions
struct <- FindClusters(object = struct, graph.name = "sct.snn",
                       resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("sct.snn_res.", res, sep = "")
  Idents(struct) <- md
  p <- DimPlot(struct, label = T, repel = T, label.size = 3, 
               reduction = 'sct.umap') + NoLegend()
  ggsave(paste("struct_dimplot_sct_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(struct, prefix = "sct.snn_res.")
ggsave("clustree_output_struct_SCT.pdf", width = 9.5, height = 11)

struct <- FindClusters(struct, graph.name = "sct.snn",
                       resolution = 1.0, verbose = T)

# > Struct Calculate DEG markers  ---------------------------
DefaultAssay(struct) <- "SCT"
struct.markers <- FindAllMarkers(struct, only.pos = FALSE, 
                                 min.pct = 0.1, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
all.markers <- struct.markers %>% group_by(cluster)
write.csv(all.markers, file = "struct_cluster_biomarkers_labeled.csv", row.names = FALSE)

## Plotting top 5 marker genes on heatmap: 
top5 <- struct.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p <- DoHeatmap(subset(struct, downsample = 100),
               features = top5$gene, size = 3, angle = 30) + NoLegend()
ggsave("struct_top5_markers_heatmap_labeled.png", plot = p, width = 10, height = 6)

# > Struct Gillich paper markers:   ---------------------------
  
  FeaturePlot(struct, features = 'doublet_score')
FeaturePlot(struct, features = 'H2-Ab1')

  # aCap -  (Same as Endothelial kdr high)
  markers.to.plot <- c("Car4", "Ednrb", "Fibin", "Tbx2", "Cdkn2b", "Rprml", "Chst1", "Apln")
  DotPlot(struct, assay = "integratedSCT_", features  = markers.to.plot)
  ggsave('Gillich_aCap_struct.png', width = 8.5, height = 7)
  
  # gCap - 2 and 3 most strongly, 10 perhaps (debatable)
  markers.to.plot <- c("Cd93", "Ptprb", "Plvap", "Gpihbp1", "H2-Ab1", "Tek", "Kit", "Aplnr")
  DotPlot(struct, assay = "integratedSCT_", features  = markers.to.plot)
  ggsave('Gillich_gCap_struct.png', width = 8.5, height = 7)
  
  # Artery - N/A
  markers.to.plot <- c("Mgp", "Cdh13", "Htra1", "Bmx", "Gja5", "Fbln2", "Sulf1")
  DotPlot(struct, assay = "integratedSCT_", features  = markers.to.plot)
  ggsave('Gillich_artery_struct.png', width = 8.5, height = 7)
  
  # Vein - Same as endothelial cells Vwf hi, already labeled
  markers.to.plot <- c("Vwf", "Slc6a2", "Bst1", "Car8", "Amigo2", "Mustin1", "Vegfc", "Csrp2", "Nr2f2")
  DotPlot(struct, assay = "integratedSCT_", features  = markers.to.plot)
  ggsave('Gillich_vein_endo_struct.png', width = 8.5, height = 7)
  
  # Lymphatics - Same as lymphatic fibroblasts, already labeled
  markers.to.plot <- c("Mmrn1", "Fxyd6", "Reln", "Pdpn", "Thy1", "Nrp2", "Tbx1", "Gja1")
  DotPlot(struct, assay = "integratedSCT_", features  = markers.to.plot)
  ggsave('Gillich_lymphatics_struct.png', width = 8.5, height = 7)
  
  # Markers to define metaclusters
  markers.to.plot <- c('Ptprc', 'Pecam1', 'Prox1', 'Hopx', 'Sftpc',
                       'Trbc2', 'Ccl5', 'Retnlg', 'Mcpt8', 'Mcpt4', 'F13a1', 
                       'Epcam', 'Cdh5', 'Col1a2', 'Scgb3a2', 'Foxj1', 'Enpp2', 
                       'Mfap4', 'Ly6c2', "Akap5", "Lamp3", "Foxp3", 'Gpihbp1',
                       'Col1a1', 'Acta2')
  DefaultAssay(struct) <- "SCT"
  for (marker in markers.to.plot){
    FeaturePlot(struct, features = marker, reduction = 'sct.umap')
    filename <- paste("FeaturePlot_", marker, "_struct.png", sep = "")
    ggsave(filename = filename, width = 5, height = 5)
  }
  DefaultAssay(struct) <- "integratedSCT_"
  
struct <- RenameIdents(struct, 
                       '0' = 'gCap',
                       '1' = 'gCap',
                       '2' = 'AT 2',
                       '3' = 'gCap',
                       '4' = 'gCap',
                       '5' = 'AT 2',
                       '6' = 'aCap',
                       '7' = 'Myofib',
                       '8' = 'gCap',
                       '9' = 'gCap',
                       '10' = 'Lipofib',
                       '11' = 'Ebf1+ Fib',
                       '12' = 'Alv Fib',
                       '13' = 'AT 1',
                       '14' = 'Vein',
                       '15' = 'DOUBLET',
                       '16' = 'AT 1',
                       '17' = 'Adv Fib',
                       '18' = 'gCap',
                       '19' = 'AT 2',
                       '20' = 'Lymph Fib')

p <- DimPlot(struct, reduction = 'sct.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_labeled_struct_subcluster.png", height = 5, width = 5)

rm(all.markers, p, p1, p2, struct.markers, top5)
# > Struct apply to parent --------------------------- 

# struct <- readRDS('struct_04_23_2021.rds')
combined <- readRDS('combined_04192021.rds')

# Generate a new column called sub_cluster in the metadata
combined$sub_cluster <- Idents(combined)
Idents(combined) <- 'sub_cluster'
combined <- SetIdent(combined, cells = Cells(struct), Idents(struct))
combined$sub_cluster <- Idents(combined)

DimPlot(combined, label = TRUE, repel = TRUE, label.size = 4,
        cells = Cells(struct), reduction = 'wnn.umap') + NoLegend()
ggsave('struct_subcluster_onUMAP.png', height = 5, width = 5)

DimPlot(combined, label = TRUE, repel = TRUE, label.size = 4,
        cells = Cells(struct), reduction = 'sct.umap') + NoLegend()
ggsave('struct_subcluster_onUMAP_SCT_sketchy.png', height = 5, width = 5)

combined$celltype <- combined$sub_cluster
Idents(combined) <- 'celltype'

saveRDS(combined, 'combined_04_22_2021.rds')
saveRDS(struct, 'struct_04_23_2021.rds')

rm(struct, p)

# Alveolar Macro Subclustering  ---------------------------
# AKA Endothelial / Epithelial / Fibroblast
Idents(combined) <- 'celltype'
am <- subset(combined, idents = c('AM 1',
                                  'AM 2',
                                  'AM 3',
                                  'AM 4',
                                  'AM 5'))

DefaultAssay(am) <- 'RNA'
am <- DietSeurat(am, assays = c("RNA", "ADT"))
dim(am)

# Confirm subsetted metadata
table(am@meta.data$seurat_clusters)

# split data to rerun workflow
seurat_list <- SplitObject(am, split.by = "orig.ident")

# Rerun analysis on subsetted data:
for (i in 1:length(x = seurat_list)){
  sample_id <- names(seurat_list)[[i]]
  print(sample_id)
  
  seurat_list[[i]] <- SCTransform(seurat_list[[i]],
                                  method = "glmGamPoi",
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = 3000, # changed to 3000
                                  verbose = TRUE)
}

# > AM SCT Integration ---------------------------
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")

# to_integrate <- SelectIntegrationFeatures(object.list = seurat_list, 
#                                           nfeatures = 6000 , verbose = TRUE,
#                                           new.assay.name = "integratedSCT_")
rm(seurat_list, combined)

# IMAGE SAVED

am <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                    new.assay.name = 'integratedSCT_', verbose = T)
saveRDS(am, file = "am_sct_integrated.rds")

# > AM Regress out S and G2M scores ---------------------------
# This will take some time. 
am <- ScaleData(am, vars.to.regress = c("S.Score", "G2M.Score"))

# > AM ADT Integration ---------------------------
seurat_list <- anchors@object.list

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- rownames(seurat_list[[i]][["ADT"]]) #select all ADT as variable features
  # DSB ALREADY NORMALIZED
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
amADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I am not sure if this is a good idea:
am[["integratedADT_"]] <- amADT[["integratedADT_"]] 

rm(amADT, anchorsADT, anchors, seurat_list)
saveRDS(am, "am_integrated.rds")
# am <- readRDS("am_integrated.rds")

# > AM Select PCs ---------------------------
DefaultAssay(am) <- "integratedSCT_"
all.genes <- rownames(am)
am <- ScaleData(am, features = all.genes)
am <- RunPCA(am, verbose = T)
ElbowPlot(am, ndims = 50, reduction = 'pca') # will use 30 because SCT is more accurate 
ggsave("am_elbow_plot_sct_dsb120clean.png")

# aPCA on ADT
DefaultAssay(am) <- 'integratedADT_'
am <- ScaleData(am)
am <- RunPCA(am, reduction.name = 'apca')
ElbowPlot(am, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("am_elbow_plot_dsb120clean.png")

#Idents(am) <- "old.ident"
am <- RunUMAP(am, reduction = 'pca', dims = 1:30, 
                  assay = 'integratedSCT_', 
                  reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
am <- RunUMAP(am, reduction = 'apca', dims = 1:10, 
                  assay = 'integratedADT_', 
                  reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

am <- FindNeighbors(am, reduction = 'pca', dims = 1:30, 
                        graph.name = 'sct.snn')
am <- FindClusters(am, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

am <- FindNeighbors(am, reduction = 'apca', dims = 1:10,
                        graph.name = 'adt.snn')
am <- FindClusters(am, graph.name = 'adt.snn', resolution = 1.0, verbose = T)

Idents(am) <- 'sct.snn_res.1'
p1 <- DimPlot(am, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
Idents(am) <- 'adt.snn_res.1'
p2 <- DimPlot(am, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_am.png", width = 10, height = 5)

# Using original labels from combined.rds
Idents(am) <- 'celltype'
p1 <- DimPlot(am, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(am, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_am_origlabel.png", width = 10, height = 5)

# > AM SCT clustering  ---------------------------

resolution.range <- seq(from = 0, to = 1.0, by = 0.1)

# SCT - Find clusters using a range of resolutions
am <- FindClusters(object = am, graph.name = "sct.snn",
                       resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("sct.snn_res.", res, sep = "")
  Idents(am) <- md
  p <- DimPlot(am, label = T, repel = T, label.size = 3, 
               reduction = 'sct.umap') + NoLegend()
  ggsave(paste("am_dimplot_sct_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(am, prefix = "sct.snn_res.")
ggsave("clustree_output_am_SCT.pdf", width = 9.5, height = 11)

am <- FindClusters(object = am, graph.name = "sct.snn",
                   resolution = 0.5, verbose = T)

# > AM Calculate DEG markers  ---------------------------
DefaultAssay(am) <- "SCT"
am.markers <- FindAllMarkers(am, only.pos = FALSE, 
                                 min.pct = 0.1, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
all.markers <- am.markers %>% group_by(cluster)
write.csv(all.markers, file = "am_cluster_biomarkers_labeled.csv", row.names = FALSE)

## Plotting top 5 marker genes on heatmap: 
top5 <- am.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p <- DoHeatmap(subset(am, downsample = 100),
               features = top5$gene, size = 3, angle = 30) + NoLegend()
ggsave("am_top5_markers_heatmap_labeled.png", plot = p, width = 10, height = 6)

saveRDS(am, 'am.rds')
rm(all.markers, am, am.markers, p, p1, p2, top5, all.genes, md, res, 
   resolution.range)

# > AM Rename original AM subclusters from main  ---------------------------
# This is so obvious, we don't need to do any returning to parent
combined <- readRDS('combined_04_22_2021.rds')
Idents(combined) <- 'sub_cluster'
unique(Idents(combined))
combined <- RenameIdents(combined, 
                         'AM 1' = 'AM',
                         'AM 2' = 'AM',
                         'AM 3' = 'AM',
                         'AM 4' = 'AM',
                         'AM 5' = 'AM')
unique(Idents(combined))

p <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_OG_am_lab.png", height = 5, width = 5)

# Renaming Neutrophils for real:
combined <- RenameIdents(combined, 
                         'N 1' = 'N B',
                         'N 2' = 'N A',
                         'N 3' = 'N C',
                         'N 4' = 'N C',
                         'N 5' = 'N A',
                         'N 6' = 'N B')
unique(Idents(combined))
p <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_OG_neutro_lab.png", height = 5, width = 5)

# B Cell No Subsetting  ---------------------------
# Want to see if we can get a Volcano Plot of DEG *markers* in B cells
# change y axis of graphs to p-val adj

b3.markers <- FindMarkers(combined, 
                          ident.1 = Cells(subset(combined, idents = 'B 3')),
                          ident.2 = Cells(subset(combined, idents = c('B 1', 
                                                                      'B 2', 
                                                                      'B 3'))),
                          logfc.threshold = 0, min.pct = 0)

b3.markers$cellType <- 'B 3'
b3.markers$diffexpressed <- "NO"
b3.markers$diffexpressed[b3.markers$avg_log2FC > 0.6 & b3.markers$p_val_adj < 0.05] <- "UP"
b3.markers$diffexpressed[b3.markers$avg_log2FC < -0.6 & b3.markers$p_val_adj < 0.05] <- "DOWN"
b3.markers <- cbind(gene_symbol = rownames(b3.markers), b3.markers)
b3.markers$delabel <- NA
b3.markers$delabel[b3.markers$diffexpressed != "NO"] <- b3.markers$gene_symbol[b3.markers$diffexpressed != "NO"]

p <- ggplot(data=b3.markers, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.4) +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  NoLegend() +
  labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(p-value)')) +
  theme(text = element_text(size=8, family = "sans"),
        axis.title=element_text(size=10, family = "sans", face="bold"))

ggsave('B3_markers.png', plot = p, width = 5, height = 5)
# B 3 appears inflammed

b2.markers <- FindMarkers(combined, 
                          ident.1 = Cells(subset(combined, idents = 'B 2')),
                          ident.2 = Cells(subset(combined, idents = c('B 1', 
                                                                      'B 2', 
                                                                      'B 3'))),
                          logfc.threshold = 0, min.pct = 0)

b2.markers$cellType <- 'B 2'
b2.markers$diffexpressed <- "NO"
b2.markers$diffexpressed[b2.markers$avg_log2FC > 0.6 & b2.markers$p_val_adj < 0.05] <- "UP"
b2.markers$diffexpressed[b2.markers$avg_log2FC < -0.6 & b2.markers$p_val_adj < 0.05] <- "DOWN"
b2.markers <- cbind(gene_symbol = rownames(b2.markers), b2.markers)
b2.markers$delabel <- NA
b2.markers$delabel[b2.markers$diffexpressed != "NO"] <- b2.markers$gene_symbol[b2.markers$diffexpressed != "NO"]

p <- ggplot(data=b2.markers, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.4) +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  NoLegend() +
  labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(p-value)')) +
  theme(text = element_text(size=8, family = "sans"),
        axis.title=element_text(size=10, family = "sans", face="bold"))

ggsave('b2_markers.png', plot = p, width = 5, height = 5)
# B 2 is likely an Ig-Producing subset

b1.markers <- FindMarkers(combined, 
                          ident.1 = Cells(subset(combined, idents = 'B 1')),
                          ident.2 = Cells(subset(combined, idents = c('B 1', 
                                                                      'B 2', 
                                                                      'B 3'))),
                          logfc.threshold = 0, min.pct = 0)

b1.markers$cellType <- 'B 1'
b1.markers$diffexpressed <- "NO"
b1.markers$diffexpressed[b1.markers$avg_log2FC > 0.6 & b1.markers$p_val_adj < 0.05] <- "UP"
b1.markers$diffexpressed[b1.markers$avg_log2FC < -0.6 & b1.markers$p_val_adj < 0.05] <- "DOWN"
b1.markers <- cbind(gene_symbol = rownames(b1.markers), b1.markers)
b1.markers$delabel <- NA
b1.markers$delabel[b1.markers$diffexpressed != "NO"] <- b1.markers$gene_symbol[b1.markers$diffexpressed != "NO"]

p <- ggplot(data=b1.markers, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.4) +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  NoLegend() +
  labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(p-value)')) +
  theme(text = element_text(size=8, family = "sans"),
        axis.title=element_text(size=10, family = "sans", face="bold"))

ggsave('b1_markers.png', plot = p, width = 5, height = 5)
# B 1 is naive?

# Cleaning up some metadata:
combined$wsnn_res.0 <- NULL
combined$wsnn_res.0.1 <- NULL
combined$wsnn_res.0.2 <- NULL
combined$wsnn_res.0.3 <- NULL
combined$wsnn_res.0.4 <- NULL
combined$wsnn_res.0.5 <- NULL
combined$wsnn_res.0.6 <- NULL
combined$wsnn_res.0.7 <- NULL
combined$wsnn_res.0.8 <- NULL
combined$wsnn_res.0.9 <- NULL
combined$wsnn_res.1 <- NULL
combined$wsnn_res.1.1 <- NULL
combined$wsnn_res.1.2 <- NULL
combined$wsnn_res.1.3 <- NULL
combined$wsnn_res.1.4 <- NULL
combined$wsnn_res.1.5 <- NULL
combined$wsnn_res.1.6 <- NULL
combined$wsnn_res.1.7 <- NULL
combined$wsnn_res.1.8 <- NULL
combined$wsnn_res.1.9 <- NULL
combined$wsnn_res.2 <- NULL

saveRDS(combined, 'combined_04_25_2021.rds')
rm(b1.markers, b2.markers, b3.markers, p)

# T & NK Cell Subclustering  ---------------------------
Idents(combined) <- 'celltype'
tcells <- subset(combined, idents = c('CD4 T 1',
                                  'CD4 T 2',
                                  'CD4 T 3',
                                  'CD4 T 4',
                                  'CD8 T',
                                  'Treg',
                                  'Nuocyte',
                                  'NK'))

DefaultAssay(tcells) <- 'RNA'
tcells <- DietSeurat(tcells, assays = c("RNA", "ADT"))
dim(tcells)

# Confirm subsetted metadata
table(tcells@meta.data$seurat_clusters)

# split data to rerun workflow
seurat_list <- SplitObject(tcells, split.by = "orig.ident")

# Rerun analysis on subsetted data:
for (i in 1:length(x = seurat_list)){
  sample_id <- names(seurat_list)[[i]]
  print(sample_id)
  
  seurat_list[[i]] <- SCTransform(seurat_list[[i]],
                                  method = "glmGamPoi",
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = 3000, 
                                  verbose = TRUE)
}

# > tcells SCT Integration ---------------------------
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")

# to_integrate <- SelectIntegrationFeatures(object.list = seurat_list, 
#                                           nfeatures = 6000 , verbose = TRUE,
#                                           new.assay.name = "integratedSCT_")
rm(seurat_list, combined)

# IMAGE SAVED

tcells <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                    new.assay.name = 'integratedSCT_', verbose = T)
saveRDS(tcells, file = "tcells_sct_integrated.rds")

# > tcells Regress out S and G2M scores ---------------------------
# This will take some time. 
tcells <- ScaleData(tcells, vars.to.regress = c("S.Score", "G2M.Score"))

# > tcells ADT Integration ---------------------------
seurat_list <- anchors@object.list

for (i in 1:length(x = seurat_list)){
  DefaultAssay(seurat_list[[i]]) = "ADT"
  VariableFeatures(seurat_list[[i]]) <- rownames(seurat_list[[i]][["ADT"]]) #select all ADT as variable features
  # DSB ALREADY NORMALIZED
  # seurat_list[[i]] <- NormalizeData(seurat_list[[i]], normalization.method = 'CLR', margin = 2) #CLR normalization for ADT
}

# Finding integration anchors for ADT
anchorsADT <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30) # may need to alter dims
tcellsADT <- IntegrateData(anchorset = anchorsADT, dims = 1:30, new.assay.name = "integratedADT_")

#I tcells not sure if this is a good idea:
tcells[["integratedADT_"]] <- tcellsADT[["integratedADT_"]] 

rm(tcellsADT, anchorsADT, anchors, seurat_list)
saveRDS(tcells, "tcells_integrated.rds")
# tcells <- readRDS("tcells_integrated.rds")

# > tcells Select PCs ---------------------------
DefaultAssay(tcells) <- "integratedSCT_"
all.genes <- rownames(tcells)
tcells <- ScaleData(tcells, features = all.genes)
tcells <- RunPCA(tcells, verbose = T)
ElbowPlot(tcells, ndims = 50, reduction = 'pca') # will use 30 because SCT is more accurate 
ggsave("tcells_elbow_plot_sct_dsb120clean.png")

# aPCA on ADT
DefaultAssay(tcells) <- 'integratedADT_'
tcells <- ScaleData(tcells)
tcells <- RunPCA(tcells, reduction.name = 'apca')
ElbowPlot(tcells, ndims = 30, reduction = "apca") # true dimensionality ~13, may be low
ggsave("tcells_elbow_plot_dsb120clean.png")

#Idents(tcells) <- "old.ident"
tcells <- RunUMAP(tcells, reduction = 'pca', dims = 1:30, 
              assay = 'integratedSCT_', 
              reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
tcells <- RunUMAP(tcells, reduction = 'apca', dims = 1:11, 
              assay = 'integratedADT_', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

tcells <- FindNeighbors(tcells, reduction = 'pca', dims = 1:30, 
                    graph.name = 'sct.snn')
tcells <- FindClusters(tcells, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

tcells <- FindNeighbors(tcells, reduction = 'apca', dims = 1:11,
                    graph.name = 'adt.snn')
tcells <- FindClusters(tcells, graph.name = 'adt.snn', resolution = 1.0, verbose = T)

Idents(tcells) <- 'sct.snn_res.1'
p1 <- DimPlot(tcells, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
Idents(tcells) <- 'adt.snn_res.1'
p2 <- DimPlot(tcells, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_tcells.png", width = 10, height = 5)

# Using original labels from combined.rds
Idents(tcells) <- 'celltype'
p1 <- DimPlot(tcells, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(tcells, reduction = 'adt.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave("SCT_or_ADT_only_clustering_tcells_origlabel.png", width = 10, height = 5)

# > tcells WNN clustering ---------------------------
# On the integrated SCT and integrated ADT data
# Identify multimodal neighbors:

tcells <- FindMultiModalNeighbors(
  tcells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:11), modality.weight.name = "SCT.weight")

tcells <- RunUMAP(tcells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

resolution.range <- seq(from = 0, to = 2.0, by = 0.1)

# WNN - Find clusters using a range of resolutions
tcells <- FindClusters(object = tcells, graph.name = "wsnn",
                       reduction.name = "wnn.umap", algorithm = 3, 
                       resolution = resolution.range, verbose = T)

for (res in resolution.range){
  md <- paste("wsnn_res.", res, sep = "")
  Idents(tcells) <- md
  p <- DimPlot(tcells, label = T, repel = T, label.size = 3, reduction = 'wnn.umap') + NoLegend()
  ggsave(paste("tcells_dimplot_", "res_", res, ".png", sep = ""), plot = p,
         width = 5, height = 5)
}

# Plotting
library(clustree)
clustree(tcells, prefix = "wsnn_res.")
ggsave("clustree_output_tcells_WNN.pdf", width = 9.5, height = 11)

Idents(tcells) <- 'celltype'
DimPlot(tcells, reduction = 'wnn.umap', label = TRUE, repel = TRUE, 
        label.size = 2.5) + NoLegend()
ggsave('tcells_wnnumap_OG_lab.png', height = 5, width = 5)

tcells <- FindClusters(object = tcells, graph.name = "wsnn",
                       reduction.name = "wnn.umap", algorithm = 3, 
                       resolution = 0.8, verbose = T)

# > tcells Calculate DEG markers  ---------------------------
DefaultAssay(tcells) <- "SCT"
tcells.markers <- FindAllMarkers(tcells, only.pos = FALSE, 
                             min.pct = 0.1, logfc.threshold = 0.6, 
                             max.cells.per.ident = Inf)
all.markers <- tcells.markers %>% group_by(cluster)
write.csv(all.markers, file = "tcells_cluster_biomarkers_labeled_cite.csv", row.names = FALSE)

## Plotting top 5 marker genes on heatmap: 
top5 <- tcells.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p <- DoHeatmap(subset(tcells, downsample = 100),
               features = top5$gene, size = 3, angle = 30) + NoLegend()
ggsave("tcells_top5_markers_heatmap_labeled_cite.png", plot = p, width = 10, height = 6)

DefaultAssay(tcells) <- 'SCT'
FeaturePlot(tcells, features = 'Foxp3', reduction = 'wnn.umap')
FeaturePlot(tcells, features = 'Cd3e', reduction = 'wnn.umap')

DefaultAssay(tcells) <- 'ADT'
FeaturePlot(tcells, features = 'CD4-TotalA', reduction = 'wnn.umap',
            min.cutoff = 0, max.cutoff = 'q99')
FeaturePlot(tcells, features = 'CD8b-TotalA', reduction = 'wnn.umap',
            min.cutoff = 0, max.cutoff = 'q99')

# wnn.umap labels on sct umap
DimPlot(tcells, reduction = 'sct.umap', label = TRUE, repel = TRUE, 
        label.size = 2.5) + NoLegend()

FeaturePlot(tcells, features = 'doublet_score', reduction = 'wnn.umap')

# Rename
tcells <- RenameIdents(tcells, 
                       '0' = 'CD4 T', 
                       '1' = 'NK',
                       '2' = 'CD4 T', 
                       '3' = 'CD8 T',
                       '4' = 'Nuocyte',
                       '5' = 'CD4 T', # activated Tbc1d4 and Tnfsf8 as exhaustion signatures, Icos as activation
                       '6' = 'Treg',
                       '7' = 'CD4 T', # unsure!, but clusters with CD4 T on sct umap
                       '8' = 'NKT',
                       '9' = 'NKT', # T and NK markers
                       '10' = 'gd T',
                       '11' = 'Treg', #proliferating
                       '12' = 'CD4 T') # activated, but grouping

p <- DimPlot(tcells, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_tcells_cd4_grouped_lab.png", height = 5, width = 5)

rm(all.markers, p, p1, p2, tcells.markers, top5, all.genes, features)

saveRDS(tcells, 'tcells_cite.rds')

# > tcells apply to parent --------------------------- 

combined <- readRDS('combined_04_25_2021.rds')

Idents(combined) <- 'sub_cluster'
unique(Idents(combined))
combined <- SetIdent(combined, cells = Cells(tcells), Idents(tcells))
combined$sub_cluster <- Idents(combined)

DimPlot(combined, label = TRUE, repel = TRUE, label.size = 4,
        cells = Cells(tcells), reduction = 'wnn.umap') + NoLegend()
ggsave('tcells_subcluster_onUMAP.png', height = 5, width = 5)

DimPlot(combined, label = TRUE, repel = TRUE, label.size = 4,
        cells = Cells(tcells), reduction = 'sct.umap') + NoLegend()
ggsave('tcells_subcluster_onUMAP_SCT.png', height = 5, width = 5)

combined$celltype <- combined$sub_cluster
Idents(combined) <- 'celltype'
unique(Idents(combined))

# Dealing with cDC 2  ---------------------------
combined <- RenameIdents(combined, 'cDC 2' = 'NC Mono')
combined <- RenameIdents(combined, 'NC Mono ' = 'NC Mono') # error in naming earlier

# Removing doublets  ---------------------------
combined <- readRDS('combined_04_25_2021b.rds')

combined <- subset(combined, idents = 'DOUBLET', invert = TRUE)

p <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 3.5) + NoLegend()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_OG_lab.png", height = 5, width = 5)

saveRDS(combined, 'combined_04_25_2021c.rds')
# combined <- readRDS('combined_04_25_2021c.rds')

# Figuring out Neutrophil / Eosinophils  ---------------------------
combined  <- readRDS('combined_04_25_2021c.rds')

DefaultAssay(combined) <- 'SCT'
FeaturePlot(combined, features = c('Ptprc', 'Itgax', 'Ly6g',    
                                   'Siglecf', 'Ccr3', 'Prg2',
                                   'Epx', 'Il13', 'Ly6c1'), 
            reduction = 'wnn.umap', ncol = 3)
ggsave('some_eosinophil_genes.png', width = 12, height = 12)

DefaultAssay(combined) <- 'ADT'
FeaturePlot(combined, features = c('Ly6g.Ly6c-TotalA', 'Ly-6G-TotalA', 'CDLy6C-TotalA'),
            cols = c("lightgrey","darkgreen"),
            reduction = 'wnn.umap', ncol = 3, min.cutoff = 0, max.cutoff = 'q99')
ggsave('some_eosinophil_adt.png', width = 12, height = 4)

# Condensed Table of celltypes vs. trt as above  ---------------------------
combined$celltype <- Idents(combined)

# number of cells per treatment
t1 <- table(combined$orig.ident)
# number of cells per type per treatment
t2 <- table(combined$celltype, col.names = combined$orig.ident)
col_order <- c('naive', 'p4', 'mp4', 'p24', 'mp24')
t2 <- t2[ ,col_order]
# number of cells per cluster
t3 <- table(combined@active.ident)

write.csv(t1, file = "t1_condensed.csv", row.names = F)
write.csv(t2, file = "t2_condensed.csv", row.names = T)
write.csv(t3, file = "t3_condensed.csv", row.names = F)

# DEGs by Treatment vs. Naive  ---------------------------
combined  <- readRDS('combined_04_25_2021c.rds')
unique(Idents(combined))
unique(combined$celltype) 
# wierd stuff here... like cDC 2 which was condensed, will rename celltype
combined$celltype <- Idents(combined)

freq_table <- read.csv('t2_condensed.csv', row.names = 'X')

# Create metadata colum for celltype + treatment condition
combined$celltype.trt <- paste(Idents(combined), combined$orig.ident, sep = "_")
Idents(combined) <- "celltype.trt"

# DEGs will be calculated relative to naive
DefaultAssay(combined) <- "SCT" # you definitely dont want to do this on integratedSCT_
cellList <- unique(combined$celltype)
cellList <- levels(cellList)

# check freq_table for any celltypes with less than 3 cells in them...
# will not include these for DEG analysis
tmp_cellList <- cellList
for (i in 1:length(cellList)){
  if ( any(freq_table[cellList[[i]], ] < 3)){
    print(cellList[[i]])
    tmp_cellList <- tmp_cellList[tmp_cellList != cellList[i]]
  }
}
cellList <- tmp_cellList; rm(tmp_cellList)

trtList <- list("p4", "mp4", "p24", "mp24")

graphList <- vector(mode = "list", length = length(cellList)*length(trtList))
degList <- vector(mode = "list", length = length(cellList)*length(trtList))

# loop to automate and store DEG lists by celltype and trt vs. naive
for (i in 1:length(x = cellList)){
  cell <- cellList[[i]] # current celltype
  base <- paste(cell, "naive", sep = "_")
  
  for (j in 1:length(x = trtList)){
    trt <- trtList[[j]] # current treatment
    ident.1 <- paste(cell, trt, sep = "_") # current celltype + treatment ident
    print(ident.1)
    
      
      degList[[(i-1)*4+j]] <- (FindMarkers(combined, ident.1 = ident.1, ident.2 = base, logfc.threshold = 0))
      degList[[(i-1)*4+j]]$cellType <- ident.1
      
      # add a column to tell whether the genes are up or down regulated
      degList[[(i-1)*4+j]]$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
      degList[[(i-1)*4+j]]$diffexpressed[degList[[(i-1)*4+j]]$avg_log2FC > 0.6 & degList[[(i-1)*4+j]]$p_val < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      degList[[(i-1)*4+j]]$diffexpressed[degList[[(i-1)*4+j]]$avg_log2FC < -0.6 & degList[[(i-1)*4+j]]$p_val < 0.05] <- "DOWN"
      # adding gene symbols in a column for easy accesss
      degList[[(i-1)*4+j]] <- cbind(gene_symbol = rownames(degList[[(i-1)*4+j]]), degList[[(i-1)*4+j]])
      # defining whether a gene is sufficiently differentially expressed or not
      degList[[(i-1)*4+j]]$delabel <- NA
      degList[[(i-1)*4+j]]$delabel[degList[[(i-1)*4+j]]$diffexpressed != "NO"] <- degList[[(i-1)*4+j]]$gene_symbol[degList[[(i-1)*4+j]]$diffexpressed != "NO"]
      # plotting
      
      graphList[[(i-1)*4+j]] <- ggplot(data=degList[[(i-1)*4+j]], aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
        geom_point(size = 0.4) +
        theme_classic() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) + 
        NoLegend() +
        labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(adj. p-value)')) +
        theme(text = element_text(size=8, family = "sans"),
              axis.title=element_text(size=10, family = "sans", face="bold"))
  }
}

saveRDS(degList, file = "04292021.degList.vs.naive.rds")
saveRDS(graphList, file = "04292021.degList.vs.naive.rds")

# Saving plots as png:
for (i in 1:length(x = graphList)){
  celltype <- degList[[i]]$cellType[1]
  filename = paste(celltype, "_v naive", ".png", sep = "")
  filename <- chartr("/", " ", filename)
  
  ggsave(filename = filename, plot = graphList[[i]],
         width = 3.75, height = 3.75)
}

degList <- do.call(rbind, degList)
write.csv(degList, file = "2021-04-29 DEGs.csv")

# GSEA  ---------------------------
library(fgsea)
library(msigdbr)
# Edited by MCK based on analysis of Hurskainen et all. 2021 Nat. Comm.
# Using version 7.4 from the MSigDB

result_table <- read.csv('2021-04-29 DEGs.csv', row.names = 'X')

# test <- result_table$cellType
# result_table$cellType <- gsub("(.*)_.*","\\1",test)

# Collect the Hallmark, KEGG, GO, and Reactome datasets
# THESE ARE FOR HUMAN
# hallmarks <- fgsea::gmtPathways("./GeneLists/h.all.v7.4.symbols.gmt")
# kegg <- fgsea::gmtPathways("./GeneLists/c2.cp.kegg.v7.4.symbols.gmt")
# go <- fgsea::gmtPathways("./GeneLists/c5.go.bp.v7.4.symbols.gmt")
# reactome <- fgsea::gmtPathways("./GeneLists/c2.cp.reactome.v7.4.symbols.gmt")

# Define a function to run GSEA for a single cluster

# Msigdbr usage TEST:
msigdbr_collections()

hallmarks <- msigdbr(species = "Mus musculus", category = "H")
hallmarks <- split(x = hallmarks$gene_symbol, f = hallmarks$gs_name)
kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
kegg <- split(x = kegg$gene_symbol, f = kegg$gs_name)
go <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
go <- split(x = go$gene_symbol, f = go$gs_name)
reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
reactome <- split(x = reactome$gene_symbol, f = reactome$gs_name)

gene_sets <- c(hallmarks, kegg, go, reactome)

runGSEA <- function(cluster){
  print(cluster)
  results <- filter(result_table, cellType == cluster)
  
  # results <- filter(results, p_val <= 0.05 &
  #                     abs(avg_log2FC) >= log2(1.5)) # +/-1.5 FC cutoff
  results <- filter(results, pct.1 >= 0.05)
  
  results <- arrange(results, desc(avg_log2FC))
  print(dim(results)[1])
  
  cluster_genes <- results$avg_log2FC
  names(cluster_genes) <- results$gene
  
  gsea <- fgsea(pathways = gene_sets,
                stats = cluster_genes,
                minSize=15,
                maxSize=500,
                nproc = 3,
                eps = 0) # no eps cutoff so we have "true" P-val
  gsea$cluster <- cluster
  
  return(gsea)
}

cluster_list <- unique(result_table$cellType)
fgsea_results <- lapply(cluster_list, runGSEA)

saveRDS(cluster_list, 'cluster_list_pos_neg_PCT_EPS.rds')
saveRDS(fgsea_results, 'fgsea_results_pos_neg_PCT_EPS.rds')

fgsea_results <- do.call("rbind", fgsea_results)
fgsea_results <- as.data.frame(fgsea_results)
fgsea_results$leadingEdge <- as.character(fgsea_results$leadingEdge)
write.csv(fgsea_results, 'fgsea_results_05012021_PCT_EPS.csv')

fgsea_results <- filter(fgsea_results, padj < 0.05)
write.csv(fgsea_results, 'fgsea_results_05012021_PCT_EPS_FILT.csv')

# DEG / GSEA Visualization ---------------------------

result_table <- read.csv('2021-04-29 DEGs.csv', row.names = 'X')

# Make bar charts of overall number of DEGs: 
# (colored regions on volcano plots)

test <- result_table$cellType
result_table$overallCellType <- gsub("(.*)_.*","\\1",test)

cellList <- unique(result_table$overallCellType)

plot_list <- vector(mode = 'list', length = length(cellList))
  
for (i in 1:length(cellList)){
  df <- filter(result_table, overallCellType == cellList[i])

  groups <- unique(df$cellType)
  
  if (length(groups) < 4){print('STOP')}
  
  num_pos <- numeric(length(groups))
  num_neg <- numeric(length(groups))
  for (j in 1:length(groups)){
    sub_df <- filter(df, cellType == groups[j])
    num_pos[j] <- sum(sub_df$p_val < 0.05 & sub_df$avg_log2FC > log2(1.5) )
    num_neg[j] <- sum(sub_df$p_val  < 0.05 & sub_df$avg_log2FC < -log2(1.5) )
  }
  
  groups <- rep(groups,2)
  deg <- c('pos', 'pos', 'pos', 'pos', 'neg', 'neg', 'neg', 'neg')
  num <- c(num_pos, num_neg)
  
  groups <- gsub(".*_", "", groups)
  
  data <- data.frame(groups, deg, num )
  
  p <- ggplot(data, aes(fill=deg, y=num, x=groups)) + 
    geom_bar(position="stack", stat="identity") + 
    labs(title=paste(cellList[i])) + 
    theme(axis.title.x=element_blank(), 
          plot.title = element_text(hjust = 0.5)) + NoLegend()
    
  p$data$groups <- factor(x = p$data$groups, levels = groups[1:4])
  
  plot_list[[i]] <- p
  
  ggsave(paste(cellList[i] , '_totalDEG.png', sep = ''), plot = p, width = 4, height = 4)
}

saveRDS(plot_list, 'GSEA_plot_list.rds')

# Define function to write a frequency table with top 10 +/- GSEA categories
top_GSEA <- function(df, n = 10){
  # Pathways that I do not want to be represented in my pseudo-heuristic method
  uninteresting_paths <- c('KEGG_PARKINSONS_DISEASE',
                           'KEGG_HUNTINGTONS_DISEASE',
                           'KEGG_ALZHEIMERS_DISEASE',
                           'KEGG_GRAFT_VERSUS_HOST_DISEASE',
                           'KEGG_PRION_DISEASES',
                           'REACTOME_NERVOUS_SYSTEM_DEVELOPMENT',
                           'GO_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION',
                           'GO_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT',
                           'GO_NERVE_DEVELOPMENT',
                           'GO_NEGATIVE_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT',
                           'GO_NERVOUS_SYSTEM_PROCESS',
                           'GO_POSITIVE_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT')
  
  df <- filter(df, !(pathway %in% uninteresting_paths))
  
  pos <- filter(df, NES > 0)
  neg <- filter(df, NES < 0)
  
  pos <- as.data.frame(table(pos$pathway))
  neg <- as.data.frame(table(neg$pathway))
  
  pos <- pos[order(pos$Freq, decreasing = T), ]
  neg <- neg[order(neg$Freq, decreasing = T), ]
  
  colnames(pos) <- c('pos_category', 'pos_freq')
  colnames(neg) <- c('neg_category', 'neg_freq')
  
  pos_neg <- data.frame(pos[1:n, ], neg[1:n, ])
  return(pos_neg)
}

# Plotting Select GSEA Terms:
library(readxl)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)

all_filt_res <- read_excel('fgsea_results_05012021_PCT_EPS_FILT.xlsx')
term_frequencies <- table(all_filt_res$pathway)
write.csv(term_frequencies, 'fgsea_results_frequencies_all_numb.csv')

df_total <- top_GSEA(all_filt_res, n = 30)
write.csv(df_total, 'fgsea_results_30_top_pos_neg.csv')

all_filt_res$cellType <- gsub("(.*)_.*","\\1",all_filt_res$cluster)

# > T / NK cell metacluster ---------------------------
# T cell, Nuocyte, and NK / NKT cell metacluster overall DEG plot grid:
cowplot::plot_grid(plotlist = plot_list[c(2,4,5,6,7,8)])
ggsave('t_NK_metacluster_DEG.png', width = 9, height = 6)

# GSEA Plot
results <- filter(all_filt_res, cellType == 'NK' | cellType == 'gd T' | 
                    cellType == 'Nuocyte' | cellType ==  'Treg' | 
                    cellType == 'NKT' | cellType == 'CD4 T' |  
                    cellType == 'CD8 T')
t_nk_df <- top_GSEA(results, n = 10)
write.csv(t_nk_df, 'fgsea_T_NK_frequencies_PCT.csv', row.names = F)

format_GSEA <- function(df_results, top_df, n = 25){
  df_results <- filter(df_results, pathway %in% top_df$pos_category | 
                      pathway %in% top_df$neg_category)
  bad_str <- c('HALLMARK_', 'REACTOME_', 'GO_', 'KEGG_')
  df_results$pathway <- str_remove(df_results$pathway, paste(bad_str, collapse = '|'))
  df_results$pathway <- str_trunc(df_results$pathway, n)
  return(df_results)
}

results <- format_GSEA(results, t_nk_df)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('T_NK_GSEA_dot_plot.png', width = 13, height = 5)

# > LESS CELLS FOR NIH ---------------------------
results <- filter(all_filt_res, cellType == 'CD4 T' |  
                    cellType == 'CD8 T' | cellType == 'B 1' | cellType == 'B 2')
t_nk_df <- top_GSEA(results, n = 10)
write.csv(t_nk_df, 'fgsea_T_NK_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, t_nk_df)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('NIH_GSEA_dot_plot.png', width = 9.75, height = 4.2)


# > B cell metacluster ---------------------------
# Same for B cells
cowplot::plot_grid(plotlist = plot_list[c(23, 27)])
ggsave('B_metacluster_DEG.png', width = 6, height = 3)

# GSEA Plot
results <- filter(all_filt_res, cellType == 'B 1' | cellType == 'B 2')
b_df <- top_GSEA(results, n = 10)
write.csv(b_df, 'fgsea_B_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, b_df, n = 35)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('B_GSEA_dot_plot.png', width = 8, height = 5)

# > Endothelial cell metacluster ---------------------------
cellList

cowplot::plot_grid(plotlist = plot_list[c(13,15,21)], ncol = 3)
ggsave('endo_metacluster_DEG.png', width = 9, height = 3)

# GSEA Plot
results <- filter(all_filt_res, cellType == 'gCap' | cellType == 'aCap' |
                    cellType == 'Vein')
endo_df <- top_GSEA(results, n = 10)
write.csv(endo_df, 'fgsea_endo_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, endo_df, n = 40)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('Endo_GSEA_dot_plot.png', width = 10, height = 5)

# > Fibroblast cell metacluster ---------------------------
cellList

cowplot::plot_grid(plotlist = plot_list[c(16,17,18,19,22)], ncol = 5)
ggsave('fibro_metacluster_DEG.png', width = 15, height = 3)

# GSEA Plot
fib <- c('Myofib', 'Lipofib', 'Ebf1+ Fib', 'Alv Fib', 'Adv Fib')
results <- filter(all_filt_res, cellType %in% fib)
fib_df <- top_GSEA(results, n = 10)
write.csv(fib_df, 'fgsea_fib_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, fib_df, n = 40)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('Fib_GSEA_dot_plot.png', width = 12, height = 5)

# > Epithelial cell metacluster ---------------------------
cellList

cowplot::plot_grid(plotlist = plot_list[c(14,20)], ncol = 2)
ggsave('epithelial_metacluster_DEG.png', width = 6, height = 3)

# GSEA Plot
results <- filter(all_filt_res, cellType == 'AT 2' | cellType == 'AT 1')
at_df <- top_GSEA(results, n = 10)
write.csv(at_df, 'fgsea_AT_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, at_df, n = 40)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('Epithelial_GSEA_dot_plot.png', width = 9.5, height = 5)

# > Monocyte cell metacluster ---------------------------
cellList

cowplot::plot_grid(plotlist = plot_list[c(1,24,25,26,29,30)], ncol = 6)
ggsave('Monocyte_metacluster_DEG.png', width = 6*3, height = 3)

# GSEA Plot
mono <- c('NC Mono', 'C Mono', 'IM', 'cDC 1', 'pDC')
results <- filter(all_filt_res, cellType %in% mono)
mono_df <- top_GSEA(results, n = 10)
write.csv(mono_df, 'fgsea_mono_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, mono_df, n = 40)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('Mono_GSEA_dot_plot.png', width = 12.5, height = 5)

# > Neutro / Alv. Macro metacluster ---------------------------
cellList

cowplot::plot_grid(plotlist = plot_list[c(9,10,11,12,28)], ncol = 5)
ggsave('AM_neutro_metacluster_DEG.png', width = 3*5, height = 3)

# GSEA Plot
n_am <- cellList[c(9,10,11,12,28)]
results <- filter(all_filt_res, cellType %in% n_am)
n_am_df <- top_GSEA(results, n = 10)
write.csv(n_am_df, 'fgsea_neutro_AM_frequencies_PCT.csv', row.names = F)

results <- format_GSEA(results, n_am_df, n = 40)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('N_AM_GSEA_dot_plot.png', width = 12.5, height = 5)

# Important Violin and Feature Plots ----
# combined <- readRDS('combined_04_25_2021c.rds')
# combined <- readRDS('combined_NIH_06222021.rds')

combined@meta.data$orig.ident <-
  factor(x = combined@meta.data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))

p1 <- VlnPlot(combined, features = c('CD200-TotalA'), idents = c('AT 2'), 
              split.by = 'orig.ident', assay = 'ADT') + NoLegend()
p2 <- VlnPlot(combined, features = c('CD200r-TotalA'), idents = c('AM'), 
              split.by = 'orig.ident', assay = 'ADT')
p1 | p2
ggsave('CD200_ADT.png', width = 10, height = 5)

VlnPlot(combined, features = 'CD19-TotalA', idents = 'B 1', 
        split.by = 'orig.ident', assay  = 'ADT')

'Ccr2' %in% rownames(combined@assays[["SCT"]]@scale.data)

combined <- ScaleData(combined, features = 'Ccr2')
VlnPlot(combined, features = 'Ccr2', idents = c('C Mono', 'NC Mono', 'IM', 'cDC 1'), 
        split.by = 'orig.ident', assay  = 'SCT', slot = 'scale.data')

FeaturePlot(combined, features = c('Il1b', 'S100a8'), reduction = 'wnn.umap')

FeaturePlot(combined, features = c('Ccr2'), reduction = 'wnn.umap')

# Figure 1 Marker ADTs:
DefaultAssay(combined) <- 'ADT'
adt.to.plot <- c('CD45','ESAM', 'CD326', 'Ly-6G', 'CD19',  
                 'CD11b', 'CD170', 'TCRB', 'CD4', 'CD8b')
gg <- list()
for (i in 1:length(adt.to.plot)){
  adt <- paste(adt.to.plot[[i]], '-TotalA', sep = '')
  gg[[i]] <- FeaturePlot(combined, features = adt, reduction = 'wnn.umap',
                         cols = c("lightgrey","darkgreen"), min.cutoff=0, 
                         max.cutoff = "q99") + 
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank())
}
ggsave('Fig1_ADT.png', plot = cowplot::plot_grid(plotlist = gg, ncol = 5),
       height = 2*3, width = 5*3)


# Neutrophil stuff:
DefaultAssay(combined) <- 'SCT'
FeaturePlot(combined, features = 'Cxcr2', reduction = 'wnn.umap')

# For Bhawana BMES:
p <- DimPlot(combined, reduction = 'wnn.umap', group.by = 'orig.ident', shuffle = T)
p <- p + theme(title = element_blank(), legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank())
ggsave('BMES_Bhawana_Cb.png', plot = p, width = 4, height = 2.75)


###############################################################################
# Generating DEGs relative for mp vs. p treatments at 4 and 24 hours

#combined <- readRDS("02-05-2020 Clustering.V2.rds")

Idents(combined) <- "celltype.trt"
DefaultAssay(combined) <- "RNA"

# cellList <- list("Lipo Fibroblasts",
#                  "gCap Endothelial",
#                  "AT1",
#                  "Neutrophils",
#                  "CD103+ DC", 
#                  "Alveolar Macro")

cellList <- c("gCap Endothelial",
              "AT1")

degList <- vector(mode = "list", length = length(cellList)*2)
timepoints <- c("4", "24")
graphList <- vector(mode = "list", length = length(cellList)*2)

for (i in 1:length(x = cellList)){
  cell <- cellList[[i]]
  print(cell)
  
  for (j in 1:2){
    time <- timepoints[[j]]
    
    ident.1 <- paste(cell, "_mp", time, sep = "")
    ident.2 <- paste(cell, "_p", time, sep = "")
    
    if (ident.1 == "Lymphatic Fibroblasts_mp24" | ident.2 == "Myeloid Cells_p24"){
      next
    } else {
      
      degList[[ (i-1)*2 + j ]] <- FindMarkers(combined, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0)
      degList[[ (i-1)*2 + j ]]$cellType <- paste(cell, time, "hr")
      
      # add a column to tell whether the genes are up or down regulated
      degList[[(i-1)*2 + j]]$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
      degList[[(i-1)*2 + j]]$diffexpressed[degList[[(i-1)*2 + j]]$avg_log2FC > 0.6 & degList[[(i-1)*2 + j]]$p_val < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      degList[[(i-1)*2 + j]]$diffexpressed[degList[[(i-1)*2 + j]]$avg_log2FC < -0.6 & degList[[(i-1)*2 + j]]$p_val < 0.05] <- "DOWN"
      # adding gene symbols in a column for easy accesss
      degList[[(i-1)*2 + j]] <- cbind(gene_symbol = rownames(degList[[(i-1)*2 + j]]), degList[[(i-1)*2 + j]])
      # defining whether a gene is sufficiently differentially expressed or not
      degList[[(i-1)*2 + j]]$delabel <- NA
      degList[[(i-1)*2 + j]]$delabel[degList[[(i-1)*2 + j]]$diffexpressed != "NO"] <- degList[[(i-1)*2 + j]]$gene_symbol[degList[[(i-1)*2 + j]]$diffexpressed != "NO"]
      
      # dealing with bad colors - DOES THIS WORK???????????
      if ( !any(degList[[(i-1)*2 + j]]$diffexpressed == "DOWN") ){
        colorValues <- c("black", "red")
      } else {
        colorValues <- c("blue", "black", "red")
      }
      
      # plotting
      graphList[[(i-1)*2 + j]] <- ggplot(data=degList[[(i-1)*2 + j]], aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
        geom_point(size = 0.4) +
        theme_classic() +
        geom_text_repel() +
        scale_color_manual(values=colorValues) + 
        NoLegend() +
        labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(p-value)')) +
        theme(text = element_text(size=8, family = "sans"),
              axis.title=element_text(size=10, family = "sans", face="bold"))
    } 
  }
}

saveRDS(degList, "2021-02-12 big DEG List mp vs p.rds")
saveRDS(graphList, "2021-02-12 big Graph List mp vs p.rds")

#degList <- readRDS("2021-02-11 DEG List mp vs p.rds")
#graphList <- readRDS("2021-02-11 Graph List mp vs p.rds")

# Saving plots as png:
for (i in 1:length(x = graphList)){
  celltype <- degList[[i]]$cellType[1]
  filename <- paste(celltype, "_mp vs p", ".png", sep = "")
  filename <- chartr("/", " ", filename)
  
  ggsave(filename = filename, plot = graphList[[i]],
         width = 3.5, height = 3.5)
}

degList <- do.call(rbind, degList)
write.csv(degList, file = "2021-02-12 mp vs. p DEGs for Epithelial, Endo, Fibroblast, and Myeloid-Derived Cells.csv")

################################################################################
# Plotting interesting features
#combined <- readRDS("02-05-2020 Clustering.V2.rds")

combined@meta.data$orig.ident <-
  factor(x = combined@meta.data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))

# this should be done above
DefaultAssay(combined) <- "RNA"
# combined <- ScaleData(combined) # not enough memory

cellList <- c("Ebf1+ Fibroblasts",
              "Myofibroblasts",
              "Lipo Fibroblasts",
              "Lymphatic Fibroblasts",
              "aCap Endothelial",
              "Vein Endothelial",
              "gCap Endothelial",
              "AT1",
              "AT2",
              "Myeloid Cells",
              "Neutrophils",
              "Eosinophils",
              "Alveolar Macro",
              "Interstitial Macro",
              "CD103+ DC",
              "Classical Mono",
              "Non-classical Mono")

# looking for IL17 signaling
il17.genes <- c("Il17a", 
                "Il17b",
                "Il17c",
                "Il17d",
                "Il25",
                "Il17f",
                "Il23a") #IL-25 is an alias for IL-17e

combined <- ScaleData(combined, features = il17.genes, assay = "RNA")

DotPlot(combined, features = il17.genes, idents = cellList.trt) # nothing much

#### looking for NF-kB genes

nfkb.genes <- c("Nfkb1", "Nfkb2", "Rela", "Nfkbia", "Ikbkb", "Cebpb", 
                "S100a7a", "S100a8", "S100a9", "Lcn2")

combined <- ScaleData(combined, features = nfkb.genes, assay = "RNA")

# seems like a lot for naive mice w/ Nuocytes
VlnPlot(combined, features = c("Nfkb1"), idents = c("Nuocytes"), split.by = "orig.ident")

Idents(combined) <- "celltype"
allCellTypes <- as.character(unique(Idents(combined)))
allCellTypes.trt <- vector(mode = "character", length = length(allCellTypes)*5)

Idents(combined) <- "celltype.trt"

cellList.trt <- vector(mode = "character", length = length(cellList)*5)
trtList <- c("naive", "p4", "mp4", "p24", "mp24")

for (i in 1:length(cellList)){
  for (j in 1:length(trtList)){
    cellList.trt[(i-1)*5 + j] <- paste(cellList[[i]], trtList[[j]], sep = "_")
  }
}

for (i in 1:length(allCellTypes)){
  for (j in 1:length(trtList)){
    allCellTypes.trt[(i-1)*5 + j] <- paste(allCellTypes[[i]], trtList[[j]], sep = "_")
  }
}

combined@active.ident <- factor(combined@active.ident, levels = allCellTypes.trt)

DotPlot(combined, features = c("Nfkb1", "Nfkb2", "Rela", "Nfkbia", "Ikbkb", "Cebpb"), idents = cellList.trt[31:35])
DotPlot(combined, features = c("S100a7a", "S100a8", "S100a9", "Lcn2"), idents = cellList.trt[31:35])
VlnPlot(combined, features = c("S100a8"), idents = cellList.trt[31:35])

# Caspases: Genes for 5, 10, 11 not present in mouse, but casp11 protein is coded by casp4 gene (ugh)

casp.genes <- c("Casp1", "Casp2", "Casp3", "Casp4", "Casp5", 
                "Casp6", "Casp7", "Casp8", "Casp9", "Casp10", 
                "Casp11", "Casp12", "Casp14", 
                "Apaf1", "Cycs", "Gsdmd", "Trpc1", "Fas", "Fadd", "Ifih1", "Calr")
combined <- ScaleData(combined, features = casp.genes)

# all cell types I want
p <- DotPlot(combined, 
             features = casp.genes, idents = cellList.trt, cols = c("blue", "red")) + 
  theme(text = element_text(size=8, family = "sans"), axis.text = element_text(size=8, family = "sans"))
ggsave(filename = "caspGenesDotPlot.tiff", plot = p, width = 8.5, height = 14)

# plotting only relevant casp and related genes
relevant.casp.genes <- casp.genes <- c("Casp1", "Casp3", "Casp4", "Casp5", "Gsdmd",
                                       "Casp6", "Casp7", "Casp8", "Fas", "Ifih1", "Calr",
                                       "Casp9", "Apaf1", "Cycs",  "Trpc1")
p <- DotPlot(combined, 
             features = casp.genes, idents = cellList.trt, cols = c("blue", "red")) + 
  theme(text = element_text(size=8, family = "sans"), axis.text = element_text(size=8, family = "sans"))
ggsave(filename = "caspGenesDotPlot.onlyExpressed.tiff", plot = p, width = 8.5, height = 14)

#### looking for IRFs and related genes
irf.genes <- c("Irf1", "Irf2", "Irf3", "Irf4", "Irf5", "Irf6", "Irf7", "Irf8", "Irf9")
combined <- ScaleData(combined, features = irf.genes)

p <- DotPlot(combined, 
             features = irf.genes, idents = cellList.trt, cols = c("blue", "red")) + 
  theme(text = element_text(size=8, family = "sans"), axis.text = element_text(size=8, family = "sans"))
ggsave(filename = "irfGenes.tiff", plot = p, width = 8.5, height = 14)

#### misc genes

# cxcl6 = cxcl5 in mice
cxcl.genes <- c("Cxcl1", "Cxcl2", "Cxcl3", "Pf4", "Cxcl5", "Cxcl9",
                "Cxcl10", "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl15",
                "Cxcr1","Cxcr2", "Cxcr3", "Cxcr4")
combined <- ScaleData(combined, cxcl.genes)
p <- DotPlot(combined, 
             features = cxcl.genes, idents = cellList.trt, cols = c("blue", "red")) + 
  theme(text = element_text(size=8, family = "sans"), axis.text = element_text(size=8, family = "sans"))
ggsave(filename = "cxclGenes.tiff", plot = p, width = 8.5, height = 14)


# S100 and other inflammatory genes
inflam.genes <- c("S100a8", "S100a9", "Lcn2", "Mt1", "Mt2", "Ngp", "Saa3")
combined <- ScaleData(combined, inflam.genes)
p <- DotPlot(combined, 
             features = inflam.genes, idents = cellList.trt, cols = c("blue", "red")) + 
  theme(text = element_text(size=8, family = "sans"), axis.text = element_text(size=8, family = "sans"))
ggsave(filename = "inflammatoryGenes.tiff", plot = p, width = 8.5, height = 14)

# Endothelial genes/adt for proliferation from Niethamer et al.
endoCells <- c("gCap Endothelial", "aCap Endothelial", "Vein Endothelial")

endoCells.trt <- vector(mode = "character", length = length(cellList)*5)
trtList <- c("naive", "p4", "mp4", "p24", "mp24")

for (i in 1:length(endoCells)){
  for (j in 1:length(trtList)){
    endoCells.trt[(i-1)*5 + j] <- paste(endoCells[[i]], trtList[[j]], sep = "_")
  }
}

Idents(combined) <- "celltype.trt"

VlnPlot(combined, features = c("CD34-TotalA", "CD309-TotalA"), 
        assay = "integratedADT_", 
        idents = endoCells.trt)

endoRNA <- c("Cd34", "Kdr", "Car4", "Mki67")
combined <- ScaleData(combined, features = endoRNA)

VlnPlot(combined, features = endoRNA, 
        assay = "RNA", 
        idents = endoCells.trt)

DotPlot(combined, assay = "RNA", features = endoRNA, idents = cellList.trt)


###############################################################################
# Trying to apply the same type of DEG analysis to DE_ADT
# who knows if this will work...

# also going to make this more clean

Idents(combined) <- "celltype.trt"
DefaultAssay(combined) <- "integratedADT_"

DefaultAssay(combined) <- "RNA"

trtList <- c("naive", "p4", "naive", "mp4", "naive", "p24", "naive", "mp24", "p4", "mp4", "p24", "mp24")

degList <- vector(mode = "list", length = length(cellList)*length(trtList)/2) # 6 comparisons for each trt and p vs. mp

for (i in 1:length(x = cellList)){
  cell <- cellList[[i]] # current celltype
  print(cell)
  
  for (j in seq(1, length(trtList), 2)){
    ident.1 <- paste(cell, trtList[[j+1]], sep = "_")
    ident.2 <- paste(cell, trtList[[j]], sep = "_")
    
    if(ident.1 == "Lymphatic Fibroblasts_p24" | ident.1 == "Lymphatic Fibroblasts_mp24" | ident.1 == "Myeloid Cells_p24" | ident.2 == "Myeloid Cells_p24"){
      next # if I knew how to catch exceptions... I would try and fix this better
    } else {
      
      idx <- (i-1)*6 + (j+1)/2 # index for degList
      
      degList[[idx]] <- (FindMarkers(combined, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0))
      
      comparison <- paste(ident.1, ident.2, sep = " vs. ")
      degList[[idx]]$cellType <- paste(comparison, DefaultAssay(combined))
      
      # add a column to tell whether the genes are up or down regulated
      degList[[idx]]$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
      degList[[idx]]$diffexpressed[degList[[idx]]$avg_log2FC > 0.6 & degList[[idx]]$p_val < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      degList[[idx]]$diffexpressed[degList[[idx]]$avg_log2FC < -0.6 & degList[[idx]]$p_val < 0.05] <- "DOWN"
      # adding gene symbols in a column for easy accesss
      degList[[idx]] <- cbind(gene_symbol = rownames(degList[[idx]]), degList[[idx]])
      # defining whether a gene is sufficiently differentially expressed or not
      degList[[idx]]$delabel <- NA
      degList[[idx]]$delabel[degList[[idx]]$diffexpressed != "NO"] <- degList[[idx]]$gene_symbol[degList[[idx]]$diffexpressed != "NO"]
      
    }
  }
}

# saving
saveRDS(degList, file = "degListForRNA-Data.AllComparisons.rds")

# plotting

graphList <- vector(mode = "list", length = length(cellList)*length(trtList)/2)

for (i in 1:length(graphList)){
  if (length(degList[[i]])==0){next}
  
  graphList[[i]] <- ggplot(data=degList[[i]], aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
    geom_point(size = 0.4) +
    theme_classic() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) + 
    NoLegend() +
    labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(p-value)')) +
    theme(text = element_text(size=8, family = "sans"),
          axis.title=element_text(size=10, family = "sans", face="bold"))
  
}

# Saving plots as png:
for (i in 1:length(x = graphList)){
  filename <- paste(degList[[i]]$cellType[1], ".png", sep = "")
  filename <- chartr("/", " ", filename)
  ggsave(filename = filename, plot = graphList[[i]],
         width = 3.5, height = 3.5)
}

degList <- do.call(rbind, degList)
write.csv(degList, file = "2021-02-23 DEGs for Epithelial, Endo, Fibroblast, and Myeloid-Derived Cells.csv")

###############################################################################
# Exporting information for CellPhoneDB:
combined <- readRDS("02-05-2020 Clustering.V2.rds")

Idents(combined) <- "orig.ident"

# Trying for naive first

# Take raw data and normalize it:
# There may be a better way to normalize, but can follow the Efremova et al. 2020 Nature Protocol

#naive <- subset(combined, idents = "naive")
#p4 <- subset(combined, idents = "p4")
mp4 <- subset(combined, idents = "mp4")

rm(combined)

library(tibble)
library(biomaRt)
library(tidyr)
library(dplyr)

naive <- mp4 #FAKE

count_raw <- naive[["RNA"]]@counts

count_norm <- NormalizeData(object = count_raw, normalization.method = "LogNormalize",
                            scale.factor = 10000, verbose = TRUE)
#count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000) #from Efremova et al.

allgenes <- rownames(naive)
count_norm <- as.data.frame(count_norm)
#unsure if this is good or bad due to null dist calculations in CellPhoneDB:
count_norm <- count_norm[rowSums(count_norm[,2:dim(count_norm)[2]])!=0,]

# CellPhoneDB can only use human genes, so need to convert mouse to human orthologs
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
# see values = if there is a problem with the next line
genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(naive@assays$RNA@data), mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
head(genesV2)

#count_norm <- count_norm[match(genesV2$MGI.symbol, rownames(naive), nomatch=F),] # old way errors
count_norm <- count_norm[match(genesV2$MGI.symbol, rownames(count_norm)),]
count_norm$gene <- genesV2$Gene.stable.ID

# remove NA rows... this might be bad
to_keep <- !is.na(count_norm[,1])
count_norm <- count_norm[to_keep, ]

# iterating through genes, non-unique expression values will be summed
# alternative: max or bitwiseor

unique_genes <- unique(count_norm$gene)
num_col <- length(count_norm[1,])

# ISSUE: NEED TO DEAL WITH NA ROWS! Not sure where they come from.

for (i in 1:length(unique_genes)){
  g <- unique_genes[i] # target gene
  
  if (sum(count_norm$gene == g) == 1){next}
  print(i)
  
  rs <- count_norm[count_norm$gene == g, ] # rows of interest
  # rs <- rs[!is.na(rs[,1]), ] # removing NA terms..
  gSum <- colSums(rs[!is.na(rs[,1]), ][, 1:(num_col-1) ] ) # NEED TO FIX bc slow
  
  idx <- match(g, count_norm$gene) # will be the only index for this gene, its first time in the matrix
  count_norm[idx, 1:(num_col-1)] <- gSum
  
  to_remove <- count_norm$gene == g
  to_remove[idx] <- FALSE # saving the one row
  
  count_norm <- count_norm[!to_remove, ]
}

saveRDS(count_norm, file = "count_mp4_maybeSketchy.rds")

#reordering column names:
col_n <- c(num_col, 1:num_col-1)
count_norm <- count_norm[, col_n]

write.table(count_norm, 'mp4_count.txt', sep='\t', row.names = F)

# Generating metadata file.
meta_data <- cbind(rownames(naive@meta.data), naive@meta.data[,'celltype', drop=F])
write.table(meta_data, 'mp4_meta.txt', sep='\t', quote=F, row.names=F)

################## Debugging p4
count_norm <- read.table(file = 'mp4_count.txt', sep = "\t")
cnames <- count_norm[1,]
meta_data <- read.table(file = 'mp4_meta.txt', sep = "\t")
rnames <- meta_data[,1]

#write.table(count_norm, 'p4_count_noQuote.txt', sep='\t', quote = F, row.names = F, col.names = F)
#write.table(meta_data, file = 'p4_meta_noQuote.txt', quote=F, sep = "\t", row.names = F, col.names = F)

# cnames and rnames are of equal length, spotchecks verified

len <- dim(count_norm)[1]

colnames(count_norm) <- cnames
count_norm <- count_norm[2:len, ]

numgenes <- dim(count_norm)[1]

# count_norm_small <- count_norm[, 1:51] # testing 50 cells
# count_norm_small <- count_norm[1:round(length(count_norm[1,])/2), ] # top half
#topLine <- colnames(count_norm)
#count_norm_small <- count_norm[round(length(count_norm[1,])/2)+1:length(count_norm[1,]), ] # bottom half
#count_norm_small <- count_norm[1:round(length(count_norm[1,])/8), ] # top eighth 

#count_norm_small <- count_norm[1:floor(numgenes/16), ] # top 16th
#count_norm_small <- count_norm[1:floor(numgenes/32), ] # top 32nd
#count_norm_small <- count_norm[floor(numgenes/16)+1:floor(numgenes/8), ] # second 16th 
#count_norm_small <- count_norm[floor(numgenes/32)+1:floor(numgenes/32)+1, ] # second 32nd 

#count_norm_small <- rbind(topLine, count_norm_small) # required to have cell names
# saveRDS(count_norm_small, file = "count_p4_topHalf_maybeSketchy.rds")
write.table(count_norm, 'mp4_count_noQuote.txt', sep='\t', quote = F, row.names = F, col.names = F)

# meta_data_small <- meta_data[1:51, ]
meta_data_small <- meta_data
# meta_data_small <- meta_data[1:round(length(count_norm[1,])/2), ]
write.table(meta_data_small, file = 'p4_meta_noQuote.txt', quote=F, sep = "\t", row.names = F, col.names = F)

###############################################################################
# Cell cycle scoring for all clusters:


