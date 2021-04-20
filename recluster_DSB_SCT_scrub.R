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

# DEGs by Treatment  ---------------------------

# Create metadata colum for celltype + treatment condition
combined$celltype.trt <- paste(Idents(combined), combined$orig.ident, sep = "_")
Idents(combined) <- "celltype.trt"

# DEGs will be calculated relative to naive
DefaultAssay(combined) <- "SCT" # you definitely dont want to do this on integratedSCT_
cellList <- unique(combined$celltype)

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
    
    if(ident.1 == "Lymphatic Fibroblasts_p24" | ident.1 == "Lymphatic Fibroblasts_mp24" | ident.1 == "Myeloid Cells_p24"){
      next # if I knew how to catch exceptions... I would try and fix this better
    } else {
      
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
        labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(p-value)')) +
        theme(text = element_text(size=8, family = "sans"),
              axis.title=element_text(size=10, family = "sans", face="bold"))
    }
  }
}

saveRDS(degList, file = "test.degList.naive.rds")
saveRDS(graphList, file = "test.degList.naive.rds")

# Saving plots as png:
for (i in 1:length(x = graphList)){
  celltype <- degList[[i]]$cellType[1]
  filename = paste(celltype, "_v naive", ".png", sep = "")
  filename <- chartr("/", " ", filename)
  
  ggsave(filename = filename, plot = graphList[[i]],
         width = 3.5, height = 3.5)
}

degList <- do.call(rbind, degList)
write.csv(degList, file = "2021-04-19 DEGs.csv")

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


