# NOTES ---------------------------
# CITE-Seq Analysis of Lung Cells w/ PUUC and MPLA+PUUC at 4 and 24 hour
# 1. Improved Scrublet to remove putative doublets conservatively
# 2. Used DSB method of normalization for ADT data incorporating isotype ctrls
# 3. Applied QC filters based on DSB tutorials
# 4. This script is for SCTransform, clustering, and generating DEGs w/ Seurat
# 5. Run this script after 3-QC.R

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

# BiocManager::install("glmGamPoi")

plan()

# Load Data  ---------------------------
setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
# setwd("C:/Users/UPDATE/Desktop/COVID Lung CITE-Seq")

# Bring in QC-filtered cells, keep DSB-normalized data
combined <- readRDS('combined_qc.rds')
DefaultAssay(combined) <- 'ADT'
combined <- DietSeurat(combined, assays = 'ADT')

rna <- readRDS('seuratObj_rna.rds')
rna <- subset(rna, cells = Cells(combined))

combined[['RNA']] <- rna@assays[["RNA"]]; rm(rna)

seurat_list <- SplitObject(combined, split.by = 'orig.ident'); rm(combined)

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
  ggsave(paste(sample_id, "_QC_v4", ".png", sep = ""), plot = qc_plot)
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
# Read previously-generated g2m- and s-phase genes mapped onto Mouse gene names
g2m_genes <- read.csv('g2m_genes.csv')$x
s_genes <- read.csv('s_genes.csv')$x

for (i in 1:length(x = seurat_list)){ # I can move this into the above for loop...
  DefaultAssay(seurat_list[[i]]) <- 'SCT'
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
# save.image("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/pre_integration_workspace_07032021.RData")

combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                          features.to.integrate = to_integrate, verbose = T,
                          new.assay.name = "integratedSCT_")

# Regress out S and G2M scores ---------------------------
# This will take some time. 
DefaultAssay(combined) <- 'integratedSCT_'
combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"))

saveRDS(combined, 'combined_integrated_SCT_07032021.rds')

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
saveRDS(combined, "combined_integrated_SCT_DSB_07032021.rds")

# Separate RNA and ADT Clustering ---------------------------
# PCA on SCT
DefaultAssay(combined) <- "integratedSCT_"
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)
combined <- RunPCA(combined, verbose = T)
ElbowPlot(combined, ndims = 50, reduction = 'pca')
ggsave("elbow_plot_sct.png")

# aPCA on ADT
DefaultAssay(combined) <- 'integratedADT_'
combined <- ScaleData(combined)
combined <- RunPCA(combined, reduction.name = 'apca')
ElbowPlot(combined, ndims = 30, reduction = "apca") # true dimensionality ~11, may be low
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
ggsave("SCT_or_ADT_only_clustering_v2.png", width = 10, height = 5)

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

saveRDS(combined, 'combined_PCA_v2.rds')
combined@assays$RNA <- NULL
saveRDS(combined, 'combined_PCA_noRNA_v2.rds')

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
ggsave("clustree_output_wnn_v2.pdf", width = 9.5, height = 11*3)

saveRDS(combined, 'combined_wnn_v2.rds')
# combined <- readRDS('combined_wnn_v2.rds')

res <- 1.5
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = FALSE)

# Data visualization:
DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
ggsave(paste("dimplot_", "res_", res, "_used_v2.png", sep = ""), width = 5, height = 5)

FeaturePlot(combined, features = c('S.Score', 'G2M.Score'), reduction = 'wnn.umap', 
            min.cutoff = 0)
ggsave('cell_cycle_scores_v2.wnn.png', width = 10, height = 5)

# Run: 4.1-cell_markers.R

# Marker Genes for Unlabeled Clsuters ---------------------------
## Exporting top 10 marker genes:
DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(combined.markers, file = "gene_biomarkers_unlab.csv", row.names = FALSE)
write.csv(top10, file = "gene_biomarkers_unlab_top10_v2.csv", row.names = FALSE)

DefaultAssay(combined) <- "ADT"
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(combined.markers, file = "adt_biomarkers_unlab.csv", row.names = FALSE)
write.csv(top10, file = "adt_biomarkers_unlab_top10_v2.csv", row.names = FALSE)

# Tables:
unlab_freq <- table(combined@active.ident, col.names = combined$orig.ident)
col_order <- c('naive', 'p4', 'mp4', 'p24', 'mp24')
unlab_freq <- unlab_freq[ ,col_order]
write.csv(unlab_freq, file = "unlab_freq_v2.csv", row.names = T)

# Need to go deeper on T/NK cell subcluster: ------
resolution.range <- seq(from = 0, to = 3.0, by = 0.1)
combined <- FindClusters(object = combined, graph.name = "wsnn",
                         reduction.name = "wnn.umap", algorithm = 3, 
                         resolution = resolution.range, verbose = T)
Idents(combined) <- 'wsnn_res.1.5'

tcells <- subset(combined, idents = c(6,8,13,14,23,26,27,38,39))
DimPlot(tcells, reduction = 'wnn.umap', label = T) + NoLegend()
# Go deeper into clustering resolution to separate CD4 and CD8 cells based on ADT
Idents(tcells) <- "wsnn_res.3.3"
DimPlot(tcells, reduction = 'wnn.umap', label = T) + NoLegend()
ggsave('tcells_res3.3_wnnUMAP_v2.png', width = 5, height = 5)

## Exporting top 10 marker genes:
DefaultAssay(tcells) <- "SCT"
tcells.markers <- FindAllMarkers(tcells, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
write.csv(tcells.markers, file = "tcells_gene_biomarkers_res3.3_v2.csv", row.names = FALSE)

## Exporting top 10 marker genes:
DefaultAssay(tcells) <- "ADT"
tcells.markers <- FindAllMarkers(tcells, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
write.csv(tcells.markers, file = "tcells_adt_biomarkers_res3.3_v2.csv", row.names = FALSE)

DefaultAssay(tcells) <- 'SCT'
FeaturePlot(tcells, features = 'Cd3e', reduction = 'wnn.umap')
FeaturePlot(tcells, features = 'Ccr7', reduction = 'wnn.umap')

# Tables:
unlab_freq <- table(tcells@active.ident, col.names = tcells$orig.ident)
col_order <- c('naive', 'p4', 'mp4', 'p24', 'mp24')
unlab_freq <- unlab_freq[ ,col_order]
write.csv(unlab_freq, file = "unlab_freq_v2_tcells.csv", row.names = T)

tcells$celltype <- Idents(tcells)
tcells <- RenameIdents(tcells, 
                       '5' = 'CD4 T', # Ccr7+
                       '8' = 'NK',
                       '15' = 'CD8 T', # Not Ccr7
                       '17' = 'CD4 T', # Ccr7+
                       '24' = 'DOUBLET',
                       '26' = 'Nuocyte',
                       '30' = 'CD4 T',
                       '37' = 'Treg',
                       '39' = 'CD8 T', # Not Ccr7
                       '48' = 'NK',
                       '50' = 'gd T',
                       '54' = 'Prolif T') # proliferating T cells
plot <- DimPlot(tcells, reduction = 'wnn.umap', label = T) + NoLegend()
ggsave('tcells_wnnUMAP_labeled_v2.png', plot = plot, width = 5, height = 5)

# Select obvious doublets:
select.cells <- CellSelector(plot = plot)
# [1] "p24_GATTCTTGTTGTGGCC-1" "mp4_AAAGTCCGTGCACAAG-1" "mp4_AAGTTCGGTCGCAGTC-1" "mp4_CTTGATTTCTAGTTCT-1"
# [5] "mp4_GAGTGTTTCGAGATAA-1" "mp4_GCATGATCACATTCTT-1"

# Return tcells idents to parents and remove several doublets found:
Idents(combined) <- 'wsnn_res.1.5'
combined <- SetIdent(combined, cells = Cells(tcells), Idents(tcells))
combined <- subset(combined, cells = select.cells, invert = TRUE)

DimPlot(combined, reduction = 'wnn.umap', label = TRUE) + NoLegend()
# No doublets in file. Will remove ident. 
combined <- subset(combined, idents = 'DOUBLET', invert = TRUE)
DimPlot(combined, reduction = 'wnn.umap', label = TRUE) + NoLegend() # Gone :)

# Rename Idents ---------------------------
# Saving old labels:
combined[["old.ident"]] <- Idents(combined)

# N = Neutrophil, T = T cell, B = B cell, gCap = gCap Endothelial, 
# AM = Alveolar Macrophage, IM = Intersitial Macro, Vein = Vein Endo, etc.

# Renaming:
combined <- RenameIdents(object = combined, 
                         '0' = "N", 
                         '1' = "gCap", 
                         '2' = "N",
                         '3' = "AM",
                         '4' = "gCap",
                         '5' = "Immature B",
                         
                         '7' = "AT 2",
                         
                         '9' = "C Mono",
                         '10' = "N",
                         '11' = "N",
                         '12' = "AM",
                         '15' = "AT 1",
                         '16' = "NC Mono",
                         '17' = "N",
                         '18' = "aCap",
                         '19' = "Myofib",
                         '20' = "IM",
                         '21' = "CD209 DC",
                         '22' = "Ebf1 Fib",

                         '24' = "Lipofib",
                         '25' = "Mature B",

                         
                         '28' = "Vein",
                         '29' = "Lipofib",
                         '30' = "AM",
                         '31' = "N",
                         '32' = "Baso",
                         '33' = "pDC",
                         '34' = "AM",
                         '35' = "N",
                         '36' = "Adv Fib",
                         '37' = "Immature B",


                         '40' = "N",
                         '41' = "AT 2",
                         '42' = "Lymph",
                         '43' = "Cilliated")
combined[["celltype"]] <- Idents(combined)
# p4 <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
# p4

p <- DimPlot(combined, reduction = 'wnn.umap', label = TRUE, repel = TRUE, 
             label.size = 4) + NoLegend()
ggsave("umap_res1.5_lab.png", height = 6, width = 6)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("umap_res1.5_lab_clean_v2.png", height = 5, width = 5)

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
p <- dittoSeq::dittoBarPlot(object = combined, var = combined@active.ident, group.by = "orig.ident",
                  scale = c('percent'), color.panel = color_df$Color, 
                  x.labels.rotate = F, ylab = 'Fraction of Cells')
p$data$grouping <- factor(x = p$data$grouping, levels = c("naive", "p4", "mp4", "p24", "mp24"))
p + theme(axis.title.x=element_blank())
ggsave("celltype_distribution_by_trt_leg_v2.png", width = 6, height = 4.2)
p + NoLegend()
ggsave("celltype_distribution_by_trt_noLeg_v2.png", width = 4.2, height = 4.2)

# Clusters with Legend
p <- DimPlot(combined, reduction = 'wnn.umap', label = FALSE, repel = TRUE, pt.size = 0.1)+ 
  theme(text = element_text(size=8, family = "sans"),  
        axis.title=element_text(size=8,family = "sans", face="bold"))
ggsave("rawClusters_v2.png", plot = p, width = 5, height = 4)

# Dimplot Split by Treatment
p <- DimPlot(combined, reduction = 'wnn.umap', split.by = "orig.ident") + NoLegend()
# change order in plot: 
p$data$orig.ident <- factor(x = p$data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))
p
ggsave("dimplot_splitby_trt_v2.png", plot = p, width = 15, height = 5)

# Cluster by Group
p5 <- DimPlot(combined, reduction = 'wnn.umap', group.by = "orig.ident") + 
  ggtitle(NULL) + NoLegend() + 
  theme(text = element_text(size=8, family = "sans"),  
        axis.title=element_text(size=8,family = "sans", face="bold"))
ggsave("clusterByGroup_v2.png", plot = p5, width = 5, height = 5)

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

saveRDS(combined, 'combined_07072021.rds')
# combined <- readRDS('combined_07072021.rds')

# Bring back in RNA for needs:

rna <- readRDS('seuratObj_rna.rds')
rna <- subset(rna, cells = Cells(combined))

combined[['RNA']] <- rna@assays[["RNA"]]; rm(rna)

saveRDS(combined, 'combined_07072021_rna.rds')
# combined <- readRDS('combined_07072021_rna.rds') #7/7/2021

# Plot nonintegrated SCT UMAP
DefaultAssay(combined) <- "SCT"
VariableFeatures(combined) <- rownames(combined)
combined <- RunPCA(combined, assay = 'SCT', verbose = T, reduction.name = 'SCTPC', reduction.key = 'SCTPC_')
ElbowPlot(combined, ndims = 50, reduction = 'SCTPC')
ggsave("elbow_plot_sct_non_integrated.png")
combined <- RunUMAP(combined, reduction = 'SCTPC', dims = 1:40, assay = 'SCT', 
                    reduction.name = 'sct.umap.raw', reduction.key = 'sct.umap.raw_')
DimPlot(combined, reduction = 'sct.umap.raw', label = T, repel = T) + NoLegend()
ggsave('sct.umap.raw.combined.png', width = 5, height = 5)
DimPlot(combined, split.by = 'orig.ident', reduction = 'sct.umap.raw', label = F) + NoLegend()
ggsave('sct.umap.raw.combined_split.png', width = 20, height = 4)
DimPlot(combined, group.by = 'orig.ident', reduction = 'sct.umap.raw', label = F) + NoLegend()
ggsave('sct.umap.raw.combined_groupedt.png', width = 5, height = 5)

# Plot nonintegrated ADT UMAP
DefaultAssay(combined) <- "ADT"
VariableFeatures(combined) <- c(rownames(combined)[1:14],
                                rownames(combined)[18:120])
combined <- ScaleData(combined)
combined <- RunPCA(combined, assay = 'ADT', verbose = T, reduction.name = 'ADTPC', reduction.key = 'ADTPC_')
ElbowPlot(combined, ndims = 50, reduction = 'ADTPC')
ggsave("elbow_plot_adt_non_integrated.png")
combined <- RunUMAP(combined, reduction = 'ADTPC', dims = 1:12, assay = 'ADT', 
                    reduction.name = 'sct.umap.raw', reduction.key = 'sct.umap.raw_')
DimPlot(combined, reduction = 'sct.umap.raw', label = T, repel = T) + NoLegend()
ggsave('adt.umap.raw.combined.png', width = 5, height = 5)
DimPlot(combined, split.by = 'orig.ident', reduction = 'sct.umap.raw', label = F) + NoLegend()
ggsave('adt.umap.raw.combined_split.png', width = 20, height = 4)
DimPlot(combined, group.by = 'orig.ident', reduction = 'sct.umap.raw', label = F) + NoLegend()
ggsave('adt.umap.raw.combined_groupedt.png', width = 5, height = 5)

saveRDS(combined, 'combined_07122021.rds')
# combined <- readRDS('combined_07122021.rds')

# Neutrophil (N) Subclustering  ---------------------------
Idents(combined) <- 'celltype'
neutro <- subset(combined, idents = 'N'); rm(combined)

DefaultAssay(neutro) <- 'RNA'
neutro <- DietSeurat(neutro, assays = c("RNA", "ADT"), 
                     dimreducs = c('sct.umap','adt.umap','wnn.umap', 'sct.umap.raw'))
dim(neutro) #[1] 32285  6533

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

# > UMAP without Integration ---------------------------
neutro <- merge(x=seurat_list[[1]], y=seurat_list[2:5])
DefaultAssay(neutro) <- "SCT"
VariableFeatures(neutro) <- rownames(neutro)
neutro <- RunPCA(neutro, assay = 'SCT', verbose = T, reduction.name = 'neutroPC', reduction.key = 'neutroPC_')
ElbowPlot(neutro, ndims = 50, reduction = 'neutroPC')
ggsave("elbow_plot_neutro_non_integrated.png")
neutro <- RunUMAP(neutro, reduction = 'neutroPC', dims = 1:40, assay = 'SCT', 
                    reduction.name = 'sct.umap.neutro', reduction.key = 'sct.umap.neutro_')
DimPlot(neutro, reduction = 'sct.umap.neutro', label = T, repel = T) + NoLegend()
ggsave('neutro_unintegrated.png', width = 5, height = 5)

DimPlot(neutro, reduction = 'sct.umap.neutro', split.by = 'orig.ident')
ggsave('neutro_unintegrated_split.png', width = 5*5, height = 5)

neutro <- FindNeighbors(neutro, reduction = 'neutroPC', dims = 1:40, 
                        graph.name = 'neutro.snn.unint')
neutro <- FindClusters(neutro, graph.name = 'neutro.snn.unint', resolution = 1.1, verbose = T)
DimPlot(neutro, label = T, reduction = 'sct.umap.neutro') + NoLegend()
ggsave('neutro_unintegrated_res1.1.png', width = 5, height = 5)
DimPlot(neutro, )

DefaultAssay(neutro) <- "SCT"
neutro.markers <- FindAllMarkers(neutro, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.6, 
                                   max.cells.per.ident = Inf)
top10 <- neutro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(neutro.markers, file = "neutro_gene_biomarkers_subcluster.csv", row.names = FALSE)

# > Remove neutro doublets
# Cluster 8: (expresses fibrillation markers, Col4a1)
dim(neutro) # [1] 11286  6533
neutro <- subset(neutro, idents = 8, invert = T)
dim(neutro) # [1] 11286  6492

DimPlot(neutro, label = T, reduction = 'sct.umap.neutro') + NoLegend()
ggsave('neutro_unintegrated_res1.1_v2.png', width = 5, height = 5)

# > N SCT Integration ---------------------------
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")

neutro <- IntegrateData(anchorset = anchors, new.assay.name = 'integratedSCT_',
                        normalization.method = "SCT", verbose = T)

# > N Regress out S and G2M scores ---------------------------
# This will take some time.
DefaultAssay(neurto) <- 'integrated'
neutro <- ScaleData(neutro, vars.to.regress = c("S.Score", "G2M.Score"))

# > N Select PCs ---------------------------
DefaultAssay(neutro) <- "integratedSCT_"
all.genes <- rownames(neutro)
neutro <- ScaleData(neutro, features = all.genes)
neutro <- RunPCA(neutro, verbose = T)
ElbowPlot(neutro, ndims = 50, reduction = 'pca') # will use 40 because SCT is more accurate 
ggsave("neutro_elbow_plot_integrated.png")

neutro <- RunUMAP(neutro, reduction = 'pca', dims = 1:30, 
                  assay = 'integratedSCT_', 
                  reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
neutro <- FindNeighbors(neutro, reduction = 'pca', dims = 1:30, 
                          graph.name = 'sct.snn')
neutro <- FindClusters(neutro, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

Idents(neutro) <- 'sct.snn_res.1'
DimPlot(neutro, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5)
ggsave("SCT_clustering_int_neutro.png", width = 5, height = 5)

DefaultAssay(neutro) <- "SCT"
neutro.markers <- FindAllMarkers(neutro, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
top10 <- neutro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(neutro.markers, file = "neutro_gene_biomarkers_subcluster_int.csv", row.names = FALSE)

# Remove cluster 10 which is doublets: (Col4a1 gene again, Ace gene)
dim(neutro) # [1] 11286  6533
neutro <- subset(neutro, idents = 10, invert = T)
dim(neutro) # [1] 11286  6490

DimPlot(neutro, reduction = 'sct.umap', label = TRUE,
        repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_int_neutro_v2.png", width = 5, height = 5)

DimPlot(neutro, reduction = 'sct.umap', group.by = "wsnn_res.1.5", label = TRUE,
        repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_int_neutro_v2_wssn_ident.png", width = 5, height = 5)

saveRDS(neutro, 'neutro_integrated_07212021.rds')

# Structural Cell Subclustering  ---------------------------
# (Endothelial / Epithelial / Fibroblast)
# combined <- readRDS('combined_07122021.rds'); rm(neutro, anchors, neutro.markers, seurat_list, top10)
Idents(combined) <- 'celltype'
struct <- subset(combined, idents = c('AT 2', 'AT 1', 'Lipofib', 'aCap',
                                      'Ebf1 Fib', 'gCap', 'Vein', 'Lymph',
                                      'Adv Fib', 'Myofib'))
rm(combined)

DefaultAssay(struct) <- 'RNA'
struct <- DietSeurat(struct, assays = c("RNA", "ADT"), 
                     dimreducs = c('sct.umap', 'adt.umap', 'wnn.umap', 'sct.umap.raw'))
dim(struct) # [1] 32285  7935

# split data to rerun workflow
seurat_list <- SplitObject(struct, split.by = "orig.ident")

# Rerun SCT on subsetted data:
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
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")
struct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                        new.assay.name = 'integratedSCT_', verbose = T)
rm(anchors, seurat_list)

# > Struct Regress out S and G2M scores ---------------------------
# This will take some time. 
struct <- ScaleData(struct, vars.to.regress = c("S.Score", "G2M.Score"))

# > Struct Select PCs ---------------------------
DefaultAssay(struct) <- "integratedSCT_"
all.genes <- rownames(struct)
struct <- ScaleData(struct, features = all.genes)
struct <- RunPCA(struct, verbose = T)
ElbowPlot(struct, ndims = 50, reduction = 'pca') 
ggsave("struct_elbow_plot_sct_dsb120clean.png")

#Idents(struct) <- "old.ident"
struct <- RunUMAP(struct, reduction = 'pca', dims = 1:40, 
                  assay = 'integratedSCT_', 
                  reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')

struct <- FindNeighbors(struct, reduction = 'pca', dims = 1:40, 
                        graph.name = 'sct.snn')
struct <- FindClusters(struct, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

Idents(struct) <- 'sct.snn_res.1'
DimPlot(struct, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_struct_int.png", width = 5, height = 5)

# Using original labels from combined.rds
Idents(struct) <- 'celltype'
DimPlot(struct, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_struct_origlabel.png", width = 5, height = 5)

# > Struct Calculate markers  ---------------------------
Idents(struct) <- 'sct.snn_res.1'
DefaultAssay(struct) <- "SCT"
struct.markers <- FindAllMarkers(struct, only.pos = FALSE, 
                                 min.pct = 0.1, logfc.threshold = 0.6, 
                                 max.cells.per.ident = Inf)
all.markers <- struct.markers %>% group_by(cluster)
write.csv(all.markers, file = "struct_cluster_biomarkers_int.csv", row.names = FALSE)

# Remove doublets: Cluster 16 is Ptprc+ Csf3r+
struct <- subset(struct, idents = 16, invert = T)
dim(struct) # [1] 15947  7785

# > Struct Gillich paper markers:   ---------------------------

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
  DefaultAssay(struct) <- "SCT"
  
saveRDS(struct, 'struct_integrated_07212021.rds')
rm(struct, all.markers, struct.markers, all.genes, features, filename, i, marker, 
   markers.to.plot, sample_id)

# Alveolar Macro Subclustering  ---------------------------
# combined <- readRDS('combined_07122021.rds')
Idents(combined) <- 'celltype'
am <- subset(combined, idents = 'AM'); rm(combined)

DefaultAssay(am) <- 'RNA'
am <- DietSeurat(am, assays = c("RNA", "ADT"), 
                 dimreducs = c('sct.umap', 'adt.umap', 'wnn.umap', 'sct.umap.raw'))
dim(am) # [1] 32285  2689

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
am <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                    new.assay.name = 'integratedSCT_', verbose = T)

# > AM Regress out S and G2M scores ---------------------------
# Skipping this to see what happens
# am <- ScaleData(am, vars.to.regress = c("S.Score", "G2M.Score"))

# > AM Select PCs ---------------------------
DefaultAssay(am) <- "integratedSCT_"
all.genes <- rownames(am)
am <- ScaleData(am, features = all.genes)
am <- RunPCA(am, verbose = T)
ElbowPlot(am, ndims = 50, reduction = 'pca') # will use 30 because SCT is more accurate 
ggsave("am_elbow_plot_sct_dsb120clean.png")

#Idents(am) <- "old.ident"
am <- RunUMAP(am, reduction = 'pca', dims = 1:30, 
                  assay = 'integratedSCT_', 
                  reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')
am <- FindNeighbors(am, reduction = 'pca', dims = 1:30, 
                        graph.name = 'sct.snn')
am <- FindClusters(am, graph.name = 'sct.snn', resolution = 1.1, verbose = T)

Idents(am) <- 'sct.snn_res.1.1'
DimPlot(am, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_am_int.png", width = 5, height = 5)

FeaturePlot(am, features = c('S.Score', 'G2M.Score'))
ggsave("am_non_regressed_G2M_S.png", width = 10, height = 5)

# Remove doublets:
# 12 - Col4a1, Tmem100 expressing (structural cells)
# 13 - Ptprcap and Cd3e and Bcl2 epxressiong (lymphocyte doublet)
dim(am) # [1] 12868  2689
am <- subset(am, idents = c(12, 13), invert = T)
dim(am) # [1] 12868  2608

saveRDS(am, 'am_int_07212021.rds')
rm(all.markers, am, am.markers, all.genes, res, resolution.range, seurat_list,
   anchors, features)

# T, B, & NK Cell Subclustering  ---------------------------
# combined <- readRDS('combined_07122021.rds')
Idents(combined) <- 'celltype'
tcells <- subset(combined, 
                 idents = c('Immature B', 'Mature B', 'CD4 T', 'CD8 T', 'NK',
                            'Nuocyte', 'Treg', 'gd T', 'Prolif T'))
rm(combined)

DefaultAssay(tcells) <- 'RNA'
tcells <- DietSeurat(tcells, assays = c("RNA", "ADT"),
                     dimreducs = c('sct.umap', 'adt.umap', 'wnn.umap', 'sct.umap.raw'))
dim(tcells) # [1] 32285  6110

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
tcells <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                    new.assay.name = 'integratedSCT_', verbose = T)

# > tcells Regress out S and G2M scores ---------------------------
# This will take some time. 
tcells <- ScaleData(tcells, vars.to.regress = c("S.Score", "G2M.Score"))

# > tcells Select PCs ---------------------------
DefaultAssay(tcells) <- "integratedSCT_"
all.genes <- rownames(tcells)
tcells <- ScaleData(tcells, features = all.genes)
tcells <- RunPCA(tcells, verbose = T)
ElbowPlot(tcells, ndims = 50, reduction = 'pca') # will use 30 because SCT is more accurate 
ggsave("tcells_elbow_plot_sct_dsb120clean.png")

#Idents(tcells) <- "old.ident"
tcells <- RunUMAP(tcells, reduction = 'pca', dims = 1:30, 
              assay = 'integratedSCT_', 
              reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')

tcells <- FindNeighbors(tcells, reduction = 'pca', dims = 1:30, 
                    graph.name = 'sct.snn')
tcells <- FindClusters(tcells, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

Idents(tcells) <- 'sct.snn_res.1'
DimPlot(tcells, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_lymphocytes_int.png", width = 5, height = 5)

# Using original labels from combined.rds
Idents(tcells) <- 'celltype'
DimPlot(tcells, reduction = 'sct.umap', label = TRUE,
              repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_lymphocytes_origlabel.png", width = 5, height = 5)

FeaturePlot(tcells, features = c('Cd8a', 'Cd4', 'Foxp3', 'S.Score', 'G2M.Score'))

# > tcells Calculate markers  ---------------------------
Idents(tcells) <- 'sct.snn_res.1'
DefaultAssay(tcells) <- "SCT"
tcells.markers <- FindAllMarkers(tcells, only.pos = FALSE, 
                             min.pct = 0.1, logfc.threshold = 0.6, 
                             max.cells.per.ident = Inf)
all.markers <- tcells.markers %>% group_by(cluster)
write.csv(all.markers, file = "tcells_cluster_biomarkers_labeled_int.csv", row.names = FALSE)

# remove doublets: ident 16 - Csf1r, Csf3r, Csf3ra (MaMoDC markers)
dim(tcells) # [1] 13397  6110
tcells <- subset(tcells, idents = 16, invert = T)
dim(tcells) # [1] 13397  6087

saveRDS(tcells, 'lymphocytes_int_07212021.rds')
rm(tcells, all.markers, anchors, seurat_list, tcells.markers, all.genes, features)

# Macrophage/Monocyte/DC Subclustering  ---------------------------
# combined <- readRDS('combined_07122021.rds')
Idents(combined) <- 'celltype'
mono <- subset(combined, idents = c('C Mono', 'NC Mono', 'IM', 'CD209 DC', 'pDC'))
rm(combined)

DefaultAssay(mono) <- 'RNA'
mono <- DietSeurat(mono, assays = c("RNA", "ADT"),
                   dimreducs = c('sct.umap', 'adt.umap', 'wnn.umap', 'sct.umap.raw'))
dim(mono) # [1] 32285  2255

# split data to rerun workflow
seurat_list <- SplitObject(mono, split.by = "orig.ident")

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

# > mono SCT Integration ---------------------------
features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                      nfeatures = 3000 , verbose = TRUE)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                  anchor.features = features, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                  anchor.features = features,
                                  normalization.method = "SCT")
mono <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                      new.assay.name = 'integratedSCT_', verbose = T)

# > mono Regress out S and G2M scores ---------------------------
# This will take some time. 
mono <- ScaleData(mono, vars.to.regress = c("S.Score", "G2M.Score"))

# > mono Select PCs ---------------------------
DefaultAssay(mono) <- "integratedSCT_"
all.genes <- rownames(mono)
mono <- ScaleData(mono, features = all.genes)
mono <- RunPCA(mono, verbose = T)
ElbowPlot(mono, ndims = 50, reduction = 'pca') # will use 30 because SCT is more accurate 
ggsave("mono_elbow_plot_sct_dsb120clean.png")

mono <- RunUMAP(mono, reduction = 'pca', dims = 1:30, 
                assay = 'integratedSCT_', 
                reduction.name = 'sct.umap', reduction.key = 'sctUMAP_')

mono <- FindNeighbors(mono, reduction = 'pca', dims = 1:30, 
                      graph.name = 'sct.snn')
mono <- FindClusters(mono, graph.name = 'sct.snn', resolution = 1.0, verbose = T)

Idents(mono) <- 'sct.snn_res.1'
DimPlot(mono, reduction = 'sct.umap', label = TRUE,
        repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_mono_int.png", width = 5, height = 5)

# Using original labels from combined.rds
Idents(mono) <- 'celltype'
DimPlot(mono, reduction = 'sct.umap', label = TRUE,
        repel = TRUE, label.size = 2.5) + NoLegend()
ggsave("SCT_clustering_mono_origlabel.png", width = 5, height = 5)

DimPlot(mono, reduction = 'sct.umap', label = TRUE,
        repel = TRUE, label.size = 3, split.by = 'orig.ident') + NoLegend()
ggsave("SCT_clustering_mono_split.png", width = 25, height = 5)

# > mono Calculate markers  ---------------------------
Idents(mono) <- 'sct.snn_res.1'
DefaultAssay(mono) <- "SCT"
mono.markers <- FindAllMarkers(mono, only.pos = FALSE, 
                               min.pct = 0.1, logfc.threshold = 0.6, 
                               max.cells.per.ident = Inf)
all.markers <- mono.markers %>% group_by(cluster)
write.csv(all.markers, file = "mono_cluster_biomarkers_labeled_int.csv", row.names = FALSE)

# Look at ADT markers:
FeaturePlot(mono, features = c('CD170-TotalA'))
FeaturePlot(mono, features = c('CD103-TotalA', 'CD11c-TotalA', 
                               'CD11b-TotalA', 'CD24-TotalA'))

# remove doublets: idents
# none

saveRDS(mono, 'mono_int_07212022.rds')
rm(mono, all.markers, anchors, seurat_list, mono.markers, all.genes, features)

# combined <- readRDS('combined_07122021.rds')
Idents(combined) <- 'celltype'
combined <- subset(combined, idents = c('C Mono', 'NC Mono', 'IM', 'CD209 DC', 'pDC'))
DimPlot(combined, reduction = 'wnn.umap')
ggsave('mono_sub_v1_label_OLD.png', width = 5, height = 5)
Idents(combined) <- 'wsnn_res.4'
DimPlot(combined, reduction = 'wnn.umap')
DefaultAssay(combined) <- 'ADT'
ggsave('mono_sub_res4.png', width = 5, height = 5)
FeaturePlot(combined, features = 'CD103-TotalA', reduction = 'wnn.umap')
ggsave('mono_sub_CD103.png', width = 5, height = 5)
FeaturePlot(combined, features = 'CD11b-TotalA', reduction = 'wnn.umap')

DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = FALSE, 
                               min.pct = 0.1, logfc.threshold = 0.6, 
                               max.cells.per.ident = Inf)
all.markers <- combined.markers %>% group_by(cluster)
write.csv(all.markers, file = "combined_sub_mono_cluster_biomarkers_wnn.csv", row.names = FALSE)

dim(combined) # [1] 18010  2255
saveRDS(combined, 'mono_clust_by_resolution_07282021.rds')

# Apply all subclustering to parent: ---------------
combined <- readRDS('combined_07122021.rds') # main

# Generate a new column called sub_cluster in the metadata
combined$sub_cluster <- Idents(combined)
Idents(combined) <- 'sub_cluster'

# Neutro:
neutro <- readRDS('neutro_integrated_07212021.rds') # subs w/o doublets
orig <- colnames( subset(combined, idents = 'N') )
length(orig) # [1] 6533
dim(neutro) # [1] 11286  6490
dub <- setdiff(orig, colnames(neutro)) # doublets
length(dub) # [1] 43

dim(combined)
combined <- subset(combined, cells = dub, invert = TRUE)
dim(combined)




###
### DO TO ALL SUBCLUSTERS SIMILARLY ###
###

struct <- readRDS('struct_integrated_07212021.rds')
am <- readRDS('am_int_07212021.rds')
lympho <- readRDS('lymphocytes_int_07212021.rds')
mono <- readRDS('mono_clust_by_resolution_07282021.rds')








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

# Struct apply to parent (old) --------------------------- 

# struct <- readRDS('struct_04_23_2021.rds')


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


