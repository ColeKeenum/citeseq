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
#plan("multiprocess", workers = 12)

# Load Data  ---------------------------
setwd('~/Documents/citeseq-code')
# setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
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
                         '5' = "Immature B", # renamed later in code
                         
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
                         '25' = "Mature B", # renamed later in code

                         
                         '28' = "Vein",
                         '29' = "Lipofib",
                         '30' = "AM",
                         '31' = "N",
                         '32' = "Baso",
                         '33' = "pDC",
                         '34' = "AM",
                         '35' = "N",
                         '36' = "Adv Fib",
                         '37' = "Immature B", # renamed later in code


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
Idents(comDbined) <- 'celltype'
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

# Define function to apply removed doublets to combined seurat
remove_dub <- function(obj, sub_obj, ident_list){
  orig <- colnames( subset(obj, idents = ident_list) )
  print(length(orig))
  print(dim(sub_obj))
  dub <- setdiff(orig, colnames(sub_obj)) 
  print(length(dub))
  
  print(dim(obj))
  obj <- subset(obj, cells = dub, invert = TRUE)
  print(dim(obj))
  return(obj)
}

# Neutro:
neutro <- readRDS('neutro_integrated_07212021.rds') # subs w/o doublets
ident_list <- c('N')
combined <- remove_dub(obj = combined, sub_obj = neutro, ident_list = ident_list)
rm(neutro)

# Struct
struct <- readRDS('struct_integrated_07212021.rds')
ident_list <- c('AT 2', 'AT 1', 'Lipofib', 'aCap', 'Ebf1 Fib', 'gCap', 'Vein', 
                'Lymph', 'Adv Fib', 'Myofib')
combined <- remove_dub(combined, struct, ident_list)
rm(struct)

# Alveolar Macro
am <- readRDS('am_int_07212021.rds')
ident_list <- c('AM')
combined <- remove_dub(obj = combined, sub_obj = am, ident_list = ident_list)
rm(am)

# Lymphocytes
lympho <- readRDS('lymphocytes_int_07212021.rds')
ident_list <- c('Immature B', 'Mature B', 'CD4 T', 'CD8 T', 'NK', 'Nuocyte', 'Treg', 'gd T', 'Prolif T')
combined <- remove_dub(obj = combined, sub_obj = lympho, ident_list = ident_list)
rm(lympho)

# Macro / Mono / DC
mono <- readRDS('mono_clust_by_resolution_07282021.rds')
ident_list <- c('C Mono', 'NC Mono', 'IM', 'CD209 DC', 'pDC')
combined <- remove_dub(obj = combined, sub_obj = mono, ident_list = ident_list)
# Apply Macro / Mono / DC Subcluster Identities
Idents(mono) <- 'wsnn_res.4'
DimPlot(mono, reduction = 'wnn.umap')
mono <- RenameIdents(mono, '12' = 'C Mono', '41' = 'C Mono', '16' = 'NC Mono', 
                     '20' = 'IM', '33' = 'H2 Hi DC', '46' = 'CD103 DC', '47' = 'pDC')
p <- DimPlot(mono, reduction = 'wnn.umap'); p
dub <- CellSelector(p)
dub # [1] "p4_CAGCAGCGTGTGTGGA-1"  "p24_GGTAACTGTAGTACGG-1" "p24_TGAGTCATCGGCATAT-1"

rm(mono)

# Return ident to parent and remove these three doublets
combined <- SetIdent(combined, cells = Cells(mono), Idents(mono))
dim(combined) # [1]   120 25442
combined <- subset(combined, cells = dub, invert = TRUE)
dim(combined) # [1]   120 25439

# Figuring out Neutrophil / Eosinophils  ---------------------------
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

saveRDS(combined, 'combined_subclustered_07292021.rds')

# DEGs by Treatment vs. Naive  ---------------------------
# combined <- readRDS('combined_subclustered_07292021.rds')
unique(Idents(combined))
unique(combined$celltype) 
freq_table <- read.csv('t2_condensed.csv', row.names = 'X')

# Create metadata colum for celltype + treatment condition
combined$celltype.trt <- paste(Idents(combined), combined$orig.ident, sep = "_")
Idents(combined) <- "celltype.trt"

# DEGs will be calculated relative to naive
DefaultAssay(combined) <- "SCT" 
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

saveRDS(degList, file = "degList.vs.naive.07292021.rds")
saveRDS(graphList, file = "degList.vs.naive.07292021.rds")

# Saving plots as png:
for (i in 1:length(x = graphList)){
  celltype <- degList[[i]]$cellType[1]
  filename = paste(celltype, "_v naive_v2", ".png", sep = "")
  filename <- chartr("/", " ", filename)
  
  ggsave(filename = filename, plot = graphList[[i]],
         width = 3.75, height = 3.75)
}

degList <- do.call(rbind, degList)
write.csv(degList, file = "2021-07-29 DEGs.csv")
saveRDS(combined, 'combined_07292021_v2.rds')

# Rename B cells: ------
# combined <- readRDS('combined_07292021_v2.rds')
Idents(combined) <- 'celltype'
combined <- RenameIdents(combined, 'Immature B' = 'IgD Hi B',
                         'Mature B' = 'IgD Lo B')
genes <- FindMarkers(combined, ident.1 = 'IgD Hi B', ident.2 = 'IgD Lo B',
                     logfc.threshold = 0, max.cells.per.ident = Inf, only.pos = F,
                     min.pct = 0)
write.csv(genes, 'bcells_markers.csv')

format_DEG <- function(df){
  df$diffexpressed <- "NO"
  df$diffexpressed[df$avg_log2FC > 0.6 & df$p_val_adj < 0.05] <- "UP"
  df$diffexpressed[df$avg_log2FC < -0.6 & df$p_val_adj < 0.05] <- "DOWN"
  df <- cbind(gene_symbol = rownames(df), df)
  df$delabel <- NA
  df$delabel[df$diffexpressed != "NO"] <- df$gene_symbol[df$diffexpressed != "NO"]
  return(df)
} 
df <- format_DEG(genes)

p <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.4) +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  NoLegend() +
  labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(adj. p-value)')) +
  theme(text = element_text(size=8, family = "sans"),
        axis.title=element_text(size=10, family = "sans", face="bold"))
ggsave('B_subtypes_DEGs.png', width = 5, height = 5)

DimPlot(combined, reduction = 'wnn.umap', label = T) + NoLegend()

saveRDS(combined, 'combined_01012022.rds')

# GSEA  ---------------------------
library(fgsea)
# Edited by MCK based on analysis of Hurskainen et all. 2021 Nat. Comm.
# Using version 7.4 from the MSigDB

result_table <- read.csv('2021-07-29 DEGs.csv', row.names = 'X')

gene_set_maker <- function(){
  library(msigdbr)
  
  # msigdbr_collections()
  # packageVersion("msigdbr")
  # [1] ‘7.4.1’
  
  hallmarks <- msigdbr(species = "Mus musculus", category = "H")
  hallmarks <- split(x = hallmarks$gene_symbol, f = hallmarks$gs_name)
  kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
  kegg <- split(x = kegg$gene_symbol, f = kegg$gs_name)
  go <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
  go <- split(x = go$gene_symbol, f = go$gs_name)
  reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
  reactome <- split(x = reactome$gene_symbol, f = reactome$gs_name)
  
  gene_sets <- c(hallmarks, kegg, go, reactome)
  return(gene_sets)
}


# Define a function to run GSEA for a single cluster
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
                maxSize=Inf,
                nproc = 12,
                eps = 0) # no eps cutoff so we have "true" P-val
  gsea$cluster <- cluster
  
  return(gsea)
}

cluster_list <- unique(result_table$cellType)
fgsea_results <- lapply(cluster_list, runGSEA)

saveRDS(cluster_list, 'cluster_list_pos_neg_PCT_EPS_12292021.rds')
saveRDS(fgsea_results, 'fgsea_results_pos_neg_PCT_EPS_12292021.rds')

fgsea_results <- do.call("rbind", fgsea_results)
fgsea_results <- as.data.frame(fgsea_results)
fgsea_results$leadingEdge <- as.character(fgsea_results$leadingEdge)
write.csv(fgsea_results, 'fgsea_results_12292021_PCT_EPS.csv')

fgsea_results <- filter(fgsea_results, padj < 0.05)
write.csv(fgsea_results, 'fgsea_results_12292021_PCT_EPS_FILT.csv')

# DEG / GSEA Visualization ---------------------------
result_table <- read.csv('2021-07-29 DEGs.csv', row.names = 'X')

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

saveRDS(plot_list, 'GSEA_plot_list_08072021.rds')

# Define function to write a frequency table with top 10 +/- GSEA categories
top_GSEA <- function(df, n = 10){
  # Pathways that I do not want to be represented
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
                           'GO_POSITIVE_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT',
                           'GOBP_RESPONSE_TO_ARSENIC_CONTAINING_SUBSTANCE')
  
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

all_filt_res <- read.csv('fgsea_results_12292021_PCT_EPS_FILT.csv')
# term_frequencies <- table(all_filt_res$pathway)
# write.csv(term_frequencies, 'fgsea_results_frequencies_all_numb_08072021.csv')

df_total <- top_GSEA(all_filt_res, n = 30)
# write.csv(df_total, 'fgsea_results_30_top_pos_neg_08072021.csv')

all_filt_res$cellType <- gsub("(.*)_.*","\\1",all_filt_res$cluster)

format_GSEA <- function(df_results, top_df, n = 25){
  df_results <- filter(df_results, pathway %in% top_df$pos_category | 
                         pathway %in% top_df$neg_category)
  # bad_str <- c('HALLMARK_', 'REACTOME_', 'GO_', 'KEGG_', 'GOBP_')
  # df_results$pathway <- str_remove(df_results$pathway, paste(bad_str, collapse = '|'))
  df_results$pathway <- str_trunc(df_results$pathway, n)
  return(df_results)
}

# > T cell subsets ---------------------------
vec <- c('CD4 T', 'CD8 T', 'gd T', 'Treg')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 3)
}

df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)

results <- format_GSEA(results, df, n=30)

# Adjust the orders of the x-axis in the plot
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('T_GSEA_08072021.png', width = 7, height = 6)

# > B cell subsets ---------------------------
vec <- c('Mature B', 'Immature B')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 3)
}

df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)

results <- format_GSEA(results, df, n=30)

# Adjust the orders of the x-axis in the plot
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('B_GSEA_08072021.png', width = 7, height = 6)

# > Lymphocyte cell subsets ---------------------------
# I like this one!
vec <- c('CD4 T', 'CD8 T', 'Mature B', 'Immature B', 'NK')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 3)
}
df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)
results <- format_GSEA(results, df, n=30)
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}
ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('lympho_GSEA_08072021_v2.png', width = 11.5, height = 6)

# > Endothelial cell subsets ---------------------------
vec <- c('gCap', 'aCap', 'Vein')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 3)
}
df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)
results <- format_GSEA(results, df, n=40)
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}
ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('endo_GSEA_08072021.png', width = 9.5, height = 5)

# > Fibroblast subsets ---------------------------
vec <- c('Lipofib', 'Myofib', 'Ebf1 Fib', 'Adv Fib')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 3)
}
df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)
results <- format_GSEA(results, df, n=40)
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}
ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('fibro_GSEA_08072021.png', width = 10.75, height = 5)

# > Epithelial subsets ---------------------------
vec <- c('AT 1', 'AT 2')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 4)
}
df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)
results <- format_GSEA(results, df, n=40)
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}
ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('epithelial_GSEA_08072021.png', width = 8.75, height = 5)

# > MaMoDC + AM + Baso (Myeloid) cell metacluster ---------------------------
vec <- c('C Mono', 'NC Mono', 'AM', 'Baso')
res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(all_filt_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 4)
}
df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)
results <- format_GSEA(results, df, n=40)
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4", "p24", "mp24")
for (i in 1:length(vec)-1){
  for (j in 1:length(trtList)){
    orders[i*4+j] <- paste(vec[i+1], trtList[j], sep = '_')
  }
}
ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('myeloid_GSEA_08072021.png', width = 10.5, height = 5)

# Done 08072021

# > NK Cells Only ---------------------------
vec <- c('NK')

test_res <- filter(all_filt_res, cluster %in% c('NK_p4', 'NK_mp4'))

res_vec <- vector(mode = 'list', length = length(vec))
df_vec <- vector(mode = 'list', length = length(vec))
for (i in 1:length(vec)){
  print(vec[[i]])
  res_vec[[i]] <- filter(test_res, cellType %in% vec[[i]])
  df_vec[[i]] <- top_GSEA(res_vec[[i]], n = 7)
}
df <- do.call(rbind, df_vec)
results <- do.call(rbind, res_vec)
results <- format_GSEA(results, df, n=85)
orders <- vector(mode = 'character', length = length(vec)*4)
trtList <- c("p4", "mp4")

orders <- c('NK_p4', 'NK_mp4')

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  scale_x_discrete(limits=orders) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('nk_GSEA_04032021.png', width = 11, height = 5)

# Important Violin and Feature Plots --------------
# combined <- readRDS('combined_07292021_v2.rds')

combined@meta.data$orig.ident <-
  factor(x = combined@meta.data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))

Idents(combined) <- 'celltype'

p1 <- VlnPlot(combined, features = c('CD200-TotalA'), idents = c('AT 2'), 
              split.by = 'orig.ident', assay = 'ADT') + NoLegend()
p2 <- VlnPlot(combined, features = c('CD200r-TotalA'), idents = c('AM'), 
              split.by = 'orig.ident', assay = 'ADT')
p1 | p2
ggsave('CD200_ADT_v2.png', width = 10, height = 5)

# Figure 1 Marker ADTs:
DefaultAssay(combined) <- 'ADT'
adt.to.plot <- c('CD45','CD326', 'CDLy6C', 'Ly-6G', 'CD19', 'IgD', 
                 'CD11b', 'CD170', 'CD4', 'CD8b')
gg <- list()
for (i in 1:length(adt.to.plot)){
  adt <- paste(adt.to.plot[[i]], '-TotalA', sep = '')
  gg[[i]] <- FeaturePlot(combined, features = adt, reduction = 'wnn.umap',
                         cols = c("lightgrey","darkgreen"), min.cutoff=0, 
                         max.cutoff = "q99") + # NoLegend() +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank()) + ggtitle(adt.to.plot[[i]])
}
ggsave('Fig1_ADT_v3.png', plot = cowplot::plot_grid(plotlist = gg, ncol = 5),
       height = 2*3, width = 5*3)

# Plotting a more clean UMAP:
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
ggsave("celltype_distribution_by_trt_leg_v3.png", width = 6, height = 4.2)
p + NoLegend()
ggsave("celltype_distribution_by_trt_noLeg_v3.png", width = 4.2, height = 4.2)

# Dimplot Split by Treatment
p <- DimPlot(combined, reduction = 'wnn.umap', split.by = "orig.ident") + NoLegend()
# change order in plot: 
p$data$orig.ident <- factor(x = p$data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24")) 
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave("dimplot_splitby_trt_v3.png", plot = p, width = 15, height = 4)

# 08/08/2021 Old neutro feature plots --------
# Will port into new neutrophil subs
load("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/pre_webGL_workspace.RData")
FeaturePlot(subset(neutro, orig.ident == 'mp24'), 'Camp') + theme(panel.grid = element_blank(),
                                                                  axis.title = element_blank(),
                                                                  axis.text = element_blank(),
                                                                  axis.ticks = element_blank())
ggsave('Camp_expression_neutro_v1.png', width = 5, height = 5)

# 08/08/2021 Violin Plots --------
combined@meta.data$orig.ident <-
  factor(x = combined@meta.data$orig.ident, levels = c("naive", "p4", "mp4", "p24", "mp24"))

# combined <- readRDS('combined_07292021_v2.rds')
Idents(combined) <- 'celltype'
VlnPlot(combined, features = c('Cxcl10'), idents = 'Lipofib', split.by = 'orig.ident')
ggsave('Cxcl10_lipofib.png', width = 4, height = 4)
VlnPlot(combined, features = c('Ccl4'), idents = 'Lipofib', split.by = 'orig.ident')
ggsave('Ccl4_lipofib.png', width = 4, height = 4)

# For Bhawana BMES:
p <- DimPlot(combined, reduction = 'wnn.umap', group.by = 'orig.ident', shuffle = T)
p <- p + theme(title = element_blank(), legend.position = "none",
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank())
ggsave('BMES_Bhawana_Cb_v2.png', plot = p, width = 4, height = 2.75)

# DEGs for mp vs. p treatments at 4 and 24 hours ------------
# Do not parallelize this section. Errors happen when calculating DEGs.
# combined <- readRDS('combined_07292021_v2.rds')

freq_table <- read.csv('t2_condensed.csv', row.names = 'X')

Idents(combined) <- "celltype.trt"

# DEGs will be calculated relative to naive
DefaultAssay(combined) <- "SCT" 
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

saveRDS(degList, "2021-08-21 DEG List mp vs p.rds")
saveRDS(graphList, "2021-08-21 Graph List mp vs p.rds")

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
write.csv(degList, file = "2021-08-21 mp vs. p DEGs.csv")

# Exporting UMAP coordinates for scVelo -----
# combined <- readRDS('combined_07292021_v2.rds')
Idents(combined) <- 'celltype'

# Cluster based off unintegrated SCT values:
DefaultAssay(combined) <- "SCT"
VariableFeatures(combined) <- rownames(combined)
combined <- RunPCA(combined, assay = 'SCT', verbose = T, reduction.name = 'sct_unintPC', reduction.key = 'sct_unintPC_')
ElbowPlot(combined, ndims = 50, reduction = 'sct_unintPC')
ggsave("elbow_plot_combined_sct_non_integrated.png")
combined <- RunUMAP(combined, reduction = 'sct_unintPC', dims = 1:40, assay = 'SCT', 
                  reduction.name = 'sct.umap.unint', reduction.key = 'sct.umap.unint_')
DimPlot(combined, reduction = 'sct.umap.unint', label = T, repel = T) + NoLegend()
ggsave('combined_unintegrated.png', width = 5, height = 5)

DimPlot(combined, reduction = 'sct.umap.unint', split.by = 'orig.ident')
ggsave('combined_unintegrated_split.png', width = 5*5, height = 5)s

Idents(combined) <- 'orig.ident'
DimPlot(combined, reduction = 'sct.umap.unint') 
ggsave('combined_sct_unint_orig_ident.png', width = 6, height = 5)

# Remove cells not in the UMAP region: -5 to +5 and less than -1 approx.
# Likely small amount of doublets
Idents(combined) <- 'celltype'
neutro <- subset(combined, idents = 'N')
plot <- DimPlot(neutro, reduction = 'sct.umap.unint')
dubs <- CellSelector(plot = plot)
dubs
# [1] "p4_GATTCTTGTCGGCCTA-1"   "p24_ACAGAAAGTAAGGCCA-1"  "p24_TGCATGATCAGTGATC-1"  "mp24_CGCGTGAGTAGCTTTG-1"
dubs2 <- CellSelector(plot = plot)
dubs2
# [1] "p24_TGTCCACCACTCGATA-1"  "mp24_GGGTGAAAGCCGATAG-1" "mp24_GTGCTTCCAGTTAAAG-1"
neutro <- subset(neutro, cells = c(dubs, dubs2), invert = TRUE)

# Plot with labeled orig.ident:
Idents(neutro) <- 'orig.ident'
DimPlot(neutro, reduction = 'sct.umap.unint')
ggsave('neutro_orig_ident.png', width = 6, height = 5)

# Find distinct phenotypic states of neutrophils:
neutro <- FindNeighbors(neutro, reduction = 'sct_unintPC')
neutro <- FindClusters(neutro, resolution = 0.1)
p <- DimPlot(neutro)
p <- p + theme(panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave('neutro_res0.1_clustering.png', width = 6, height = 5)

neutro <- RenameIdents(neutro,
                       '0' = 'N0',
                       '1' = 'N1',
                       '2' = 'N2',
                       '3' = 'N3')
neutro$neutro_idents <- Idents(neutro)

DimPlot(neutro, split.by = 'orig.ident')
ggsave('neutro_split_by.png', width = 25, height = 5)

# Extracting metadata including unintegrated UMAP coordinates
# From: https://github.com/basilkhuder/Seurat-to-RNA-Velocity
write.csv(Cells(neutro), file = "cellID_obs.csv", row.names = FALSE)
write.csv(Embeddings(neutro, reduction = "sct.umap.unint"), file = "cell_embeddings.csv")
write.csv(neutro@meta.data$neutro_idents, file = "clusters.csv")

# Neutrophil gene expression scores ----
library(readxl)
library(ggpubr)
xie_2020 <- read_excel("neutrophil_scores_Xie_2020.xlsx")
# compare only features in the data
neutro_features <- list(xie_2020$`Positive regulation of apoptotic process (GO:0043065)`, 
                        xie_2020$`Necroptosis (GO:0070266)`, 
                        xie_2020$`Chemotaxis (GO:0030593)`,
                        xie_2020$`Phagocytosis (GO:0006911)`,
                        xie_2020$`NADPH oxidase (Henderson and Chappel, 1996)`)
neutro <- AddModuleScore(neutro, features = neutro_features, name = 'xie')

FeaturePlot(neutro, features = c('xie1', 'xie2', 'xie3'))
FeaturePlot(neutro, features = c('xie4', 'xie5'))

compare <- list(c('N0', 'N1'), c('N0', 'N1'), c('N0', 'N2'), c('N0', 'N3'),
                c('N1', 'N2'), c('N1', 'N3'), c('N2', 'N3'))
VlnPlot(neutro, features = 'xie1', pt.size = 0) + ylim(0.05, 0.3) + NoLegend() + geom_boxplot() + ggtitle('apoptosis score') +
  stat_compare_means(comparisons = compare, label = "p.signif")
ggsave(filename = 'apoptosis_neutro_sig.png')
VlnPlot(neutro, features = 'xie2', pt.size = 0) + ylim(-0.2, 0.7) + NoLegend() + geom_boxplot() + ggtitle('necroptosis score') +
  stat_compare_means(comparisons = compare, label = "p.signif")
ggsave(filename = 'necroptosis_neutro_sig.png')
VlnPlot(neutro, features = 'xie3', pt.size = 0) + ylim(0, 1.3) + NoLegend() + geom_boxplot() + ggtitle('chemotaxis score') +
  stat_compare_means(comparisons = compare, label = "p.signif")
ggsave(filename = 'chemotaxis_neutro_sig.png')
VlnPlot(neutro, features = 'xie4', pt.size = 0) + ylim(0, 0.85) + NoLegend() + geom_boxplot() + ggtitle('phagocytosis score') +
  stat_compare_means(comparisons = compare, label = "p.signif")
ggsave(filename = 'phagocytosis_neutro_sig.png')
VlnPlot(neutro, features = 'xie5', pt.size = 0) + ylim(0, 2.4) + NoLegend() + geom_boxplot() + ggtitle('NADPH oxidase score') +
  stat_compare_means(comparisons = compare, label = "p.signif")
ggsave(filename = 'NADPH_oxidase_neutro_sig.png')

# Plot specific early neutrophil genes:
DefaultAssay(neutro) <- 'SCT'
p <- FeaturePlot(neutro, features = 'Camp')
p <- p + theme(panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave('camp_expression.png', width = 6, height = 5)

p <- FeaturePlot(neutro, features = 'Mmp8')
p <- p + theme(panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave('Mmp8_expression.png', width = 6, height = 5)

p <- FeaturePlot(neutro, features = 'Ngp')
p <- p + theme(panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave('Ngp_expression.png', width = 6, height = 5)

p <- FeaturePlot(neutro, features = 'Ltf')
p <- p + theme(panel.grid = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())
ggsave('Ltf_expression.png', width = 6, height = 5)

saveRDS(neutro, 'neutro_12252021.rds')
# neutro <- readRDS('neutro_12252021.rds')

# Find markers for different neutrophil clusters: ----
DefaultAssay(neutro) <- "SCT"
neutro.markers <- FindAllMarkers(neutro, only.pos = FALSE, 
                                 min.pct = 0, logfc.threshold = 0, 
                                 max.cells.per.ident = Inf)
write.csv(neutro.markers, file = "neutro_unint_all_markers.csv", row.names = TRUE)

neutro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(neutro, features = top10$gene) + NoLegend()
neutro.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC) -> top5
DoHeatmap(neutro, features = top5$gene) + NoLegend()
ggsave('neutro_heatmap_markers.png', width = 5, height = 5)

markers_24h_neutro <- FindMarkers(neutro, only.pos = FALSE, 
                                     ident.1 = 'N0', ident.2 = 'N2',
                                     min.pct = 0, logfc.threshold = 0, 
                                     max.cells.per.ident = Inf)
write.csv(markers_24h_neutro, file = "neutro_unint_N0_vs_N2_markers.csv", row.names = TRUE)
# markers_24h_neutro <- read.csv(file = "neutro_unint_N0_vs_N2_markers.csv", row.names = 'X')

markers_4h_neutro <- FindMarkers(neutro, only.pos = FALSE, 
                                  ident.1 = 'N1', ident.2 = 'N3',
                                  min.pct = 0, logfc.threshold = 0, 
                                  max.cells.per.ident = Inf)
write.csv(markers_4h_neutro, file = "neutro_unint_N1_vs_N3_markers.csv", row.names = TRUE)
# markers_4h_neutro <- read.csv(file = "neutro_unint_N1_vs_N3_markers.csv", row.names = 'X')

# plotting neutro DEG results: ----
format_DEG <- function(df){
  df$diffexpressed <- "NO"
  df$diffexpressed[df$avg_log2FC > 0.6 & df$p_val_adj < 0.05] <- "UP"
  df$diffexpressed[df$avg_log2FC < -0.6 & df$p_val_adj < 0.05] <- "DOWN"
  df <- cbind(gene_symbol = rownames(df), df)
  df$delabel <- NA
  df$delabel[df$diffexpressed != "NO"] <- df$gene_symbol[df$diffexpressed != "NO"]
  return(df)
} 
markers_24h_neutro <- format_DEG(markers_24h_neutro)
p <- ggplot(data=markers_24h_neutro, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.4) +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("green4", "black", "blue")) + 
  NoLegend() +
  labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(adj. p-value)')) +
  theme(text = element_text(size=8, family = "sans"),
        axis.title=element_text(size=10, family = "sans", face="bold"))
ggsave('neutro_N0_vs_N2.png', width = 5, height = 5)

markers_4h_neutro <- format_DEG(markers_4h_neutro)
p <- ggplot(data=markers_4h_neutro, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.4) +
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("red3", "black", "orangered")) + 
  NoLegend() +
  labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(adj. p-value)')) +
  theme(text = element_text(size=8, family = "sans"),
        axis.title=element_text(size=10, family = "sans", face="bold"))
p
ggsave('neutro_N1_vs_N3.png', width = 5, height = 5)

# Run GSEA on neutro categories: ----
library(fgsea)
library(msigdbr)
msigdbr_collections()

hallmarks <- msigdbr(species = "Mus musculus", category = "H")
hallmarks <- split(x = hallmarks$gene_symbol, f = hallmarks$gs_name)
kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
kegg <- split(x = kegg$gene_symbol, f = kegg$gs_name)
go <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
go <- split(x = go$gene_symbol, f = go$gs_name)
reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
reactome <- split(x = reactome$gene_symbol, f = reactome$gs_name)

# 24h clusters
stats_list <- markers_24h_neutro$avg_log2FC
names(stats_list) <- markers_24h_neutro$gene_symbol
gene_sets <- c(hallmarks, kegg, go, reactome)

gsea <- fgsea(pathways = gene_sets,
              stats = stats_list,
              minSize=15,
              maxSize=Inf,
              nproc = 12,
              eps = 0) # no eps cutoff so we have "true" P-val
gsea$leadingEdge <- as.character(gsea$leadingEdge)
write.csv(gsea, file = 'neutro_N0_v_N2_gsea.csv')
gsea <- filter(gsea, padj < 0.05)
write.csv(gsea, file = 'neutro_N0_v_N2_gsea_sig.csv')

# 4h clusters
stats_list <- markers_4h_neutro$avg_log2FC
names(stats_list) <- markers_4h_neutro$gene_symbol

gsea <- fgsea(pathways = gene_sets,
              stats = stats_list,
              minSize=15,
              maxSize=Inf,
              nproc = 12,
              eps = 0) # no eps cutoff so we have "true" P-val
gsea$leadingEdge <- as.character(gsea$leadingEdge)
write.csv(gsea, file = 'neutro_N1_v_N3_gsea.csv')
gsea <- filter(gsea, padj < 0.05)
write.csv(gsea, file = 'neutro_N1_v_N3_gsea_sig.csv')

# Due to this, plot IFNgamma responsive scores on the feature plot:
neutro <- AddModuleScore(neutro, features = list(hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE), name = 'ifng')
FeaturePlot(neutro, features = 'ifng1') + theme(panel.grid = element_blank(),
                                                 axis.title = element_blank(),
                                                 axis.text = element_blank(),
                                                 axis.ticks = element_blank()) +
  ggtitle(expression('IFN'~gamma ~ 'score'))
ggsave('neutro_ifng_score.png', width = 6, height = 5)

# Plot as vln plot
compare <- list(c('N0', 'N1'), c('N0', 'N1'), c('N0', 'N2'), c('N0', 'N3'),
                c('N1', 'N2'), c('N1', 'N3'), c('N2', 'N3'))
VlnPlot(neutro, features = 'ifng1', pt.size = 0) + ylim(-0.05, 1.5) + NoLegend() + geom_boxplot() + ggtitle(expression('IFN'~gamma ~ 'score')) +
  stat_compare_means(comparisons = compare, label = "p.signif")
ggsave(filename = 'IFNg_neutro_sig.png')

# Root cell identification ------
# Biologically known to be cells that are synthesizing granules, i.e. Camp transcripts
Idents(neutro) <- 'orig.ident'
p <- DimPlot(subset(neutro, idents = 'p4'))
CellSelector(p)
# [1] "p4_CTACTATAGGAACTCG-1"

p <- DimPlot(subset(neutro, idents = 'mp4'))
CellSelector(p)
# [1] "mp4_ACAGAAAAGATGTTAG-1"

# Plotting IFN stuff with Seurat -----
library(ggpubr)
# combined <- readRDS('combined_07292021_v2.rds')
Idents(combined) <- 'celltype'
combined@meta.data$orig.ident <- factor(x = combined@meta.data$orig.ident, 
                                        levels = c('naive', 'p4', 'mp4', 'p24', 'mp24'))

compare <- list(c('naive', 'p4'), c('naive', 'mp4'), c('naive', 'p24'), c('naive', 'mp24'))
VlnPlot(combined, idents = 'NK', features = 'Ifng', split.by = 'orig.ident') + stat_compare_means(comparisons = compare, label = 'p.signif')
ggsave('Ifng_NK.png')

# GSEA Interferon Results with all celltypes ----
gsea <- read.csv('fgsea_results_12292021_PCT_EPS_FILT.csv')
gsea$trt <- str_split_fixed(gsea$cluster, '_', 2)[,2]

paths <- c('HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'GOBP_RESPONSE_TO_INTERFERON_GAMMA',
           'HALLMARK_INTERFERON_ALPHA_RESPONSE',
           'REACTOME_INTERFERON_SIGNALING',
           'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING')

trts <- c('p24', 'mp24')
results <- gsea %>% filter(pathway %in% paths) %>% filter(trt %in% trts)

ggplot(results, aes(x = pathway, y = cluster, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_interferons_24.png', height = 12)

paths <- c('HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'GOBP_RESPONSE_TO_INTERFERON_GAMMA',
           'HALLMARK_INTERFERON_ALPHA_RESPONSE',
           'REACTOME_INTERFERON_SIGNALING',
           'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING')

trts <- c('p4', 'mp4')
results <- gsea %>% filter(pathway %in% paths) %>% filter(trt %in% trts)

ggplot(results, aes(x = pathway, y = cluster, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_interferons_4.png', height = 12)

# Bar plots:
gsea <- read.csv('fgsea_results_12292021_PCT_EPS.csv')
gsea$trt <- str_split_fixed(gsea$cluster, '_', 2)[,2]
trts <- c('p24', 'mp24')
results <- gsea %>% filter(pathway %in% c('HALLMARK_INTERFERON_GAMMA_RESPONSE')) %>% filter(trt %in% trts)
ggplot(results, aes(x = cluster, y = -log10(padj))) + 
  geom_bar(stat = 'identity', fill = 'red3') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  coord_flip() + theme(axis.title.y = element_blank())
ggsave('hallmark_Ifng_response_24h.png', height = 12)

results <- gsea %>% filter(pathway %in% c('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE')) %>% filter(trt %in% trts)
results$pathway <- str_replace(results$pathway, 'HALLMARK_INTERFERON', 'IFN')
ggplot(results, aes(fill = pathway, x = cluster, y = -log10(padj))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  coord_flip() + theme(axis.title.y = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave('hallmark_dual_IFN_response_24h.png', height = 12, width = 5)

results <- gsea %>% filter(pathway %in% c('HALLMARK_INTERFERON_ALPHA_RESPONSE')) %>% filter(trt %in% trts)
ggplot(results, aes(x = cluster, y = -log10(padj))) + 
  geom_bar(stat = 'identity', fill = 'red4') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  coord_flip() + theme(axis.title.y = element_blank())
ggsave('hallmark_Ifna_response_24h.png', height = 12)

trts <- c('p4', 'mp4')
results <- gsea %>% filter(pathway %in% c('HALLMARK_INTERFERON_GAMMA_RESPONSE')) %>% filter(trt %in% trts)
ggplot(results, aes(x = cluster, y = -log10(padj))) + 
  geom_bar(stat = 'identity', fill = 'red3') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  coord_flip() + theme(axis.title.y = element_blank())

# Heatmap with interferon-related genes: ----
library(ComplexHeatmap)

# Get a list of IFN-related genes: use both hallmarks
gene_sets <- gene_set_maker()
genes <- Reduce(union, gene_sets[paths])

degs <- read.csv('2021-07-29 DEGs.csv')
degs$trt <- str_split_fixed(degs$cellType, '_', 2)[,2]

# note: hallmark IFN gene sets are fairly similar
genes <- intersect(genes, degs$gene_symbol)
quantile(degs$avg_log2FC)

degs <- filter(degs, trt %in% c('mp24', 'p24'))
unique(intersect(degs$gene_symbol, genes))
degs <- filter(degs, gene_symbol %in% genes)

degs <- degs[, c("gene_symbol", "avg_log2FC", 'cellType')]
degs <- pivot_wider(degs, 
                    names_from = cellType, 
                    values_from = avg_log2FC,
                    values_fill = 0) %>% as.data.frame()
rownames(degs) <- degs$gene_symbol
degs <- select(degs, -c('gene_symbol'))

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Massive heatmap for data exploration:
Heatmap(degs, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

# Get top 25 genes by average expression from these pathways:
top25 <- sort(rowMeans(degs), decreasing = TRUE)[1:25]
degs <- degs[names(top25), ]

col_fun = circlize::colorRamp2(c(-4, 0, 7), c("blue3", "white", "red3"))
Heatmap(degs, name = 'log2FC',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))
# save as 5.5 x 11 in pdf

# GSEA Ribosome Results on all celltypes -----
gsea <- read.csv('fgsea_results_12292021_PCT_EPS_FILT.csv')
gsea$trt <- str_split_fixed(gsea$cluster, '_', 2)[,2]

paths <- c('REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION',
           'REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE',
           'REACTOME_EUKARYOTIC_TRANSLATION_INITIATION',
           'REACTOME_TRANSLATION',
           'GOBP_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE',
           'GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM',
           'GOBP_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM',
           'GOBP_TRANSLATIONAL_INITIATION',
           'KEGG_RIBOSOME')

trts <- c('p24', 'mp24')
results <- gsea %>% filter(pathway %in% paths) %>% filter(trt %in% trts)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_ribosome_24.png', width = 16)

trts <- c('p4', 'mp4')
results <- gsea %>% filter(pathway %in% paths) %>% filter(trt %in% trts)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_ribosome_4.png', width = 16.5)

# Bar plots:
# Pull all GSEA including non significant
gsea <- read.csv('fgsea_results_12292021_PCT_EPS.csv')
gsea$trt <- str_split_fixed(gsea$cluster, '_', 2)[,2]
trts <- c('p4', 'mp4')

results <- gsea %>% filter(pathway %in% c('KEGG_RIBOSOME', 'REACTOME_TRANSLATION', 'GOBP_TRANSLATIONAL_INITIATION')) %>% filter(trt %in% trts)
ggplot(results, aes(fill = pathway, x = cluster, y = -log10(padj))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  coord_flip() + theme(axis.title.y = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave('ribosome_translation_4hr_width5.png', height = 12, width = 5)
ggsave('ribosome_translation_4hr_width6.png', height = 12, width = 6)
ggsave('ribosome_translation_4hr_width7.png', height = 12, width = 7)

trts <- c('p24', 'mp24')
results <- gsea %>% filter(pathway %in% c('KEGG_RIBOSOME', 'REACTOME_TRANSLATION', 'GOBP_TRANSLATIONAL_INITIATION')) %>% filter(trt %in% trts)
ggplot(results, aes(fill = pathway, x = cluster, y = -log10(padj))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + 
  coord_flip() + theme(axis.title.y = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave('ribosome_translation_24hr_width5.png', height = 12, width = 5)
ggsave('ribosome_translation_24hr_width6.png', height = 12, width = 6)
ggsave('ribosome_translation_24hr_width7.png', height = 12, width = 7)

# Heatmap with ribosome-realted genes: ----
library(ComplexHeatmap)
library(stringr)

# Get a list of ribosome-related genes:
gene_sets <- gene_set_maker()
genes <- Reduce(union, gene_sets[paths])

degs <- read.csv('2021-07-29 DEGs.csv')
degs$trt <- str_split_fixed(degs$cellType, '_', 2)[,2]

genes <- intersect(genes, degs$gene_symbol)
quantile(degs$avg_log2FC)

degs <- filter(degs, trt %in% c('p4', 'mp4'))
unique(intersect(degs$gene_symbol, genes))
degs <- filter(degs, gene_symbol %in% genes)

degs <- degs[, c("gene_symbol", "avg_log2FC", 'cellType')]
degs <- tidyr::pivot_wider(degs, 
                    names_from = cellType, 
                    values_from = avg_log2FC,
                    values_fill = 0) %>% as.data.frame()
rownames(degs) <- degs$gene_symbol
degs <- select(degs, -c('gene_symbol'))

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Massive heatmap for data exploration:
Heatmap(degs, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

# Get top (bottom) 25 genes by average expression from these pathways:
top50 <- sort(rowMeans(degs), decreasing = FALSE)[1:50]
degs50 <- degs[names(top50), ]

col_fun = circlize::colorRamp2(c(-1.5, 0, 1), c("blue3", "white", "red3"))
Heatmap(degs50, name = 'log2FC',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))
# save as 5.5 x 11 in pdf

# Plot ribosome-related genes on KEGG pathway -----
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pathview)

genes <- rowMeans(degs)
ids <- mapIds(org.Mm.eg.db, keys=names(genes), column="ENTREZID", keytype="SYMBOL", 
              multiVals="first")
names(genes) <- ids

pv.out <- pathview(gene.data = genes, pathway.id = "03010",
                   species = "mmu", limit = list(gene=1, cpd=1),
                   out.suffix = "ribosome_green")

# Export markers for all celltypes: --------
# combined <- readRDS('combined_07292021_v2.rds')
Idents(combined) <- 'celltype'
DefaultAssay(combined) <- 'SCT'
biomarkers <- FindAllMarkers(combined, logfc.threshold = 0, only.pos = FALSE,
                             max.cells.per.ident = Inf)
top10 <- biomarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(biomarkers, file = "gene_biomarkers_all.csv", row.names = FALSE)
write.csv(top10, file = "gene_biomarkers_top10.csv", row.names = FALSE)

# Make feature plots for the most important marker genes: -----
DefaultAssay(combined) <- 'SCT'
rna.to.plot <- c('Ly6c2', 'Gngt2', 'Gpihbp1', 'Kdr', 'Vwf', 'Sftpc', 
                 'Hopx', 'Tagln', 'Postn', 'Inmt', 'Gzma', 'Rora')
gg <- list()
for (i in 1:length(rna.to.plot)){
  gg[[i]] <- FeaturePlot(combined, features = rna.to.plot[[i]], reduction = 'wnn.umap') + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) + ggtitle(rna.to.plot[[i]])
}
ggsave('Fig1_RNA.png', plot = cowplot::plot_grid(plotlist = gg, ncol = 6),
       height = 2*3, width = 6*3)

# Plot all ADTs
DefaultAssay(combined) <- 'ADT'
protein.to.plot <- rownames(combined@assays$ADT)
for (i in 1:length(protein.to.plot)){
  p <- FeaturePlot(combined, features = protein.to.plot[[i]], cols = c("lightgrey","darkgreen"),
                   min.cutoff = 0, max.cutoff = 'q99', reduction = 'wnn.umap')
  ggsave(gsub('/','',paste0(protein.to.plot[[i]],'.png', sep = '')), 
         plot = p)
}

# Calculate all DEGs for orig.ident comparion: ------
# combined <- readRDS('combined_07292021_v2.rds')
Idents(combined) <- 'orig.ident'
DefaultAssay(combined) <- 'SCT'
degs.orig.ident <- FindAllMarkers(combined, logfc.threshold = 0, only.pos = FALSE,
                                  max.cells.per.ident = Inf)
write.csv(degs.orig.ident, file = "degs_orig_ident.csv", row.names = FALSE)

# make comparison to naive
degList <- vector(mode = 'list', length = 4)
trtList <- c('p4', 'mp4', 'p24', 'mp24')
for (i in 1:length(degList)){
  print(trtList[[i]])
  degList[[i]] <- FindMarkers(combined, ident.1 = trtList[[i]], ident.2 = 'naive',
                              logfc.threshold = 0, only.pos = FALSE,
                              max.cells.per.ident = Inf)
  degList[[i]]$trt <- trtList[[i]]
}
write.csv(do.call(rbind, degList), 'overall_DEGs_vs_naive.csv')

# Export total nCountRNA values: -----
library(ggpubr)
# combined <- readRDS('combined_07292021_v2.rds')
combined@meta.data$orig.ident <- factor(x = combined@meta.data$orig.ident, 
                                        levels = c('naive', 'p4', 'mp4', 'p24', 'mp24'))
Idents(combined) <- 'orig.ident'
compare <- list(c('naive', 'p4'), c('naive', 'mp4'), c('naive', 'p24'), c('naive', 'mp24'))
p <- VlnPlot(combined, features = 'nCount_RNA', split.by = 'orig.ident', pt.size = 0, y.max = 60000)  + 
  stat_compare_means(comparisons = compare, label = 'p.signif') +
  stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
p
ggsave(plot = p, filename = 'nCount_RNA_global.png', width = 5, height = 5)

# Way to get colors on UMAP with scanpy:
# old method pre-christmas 2021
# save metadata table:
neutro$barcode <- colnames(neutro)
neutro$UMAP_1 <- neutro@reductions$sct.umap.unint@cell.embeddings[,1]
neutro$UMAP_2 <- neutro@reductions$sct.umap.unint@cell.embeddings[,2]
write.csv(neutro@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(neutro, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(neutro@reductions$sct_unintPC@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

# Compute GSEA for MP vs. P categories: -----------
library(fgsea)

result_table <- read.csv('2021-08-21 mp vs. p DEGs.csv', row.names = 'X')

gene_set_maker <- function(){
  library(msigdbr)
  
  # msigdbr_collections()
  
  hallmarks <- msigdbr(species = "Mus musculus", category = "H")
  hallmarks <- split(x = hallmarks$gene_symbol, f = hallmarks$gs_name)
  kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
  kegg <- split(x = kegg$gene_symbol, f = kegg$gs_name)
  go <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
  go <- split(x = go$gene_symbol, f = go$gs_name)
  reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
  reactome <- split(x = reactome$gene_symbol, f = reactome$gs_name)
  
  gene_sets <- c(hallmarks, kegg, go, reactome)
  return(gene_sets)
}

# Define a function to run GSEA for a single cluster
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
                maxSize=Inf,
                nproc = 12,
                eps = 0) # no eps cutoff so we have "true" P-val
  gsea$cluster <- cluster
  
  return(gsea)
}

gene_sets <- gene_set_maker()

cluster_list <- unique(result_table$cellType)
fgsea_results <- lapply(cluster_list, runGSEA)

saveRDS(cluster_list, 'cluster_list_pos_neg_mp_vs_p_02132022.rds')
saveRDS(fgsea_results, 'fgsea_results_pos_neg_mp_vs_p_02132022.rds')

fgsea_results <- do.call("rbind", fgsea_results)
fgsea_results <- as.data.frame(fgsea_results)
fgsea_results$leadingEdge <- as.character(fgsea_results$leadingEdge)
write.csv(fgsea_results, 'fgsea_results_mp_vs_p_02132022.csv')

fgsea_results <- filter(fgsea_results, padj < 0.05)
write.csv(fgsea_results, 'fgsea_results_mp_vs_p_FILT_02132022.csv')

# Making table of GSEA frequencies:

all_filt_res <- read.csv('fgsea_results_mp_vs_p_FILT_02132022.csv')
term_frequencies <- table(all_filt_res$pathway)
write.csv(term_frequencies, 'fgsea_results_mp_vs_p_FILT_all_numb_02132022.csv')

df_total <- top_GSEA(all_filt_res, n = 30)
write.csv(df_total, 'fgsea_results_30_top_pos_neg_mp_vs_p_02132022.csv')

# Visualizing most represented GSEA results across all celltypes mp vs p: -------
library(stringr)
gsea <- read.csv('fgsea_results_mp_vs_p_FILT_02132022.csv')

# Neutrophil chemotaxis
paths <- c('GOBP_NEUTROPHIL_MIGRATION',
           'GOBP_GRANULOCYTE_CHEMOTAXIS',
           'GOBP_NEUTROPHIL_CHEMOTAXIS',
           'GOBP_GRANULOCYTE_MIGRATION',
           'GOBP_LEUKOCYTE_CHEMOTAXIS')

results <- gsea %>% filter(pathway %in% paths)

ggplot(results, aes(x = pathway, y = cluster, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_mp_vs_p_neutrophil_chemotaxis.png', height = 12)

# Bar plots:
gsea <- read.csv('fgsea_results_mp_vs_p_02132022.csv')
results <- gsea %>% filter(pathway %in% c('GOBP_NEUTROPHIL_CHEMOTAXIS', 'GOBP_NEUTROPHIL_MIGRATION'))
p <- ggplot(results, aes(fill = pathway, x = cluster, y = -sign(NES)*log10(padj))) + 
            geom_bar(position = 'dodge', stat = 'identity') + 
            geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + geom_hline(yintercept = log10(0.05), linetype = 'dotted') +
            coord_flip() + theme(axis.title.y = element_blank()) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")); p
ggsave('gobp_neutrophil_chemotaxis_mp_vs_p.png', height = 12)
p + NoLegend()
ggsave('gobp_neutrophil_chemotaxis_mp_vs_p_nolegend.png', height = 12)

# Chemotaxis Heatmap: ----
library(ComplexHeatmap)
library(stringr)

# Get a list of ribosome-related genes:
gene_sets <- gene_set_maker()
genes <- Reduce(union, gene_sets[paths])

degs <- read.csv('2021-08-21 mp vs. p DEGs.csv')
degs$trt <- str_split_fixed(degs$cellType, '_', 2)[,2]

genes <- intersect(genes, degs$gene_symbol)
quantile(degs$avg_log2FC)

#degs <- filter(degs, trt %in% c('p4', 'mp4'))
unique(intersect(degs$gene_symbol, genes))
degs <- filter(degs, gene_symbol %in% genes)

# Include only cells that have significant increases in neutrophil chemotaxis
# And include genes that are > 0
sig_sets <- gsea %>% filter(padj < 0.05 & NES > 0 & pathway %in% c('GOBP_NEUTROPHIL_CHEMOTAXIS', 'GOBP_NEUTROPHIL_MIGRATION'))
sig_celltypes <- sig_sets$cluster %>% unique()
degs <- filter(degs, cellType %in% sig_celltypes)

degs <- degs[, c("gene_symbol", "avg_log2FC", 'cellType')]
degs <- tidyr::pivot_wider(degs, 
                           names_from = cellType, 
                           values_from = avg_log2FC,
                           values_fill = 0) %>% as.data.frame()
rownames(degs) <- degs$gene_symbol
degs <- select(degs, -c('gene_symbol'))

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Massive heatmap for data exploration:
Heatmap(degs, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

# Get top (bottom) 25 genes by average expression from these pathways:
top25 <- sort(rowMeans(degs), decreasing = TRUE)[1:25]
degs25 <- degs[names(top25), ]

col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue3", "white", "red3"))
Heatmap(degs25, name = 'log2FC',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))
# save as 5.5 x 6.5 in pdf
# neutrophil_chemotaxis_heatmap

# Alveolar macrophage 4 hr volcano plot ------

# define function to label DEGs on volcano plot for a given table:
volcano_plotter <-  function(tab){
  # add a column to tell whether the genes are up or down regulated
  tab$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue-adj < 0.05, set as "UP"
  tab$diffexpressed[tab$avg_log2FC > 0.5 & tab$p_val_adj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue-adj < 0.05, set as "DOWN"
  tab$diffexpressed[tab$avg_log2FC < -0.5 & tab$p_val_adj < 0.05] <- "DOWN"
  # adding gene symbols in a column for easy accesss
  # tab <- cbind(gene_symbol = rownames(tab), tab)
  # defining whether a gene is sufficiently differentially expressed or not
  tab$delabel <- NA
  tab$delabel[tab$diffexpressed != "NO"] <- tab$gene_symbol[tab$diffexpressed != "NO"]
  
  # dealing with bad colors - DOES THIS WORK???????????
  if ( !any(tab$diffexpressed == "DOWN") ){
    colorValues <- c("black", "red")
  } else {
    colorValues <- c("blue", "black", "red")
  }
  
  # plotting
  ggplot(data=tab, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
    geom_point(size = 0.4) +
    theme_classic() +
    geom_text_repel() +
    scale_color_manual(values=colorValues) + 
    NoLegend() +
    labs(x = expression('log' [2] * '(FC)'), y = expression('log' [10] * '(adj. p-value)')) +
    theme(text = element_text(size=8, family = "sans"),
          axis.title=element_text(size=10, family = "sans", face="bold"))
}

deg_AM <- read.csv('2021-08-21 mp vs. p DEGs.csv', row.names = 'X') %>% filter(cellType == 'AM 4 hr') 
rownames(deg_AM) <- deg_AM$gene_symbol
volcano_plotter(deg_AM); ggsave('AM_4_mp_vs_p_pub.png', width = 4.25, height = 4.25)

deg_AM <- read.csv('2021-08-21 mp vs. p DEGs.csv', row.names = 'X') %>% filter(cellType == 'AM 24 hr') 
rownames(deg_AM) <- deg_AM$gene_symbol
volcano_plotter(deg_AM); ggsave('AM_24_mp_vs_p_pub.png', width = 4.25, height = 4.25)

# TNF-a signaling: ----
paths <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
           'REACTOME_INTERLEUKIN_10_SIGNALING',
           'GOBP_INFLAMMATORY_RESPONSE',
           'GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN')

results <- gsea %>% filter(pathway %in% paths)

ggplot(results, aes(x = pathway, y = cluster, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_mp_vs_p_nfkb.png', height = 12)

# Bar plots:
gsea <- read.csv('fgsea_results_mp_vs_p_02132022.csv')
results <- gsea %>% filter(pathway %in% c('HALLMARK_TNFA_SIGNALING_VIA_NFKB'))
p <- ggplot(results, aes(fill = pathway, x = cluster, y = -sign(NES)*log10(padj))) + 
            geom_bar(position = 'dodge', stat = 'identity') + 
            geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + geom_hline(yintercept = log10(0.05), linetype = 'dotted') +
            coord_flip() + theme(axis.title.y = element_blank()) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")); p
ggsave('gobp_nfkb_mp_vs_p.png', height = 12)
p + NoLegend()
ggsave('gobp_nfkb_mp_vs_p_nolegend.png', height = 12)

# Ion homeostasis gene ontology -----
paths <- c('GOBP_CELLULAR_TRANSITION_METAL_ION_HOMEOSTASIS',
           'GOBP_TRANSITION_METAL_ION_HOMEOSTASIS')

results <- gsea %>% filter(pathway %in% paths)

ggplot(results, aes(x = pathway, y = cluster, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_mp_vs_p_metal_ions.png', height = 12)

# Bar plots:
gsea <- read.csv('fgsea_results_mp_vs_p_02132022.csv')
results <- gsea %>% filter(pathway %in% c('GOBP_CELLULAR_TRANSITION_METAL_ION_HOMEOSTASIS', 'GOBP_TRANSITION_METAL_ION_HOMEOSTASIS'))
p <- ggplot(results, aes(fill = pathway, x = cluster, y = -sign(NES)*log10(padj))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted') + geom_hline(yintercept = log10(0.05), linetype = 'dotted') +
  coord_flip() + theme(axis.title.y = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")); p
ggsave('gobp_ion_mp_vs_p_homeostasis.png', height = 12)
p + NoLegend()
ggsave('gobp_ion_mp_vs_p_homneostasis_nolegend.png', height = 12)

# Ion homeostasis heatmap ----
library(ComplexHeatmap)
library(stringr)

# Get a list of ribosome-related genes:
gene_sets <- gene_set_maker()
genes <- Reduce(union, gene_sets[paths])

degs <- read.csv('2021-08-21 mp vs. p DEGs.csv')
degs$trt <- str_split_fixed(degs$cellType, '_', 2)[,2]

genes <- intersect(genes, degs$gene_symbol)
quantile(degs$avg_log2FC)

#degs <- filter(degs, trt %in% c('p4', 'mp4'))
unique(intersect(degs$gene_symbol, genes))
degs <- filter(degs, gene_symbol %in% genes)

# Include only cells that have significant increases in iron homeostasis
# And include genes that are > 0
sig_sets <- gsea %>% filter(padj < 0.05 & NES > 0 & pathway %in% c('GOBP_CELLULAR_TRANSITION_METAL_ION_HOMEOSTASIS', 'GOBP_TRANSITION_METAL_ION_HOMEOSTASIS'))
sig_celltypes <- sig_sets$cluster %>% unique()
degs <- filter(degs, cellType %in% sig_celltypes)

degs <- degs[, c("gene_symbol", "avg_log2FC", 'cellType')]
degs <- tidyr::pivot_wider(degs, 
                           names_from = cellType, 
                           values_from = avg_log2FC,
                           values_fill = 0) %>% as.data.frame()
rownames(degs) <- degs$gene_symbol
degs <- select(degs, -c('gene_symbol'))

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Massive heatmap for data exploration:
Heatmap(degs, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

# Get top (bottom) 10 genes by average expression from these pathways:
top10 <- sort(rowMeans(degs), decreasing = TRUE)[1:25]
degs10 <- degs[names(top10), ]

col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue3", "white", "red3"))
Heatmap(degs10, name = 'log2FC',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))
# save as 3.5 x 9.5 in pdf
# ion_homeostasis_mp_vs_p

# Combined heatmap -----
paths <- c('GOBP_NEUTROPHIL_MIGRATION',
           'GOBP_GRANULOCYTE_CHEMOTAXIS',
           'GOBP_NEUTROPHIL_CHEMOTAXIS',
           'GOBP_GRANULOCYTE_MIGRATION',
           'GOBP_LEUKOCYTE_CHEMOTAXIS',
           'GOBP_CELLULAR_TRANSITION_METAL_ION_HOMEOSTASIS', 
           'GOBP_TRANSITION_METAL_ION_HOMEOSTASIS',
           'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
           'REACTOME_INTERLEUKIN_10_SIGNALING',
           'GOBP_INFLAMMATORY_RESPONSE',
           'GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN')

library(ComplexHeatmap)
library(stringr)

# Get a list of ribosome-related genes:
gene_sets <- gene_set_maker()
genes <- Reduce(union, gene_sets[paths])

degs <- read.csv('2021-08-21 mp vs. p DEGs.csv')
degs$trt <- str_split_fixed(degs$cellType, '_', 2)[,2]

genes <- intersect(genes, degs$gene_symbol)
quantile(degs$avg_log2FC)

#degs <- filter(degs, trt %in% c('p4', 'mp4'))
unique(intersect(degs$gene_symbol, genes))
degs <- filter(degs, gene_symbol %in% genes)

sig_sets <- gsea %>% filter(padj < 0.05 & NES > 0 & pathway %in% paths)
sig_celltypes <- sig_sets$cluster %>% unique()
degs <- filter(degs, cellType %in% sig_celltypes)

degs <- degs[, c("gene_symbol", "avg_log2FC", 'cellType')]
degs <- tidyr::pivot_wider(degs, 
                           names_from = cellType, 
                           values_from = avg_log2FC,
                           values_fill = 0) %>% as.data.frame()
rownames(degs) <- degs$gene_symbol
degs <- select(degs, -c('gene_symbol'))

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Massive heatmap for data exploration:
Heatmap(degs, cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

# Get top (bottom) 10 genes by average expression from these pathways:
top25 <- sort(rowMeans(degs), decreasing = TRUE)[1:25]
degs25 <- degs[names(top25), ]

col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue3", "white", "red3"))
Heatmap(degs25, name = 'log2FC',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))
# save as 5.5 x 8 in pdf
# ion_homeostasis_mp_vs_p

# Combined dot plot mp vs p: -----
# Make dot plot of all the categories in paths above
library(stringr)

gsea <- read.csv('fgsea_results_mp_vs_p_FILT_02132022.csv')

results <- gsea %>% filter(pathway %in% paths)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_ribosome_24.png', width = 16)

trts <- c('p4', 'mp4')
results <- gsea %>% filter(pathway %in% paths) %>% filter(trt %in% trts)

ggplot(results, aes(x = cluster, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() + cowplot::theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_color_gradient2(low = 'blue', mid = 'white',  high = 'red')
ggsave('GSEA_mp_vs_p_big.png', width = 15, height = 4)

# Plotting NK cell DEGs and GSEA only: ----
# These cells did not downregulate ribosomes at 4 hours
# GSEA plotted above
# DEG:
deg_NK <- read.csv('2021-07-29 DEGs.csv', row.names = 'X') %>% filter(cellType == 'NK_p4') 
rownames(deg_AM) <- deg_AM$gene_symbol
volcano_plotter(deg_NK); ggsave('NK_p4_deg.png', width = 4.25, height = 4.25)

deg_NK <- read.csv('2021-07-29 DEGs.csv', row.names = 'X') %>% filter(cellType == 'NK_mp4') 
rownames(deg_AM) <- deg_AM$gene_symbol
volcano_plotter(deg_NK); ggsave('NK_mp4_deg.png', width = 4.25, height = 4.25)

deg_NK <- read.csv('2021-07-29 DEGs.csv', row.names = 'X') %>% filter(cellType == 'NK_p24') 
rownames(deg_AM) <- deg_AM$gene_symbol
volcano_plotter(deg_NK); ggsave('NK_p24_deg.png', width = 4.25, height = 4.25)

deg_NK <- read.csv('2021-07-29 DEGs.csv', row.names = 'X') %>% filter(cellType == 'NK_mp24') 
rownames(deg_AM) <- deg_AM$gene_symbol
volcano_plotter(deg_NK); ggsave('NK_mp24_deg.png', width = 4.25, height = 4.25)

# November 2022 below

# Plotting Table 1 with Cell Counts per Read for publication: -----
# combined <- readRDS('combined_07292021_v2.rds')
# number of cells per treatment
t1 <- table(combined$orig.ident)
write.csv(t1, 'Table1_Cell_Counts.csv')

# Plotting full biomarker list for publication: -----
# combined <- readRDS('combined_07292021_v2.rds')

Idents(combined) <- 'celltype'
DefaultAssay(combined) <- "SCT"
combined.markers <- FindAllMarkers(combined, only.pos = FALSE, 
                                   min.pct = 0, logfc.threshold = 0, 
                                   max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster)
write.csv(combined.markers, file = "TableS3_gene_biomarkers.csv", row.names = FALSE)

# Plotting IgD and IgM for publication of B cells: ----
# combined <- readRDS('combined_07292021_v2.rds')
# Idents(combined) <- 'celltype'

FeaturePlot(combined, features = c('IgM-TotalA', 'IgD-TotalA'))
ggsave('B_cell_markers.png')

