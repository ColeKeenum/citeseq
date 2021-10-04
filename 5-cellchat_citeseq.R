# Load packages / setup ---------------------------
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

#setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
setwd('~/Documents/citeseq-code')

# Split up treatment matrices ---------------------------
combined <- readRDS('combined_07292021_v2.rds')

Idents(combined) <- 'celltype'
# combined[['RNA']] <- NULL # save memory
# combined[['integratedSCT_']] <- NULL
# combined[['integratedADT_']] <- NULL

freq_table <- read.csv('t2_condensed.csv', row.names = 'X')

Idents(combined) <- "celltype.trt"

# DEGs will be calculated relative to naive
DefaultAssay(combined) <- "SCT" # you definitely dont want to do this on integratedSCT_
cellList <- unique(combined$celltype)
cellList <- levels(cellList)

# check freq_table for any celltypes with less than 10 cells in them...
# will not include these for DEG analysis
tmp_cellList <- cellList
for (i in 1:length(cellList)){
  if ( any(freq_table[cellList[[i]], ] < 10)){
    print(cellList[[i]])
    tmp_cellList <- tmp_cellList[tmp_cellList != cellList[i]]
  }
}
cellList <- tmp_cellList; rm(tmp_cellList)

Idents(combined) <- 'celltype'
combined <- subset(combined, idents = cellList)

seurat_list <- SplitObject(combined, split.by = 'orig.ident')

saveRDS(seurat_list, 'split_slim_combined_08212021.rds'); rm(combined)

# Preprocessing function  ---------------------------
preprocess <- function(dgCMatrix, meta){
  meta$celltype <- paste(meta$celltype) # sketchy, doing this to remove the levels
  
  # Create CellChat object
  cellchat <- createCellChat(object = dgCMatrix, meta = meta, group.by = "celltype")
  
  #cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)
  
  # > ALL signaling ---------------------------
  
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  #CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Secreted Signaling
  
  CellChatDB.use <- CellChatDB # full CellChatDB
  
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  # Issue identified!! Please check the official Gene Symbol of the following genes:  
  # H2-Q8 H2-T9 H2-T18 H2-Q9 H2-L H2-BI H2-D H60a H2-Ea-ps 
  
  # this will take a little time
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse) # used to be be PPI.human
  
  # Setting population.size = TRUE because this is on unsorted single cells
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 8)
  
  # infer cell-cell communication at a signaling pathway level:
  cellchat <- computeCommunProbPathway(cellchat)
  
  # calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

# Naive ---------------------------
# Export data from SCTransformed data slot
data.naive <- GetAssayData(object = seurat_list$naive, assay = 'SCT', slot = "data")
meta.naive <- seurat_list$naive@meta.data[c('celltype')]
cc_naive <- preprocess(data.naive, meta.naive)

rm(data.naive)
saveRDS(cc_naive, 'cc_naive_08212021.rds')

# > Visualize  ---------------------------
groupSize <- as.numeric(table(cc_naive@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_naive@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_naive@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# top = 0.1
groupSize <- as.numeric(table(cc_naive@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_naive@net$count, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_naive@net$weight, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Subsetted big networks: Secteted, cell-cell, and ECM
# Define a function to subset the pathways for visualization purposes
path_sub <- function(ref, total = pathways.show){
  CellChatDB <- CellChatDB.mouse
  db_sub <- subsetDB(CellChatDB, search = ref)
  path_names <- unique(db_sub[["interaction"]][["pathway_name"]])
  path_intersect <- intersect(path_names, total)
}

pathways.show <- cc_naive@netP$pathways
sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_naive, signaling = sec_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_naive, signaling = ecm_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_naive, signaling = cc_paths, top = 0.1, layout = "circle")

# analyze signaling from each cell group
mat <- cc_naive@net$weight
par(mfrow = c(3,4), xpd=TRUE) # broke from me
#par(mar=c(1,1,1,1), xpd=TRUE) # worked for me
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# > Emphasis on CD200 pathway  ---------------------------
# pathways.show <- c("CXCL") 
pathways.show <- c('CD200')

# Hierarchy plot
# Alveolar-side cells: AT1, AT2 (two types?), Alv. Macrophages.
#vertex.receiver = seq(2,9) # indices of intersted cells

vertex.receiver = seq(4,6)
netVisual_aggregate(cc_naive, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_naive, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cc_naive, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_naive, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# > Identify Signaling Roles  ---------------------------
# cc_naive <- readRDS('cc_naive_08212021.rds')

# Compute the network centrality scores
cc_naive <- netAnalysis_computeCentrality(cc_naive, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cc_naive, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

saveRDS(cc_naive, 'cc_naive_network_08212021.rds')

rm(cc_naive)

# P4 ---------------------------
data.p4 <- GetAssayData(object = seurat_list$p4, assay = 'SCT', slot = "data")
meta.p4 <- seurat_list$p4@meta.data[c('celltype')]
cc_p4 <- preprocess(data.p4, meta.p4)

rm(data.p4, meta.p4); saveRDS(cc_p4, 'cc_p4_08212021.rds')

# Visualize top = 0.1
groupSize <- as.numeric(table(cc_p4@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_p4@net$count, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_p4@net$weight, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- cc_p4@netP$pathways
sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_p4, signaling = sec_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_p4, signaling = ecm_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_p4, signaling = cc_paths, top = 0.1, layout = "circle")

# > Emphasis on CD200 pathway  ---------------------------
pathways.show <- c('CD200')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_p4, signaling = pathways.show, layout = "circle")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_p4, signaling = pathways.show, color.heatmap = "Reds")

# > Identify Signaling Roles  ---------------------------
# cc_p4 <- readRDS('cc_p4_08212021.rds')

# Compute the network centrality scores
cc_p4 <- netAnalysis_computeCentrality(cc_p4, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cc_p4, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

saveRDS(cc_p4, 'cc_p4_network_08212021.rds')

rm(cc_p4)

# MP4 ---------------------------
data.mp4 <- GetAssayData(object = seurat_list$mp4, assay = 'SCT', slot = "data")
meta.mp4 <- seurat_list$mp4@meta.data[c('celltype')]
cc_mp4 <- preprocess(data.mp4, meta.mp4)

rm(data.mp4, meta.mp4); saveRDS(cc_mp4, 'cc_mp4_08212021.rds')

# Visualize top = 0.1
groupSize <- as.numeric(table(cc_mp4@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_mp4@net$count, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_mp4@net$weight, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- cc_mp4@netP$pathways
sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp4, signaling = sec_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp4, signaling = ecm_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp4, signaling = cc_paths, top = 0.1, layout = "circle")

# > Emphasis on CD200 pathway  ---------------------------
pathways.show <- c('CD200')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp4, signaling = pathways.show, layout = "circle")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_mp4, signaling = pathways.show, color.heatmap = "Reds")

# > Identify Signaling Roles  ---------------------------
# cc_mp4 <- readRDS('cc_mp4_08212021.rds')

# Compute the network centrality scores
cc_mp4 <- netAnalysis_computeCentrality(cc_mp4, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cc_mp4, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

saveRDS(cc_mp4, 'cc_mp4_network_08212021.rds')

rm(cc_mp4)

# P24 ---------------------------
data.p24 <- GetAssayData(object = seurat_list$p24, assay = 'SCT', slot = "data")
meta.p24 <- seurat_list$p24@meta.data[c('celltype')]
cc_p24 <- preprocess(data.p24, meta.p24)

rm(data.p24, meta.p24); saveRDS(cc_p24, 'cc_p24_08212021.rds')

# Visualize top = 0.1
groupSize <- as.numeric(table(cc_p24@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_p24@net$count, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_p24@net$weight, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- cc_p24@netP$pathways
sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_p24, signaling = sec_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_p24, signaling = ecm_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_p24, signaling = cc_paths, top = 0.1, layout = "circle")

# > Emphasis on CD200 pathway  ---------------------------
pathways.show <- c('CD200')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_p24, signaling = pathways.show, layout = "circle")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_p24, signaling = pathways.show, color.heatmap = "Reds")

# > Identify Signaling Roles  ---------------------------
# cc_p24 <- readRDS('cc_p24_08212021.rds')

# Compute the network centrality scores
cc_p24 <- netAnalysis_computeCentrality(cc_p24, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cc_p24, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

saveRDS(cc_p24, 'cc_p24_network_08212021.rds')

rm(cc_p24)

# MP24 ---------------------------
data.mp24 <- GetAssayData(object = seurat_list$mp24, assay = 'SCT', slot = "data")
meta.mp24 <- seurat_list$mp24@meta.data[c('celltype')]
cc_mp24 <- preprocess(data.mp24, meta.mp24)

rm(data.mp24, meta.mp24); saveRDS(cc_mp24, 'cc_mp24_08212021.rds')

# Visualize top = 0.1
groupSize <- as.numeric(table(cc_mp24@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_mp24@net$count, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_mp24@net$weight, vertex.weight = groupSize, top = 0.1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- cc_mp24@netP$pathways
sec_paths<- path_sub('Secreted Signaling', total = pathways.show)
ecm_paths<- path_sub('ECM-Receptor', total = pathways.show)
cc_paths<- path_sub('Cell-Cell Contact', total = pathways.show)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp24, signaling = sec_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp24, signaling = ecm_paths, top = 0.1, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp24, signaling = cc_paths, top = 0.1, layout = "circle")

# > Emphasis on CD200 pathway  ---------------------------
pathways.show <- c('CD200')

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cc_mp24, signaling = pathways.show, layout = "circle")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_mp24, signaling = pathways.show, color.heatmap = "Reds")

# > Identify Signaling Roles  ---------------------------
# cc_mp24 <- readRDS('cc_mp24_08212021.rds')

# Compute the network centrality scores
cc_mp24 <- netAnalysis_computeCentrality(cc_mp24, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cc_mp24, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

saveRDS(cc_mp24, 'cc_mp24_network_08212021.rds')

rm(cc_mp24)

# Combined Signaling  ---------------------------
# Generate combined CellChat object
cc_naive <- readRDS('cc_naive_network_08212021.rds')
cc_p4 <- readRDS('cc_p4_network_08212021.rds')
cc_mp4 <- readRDS('cc_mp4_network_08212021.rds')
cc_p24 <- readRDS('cc_p24_network_08212021.rds')
cc_mp24 <- readRDS('cc_mp24_network_08212021.rds')

trtList <- c("naive", "p4", "mp4", "p24", "mp24")

object.list <- list(cc_naive,
                    cc_p4,
                    cc_mp4,
                    cc_p24,
                    cc_mp24)
cellchat <- mergeCellChat(object.list, add.names = trtList)

# Calculate basic interaction metrics:
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5), measure = "weight")
gg1 + gg2
ggsave('interactions_metrics_compared.png', width = 8, height = 4)

# Will compare 1,2 and 1,3 for 4 hour treatments
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", top = 0.1)
#netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", top = 0.1)
netVisual_diffInteraction(cellchat, comparison = c(1,3), weight.scale = T, measure = "weight", top = 0.1)

# Will compare 1,4 and 1,5 for 24 hour treatments
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1,4), weight.scale = T, measure = "weight", top = 0.1)
#netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", top = 0.1)
netVisual_diffInteraction(cellchat, comparison = c(1,5), weight.scale = T, measure = "weight", top = 0.1)

# Will compare 1,4 and 4,5 for 24 hour treatments (figure 1/2)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1,4), weight.scale = T, measure = "weight", top = 0.2)
#netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", top = 0.1)
netVisual_diffInteraction(cellchat, comparison = c(4,5), weight.scale = T, measure = "weight", top = 0.2)


# Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + NoLegend()
}
ggsave('sources_and_targets.png', plot = gg[[1]] | gg[[2]] | gg[[3]] | gg[[4]] | gg[[5]],
       width = 4*5, height = 4)





# old ----
interest <- c('aCap', 'N A', 'N B', 'N C', 'AM', 'Adv Fib')

for (i in 1:length(interest)){
  gg <- list()
  for (j in c(2,3,4,5)){
    gg[[j-1]] <- 
      netAnalysis_signalingChanges_scatter(cellchat, 
                                           idents.use = interest[[i]],
                                           title = trtList[[j]], 
                                           comparison = c(1, j))
  }
  p <- gg[[1]] + NoLegend() | gg[[2]] + NoLegend()| gg[[3]] + NoLegend()| gg[[4]]
  ggsave(paste(interest[[i]], '_signaling_changes.png', sep = ''), width = 4*4, height = 4)
}

# Identify conserved and context-specific signaling pathways:

# Functional Similarity:
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
ggsave('functionally_conserved_signaling.png', width = 5, height = 4)

netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 1)
  # giving errors
  
# Structural Similarity:
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
  # THIS LINE TAKES FOREVER... THEN ERROR:
  # Error in base::colSums(x, na.rm = na.rm, dims = dims, ...) : 
  # 'x' must be an array of at least two dimensions

  cellchat <- netEmbedding(cellchat, type = "structural")
  cellchat <- netClustering(cellchat, type = "structural")
  netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
  ggsave('structurally_conserved_signaling.png', width = 5, height = 4)

# All Compared
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3,4,5), stacked = F, do.stat = TRUE)
gg1 + gg2

# 4 hour timepoints
gg1 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", comparison = c(1,2,3), stacked = F, do.stat = TRUE)
gg1 + gg2

####



# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cc_naive@idents)
netVisual_chord_cell(cc_naive, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cc_naive, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cc_naive, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(2,9) # a numeric vector
netVisual_individual(cc_naive, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cc_naive, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#> [[1]]
# Chord diagram
netVisual_individual(cc_naive, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# MANY OUTPUTS
# Access all the signaling pathways showing significant communications
pathways.show.all <- cc_naive@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cc_naive@idents)
vertex.receiver = seq(2,9)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cc_naive, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cc_naive, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cc_naive, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cc_naive, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cc_naive, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cc_naive, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cc_naive, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
#> Note: The first link end is drawn out of sector 'MIF'.

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cc_naive, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cc_naive, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cc_naive, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

plotGeneExpression(cc_naive, signaling = "CXCL")

plotGeneExpression(cc_naive, signaling = "CXCL", enriched.only = FALSE)

## SYSTEMS ANALYSIS OF CELL-CELL COMMUNICATION NETWORK
# Compute the network centrality scores
cc_naive <- netAnalysis_computeCentrality(cc_naive, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cc_naive, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cc_naive)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cc_naive, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cc_naive, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cc_naive, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cc_naive, signaling = c("CXCL", "CCL"))

### Identify global communication patterns to explore how multiple cell types 
### and signaling pathways coordinate together
library(NMF)
library(ggalluvial)
selectK(cc_naive, pattern = "outgoing") # find where this starts to drop

nPatterns = 3
cc_naive <- identifyCommunicationPatterns(cc_naive, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cc_naive, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cc_naive, pattern = "outgoing")

# INCOMMING COMMUNICATION PATTERN
selectK(cc_naive, pattern = "incoming")

nPatterns = 4
cc_naive <- identifyCommunicationPatterns(cc_naive, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cc_naive, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cc_naive, pattern = "incoming")

## MANIFOLD AND CLASSIFICATION LEARNING ANALYSIS OF SIGNALING NETWORKS
cc_naive <- computeNetSimilarity(cc_naive, type = "functional")
cc_naive <- netEmbedding(cc_naive, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cc_naive <- netClustering(cc_naive, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cc_naive, type = "functional", label.size = 3.5)

