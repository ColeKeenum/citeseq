# Load packages / setup ---------------------------
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")

# Split up treatment matrices ---------------------------
combined <- readRDS('combined_04_25_2021c.rds')
Idents(combined) <- 'celltype'
combined[['RNA']] <- NULL # save memory
combined[['integratedSCT_']] <- NULL
combined[['integratedADT_']] <- NULL

combined <- RenameIdents(combined, 'cDC 2' = 'NC Mono')
combined <- RenameIdents(combined, 'NC Mono ' = 'NC Mono') # error in naming earlier

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

Idents(combined) <- 'celltype'
combined <- subset(combined, idents = cellList)

seurat_list <- SplitObject(combined, split.by = 'orig.ident')

saveRDS(seurat_list, 'split_slim_combined.rds')

rm(combined)

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
  
  # For later: consider using population.size = TRUE for probability computation
  # Consider setting raw.use = FALSE to alleviate the dropout issue
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
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
saveRDS(cc_naive, 'cc_naive_05042021.rds')

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

rm(cc_naive)

# P4 ---------------------------
data.p4 <- GetAssayData(object = seurat_list$p4, assay = 'SCT', slot = "data")
meta.p4 <- seurat_list$p4@meta.data[c('celltype')]
cc_p4 <- preprocess(data.p4, meta.p4)

rm(data.p4, meta.p4); saveRDS(cc_p4, 'cc_p4_05042021.rds')

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

rm(cc_p4)

# MP4 ---------------------------
data.mp4 <- GetAssayData(object = seurat_list$mp4, assay = 'SCT', slot = "data")
meta.mp4 <- seurat_list$mp4@meta.data[c('celltype')]
cc_mp4 <- preprocess(data.mp4, meta.mp4)

rm(data.mp4, meta.mp4); saveRDS(cc_mp4, 'cc_mp4_05042021.rds')

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

rm(cc_mp4)

# P24 ---------------------------
data.p24 <- GetAssayData(object = seurat_list$p24, assay = 'SCT', slot = "data")
meta.p24 <- seurat_list$p24@meta.data[c('celltype')]
cc_p24 <- preprocess(data.p24, meta.p24)

rm(data.p24, meta.p24); saveRDS(cc_p24, 'cc_p24_05042021.rds')

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

rm(cc_p24)

# MP24 ---------------------------
data.mp24 <- GetAssayData(object = seurat_list$mp24, assay = 'SCT', slot = "data")
meta.mp24 <- seurat_list$mp24@meta.data[c('celltype')]
cc_mp24 <- preprocess(data.mp24, meta.mp24)

rm(data.mp24, meta.mp24); saveRDS(cc_mp24, 'cc_mp24_05042021.rds')

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

rm(cc_mp24)




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

# MP4-hr OLD ---------------------------
data.input <- GetAssayData(object = seurat_list$mp4, assay = 'SCT', slot = "data")

meta <- seurat_list$mp4@meta.data[c('celltype')]
meta$celltype <- paste(meta$celltype) # sketchy, doing this to remove the levels

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

#cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# > ALL signaling ---------------------------

#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # full CellChatDB

cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

# this will take a little time
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse) # used to be be PPI.human

# For later: consider using population.size = TRUE for probability computation
# Consider setting raw.use = FALSE to alleviate the dropout issue
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# infer cell-cell communication at a signaling pathway level:
cellchat <- computeCommunProbPathway(cellchat)

# calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# visualize
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
ggsave('net_test.pdf', width = 9.5, height = 11)

# analyze signaling from each cell group
mat <- cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE) # broke from me
par(mar=c(1,1,1,1), xpd=TRUE) # worked for me
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# > Emphasis on CXCL pathway  ---------------------------
pathways.show <- c("CXCL") 

pathways.show <- c("CCL") 

pathways.show <- cellchat@netP$pathways

# Hierarchy plot
# Alveolar-side cells: AT1, AT2 (two types?), Alv. Macrophages.
vertex.receiver = seq(2,9) # indices of intersted cells

vertex.receiver = seq(33,38)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(2,9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#> [[1]]
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.



# MP24 OLD---------------------------
data.input <- GetAssayData(object = seurat_list$mp24, assay = 'SCT', slot = "data")

meta <- seurat_list$mp24@meta.data[c('celltype')]
meta$celltype <- paste(meta$celltype) # sketchy, doing this to remove the levels

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

#cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# > ALL signaling ---------------------------

#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # full CellChatDB

cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

# this will take a little time
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse) # used to be be PPI.human

# For later: consider using population.size = TRUE for probability computation
# Consider setting raw.use = FALSE to alleviate the dropout issue
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# infer cell-cell communication at a signaling pathway level:
cellchat <- computeCommunProbPathway(cellchat)

# calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# visualize
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
ggsave('net_test.pdf', width = 9.5, height = 11)

# analyze signaling from each cell group
mat <- cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE) # broke from me
par(mar=c(1,1,1,1), xpd=TRUE) # worked for me
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# > Emphasis on CXCL pathway  ---------------------------
pathways.show <- c("CXCL") 

pathways.show <- c("CCL") 

pathways.show <- cellchat@netP$pathways

# Hierarchy plot
# Alveolar-side cells: AT1, AT2 (two types?), Alv. Macrophages.
vertex.receiver = seq(2,9) # indices of intersted cells

vertex.receiver = seq(33,38)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(2,9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#> [[1]]
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

