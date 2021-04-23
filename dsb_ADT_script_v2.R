# Load Packages  ---------------------------
# Updated version of 4/15/2021
library(dsb)
library(Seurat)
library(tidyverse) # used for ggplot and %>% (tidyverse is not a dsb dependency) 

path_data<-"C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data"

# Read raw ADT data  ---------------------------
# Had to split up ADT and RNA separately due to memory issues
mp4 <- Read10X(paste0(path_data,"/mp4hr10xdata/outs/raw_feature_bc_matrix/"))[["Antibody Capture"]]
mp24 <- Read10X(paste0(path_data,"/mp24hr10xdata/mp24hr10xdata/outs/raw_feature_bc_matrix/"))[["Antibody Capture"]]
p4 <- Read10X(paste0(path_data,"/p4hr10xdata/p4hr10xdata/outs/raw_feature_bc_matrix/"))[["Antibody Capture"]]
p24 <- Read10X(paste0(path_data,"/p24hr10xdata/p24hr10xdata/outs/raw_feature_bc_matrix/"))[["Antibody Capture"]]
naive <- Read10X(paste0(path_data,"/naive10xdata/naive10xdata/outs/raw_feature_bc_matrix/"))[["Antibody Capture"]]

mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "ADT")
mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "ADT")
p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "ADT")
p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "ADT")
naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "ADT")

raw <- merge(naive, c(p4, p24, mp4, mp24),
                 add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))
rm(mp24, mp4, naive, p4, p24)

# Read10X formats the output as a list for the ADT and RNA assays: split this list 
prot <- raw[["ADT"]]@counts
saveRDS(prot, "prot_merged.rds")
rm(prot, raw)

# Read raw RNA data  ---------------------------
mp4 <- Read10X(paste0(path_data,"/mp4hr10xdata/outs/raw_feature_bc_matrix/"))[["Gene Expression"]]
mp24 <- Read10X(paste0(path_data,"/mp24hr10xdata/mp24hr10xdata/outs/raw_feature_bc_matrix/"))[["Gene Expression"]]
p4 <- Read10X(paste0(path_data,"/p4hr10xdata/p4hr10xdata/outs/raw_feature_bc_matrix/"))[["Gene Expression"]]
p24 <- Read10X(paste0(path_data,"/p24hr10xdata/p24hr10xdata/outs/raw_feature_bc_matrix/"))[["Gene Expression"]]
naive <- Read10X(paste0(path_data,"/naive10xdata/naive10xdata/outs/raw_feature_bc_matrix/"))[["Gene Expression"]]

mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "RNA")
mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "RNA")
p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "RNA")
p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "RNA")
naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "RNA")

raw <- merge(naive, c(p4, p24, mp4, mp24),
             add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))
rm(mp24, mp4, naive, p4, p24)

# Read10X formats the output as a list for the ADT and RNA assays: split this list 
rna <- raw[["RNA"]]@counts
saveRDS(rna, "rna_merged.rds")

# START HERE Read filtered ADT data  ---------------------------
path_data<-"C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data"

mp4 <- Read10X(paste0(path_data,"/mp4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
mp24 <- Read10X(paste0(path_data,"/mp24hr10xdata/mp24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
p4 <- Read10X(paste0(path_data,"/p4hr10xdata/p4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
p24 <- Read10X(paste0(path_data,"/p24hr10xdata/p24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
naive <- Read10X(paste0(path_data,"/naive10xdata/naive10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]

mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "ADT")
mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "ADT")
p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "ADT")
p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "ADT")
naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "ADT")

filt <- merge(naive, c(p4, p24, mp4, mp24),
             add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))
rm(mp24, mp4, naive, p4, p24)
prot_filt <- filt[['ADT']]@counts

filt_cells <- colnames(prot_filt)

# saveRDS(filt_cells, 'filt_cells.rds')
# filt_cells <- readRDS('filt_cells.rds')

# Input from merged raw count matrices ---------------------------
# 4/15/2021 note: using the same prot_merged and rna_merged files

prot <- readRDS("prot_merged.rds") #prot = raw$`Antibody Capture`
rna <- readRDS("rna_merged.rds") #rna = raw$`Gene Expression`

length(intersect(filt_cells, colnames(prot))) # = to length of filt_cells?
idx <- sort(match(filt_cells, colnames(prot)))
prot_background <- prot[, -idx]
all(colnames(prot_filt) == colnames(prot[ ,idx])) # check concordance

length(intersect(filt_cells, colnames(rna))) # = to length of filt_cells?
idx <- sort(match(filt_cells, colnames(rna)))
rna_background <- rna[, -idx]
all(colnames(prot_filt) == colnames(rna[ ,idx])) # check concordance

rna_filt <- rna[, idx]

# Histograms ---------------------------
# define a function that plots the histograms for rna and prot data:
plot_hists <- function(rna_mtx, prot_mtx){
  rna_size = log10(Matrix::colSums(rna_mtx))
  prot_size = log10(Matrix::colSums(prot_mtx))
  ngene = Matrix::colSums(rna_mtx > 0)
  mtgene = grep(pattern = "^mt-", rownames(rna_mtx), value = TRUE)
  propmt = Matrix::colSums(rna_mtx[mtgene, ]) / Matrix::colSums(rna_mtx)
  md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
  md$bc = rownames(md)

  p1 = ggplot(md[md$rna_size > 0, ], aes(x = rna_size)) + geom_histogram(fill = "dodgerblue") + ggtitle("RNA library size \n distribution")
  p2 = ggplot(md[md$prot_size> 0, ], aes(x = prot_size)) + geom_histogram(fill = "firebrick2") + ggtitle("Protein library size \n distribution")
  cowplot::plot_grid(p1, p2, nrow = 1)
}

# # plot background
# plot_hists(rna_filt, prot_filt)
# ggsave("filt_histogram.png", width = 10, height = 5)
# plot_hists(rna_background, prot_background)
# ggsave("background_histogram.png", width = 10, height = 5)

# Define Background ---------------------------
rna_size = log10(Matrix::colSums(rna_background))
prot_size = log10(Matrix::colSums(prot_background))
ngene = Matrix::colSums(rna_background > 0)
mtgene = grep(pattern = "^mt-", rownames(rna_background), value = TRUE)
propmt = Matrix::colSums(rna_background[mtgene, ]) / Matrix::colSums(rna_background)
md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
md$bc = rownames(md)

# define a vector of background / empty droplet barcodes based on protein library size and mRNA content

max(md$ngene) # this is interesting, impllies some cells were killed by cellranger

# background_drops = md[md$prot_size > 1.4 & md$prot_size < 2.5 & md$ngene < 80, ]$bc
rm(filt, prot, rna, propmt)
#background_drops = md[md$prot_size > 0.25, ]$bc
# 4/15/21 UP FROM 0.25 CUTOFF! Using md$ngene < 80 as well
background_drops = md[md$prot_size > 1.5 & md$ngene < 80, ]$bc 
length(background_drops)
# rm(md)
negative_mtx_rawprot = prot_background[ , background_drops] %>% as.matrix()

# define a vector of cell-containing droplet barcodes based on protein library size and mRNA content 
# positive_cells = md[md$prot_size > 2.8 & md$ngene < 3000 & md$ngene > 200 & propmt <0.2, ]$bc
# cells_mtx_rawprot = prot[ , positive_cells] %>% as.matrix()
cells_mtx_rawprot = as.matrix(prot_filt)
length(colnames(cells_mtx_rawprot))

rm(md, prot_filt, rna_filt, prot_background, rna_background, background_drops,
   filt_cells, idx, mtgene, ngene, rna_size)

# Apply DSB method ---------------------------
#normalize protein data for the cell containing droplets with the dsb method. 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot, # cell containing droplets
  empty_drop_matrix = negative_mtx_rawprot, # estimate ambient noise with the background drops 
  denoise.counts = FALSE, 
  use.isotype.control = FALSE, 
)

saveRDS(dsb_norm_prot, file = "dsb_norm_prot_V2.rds")
# dsb_norm_prot <- readRDS("dsb_norm_prot_V2.rds")
# bring in other files
# load("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/pre_dsb_norm_04052021_workspace.RData")

# create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
dsb_norm_prot <- readRDS("dsb_norm_prot_V2.rds")
s = Seurat::CreateSeuratObject(counts = rna_filt, assay = "RNA")

# add DSB normalized "dsb_norm_prot" protein data to a assay called "CITE" created in step II 
s[["CITE"]] = Seurat::CreateAssayObject(data = dsb_norm_prot)

saveRDS(s, 'dsb_seurat_obj_v2.rds')

# define Euclidean distance matrix on dsb normalized protein data (without isotype controls)
# SAVED WORKSPACE FOR BIG PC
dsb = s@assays$CITE@data
p_dist = dist(t(dsb))
p_dist = as.matrix(p_dist)

# Cluster using Seurat 
s[["p_dist"]] = Seurat::FindNeighbors(p_dist)$snn
s = Seurat::FindClusters(s, resolution = 0.5, graph.name = "p_dist")


# calculate the average of each protein separately for each cluster 
library(pheatmap)
prots = rownames(s@assays$CITE@data)
adt_plot = adt_data %>% 
  group_by(seurat_clusters) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  column_to_rownames("seurat_clusters") 
# plot a heatmap of the average dsb normalized values for each cluster
pheatmap::pheatmap(t(adt_plot), color = viridis::viridis(25, option = "B"), fontsize_row = 8, border_color = NA)



# OLD STUFF ---------------------------




# plot background
rna_size = log10(Matrix::colSums(rna_background))
prot_size = log10(Matrix::colSums(prot_background))
ngene = Matrix::colSums(rna_background > 0)
mtgene = grep(pattern = "^mt-", rownames(rna_background), value = TRUE)
propmt = Matrix::colSums(rna_background[mtgene, ]) / Matrix::colSums(rna_background)
md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
md$bc = rownames(md)

p1 = ggplot(md[md$rna_size > 0, ], aes(x = rna_size)) + geom_histogram(fill = "dodgerblue") + ggtitle("RNA library size \n distribution (Background)")
p2 = ggplot(md[md$prot_size> 0, ], aes(x = prot_size)) + geom_histogram(fill = "firebrick2") + ggtitle("Protein library size \n distribution (Background)")
cowplot::plot_grid(p1, p2, nrow = 1)

# create a metadata dataframe of simple qc stats for each droplet 
rna_size = log10(Matrix::colSums(rna))
prot_size = log10(Matrix::colSums(prot))
ngene = Matrix::colSums(rna > 0)
mtgene = grep(pattern = "^mt-", rownames(rna), value = TRUE)
propmt = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
md$bc = rownames(md)

p1 = ggplot(md[md$rna_size > 0, ], aes(x = rna_size)) + geom_histogram(fill = "dodgerblue") + ggtitle("RNA library size \n distribution")
p2 = ggplot(md[md$prot_size> 0, ], aes(x = prot_size)) + geom_histogram(fill = "firebrick2") + ggtitle("Protein library size \n distribution")
cowplot::plot_grid(p1, p2, nrow = 1)
ggsave("raw_histogram.png", width = 10, height = 5)

p1 = ggplot(md[md$rna_size > 2.5, ], aes(x = rna_size)) + geom_histogram(fill = "dodgerblue") + ggtitle("RNA library size \n distribution")
p2 = ggplot(md[md$prot_size> 2.5, ], aes(x = prot_size)) + geom_histogram(fill = "firebrick2") + ggtitle("Protein library size \n distribution")
cowplot::plot_grid(p1, p2, nrow = 1)
ggsave("zoomed_raw_histogram.png", width = 10, height = 5)
