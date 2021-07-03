# First determine doublet scores with Scrublet:
# 1-scrublet_and_QC_citeseq.R
# Then run this script.
library(dsb)
library(Seurat)
library(tidyverse)

# Test together dsb normalization:
# setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
setwd("C:/Users/UPDATE/Desktop/COVID Lung CITE-Seq")
path_data <- "Y:/Cole Keenum/CITE-seq files" # will be slow

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

# Read filtered ADT data  ---------------------------
path_data <- "Y:/Cole Keenum/CITE-seq files" # will be slow

mp4 <- Read10X(paste0(path_data,"/mp4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
mp24 <- Read10X(paste0(path_data,"/mp24hr10xdata/mp24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
p4 <- Read10X(paste0(path_data,"/p4hr10xdata/p4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
p24 <- Read10X(paste0(path_data,"/p24hr10xdata/p24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]
naive <- Read10X(paste0(path_data,"/naive10xdata/outs/filtered_feature_bc_matrix/"))[["Antibody Capture"]]

mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "ADT")
mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "ADT")
p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "ADT")
p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "ADT")
naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "ADT")

filt <- merge(naive, c(p4, p24, mp4, mp24),
              add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))
rm(mp24, mp4, naive, p4, p24)
prot_filt <- filt[['ADT']]@counts
saveRDS(prot_filt, 'prot_filt.rds')

# Apply DSB normalization to merged matrices: -----
rna <- readRDS('rna_merged.rds')
prot <- readRDS('prot_merged.rds')
prot_filt <- readRDS('prot_filt.rds')

stained_cells <- colnames(prot_filt) # length = 29389
background <- setdiff(colnames(rna), stained_cells) # length = 33945011

# md dataframe from DSB vignette:
rna_size = log10(Matrix::colSums(rna))
prot_size = log10(Matrix::colSums(prot))
ngene = Matrix::colSums(rna > 0)
mtgene = grep(pattern = "^mt-", rownames(rna), value = TRUE)
propmt = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
md$bc = rownames(md)
md$droplet_class = ifelse(test = md$bc %in% stained_cells, yes = 'cell', no = 'background')

md = md %>% dplyr::filter(rna_size > 0 & prot_size > 0 )

ggplot(md, aes(x = log10(ngene), y = prot_size )) +
  theme_bw() + 
  geom_bin2d(bins = 300) + 
  scale_fill_viridis_c(option = "C") + 
  facet_wrap(~droplet_class) 
ggsave('ADT_metrics.png', height = 5, width = 10)

cellmd = md %>% filter(droplet_class == 'cell')
plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
p1 = ggplot(cellmd, aes(x = rna_size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
p2 = ggplot(cellmd, aes(x = propmt)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
p3 = ggplot(cellmd, aes(x = log10(ngene), y = rna_size, fill = propmt )) + plot_aes
p4 = ggplot(cellmd, aes(x = ngene, y = prot_size, fill = propmt )) + plot_aes
p1+p2+p3+p4
ggsave('filtered_cells_qc.png', height = 2*4, width = 2*4)

# calculate statistical thresholds for droplet filtering. 
rna_size_min = median(cellmd$rna_size) - (3*mad(cellmd$rna_size))
rna_size_max = median(cellmd$rna_size) + (3*mad(cellmd$rna_size))
prot_size_min = median(cellmd$prot_size) - (3*mad(cellmd$prot_size))
prot_size_max = median(cellmd$prot_size) + (3*mad(cellmd$prot_size))

# filter rows based on droplet quality control metrics
positive_cells = cellmd[
  cellmd$prot_size > prot_size_min & 
    cellmd$prot_size < prot_size_max & 
    cellmd$propmt < 0.25 &  
    cellmd$rna_size > rna_size_min & 
    cellmd$rna_size < rna_size_max, ]$bc
cells_mtx_rawprot = as.matrix(prot[ , positive_cells])

length(positive_cells) # 27700

# define a vector of background droplet barcodes based on protein library size and mRNA content
background_drops = md[md$prot_size > 1.5 & md$prot_size < 2.5 & md$ngene < 100, ]$bc
negative_mtx_rawprot = as.matrix(prot[ , background_drops])

# calculate quantiles of the raw protein matrix 
d1 = data.frame(pmax = apply(cells_mtx_rawprot, 1, max)) %>% 
  rownames_to_column('prot') %>% arrange(pmax) %>% head() 

# remove non staining CD34 protein 
prot_names = rownames(cells_mtx_rawprot)
cells_mtx_rawprot = 
  cells_mtx_rawprot[!(prot_names == 'DR3-TotalA' | prot_names == 'CD22-TotalA' | 
                      prot_names == "Ly49d-TotalA" | prot_names == "IL21-TotalA" |
                      prot_names == "CD115-TotalA"), ]
dim(cells_mtx_rawprot)
negative_mtx_rawprot = negative_mtx_rawprot[!(prot_names == 'DR3-TotalA' | prot_names == 'CD22-TotalA' | 
                                                prot_names == "Ly49d-TotalA" | prot_names == "IL21-TotalA" |
                                                prot_names == "CD115-TotalA"), ]
dim(negative_mtx_rawprot)

# Adding Scrublet doublets to md:
scrub_naive <- read.csv('./Scrublet outputs/scrublet_output_naive_colnames.csv', row.names = 'X')
scrub_p4 <- read.csv('./Scrublet outputs/scrublet_output_p4_colnames.csv', row.names = 'X')
scrub_p24 <- read.csv('./Scrublet outputs/scrublet_output_p24_colnames.csv', row.names = 'X')
scrub_mp4 <- read.csv('./Scrublet outputs/scrublet_output_mp4_colnames.csv', row.names = 'X')
scrub_mp24 <- read.csv('./Scrublet outputs/scrublet_output_mp24_colnames.csv', row.names = 'X')

rownames(scrub_naive) <- paste('naive_', rownames(scrub_naive), sep = '')
rownames(scrub_p4) <- paste('p4_', rownames(scrub_p4), sep = '')
rownames(scrub_p24) <- paste('p24_', rownames(scrub_p24), sep = '')
rownames(scrub_mp4) <- paste('mp4_', rownames(scrub_mp4), sep = '')
rownames(scrub_mp24) <- paste('mp24_', rownames(scrub_mp24), sep = '')

scrub_md <- rbind(scrub_naive, scrub_p4, scrub_p24, scrub_mp4, scrub_mp24)
write.csv(scrub_md, 'scrublet_metadata.csv')

scrub_md <- filter(scrub_md, predicted_doublet == 'False')

singlets <- intersect(colnames(cells_mtx_rawprot), rownames(scrub_md))
length(singlets)

dim(cells_mtx_rawprot)
cells_mtx_rawprot <- cells_mtx_rawprot[, singlets]
dim(cells_mtx_rawprot)

# dim(cells_mtx_rawprot)
# cells_mtx_rawprot <- select(as.data.frame(cells_mtx_rawprot), all_of(singlets))
# dim(cells_mtx_rawprot)

dim(negative_mtx_rawprot)

#normalize protein data for the cell containing droplets with the dsb method. 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot, 
  empty_drop_matrix = negative_mtx_rawprot, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = rownames(cells_mtx_rawprot)[17:19] 
)

dsb_norm_prot = apply(dsb_norm_prot, 2, function(x){ ifelse(test = x < -10, yes = 0, no = x)}) 

# Find doublets via ADT expression: ----
# filter raw protein, RNA and metadata to only include cell-containing droplets 
cells_rna = rna[ ,colnames(cells_mtx_rawprot)]
md2 = md[colnames(cells_mtx_rawprot), ]

# create Seurat object !note: min.cells is a gene filter, not a cell filter
s = Seurat::CreateSeuratObject(counts = cells_rna, meta.data = md2, 
                               assay = "RNA", min.cells = 20)

# add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
s[["CITE"]] = Seurat::CreateAssayObject(data = dsb_norm_prot)

# cluster and run umap (based directly on dsb normalized values without istype controls )
prots = rownames(s@assays$CITE@data)[c(1:16, 20:183)]

s = FindNeighbors(object = s, dims = NULL, assay = 'CITE', 
                  features = prots, k.param = 30, verbose = FALSE)

# direct graph clustering: overcluster to find doublets
s = FindClusters(object = s, resolution = 2.5, algorithm = 3, graph.name = 'CITE_snn', verbose = FALSE)

# umap for visualization only; (this is optional)
s = RunUMAP(object = s, assay = "CITE", features = prots, seed.use = 1990,
            min.dist = 0.2, n.neighbors = 30, verbose = FALSE)

# make results dataframe 
d = cbind(s@meta.data, as.data.frame(t(s@assays$CITE@data)), s@reductions$umap@cell.embeddings)

# calculate the median protein expression separately for each cluster 
adt_plot = d %>% 
  dplyr::group_by(CITE_snn_res.1) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("CITE_snn_res.1") 

# plot a heatmap of the average dsb normalized values for each cluster
p <- pheatmap::pheatmap(t(adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                   scale = 'row',
                   fontsize_row = 8, border_color = NA)
ggsave('pheatmap_ADT_cluster.pdf', plot = p, width = 8.5, height = 18)

p <- Seurat::DimPlot(s, reduction = 'umap')
ggsave('ADT_dimplot.png', plot = p, width = 5, height = 5)

Seurat::FeaturePlot(s, features = c('CD3-TotalA', 'CD8a-TotalA', 'CD4-TotalA'), min.cutoff = 0,
                    max.cutoff = 'q99')

# Quick normalization just for visualization
DefaultAssay(s) = "RNA"
s = NormalizeData(s, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = 'vst', verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# Plot all ADTs:
# Will select only the good ones and then recluster
DefaultAssay(s) <- 'CITE'
protein.to.plot <- rownames(s)
library(readxl)
ADT_with_associated_genes <- read_excel("Y:/Cole Keenum/CITE-seq files/ADT_with_associated_genes.xlsx")

bad_adt <- setdiff(ADT_with_associated_genes$Protein, protein.to.plot)
length(bad_adt)

idx <-  ADT_with_associated_genes$Protein %in% bad_adt
ADT_with_associated_genes <- ADT_with_associated_genes[!idx, ]

genes.to.plot <- ADT_with_associated_genes$Gene1

# with a 99% quartile cutoff for ADT
for (i in 1:length(x = protein.to.plot)){
  protein <- protein.to.plot[[i]]
  gene <- genes.to.plot[[i]]
  DefaultAssay(s) <- "CITE"
  p1 <- FeaturePlot(s, features = protein, cols = c("lightgrey","darkgreen"),
                    min.cutoff = 0, max.cutoff = 'q99')
  DefaultAssay(s) <- "RNA"
  p2 <- FeaturePlot(s, features = gene, reduction = 'umap')
  filename <- paste("adtFeaturePlot_", protein, ".png", sep = "")
  filename <- gsub("/",  "", filename)
  filename <- gsub("-TotalA",  "", filename)
  ggsave(filename = filename, plot = p1 | p2,  width = 10, height = 5)
}

DimPlot(s, reduction = 'umap', label = T) + NoLegend()

DefaultAssay(s) <- "CITE"
combined.markers <- FindAllMarkers(s, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0, 
                                   max.cells.per.ident = Inf)
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Note: a p-value = 0 in the following is rounded due to truncation after 1E-300
write.csv(combined.markers, file = "CITE_ADT_markers.csv", row.names = FALSE)


# Remove obvious doublets: 41, 42
dim(s)
s <- subset(s, idents = c('41', '42'), invert = T)
dim(s)

# Remove low-staining or non-staining ADT ----
adt_data <- s[['CITE']]@data

library(readxl)
bad_adts <- read_excel("Y:/Cole Keenum/CITE-seq files/2021-06-30 Bad ADT List V5.xlsx")
bad_adts <- bad_adts$`DSB-Final`
for (i in 1:length(bad_adts)){bad_adts[[i]] <- paste(bad_adts[[i]], "-TotalA", sep = "")}

setdiff(bad_adts, rownames(adt_data))
intersect(bad_adts, rownames(adt_data))

idx <- match(rownames(adt_data), bad_adts)
idx <- which(is.na(idx))

dim(adt_data)
adt_data <- adt_data[-which(rownames(adt_data) %in% bad_adts), ]
dim(adt_data)

# Apply clustering with removed poor ADTs ----
# re-creating ADT assay
s[["CITE"]] = Seurat::CreateAssayObject(data = adt_data)

# cluster and run umap (based directly on dsb normalized values without istype controls )
prots = rownames(s@assays$CITE@data)[c(1:14, 18:120)]

s = FindNeighbors(object = s, dims = NULL, assay = 'CITE', 
                  features = prots, k.param = 30, verbose = FALSE)

# direct graph clustering: overcluster to find doublets
s = FindClusters(object = s, resolution = 1, algorithm = 3, graph.name = 'CITE_snn', verbose = FALSE)

# umap for visualization only; (this is optional)
s = RunUMAP(object = s, assay = "CITE", features = prots, seed.use = 1990,
            min.dist = 0.2, n.neighbors = 30, verbose = FALSE)

# make results dataframe 
d = cbind(s@meta.data, as.data.frame(t(s@assays$CITE@data)), s@reductions$umap@cell.embeddings)

# calculate the median protein expression separately for each cluster 
adt_plot = d %>% 
  dplyr::group_by(CITE_snn_res.1) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("CITE_snn_res.1") 

# plot a heatmap of the average dsb normalized values for each cluster
p <- pheatmap::pheatmap(t(adt_plot), 
                        color = viridis::viridis(25, option = "B"), 
                        scale = 'row',
                        fontsize_row = 8, border_color = NA)
ggsave('pheatmap_ADT_cluster_120.pdf', plot = p, width = 8.5, height = 18)

p <- Seurat::DimPlot(s, reduction = 'umap')
ggsave('ADT_dimplot_120.png', plot = p, width = 5, height = 5)

dsb_singlets <- s[['CITE']]@data

saveRDS(dsb_singlets, 'dsb_singlets.rds')
