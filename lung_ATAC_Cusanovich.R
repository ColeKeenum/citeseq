library(Signac)
library(Seurat)
library(S4Vectors)
library(patchwork)
library(ggplot2)
set.seed(1234)

setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/ATAC")
# setwd("C:/Users/UPDATE/Desktop/COVID Lung CITE-Seq/ATAC")

sci.metadata <- read.table(
  file = "cell_metadata.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# subset to include only the LUNG data
sci.metadata <- sci.metadata[sci.metadata$tissue == 'Lung', ]
sci.counts <- readRDS(file = "atac_matrix.binary.qc_filtered.rds")
sci.counts <- sci.counts[, rownames(x = sci.metadata)] 

# checking formatting
head(sci.counts)
head(sci.metadata)

# SKIPPING CONVERSION TO MM10 COORDINATES! (MAY BE A PROBLEM?)

sci_peaks_mm9 <- StringToGRanges(regions = rownames(sci.counts), sep = c("_", "_"))

# create object and perform some basic QC filtering
sci_assay <- CreateChromatinAssay(
  counts = sci.counts,
  assay = 'sci',
  genome = "mm9",
  ranges = sci_peaks_mm9,
  project = 'sci'
)

sci <- CreateSeuratObject(
  counts = sci_assay,
  meta.data = sci.metadata,
  project = "sci",
  assay = "sci"
)
dim(sci)

all(sci$nFeature_sci == sci$nCount_sci)
# [1] TRUE
# Why? Because the matrix is binarized according to GitHub Signac devs

# NOT removing the "unknown" cell types, as in the vignette
# not doing this:  
# sci <- sci[, sci$nFeature_sci > 2000 &
#              sci$nCount_sci > 5000 &
#              !(sci$cell_label %in% c("Collisions", "Unknown"))]

sci <- sci[, sci$nFeature_sci > 2000 & !(sci$cell_label %in% c("Collisions"))]
dim(sci)

sci$tech <- 'sciATAC'
sci <- RunTFIDF(sci)
sci <- FindTopFeatures(sci, min.cutoff = 50)
sci <- RunSVD(sci)
sci <- RunUMAP(sci, reduction = 'lsi', dims = 2:30)
sci

# Bring in naive labeled cells from our CITE-Seq Experiment
setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
combined <- readRDS('combined_04_25_2021c.rds')

naive <- subset(combined, orig.ident == 'naive')
naive <- RenameIdents(naive, 'N A' = 'N', 'N B' = 'N', 'N C' = 'N')
naive$celltype <- Idents(naive)

# Try to integrate naive mouse and ATAC

p1 <- DimPlot(naive, group.by = "celltype", label = T) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(sci, group.by = "cell_label", label = T) + NoLegend() + ggtitle("ATAC")
p1 | p2


