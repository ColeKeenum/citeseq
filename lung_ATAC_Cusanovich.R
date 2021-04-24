library(Signac)
library(Seurat)
library(S4Vectors)
library(patchwork)
set.seed(1234)

setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/ATAC")

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

# cannot allocate vector

