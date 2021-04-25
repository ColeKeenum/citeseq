library(Signac)
library(Seurat)
library(S4Vectors)
library(patchwork)
set.seed(1234)

#setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/ATAC")
setwd("C:/Users/UPDATE/Desktop/COVID Lung CITE-Seq/ATAC")

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

hist(sci$nFeature_sci)
hist(sci$nCount_sci)

all(sci$nFeature_sci == sci$nCount_sci)
# [1] TRUE
# WHYY?????

# NOT removing the "unknown" cell types, as in the vignette
# not doing this:  
# sci <- sci[, sci$nFeature_sci > 2000 &
#              sci$nCount_sci > 5000 &
#              !(sci$cell_label %in% c("Collisions", "Unknown"))]

sci <- sci[, !(sci$cell_label %in% c("Collisions"))]
dim(sci)

####### trying to figure out nCount and nFeature -----
  
sci.metadata <- read.table(
  file = "cell_metadata.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# subset to include only the brain data
sci.metadata <- sci.metadata[sci.metadata$tissue == 'PreFrontalCortex', ]
sci.counts <- readRDS(file = "atac_matrix.binary.qc_filtered.rds")
sci.counts <- sci.counts[, rownames(x = sci.metadata)]

sci_peaks_mm9 <- StringToGRanges(regions = rownames(sci.counts), sep = c("_", "_"))
mm9_mm10 <- rtracklayer::import.chain("mm9ToMm10.over.chain")
sci_peaks_mm10 <- rtracklayer::liftOver(x = sci_peaks_mm9, chain = mm9_mm10)
names(sci_peaks_mm10) <- rownames(sci.counts)

# discard any peaks that were mapped to >1 region in mm10
correspondence <- elementNROWS(sci_peaks_mm10)
sci_peaks_mm10 <- sci_peaks_mm10[correspondence == 1]
sci_peaks_mm10 <- unlist(sci_peaks_mm10)
sci.counts <- sci.counts[names(sci_peaks_mm10), ]

# rename peaks with mm10 coordinates
rownames(sci.counts) <- GRangesToString(grange = sci_peaks_mm10)

# create object and perform some basic QC filtering
sci_assay <- CreateChromatinAssay(
  counts = sci.counts,
  assay = 'sci',
  genome = "mm10",
  ranges = sci_peaks_mm10,
  project = 'sci'
)

sci <- CreateSeuratObject(
  counts = sci_assay,
  meta.data = sci.metadata,
  project = "sci",
  assay = "sci"
)

all(sci$nCount_sci == sci$nFeature_sci)
# [1] TRUE

sci <- sci[, sci$nFeature_sci > 2000 &
             sci$nCount_sci > 5000 &
             !(sci$cell_label %in% c("Collisions", "Unknown"))]

sci$tech <- 'sciATAC'
sci <- RunTFIDF(sci)
sci <- FindTopFeatures(sci, min.cutoff = 50)
sci <- RunSVD(sci)
sci <- RunUMAP(sci, reduction = 'lsi', dims = 2:30)
sci
  