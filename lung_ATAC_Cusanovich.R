library(Signac)
library(Seurat)
library(S4Vectors)
library(patchwork)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
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

all(sci$nFeature_sci == sci$nCount_sci)
# [1] TRUE
# Why? Because the matrix is binarized according to GitHub Signac devs

# NOT removing the "unknown" cell types, as in the vignette
# not doing this:  
# sci <- sci[, sci$nFeature_sci > 2000 &
#              sci$nCount_sci > 5000 &
#              !(sci$cell_label %in% c("Collisions", "Unknown"))]

# NOT SURE IF FILTERING IS NEEDED BECAUSE THIS IS QC FILTERED ALREADY?
sci <- sci[, sci$nFeature_sci > 2000 & !(sci$cell_label %in% c("Collisions"))]

dim(sci)

sci$tech <- 'sciATAC'
sci <- RunTFIDF(sci)
sci <- FindTopFeatures(sci, min.cutoff = 50)
sci <- RunSVD(sci)
sci <- RunUMAP(sci, reduction = 'lsi', dims = 2:30)
sci

# sci <- subset(sci, features = rownames(sci)[rowSums(sci) == 0], invert = T)

# Annotating to get genes from ATAC
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm9"
Annotation(sci) <- annotations

gene.activities <- GeneActivity(sci)

genes.use <- intersect(x = VariableFeatures(object = naive), y = rownames(x = sci))

DefaultAssay(sci) <- 'RNA'

# Bring in naive labeled cells from our CITE-Seq Experiment
setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
combined <- readRDS('combined_04_25_2021c.rds')

naive <- subset(combined, orig.ident == 'naive')
naive <- RenameIdents(naive, 'N A' = 'N', 'N B' = 'N', 'N C' = 'N')
naive$celltype <- Idents(naive)

DefaultAssay(naive) <- 'SCT'

# Try to integrate naive mouse and ATAC
# Unsure if gene activities are named the way my features are...
# Will try to do an analysis with ATAC annotations


# Identify anchors
transfer.anchors <- FindTransferAnchors()


# Visualizations:
p1 <- DimPlot(naive, group.by = "celltype", label = T) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(sci, group.by = "cell_label", label = T) + NoLegend() + ggtitle("ATAC")
p1 | p2






# feat <- rownames(naive@assays[["SCT"]]@scale.data)
# class(feat)


# Bring in activity scores as RNA assay:
activity.scores <- readRDS(file = "activity_scores.quantitative.rds")
dim(activity.scores)

activity.scores <- activity.scores[, colnames(x = sci)]
dim(activity.scores)

sci[['RNA']] <- CreateAssayObject(counts = activity.scores)

# Back to doing ATAC peak stuff:
DefaultAssay(sci)
