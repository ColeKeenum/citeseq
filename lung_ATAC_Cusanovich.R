# NOTE: IT MAY BE MORE ACCURATE TO MAKE THE GENES UPPERCASE IN THE CITE-SEQ 
# SEURAT OBJECT, WILL TRY THIS INSTEAD. 

library(Signac)
library(Seurat)
library(S4Vectors)
library(patchwork)
library(ggplot2)
library(EnsDb.Mmusculus.v79)

library(stringr)
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

# Bring in activity scores as RNA assay:
activity.scores <- readRDS(file = "activity_scores.quantitative.rds")
dim(activity.scores)

activity.scores <- activity.scores[, colnames(x = sci)]
dim(activity.scores)
# 
# gene_name_convert <- function(names){
#   rik <- str_detect(names, 'RIK')
#   names[!rik] <- str_to_title(names[!rik])
#   names <- str_replace_all(names, 'RIK', 'Rik')
#   return(names)
# }
# 
# rownames(activity.scores) <- gene_name_convert(names = rownames(activity.scores))

sci[['RNA']] <- CreateAssayObject(counts = activity.scores)

# Back to doing ATAC peak stuff:
DefaultAssay(sci)

all(sci$nFeature_sci == sci$nCount_sci)
# [1] TRUE
# Why? Because the matrix is binarized according to GitHub Signac devs

# NOT removing the "unknown" cell types, as in the vignette
# not doing this:  
# sci <- sci[, sci$nFeature_sci > 2000 &
#              sci$nCount_sci > 5000 &
#              !(sci$cell_label %in% c("Collisions", "Unknown"))]

# NOT SURE IF FILTERING IS NEEDED BECAUSE THIS IS QC FILTERED ALREADY?
# sci <- sci[, sci$nFeature_sci > 2000 & !(sci$cell_label %in% c("Collisions"))]
# sci <- sci[, sci$nFeature_sci > 5000 & !(sci$cell_label %in% c("Collisions"))]
sci <- sci[, !(sci$cell_label %in% c("Collisions"))]

dim(sci)

sci$tech <- 'sciATAC'
sci <- RunTFIDF(sci)
sci <- FindTopFeatures(sci, min.cutoff = 50)
sci <- RunSVD(sci)
sci <- RunUMAP(sci, reduction = 'lsi', dims = 2:30)
sci

DefaultAssay(sci) <- 'RNA'

# Bring in naive labeled cells from our CITE-Seq Experiment
setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
combined <- readRDS('combined_04_25_2021c.rds')

naive <- subset(combined, orig.ident == 'naive')
naive <- RenameIdents(naive, 'N A' = 'N', 'N B' = 'N', 'N C' = 'N')
naive$celltype <- Idents(naive)

rna_counts <- naive@assays[["RNA"]]@counts
rownames(rna_counts) <- toupper(rownames(rna_counts))
naive[['upRNA']] <- CreateAssayObject(counts = rna_counts)

DefaultAssay(naive) <- 'upRNA'
naive <- NormalizeData(naive)
naive <- FindVariableFeatures(naive, nfeatures = 5000)
naive <- ScaleData(naive)

# Try to integrate naive mouse and ATAC
# Unsure if gene activities are named the way my features are...
# Will try to do an analysis with ATAC annotations

DefaultAssay(sci) <- 'RNA'
sci <- NormalizeData(sci, scale.factor = 1e6)

# Visualizations before I integrate:
p1 <- DimPlot(naive, group.by = "celltype", label = T) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(sci, group.by = "cell_label", label = T) + NoLegend() + ggtitle("ATAC")
p1 | p2

# Checking to see if this makes any sense
FeaturePlot(sci, features = "CD3E") # it does

genes.use <- intersect(x = VariableFeatures(object = naive)[1:5000], y = rownames(x = sci))
# need to uppercase this stuff...

# Identify anchors
transfer.anchors <- FindTransferAnchors(
  reference = naive,
  query = sci,
  reduction = 'cca',
  npcs = 40, # will give warning otherwise
  dims = 1:40 # they used 40 in the other set
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = naive$celltype,
  weight.reduction = sci[['lsi']],
  dims = 2:30
)

sci <- AddMetaData(object = sci, metadata = predicted.labels)

p1 <- DimPlot(naive, group.by = "celltype", label = T) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(sci, group.by = "predicted.id", label = T, repel = T) + NoLegend() + ggtitle("ATAC")
p3 <- DimPlot(sci, group.by = "cell_label", label = T, repel = T) + NoLegend() + ggtitle("ATAC")
p1 | p2 | p3
ggsave('5000_cutoff_ATAC_naive_integration_repel.png', width = 15, height = 5)

FeaturePlot(sci, features = 'nCount_sci')

hist(sci$nCount_sci)

dim(sci)






# PROBABLY CANNOT DO THIS BECAUSE DO NOT HAVE A FRAGMENT FILE!
# Annotating to get genes from ATAC
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
Annotation(sci) <- annotations

gene.activities <- GeneActivity(sci)

genes.use <- intersect(x = VariableFeatures(object = naive), y = rownames(x = sci))
