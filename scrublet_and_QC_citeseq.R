#  Setup ---------------------------
library(Seurat)
# For paralelizing
library(future)
plan("multicore")

setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")

#  SCRUBLET CONFIRMATION ---------------------------
#  > naive ---------------------------
naive <- Read10X(paste0(getwd(),"/naive10xdata/naive10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('scrublet_output_table_naive.csv')

dim(naive)
dim(scrublet_output)
test <- as.matrix(head(naive, n = 20))
View(test)

rownames(scrublet_output) <- colnames(naive)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_naive_colnames.csv', quote = FALSE)

naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "RNA")
naive <- AddMetaData(naive, scrublet_output)

#  > p4 ---------------------------
p4 <- Read10X(paste0(getwd(),"/p4hr10xdata/p4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('scrublet_output_table_p4.csv')

dim(p4)
dim(scrublet_output)
test <- as.matrix(head(p4, n = 20))
View(test)

rownames(scrublet_output) <- colnames(p4)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_p4_colnames.csv', quote = FALSE)

p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "RNA")
p4 <- AddMetaData(p4, scrublet_output)

#  > mp4 ---------------------------
mp4 <- Read10X(paste0(getwd(),"/mp4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('scrublet_output_table_mp4.csv')

dim(mp4)
dim(scrublet_output)
test <- as.matrix(head(mp4, n = 20))
#View(test)

rownames(scrublet_output) <- colnames(mp4)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_mp4_colnames.csv', quote = FALSE)

mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "RNA")
mp4 <- AddMetaData(mp4, scrublet_output)

#  > p24 ---------------------------
p24 <- Read10X(paste0(getwd(),"/p24hr10xdata/p24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('scrublet_output_table_p24.csv')

dim(p24)
dim(scrublet_output)
test <- as.matrix(head(p24, n = 20))
#View(test)

rownames(scrublet_output) <- colnames(p24)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_p24_colnames.csv', quote = FALSE)

p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "RNA")
p24 <- AddMetaData(p24, scrublet_output)

# > mp24 ---------------------------
mp24 <- Read10X(paste0(getwd(),"/mp24hr10xdata/mp24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('scrublet_output_table_mp24.csv')

dim(mp24)
dim(scrublet_output)
test <- as.matrix(head(mp24, n = 20))
#View(test)

rownames(scrublet_output) <- colnames(mp24)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_mp24_colnames.csv', quote = FALSE)

mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "RNA")
mp24 <- AddMetaData(mp24, scrublet_output)

# > Merge seurat objects  ---------------------------
rna_merged <- merge(naive, c(p4, p24, mp4, mp24),
                    add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))

saveRDS(rna_merged, 'rna_merged_scrub_metadata_preQC.rds')

rm(naive, p4, mp4, p24, mp24)

#  mt-Gene Outlier Analysis (MAD = 2.25) ---------------------------
library(scater)
library(mvoutlier)
library(VennDiagram)
library(ggpubr)

sce <- as.SingleCellExperiment(rna_merged)

# > Calcualte mt gene outlier ---------------------------

##### Get mitochondrial genes
mito.genes <- grep("^mt-", rownames(sce))
length(mito.genes) # 13

##### Calculate QC metrics per cell
## percentage expressed mitochondrial genes 
sce$percent.mito <- (Matrix::colSums(counts(sce)[mito.genes, ])*100)/Matrix::colSums(counts(sce))
# number of expressed genes
sce$nGene <- apply(counts(sce),  2,  function(x) length(x[x > 0]))
# total UMI counts (library size)
sce$nUMI<-apply(counts(sce),  2,  sum)
sce$staticNr <- 1
dim(colData(sce))
colnames(colData(sce))

# > Get mt gene outliers ---------------------------
## Outliers are determined based on the median absolute deviation (MAD). It is possible to modify the nmads parameter 
## (minimum number of MADs away from median required for a value to be called an outlier)

sample <- "merged"

## Mitochondrial gene proportions 
sce$mito.outlier.high <- isOutlier(sce$percent.mito, nmads=2.25, type="higher", log=TRUE)
sum(sce$mito.outlier.high)


png(paste(sample, "_percent_mito_hist_MAD2.25.png"))
hist(sce$percent.mito,
     breaks = 100,
     main="Histogram of % mito genes per cell")
abline(v = min(sce$percent.mito[sce$mito.outlier.high]), col = "red")
dev.off()

###### Create violin plots 

### Before filtering
metaData=as.data.frame(colData(sce))
p3<-ggplot(metaData, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell")

### After filtering
metaData.filtered<-metaData[! metaData$mito.outlier.high ,]
#metaData.filtered<-metaData[! (metaData$nUMI.outlier.low | metaData$nGene.outlier.low | metaData$mito.outlier.high) ,]

p6<-ggplot(metaData.filtered, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell after filtering")

ggarrange(p3, p6, ncol = 2, nrow = 1)
ggsave(paste(sample, "filtering_vln_plots_v2.png", sep = ""), width = 8, height = 8)

###### Remove outliers
dim(sce)
sce_clean <- sce[,!sce$mito.outlier.high]
dim(sce_clean)

### Number of cells removed
ncol(sce_clean)-ncol(sce)

# > Seurat obj prep ---------------------------

seuratObj <- as.Seurat(sce_clean)
seuratObj$staticNr <- NULL
seuratObj$mito.outlier.high <- NULL
seuratObj$percent.mito <- NULL

saveRDS(seuratObj, "seuratObj_scrublet_mtQC.rds")

rm(p3, p6, rna_merged, sce, sce_clean, scrublet_output, test, mito.genes, 
   sample, metaData, metaData.filtered)

# ADT Library Outlier Analysis ---------------------------
# > Reading in ADT reads ---------------------------
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

adt_merged <- merge(naive, c(p4, p24, mp4, mp24),
                    add.cell.ids = c("naive", "p4", "p24", "mp4", "mp24"))
dim(adt_merged)

expr.mat <- adt_merged@assays$ADT@counts
# make sce object
sce <- SingleCellExperiment(assays=list(counts=expr.mat))
rm(expr.mat)

# > Calculate QC metrics ---------------------------

##### Calculate QC metrics per cell
sce$nADT<-apply(counts(sce),  2,  function(x) length(x[x > 0])) # number of expressed genes
sce$nUMI<-apply(counts(sce),  2,  sum) # total UMI counts (library size) # ALREADY DONE WITH RNA
sce$staticNr<-1
dim(colData(sce))
colnames(colData(sce))

########## 3. Get outliers ##########
## Outliers are determined based on the median absolute deviation (MAD). It is possible to modify the nmads parameter 
## (minimum number of MADs away from median required for a value to be called an outlier)

sample <- "merged"

##### Aim: remove outlier cells for low library size, low number of expressed genes and high mitochondrial gene proportions
## Number of expressed adt
sce$nADT.outlier.low <- isOutlier(sce$nADT, nmads=3, type="lower", log=TRUE)
sum(sce$nADT.outlier.low)

sce$nADT.outlier.high <- isOutlier(sce$nADT, nmads=3, type="higher", log=TRUE)
sum(sce$nADT.outlier.high)

## UMI counts per cell
sce$nUMI.outlier.low <- isOutlier(sce$nUMI, nmads=3, type="lower", log=TRUE)
sum(sce$nUMI.outlier.low )

sce$nUMI.outlier.high <- isOutlier(sce$nUMI, nmads=3, type="higher", log=TRUE)
sum(sce$nUMI.outlier.high )

##### Create histograms
png(paste(sample, "_num_ADT_per_cell_hist_v2.png"))
hist(sce$nADT,
     breaks = 100, 
     main="Histogram of number of ADT per cell")
abline(v = max(sce$nADT[sce$nADT.outlier.low], na.rm=T), col = "red")
#abline(v = max(sce$nADT[sce$nADT.outlier.high], na.rm=T), col = "red")
dev.off()

png(paste(sample, "_nUMI_ADT_hist_v2.png"))
hist(sce$nUMI,
     breaks = 50,
     main="Histogram of total UMI counts per cell")
abline(v = max(sce$nUMI[sce$nUMI.outlier.low]), col = "red")
abline(v = max(sce$nUMI[sce$nUMI.outlier.high]), col = "red")
dev.off()

###### Create violin plots 

### Before filtering
metaData=as.data.frame(colData(sce))
p1<-ggplot(metaData, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell")
p2<-ggplot(metaData, aes(staticNr, nADT)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nADT.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of ADT per cell")

### After filtering
metaData.filtered<-metaData[! (metaData$nUMI.outlier.high | metaData$nADT.outlier.high) ,]
p4<-ggplot(metaData.filtered, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell after filtering")
p5<-ggplot(metaData.filtered, aes(staticNr, nADT)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nADT.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of ADT per cell after filtering")

ggarrange(p1, p4, p2, p5, ncol = 2, nrow = 2)
ggsave(paste(sample, "filtering_vln_plots_ADT_v2.png", sep = ""), width = 8, height = 8)

########## 4. Remove outliers and Filter Orig. Seurat Obj ##########
dim(sce)
sce_clean <- sce[,!(sce$nUMI.outlier.high | sce$nADT.outlier.high)]
dim(sce_clean)

### Number of cells removed
ncol(sce_clean)-ncol(sce)

# save column of cellnames to apply to SeuratObj
filt_cells <- colnames(sce_clean)
length(filt_cells)

#seuratObj <- readRDS('seuratObj_scrublet_mtQC.rds')
dim(seuratObj)

length(intersect(colnames(seuratObj), filt_cells))

seuratObj <- subset(seuratObj, cells = intersect(colnames(seuratObj), filt_cells))
dim(seuratObj)

sce_clean <- sce_clean[, intersect(colnames(seuratObj), filt_cells)]
dim(sce_clean)

seuratObj[['ADT']] <- CreateAssayObject(counts = sce_clean@assays@data@listData[["counts"]])
dim(seuratObj[['ADT']])



# save for analysis
saveRDS(seuratObj, "seuratObj_scrub_mtQC_adtQC.rds")
