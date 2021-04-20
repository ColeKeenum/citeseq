suppressPackageStartupMessages({
  library(Seurat)
  library(scater)
  library(dplyr)
  library(mvoutlier)
  library(VennDiagram)
  library(harmony)
  library(clustree)
  
  library(ggpubr)
})

path_data<-"C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data"
path_output<-"C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data/QC_out"

##########  1 .Read the raw filtered ADT count data ##########

# import RDS
seuratObj <- readRDS('seuratObj_RNA_ADT.rds')
expr.mat <- seuratObj@assays$ADT@counts
# make sce object
sce <- SingleCellExperiment(assays=list(counts=expr.mat))
rm(expr.mat)

##########  2. Calculate QC metrics ##########  

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
png(paste(sample, "_num_ADT_per_cell_hist.png"))
hist(sce$nADT,
     breaks = 100, 
     main="Histogram of number of ADT per cell")
abline(v = max(sce$nADT[sce$nADT.outlier.low], na.rm=T), col = "red")
#abline(v = max(sce$nADT[sce$nADT.outlier.high], na.rm=T), col = "red")
dev.off()

png(paste(sample, "_nUMI_ADT_hist.png"))
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

ggarrange(p1, p4, p2, p5, ncol = 2, nrow = 3)
ggsave(paste(sample, "filtering_vln_plots_ADT.png", sep = ""), width = 8, height = 8)

########## 4. Remove outliers and Filter Orig. Seurat Obj ##########
dim(sce)
sce_clean <- sce[,!(sce$nUMI.outlier.high | sce$nADT.outlier.high)]
dim(sce_clean)

### Number of cells removed
ncol(sce_clean)-ncol(sce)

# save column of cellnames to apply to SeuratObj
filt_cells <- colnames(sce_clean)
seuratObj <- subset(seuratObj, cells = filt_cells)
dim(seuratObj[['ADT']])

# save for analysis
saveRDS(seuratObj, "seuratObj_RNA_ADT_V2.rds") # removed of + ADT and mt_gene outliers

rm(p1, p2, p24, p4, p5, sce, sce_clean, sce_sub, naive, p4, mp4, p24, mp24, metaData,
   metaData.filtered, adt_merged, ath_data, path_data, sample, filt_cells)
