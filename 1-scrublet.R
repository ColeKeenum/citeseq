

#  Setup ---------------------------
library(Seurat)
# For paralelizing
library(future)
plan("multicore")

setwd("C:/Users/colek/Desktop/Roy Lab/CITE-Seq Data")
path_data <- "Y:/Cole Keenum/CITE-seq files" # will be slow

# Validating Scrublet Cell Identities by Gene Counts ---------------------------
# Note 
#  > naive ---------------------------
naive <- Read10X(paste0(path_data,"/naive10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('./Scrublet outputs/scrublet_output_table_naive.csv')

dim(naive)
dim(scrublet_output)
test <- t(as.matrix(head(naive, n = 10)))[1:10,]
test2 <- as.matrix(read.table('./Scrublet outputs/test_naive.txt', sep = ' ')[,1:10])

all(test == test2) # TRUE

rownames(scrublet_output) <- colnames(naive)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_naive_colnames.csv', quote = FALSE)

# naive <- CreateSeuratObject(counts = naive, project = "naive", assay = "RNA")
# naive <- AddMetaData(naive, scrublet_output)

rm(naive, scrublet_output, test, test2)

#  > p4 ---------------------------
p4 <- Read10X(paste0(path_data,"/p4hr10xdata/p4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('./Scrublet outputs/scrublet_output_table_p4.csv')

dim(p4)
dim(scrublet_output)
test <- t(as.matrix(head(p4, n = 10)))[1:10,]
test2 <- as.matrix(read.table('./Scrublet outputs/test_p4.txt', sep = ' ')[,1:10])

all(test == test2) # TRUE

rownames(scrublet_output) <- colnames(p4)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_p4_colnames.csv', quote = FALSE)

# p4 <- CreateSeuratObject(counts = p4, project = "p4", assay = "RNA")
# p4 <- AddMetaData(p4, scrublet_output)

rm(p4, scrublet_output, test, test2)

#  > mp4 ---------------------------
mp4 <- Read10X(paste0(path_data,"/mp4hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('./Scrublet outputs/scrublet_output_table_mp4.csv')

dim(mp4)
dim(scrublet_output)
test <- t(as.matrix(head(mp4, n = 10)))[1:10,]
test2 <- as.matrix(read.table('./Scrublet outputs/test_mp4.txt', sep = ' ')[,1:10])

all(test == test2) # TRUE

rownames(scrublet_output) <- colnames(mp4)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_mp4_colnames.csv', quote = FALSE)

# mp4 <- CreateSeuratObject(counts = mp4, project = "mp4", assay = "RNA")
# mp4 <- AddMetaData(mp4, scrublet_output)

rm(mp4, scrublet_output, test, test2)

#  > p24 ---------------------------
p24 <- Read10X(paste0(path_data,"/p24hr10xdata/p24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('./Scrublet outputs/scrublet_output_table_p24.csv')

dim(p24)
dim(scrublet_output)
test <- t(as.matrix(head(p24, n = 10)))[1:10,]
test2 <- as.matrix(read.table('./Scrublet outputs/test_p24.txt', sep = ' ')[,1:10])

all(test == test2) # TRUE

rownames(scrublet_output) <- colnames(p24)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_p24_colnames.csv', quote = FALSE)

# p24 <- CreateSeuratObject(counts = p24, project = "p24", assay = "RNA")
# p24 <- AddMetaData(p24, scrublet_output)

rm(p24, scrublet_output, test, test2)

# > mp24 ---------------------------
mp24 <- Read10X(paste0(path_data,"/mp24hr10xdata/mp24hr10xdata/outs/filtered_feature_bc_matrix/"))[["Gene Expression"]]

scrublet_output <- read.csv('./Scrublet outputs/scrublet_output_table_mp24.csv')

dim(mp24)
dim(scrublet_output)
test <- t(as.matrix(head(mp24, n = 10)))[1:10,]
test2 <- as.matrix(read.table('./Scrublet outputs/test_mp24.txt', sep = ' ')[,1:10])

all(test == test2) # TRUE

rownames(scrublet_output) <- colnames(mp24)
write.csv(scrublet_output, file = './Scrublet outputs/scrublet_output_mp24_colnames.csv', quote = FALSE)

# mp24 <- CreateSeuratObject(counts = mp24, project = "mp24", assay = "RNA")
# mp24 <- AddMetaData(mp24, scrublet_output)

rm(mp24, scrublet_output, test, test2)