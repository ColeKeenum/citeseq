# RT2 PCR Analysis of Adjuvants on Mouse Lung at 4 and 24 h
# By: Cole Keenum; mkeenum3@gatech.edu

library(dplyr)

# Formatting RT2 Profiler Data for STRING ----
# Upregulated data only

# Read in fold regulation data
df4 <- read.csv('4h_fold_regulation.csv') %>% as_tibble()
hist(rbind(df4$Group.1.Fold.Regulation,
           df4$Group.2.Fold.Regulation,
           df4$Group.3.Fold.Regulation,
           df4$Group.4.Fold.Regulation))

df24 <- read.csv('24h_fold_regulation.csv' ) %>% as_tibble()
hist(rbind(df24$Group.1.Fold.Regulation,
           df24$Group.2.Fold.Regulation,
           df24$Group.3.Fold.Regulation,
           df24$Group.4.Fold.Regulation))

# Split the dataframes and place into a named list:
df_list <- vector(mode = 'list', length = 8)
j = 1
# Add in 4h data
for (i in seq(from = 2, to = dim(df4)[2], by = 2)){
  df_list[[j]] <- df4[,c(1,i,i+1)]
  j = j + 1
}

# Add in 24h data
for (i in seq(from = 2, to = dim(df4)[2], by = 2)){
  df_list[[j]] <- df24[,c(1,i,i+1)]
  j = j + 1
}

# Duplicate for downregulated genes
df_list_down <- df_list

# QC out 'C' quality genes and filter by fold regulation >1.5 or <1.5:
# Upregulated
for (i in 1:length(df_list)){
  df_list[[i]] <- 
    filter(df_list[[i]], across(3, ~ . != 'C') & across(2, ~ . > 1.5))
  colnames(df_list[[i]]) <- c('gene', 'Fold Regulation', 'Comments')
  df_list[[i]]$Comments <- NULL
}

# Downregulated
for (i in 1:length(df_list_down)){
  df_list_down[[i]] <- 
    filter(df_list_down[[i]], across(3, ~ . != 'C') & across(2, ~ . < 1.5))
  colnames(df_list_down[[i]]) <- c('gene', 'Fold Regulation', 'Comments')
  df_list_down[[i]]$Comments <- NULL
}

# Put gene lists into a dataframe and export as csv
rem <- function(a_df){a_df$`Fold Regulation` <- NULL; return(a_df)}
test <- lapply(df_list, as.data.frame) %>% lapply(rem)
nrow = 0
for (i in 1:length(df_list)){
  if(dim(test[[i]])[1] > nrow){
    nrow <- dim(test[[i]])[1]
  }
}

df <- data.frame(matrix(nrow = nrow, ncol = 8))
for (i in 1:length(df_list)){
  df[1:length(df_list[[i]]$gene), i] <- df_list[[i]]$gene
}

colnames(df) <- c('mp4', 'p4', 'rp4', 'cp4', 'mp24', 'p24', 'rp24', 'cp24')

# Export data:
write.csv

# Input into stringdb at: https://string-db.org/

# Formatting full RT PCR data for heatmap ----
# Read in fold regulation data
df4 <- read.csv('4h_fold_regulation.csv') %>% as_tibble()
df24 <- read.csv('24h_fold_regulation.csv' ) %>% as_tibble()

fc_func <- function(df){
  # function to convert data to fold change and to correct for 'C' quality genes
  for (i in seq(1, dim(df)[1])){
    for (j in seq(2, dim(df)[2], 2)){
      # Turn fold regulation to fold change
      if (df[i,j] < 0){
        df[i,j] <- -1/df[i,j]
      }
      
      # Put C quality columns as FC = 1
      if (df[i, j+1] == 'C'){
        df[i, j] == 1
      }
    }
  }
  # Remove comments columns
  df[, c(3,5,7,9)] <- NULL
  return(df)
}

df4 <- fc_func(df4)
colnames(df4) <- c('gene', 'MP4', 'P4', 'RP4', 'CP4')
df24 <- fc_func(df24)
colnames(df24) <- c('gene', 'MP24', 'P24', 'RP24', 'CP24')

# Export
table_fc <- merge(df4, df24)
write.csv(table_fc, file = 'RT2_Profiler_FC.csv', na = '', row.names = F)
table_fc[,2:9] <- log2(table_fc[,2:9])
write.csv(table_fc, file = 'RT2_Profiler_log2FC.csv', na = '', row.names = F)

# Plotting heatmap ----
library(ComplexHeatmap)
library(circlize)

# Formatting data for heatmap
rtpcr <- as.matrix(table_fc[,2:9])
rownames(rtpcr) <- table_fc[,1]

# Taking only p and mp groups
rtpcr <- rtpcr[, c(2,1,6,5)]

# Filtering genes with at least abs(log2FC) > 2.5 in one sample
filt <- abs(rtpcr) < 2.5
temp <- matrix(nrow = dim(filt)[1], ncol = 1)
for (i in 1:dim(filt)[1]){
  temp[i] <- all(filt[i,])
}
rownames(rtpcr)[temp]
rtpcr <- rtpcr[-which(temp), ]

# Color with log2FC cutoff of 5 for visualization
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
col_fun(seq(-5, 5))

Heatmap(rtpcr, name = "mat", col = col_fun, row_names_gp = gpar(fontsize = 5))
# Exported manually in R with .pdf width = 3 in, height = 8.5 in

# Formatting DAVID results .tsv files ----
p4 <- read.table(file = 'p4_DAVID.tsv', sep = '\t', header = TRUE)
mp4 <- read.table(file = 'mp4_DAVID.tsv', sep = '\t', header = TRUE)
p24 <- read.table(file = 'p24_DAVID.tsv', sep = '\t', header = TRUE)
mp24 <- read.table(file = 'mp24_DAVID.tsv', sep = '\t', header = TRUE)

p4$trt <- 'p4'
mp4$trt <- 'mp4'
p24$trt <- 'p24'
mp24$trt <- 'mp24'

david <- rbind(p4, mp4, p24, mp24)

# Exporting to one csv file:
write.csv(david, file = 'david_resutls.csv', row.names = FALSE)

# Plotting on KEGG pathways ----
# Load in packages
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pathview)

# Convert gene symbols to ENTREZ IDs for use in Pathview
df <- read.csv('RT2_Profiler_log2FC.csv')
genes <- df$gene
ids <- mapIds(org.Mm.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", 
              multiVals="first")

# Coerce P24 log2FC data into a named numeric vector
p24 <- df$P24
names(p24) <- ids

# Plot on KEGG pathways
pv.out <- pathview(gene.data = p24, pathway.id = "04620",
                   species = "mmu", limit = list(gene=1, cpd=1),
                   out.suffix = "p24_limit1")
pv.out <- pathview(gene.data = p24, pathway.id = "04620",
                   species = "mmu", limit = list(gene=5, cpd=5),
                   out.suffix = "p24_limit5")





