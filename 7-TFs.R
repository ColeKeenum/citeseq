# NOTES ---------------------------
# CITE-Seq Analysis of Lung Cells w/ PUUC and MPLA+PUUC at 4 and 24 hour
# By Cole Keenum, Medical College of Georgia, mkeenum@augusta.edu
# Run after DEGs are calculated. Can be run in parallel with other analysis.
# Running ChEA3 API in R
library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd('~/Documents/citeseq-code')

# get DEGs vs. naive ----
degs <- read.csv('2021-07-29 DEGs.csv')

# convert to human orthologs for ChEA3 from mouse genes
# Basic function to convert mouse to human gene names
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  # Last modified by JAX: 2023-01-16 08:00	version

convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}

# prepare for for loop to populate with ChEA3 queries
celltype_trt <- unique(degs$cellType)
tf_list <- vector(mode = "list", length = length(celltype_trt))

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"

for (i in 1:length(celltype_trt)){
  print(i)
  print(celltype_trt[[i]])
  deg_subset <- filter(degs, cellType == celltype_trt[[i]], p_val < 0.05, avg_log2FC < -0.06) # alter the > or < as needed
  genes <- deg_subset$gene_symbol

  # convert to human ortholog gene names (vast majority are same name capitalized, with a few exceptions)
  genes <- convert_mouse_to_human(genes)
  
  payload = list(query_name = celltype_trt[[i]], gene_set = genes)
  
  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  
  #results as list of R dataframes
  tf_list[[i]] <- fromJSON(json)
}

# save the full information in case I need it later:
saveRDS(tf_list, file = 'tf_list_full_down.rds')

# trim the data
tf_meanRank <- vector(mode = "list", length = length(tf_list))
for (i in 1:length(tf_meanRank)){
  tf_meanRank[[i]] <- tf_list[[i]]$`Integrated--meanRank`
}
tf_meanRank_master <- do.call(rbind, tf_meanRank)

saveRDS(tf_meanRank_master, 'tf_meanRank_master_downreg.rds')

# bring it together
upreg <- readRDS('tf_meanRank_master_upreg.rds')
downreg <- readRDS('tf_meanRank_master_downreg.rds')

upreg$Score <- as.numeric(upreg$Score) * (1)
downreg$Score <- as.numeric(downreg$Score)  * (-1)

tfs <- combine(upreg, downreg)
tfs$Rank <- as.numeric(tfs$Rank)

write.csv(tfs, 'tfs_meanRank.csv')

# get just top 100 tfs for each
write.csv(filter(tfs, Rank <= 100), 'tfs_meanRank_top100.csv')

# visualize mean rank scores for up and downregulated genes ----
tfs <- read.csv('tfs_meanRank.csv') %>% select(-c(Library, Overlapping_Genes))
tfs <- tfs[,2:5]
celltype_trt <- unique(tfs$Query.Name)

for (c in celltype_trt){
  tf_sub <- filter(tfs, Query.Name == c) %>% select(TF, Score)
  pos_tfs <- filter(tf_sub, Score > 0); colnames(pos_tfs) <- c('TF', 'pos_score')
  neg_tfs <- filter(tf_sub, Score < 0); colnames(neg_tfs) <- c('TF', 'neg_score')
  tf_sub <- left_join(pos_tfs, neg_tfs)
  
  top10_up <- arrange(pos_tfs, pos_score)$TF[1:10]
  top10_down <- arrange(neg_tfs, desc(neg_score))$TF[1:10]
  
  # math transformations based on panel C in Bouyer et al. Nat. Comm., 2020. 
  # https://www.nature.com/articles/s41467-021-22973-9/figures/6 
  tf_sub$pos_scaled <- log(1/tf_sub$pos_score)
  tf_sub$neg_scaled <- log(1/-tf_sub$neg_score)
  rownames(tf_sub) <- tf_sub$TF
  
  ggplot(tf_sub, aes(x=neg_scaled, y=pos_scaled, label = TF)) + 
    geom_point(size = 0.4, colour = 'cornflowerblue') + theme_classic() + geom_text_repel() +
    xlim(-8, 0) + ylim(-8, 0) + 
    xlab(expression(ln('ChEA3 MeanRank Neg Enriched'^ {-1}))) +
    ylab(expression(ln('ChEA3 MeanRank Pos Enriched'^ {-1})))
  ggsave(paste(c, 'TFs_v2.png'), width = 4.5, height = 4.5)
}

# plot heatmap of main TFs-----
library(ComplexHeatmap)
library(stringr)

# load TFs
tfs <- read.csv('tfs_meanRank.csv') %>% select(-c(Library, Overlapping_Genes))
tfs <- tfs[,2:5]
celltype_trt <- unique(tfs$Query.Name)

# get top TFs, anything ranked in the top 5 for any cell type **4 h**
tfs$trt  <- str_split_fixed(tfs$Query.Name, '_', 2)[,2]
tfs <- filter(tfs, trt %in% c('p4', 'mp4'))

tops <- unique( filter(tfs, Rank <= 3, Score > 0)$TF )
tfs5 <- filter(tfs, TF %in% tops, Score > 0) %>% select(-c(Rank, trt))
tfs5$Score <- log(1/tfs5$Score)

df <- tidyr::pivot_wider(tfs5, 
                           names_from = Query.Name, 
                           values_from = Score) %>% as.data.frame()
rownames(df) <- df$TF
df <- select(df, -c('TF'))

Heatmap(df)

col_fun = circlize::colorRamp2(c(-8, -5, 0), c("black", "black", "yellow"))
Heatmap(df, name = 'Score',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))


# Same for **24 h**
# load TFs
tfs <- read.csv('tfs_meanRank.csv') %>% select(-c(Library, Overlapping_Genes))
tfs <- tfs[,2:5]
celltype_trt <- unique(tfs$Query.Name)

# get top TFs, anything ranked in the top 5 for any cell type **4 h**
tfs$trt  <- str_split_fixed(tfs$Query.Name, '_', 2)[,2]
tfs <- filter(tfs, trt %in% c('p24', 'mp24'))

tops <- unique( filter(tfs, Rank <= 7, Score > 0)$TF )
tfs5 <- filter(tfs, TF %in% tops, Score > 0) %>% select(-c(Rank, trt))
tfs5$Score <- log(1/tfs5$Score)

df <- tidyr::pivot_wider(tfs5, 
                         names_from = Query.Name, 
                         values_from = Score) %>% as.data.frame()
rownames(df) <- df$TF
df <- select(df, -c('TF'))

col_fun = circlize::colorRamp2(c(-8, -5, 0), c("black", "black", "yellow"))
Heatmap(df, name = 'Score',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))

# Calculate TFs for MP vs. P comparison: ---------
degs <- read.csv('2021-08-21 mp vs. p DEGs.csv', row.names = 'X')

# convert to human orthologs for ChEA3 from mouse genes
# Basic function to convert mouse to human gene names
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
# Last modified by JAX: 2023-01-16 08:00	version

convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}

# prepare for for loop to populate with ChEA3 queries
celltype_trt <- unique(degs$cellType)
tf_list <- vector(mode = "list", length = length(celltype_trt))

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"

for (i in 1:length(celltype_trt)){
  print(i)
  print(celltype_trt[[i]])
  deg_subset <- filter(degs, cellType == celltype_trt[[i]], p_val < 0.05, avg_log2FC > -0.06) # alter the > or < as needed
  genes <- deg_subset$gene_symbol
  
  # convert to human ortholog gene names (vast majority are same name capitalized, with a few exceptions)
  genes <- convert_mouse_to_human(genes)
  
  payload = list(query_name = celltype_trt[[i]], gene_set = genes)
  
  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  
  #results as list of R dataframes
  tf_list[[i]] <- fromJSON(json)
}

# save the full information in case I need it later:
saveRDS(tf_list, file = 'tf_list_full_up_MP_vs_P.rds')

# trim the data
tf_meanRank <- vector(mode = "list", length = length(tf_list))
for (i in 1:length(tf_meanRank)){
  tf_meanRank[[i]] <- tf_list[[i]]$`Integrated--meanRank`
}
tf_meanRank_master <- do.call(rbind, tf_meanRank)

saveRDS(tf_meanRank_master, 'tf_meanRank_master_up_MP_vs_P.rds')

# bring it together
upreg <- readRDS('tf_meanRank_master_up_MP_vs_P.rds')
downreg <- readRDS('tf_meanRank_master_down_MP_vs_P.rds')

upreg$Score <- as.numeric(upreg$Score) * (1)
downreg$Score <- as.numeric(downreg$Score)  * (-1)

tfs <- combine(upreg, downreg)
tfs$Rank <- as.numeric(tfs$Rank)

write.csv(tfs, 'tfs_meanRank_MP_vs_P.csv')

# get just top 100 tfs for each
write.csv(filter(tfs, Rank <= 100), 'tfs_meanRank_top100_MP_vs_P.csv')

# mp vs p visualize mean rank scores for up and downregulated genes MP vs P ----
tfs <- read.csv('tfs_meanRank_MP_vs_P.csv') %>% select(-c(Library, Overlapping_Genes))
tfs <- tfs[,2:5]
celltype_trt <- unique(tfs$Query.Name)

for (c in celltype_trt){
  tf_sub <- filter(tfs, Query.Name == c) %>% select(TF, Score)
  pos_tfs <- filter(tf_sub, Score > 0); colnames(pos_tfs) <- c('TF', 'pos_score')
  neg_tfs <- filter(tf_sub, Score < 0); colnames(neg_tfs) <- c('TF', 'neg_score')
  tf_sub <- left_join(pos_tfs, neg_tfs)
  
  top10_up <- arrange(pos_tfs, pos_score)$TF[1:10]
  top10_down <- arrange(neg_tfs, desc(neg_score))$TF[1:10]
  
  # math transformations based on panel C in Bouyer et al. Nat. Comm., 2020. 
  # https://www.nature.com/articles/s41467-021-22973-9/figures/6 
  tf_sub$pos_scaled <- log(1/tf_sub$pos_score)
  tf_sub$neg_scaled <- log(1/-tf_sub$neg_score)
  rownames(tf_sub) <- tf_sub$TF
  
  ggplot(tf_sub, aes(x=neg_scaled, y=pos_scaled, label = TF)) + 
    geom_point(size = 0.4, colour = 'cornflowerblue') + theme_classic() + geom_text_repel() +
    xlim(-8, 0) + ylim(-8, 0) + 
    xlab(expression(ln('ChEA3 MeanRank Neg Enriched'^ {-1}))) +
    ylab(expression(ln('ChEA3 MeanRank Pos Enriched'^ {-1})))
  ggsave(paste(c, 'TFs_mp_vs_p.png'), width = 4.5, height = 4.5)
}


# mp vs p plot heatmap of main TFs-----
library(ComplexHeatmap)
library(stringr)

# load TFs
tfs <- read.csv('tfs_meanRank_MP_vs_P.csv') %>% select(-c(Library, Overlapping_Genes))
celltype_trt <- unique(tfs$Query.Name)

# get top TFs, anything ranked in the top 5 for any cell type **4 h**
has_4h <- grepl(" 4 h", tfs$Query.Name)
tfs <- tfs[has_4h,]

tops <- unique( filter(tfs, Rank <= 5, Score > 0)$TF )
tfs5 <- filter(tfs, TF %in% tops, Score > 0) %>% select(-c(Rank))
tfs5$Score <- log(1/tfs5$Score)

df <- tidyr::pivot_wider(tfs5, 
                         names_from = Query.Name, 
                         values_from = Score) %>% as.data.frame()
rownames(df) <- df$TF
df <- select(df, -c('TF'))

col_fun = circlize::colorRamp2(c(-8, -5, 0), c("black", "black", "yellow"))
Heatmap(df, name = 'Score',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))


# Same for **24 h**
# load TFs
tfs <- read.csv('tfs_meanRank_MP_vs_P.csv') %>% select(-c(Library, Overlapping_Genes))
celltype_trt <- unique(tfs$Query.Name)

# get top TFs, anything ranked in the top 5 for any cell type **24 h**
has_4h <- grepl(" 24 h", tfs$Query.Name)
tfs <- tfs[has_4h,]

tops <- unique( filter(tfs, Rank <= 7, Score > 0)$TF )
tfs5 <- filter(tfs, TF %in% tops, Score > 0) %>% select(-c(Rank))
tfs5$Score <- log(1/tfs5$Score)

df <- tidyr::pivot_wider(tfs5, 
                         names_from = Query.Name, 
                         values_from = Score) %>% as.data.frame()
rownames(df) <- df$TF
df <- select(df, -c('TF'))

col_fun = circlize::colorRamp2(c(-8, -5, 0), c("black", "black", "yellow"))
Heatmap(df, name = 'Score',
        cluster_rows = T, cluster_columns = T, col = col_fun,
        border_gp = gpar(col = "black", lty = 2))
# MP_vs_P_TF_24h_top





