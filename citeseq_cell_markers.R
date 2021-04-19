## Identifying celltypes  ---------------------------
format <- theme(axis.text.x = element_text(angle = 90))

# Check this and compare to dsb
# DefaultAssay(combined) <- "ADT"
# VlnPlot(combined, features = c('IgG2-TotalA'))+NoLegend()
# ggsave("IgG-TotalA_ADT_CLR_V2.png", width = 10, height=5)
# VlnPlot(combined, features = c('IgG2-TotalA'), assay = "integratedADT_")+NoLegend()
# ggsave("IgG-TotalA_ADTintegrated_CLR_V2.png", width = 10, height=5)

# # Density Plots (Glioblastoma Paper Inspired) ---------------------------
# # get dsb-normalized values
# ADT_data_plotting<-as.data.frame(t( as.matrix(GetAssayData(combined, slot="data", assay="ADT"))))
# ADT_data_plotting$orig.ident<-combined$orig.ident[rownames(ADT_data_plotting)]
# ADT_data_plotting<-reshape2::melt(ADT_data_plotting)
# #ADT_data_plotting <- ADT_data_plotting[sample(nrow(ADT_data_plotting), 5000), ]
# colnames(ADT_data_plotting)[colnames(ADT_data_plotting)=="variable"]<-"protein"
# 
# # create list of individual dataframes for each ab:
# ab_list <- as.character(unique(ADT_data_plotting$protein))
# df_list <- vector(mode = "list", length = length(ab_list))
# prot_list <- ADT_data_plotting$protein
# for (i in 1:length(ab_list)){
#   protein <- ab_list[[i]]
#   df <- ADT_data_plotting[which(prot_list == protein), ]
#   df_list[[i]] <- df
#   ggplot(df_list[[i]],  aes(x=value,color=orig.ident, fill=orig.ident))+
#     geom_density(alpha=0.2)+
#     facet_wrap(~protein)+xlab("Arsinh geom counts")+
#     geom_rug()
#   filename <- gsub("/", "", paste(protein, "_density_plot_dsb.png"))
#   ggsave(filename, width = 5, height = 5)
# }
# 
## Plot all antibodies on wnnUMAP  ---------------------------
# Will select only the good ones and then recluster
DefaultAssay(combined) <- "integratedADT_"
protein.to.plot <- rownames(combined)
library(readxl)
ADT_with_associated_genes <- read_excel("Y:/Cole Keenum/CITE-seq files/ADT_with_associated_genes.xlsx")

bad_adt <- setdiff(ADT_with_associated_genes$Protein, protein.to.plot)
length(bad_adt)

idx <-  ADT_with_associated_genes$Protein %in% bad_adt
ADT_with_associated_genes <- ADT_with_associated_genes[!idx, ]

genes.to.plot <- ADT_with_associated_genes$Gene1

# # WITHOUT a 95% quartile cutoff for ADT
# for (i in 1:length(x = protein.to.plot)){
#   protein <- protein.to.plot[[i]]
#   gene <- genes.to.plot[[i]]
#   DefaultAssay(combined) <- "ADT"
#   p1 <- FeaturePlot(combined, features = protein, cols = c("lightgrey","darkgreen"))
#   DefaultAssay(combined) <- "SCT"
#   p2 <- FeaturePlot(combined, features = gene)
#   filename <- paste("FeaturePlot_", protein, ".png", sep = "")
#   filename <- gsub("/",  "", filename)
#   filename <- gsub("-TotalA",  "", filename)
#   ggsave(filename = filename, plot = p1 | p2,  width = 10, height = 5)
# }
# DefaultAssay(combined) <- "integratedSCT_"

# with a 99% quartile cutoff for ADT
for (i in 1:length(x = protein.to.plot)){
  protein <- protein.to.plot[[i]]
  gene <- genes.to.plot[[i]]
  DefaultAssay(combined) <- "ADT"
  p1 <- FeaturePlot(combined, features = protein, cols = c("lightgrey","darkgreen"),
                    min.cutoff = 0, max.cutoff = 'q99', reduction = 'wnn.umap')
  DefaultAssay(combined) <- "SCT"
  p2 <- FeaturePlot(combined, features = gene, reduction = 'wnn.umap')
  filename <- paste("wnnFeaturePlot_", protein, ".png", sep = "")
  filename <- gsub("/",  "", filename)
  filename <- gsub("-TotalA",  "", filename)
  ggsave(filename = filename, plot = p1 | p2,  width = 10, height = 5)
  #ggsave(filename = filename, plot = p1,  width = 10, height = 5)
}
DefaultAssay(combined) <- "integratedSCT_"

## Mouse cell atlas top10 markers  ---------------------------

# AT2 - 16
markers.to.plot <-  c('Sftpc','Sftpa1','Sftpb','Lyz2','Cxcl15','Slc34a2','Lamp3',
                      'Chil1','Sftpd','Napsa','Scd1','Sfta2','Dram1','S100g','Lpcat1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_AT2.png", width = 8.5, height = 7)


# B Cells - 5 + others
markers.to.plot <- c('Ly6d',
                     'Cd79a',
                     'Ms4a1',
                     'Cd79b',
                     'Ighd',
                     'Iglc3',
                     'Igkc',
                     'Iglc2',
                     'Fcmr',
                     'Iglc1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_B_Cell.png", width = 8.5, height = 7)
VlnPlot(combined, features = c("Cd19"), pt.size = 0) + NoLegend()
ggsave("Cd19_cells.png")

# Alveolar macrophage - 3
markers.to.plot <- c('Chil3',
                     'Ear2',
                     'Ear1',
                     'Plet1',
                     'Lpl',
                     'Ccl6',
                     'Atp6v0d2',
                     'Fabp1',
                     'Ctsd',
                     'Krt79')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Alveolar_Macrophage.png", width = 8.5, height = 7)
# check:
VlnPlot(combined, features = c("Marco"))

# CD8+ T - 4, 11
markers.to.plot <- c('Trbc2',
                     'Cd8b1',
                     'Ms4a4b',
                     'Cd3d',
                     'Trbc1',
                     'Dapl1',
                     'Cd3g',
                     'Thy1',
                     'Lck',
                     'Trac')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot)+ format
ggsave("MCA_CD8_T.png", width = 8.5, height = 7)
VlnPlot(combined, features = c("Cd3e"), pt.size = 0) + NoLegend()
ggsave("Cd3e_cells.png")
VlnPlot(combined, features = c("Cd8a"), pt.size = 0) + NoLegend()
ggsave("Cd8a_cells.png")
VlnPlot(combined, features = c("Cd8b1"), pt.size = 0) + NoLegend()
ggsave("Cd8b1_cells.png")
VlnPlot(combined, features = c("Cd4"), pt.size = 0) + NoLegend()
ggsave("Cd4_cells.png")

# pDC - 9 probably
markers.to.plot <- c('Ms4a6c',
                     'Plac8',
                     'Ly6c2',
                     'S100a4',
                     'Ifitm3',
                     'Ms4a4c',
                     'Ccl9',
                     'F13a1',
                     'Wfdc17',
                     'Cybb')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_pDC.png", width = 8.5, height = 7)

# NK cell - 7
markers.to.plot <- c('Gzma',
                     'Ccl5',
                     'Nkg7',
                     'Klra8',
                     'AW112010',
                     'Klra4',
                     'Gzmb',
                     'Ncr1',
                     'Klrb1c',
                     'Klra13-ps')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_NK.png", width = 8.5, height = 7)

# Eosinophil - 10 (neutrophil), or 17
markers.to.plot <- c('G0s2',
                     'Clec4d',
                     'Il1r2',
                     'Cxcl2',
                     'Il1b',
                     'S100a9',
                     'Acod1',
                     'S100a8',
                     'Trem1',
                     'Il1rn')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Eosinophil.png", width = 8.5, height = 7)

# Interstitial macrophage - N/A
markers.to.plot <- c('C1qc',
                     'C1qa',
                     'C1qb',
                     'Pf4',
                     'Apoe',
                     'Ccl8',
                     'Mgl2',
                     'Mmp12',
                     'Csf1r',
                     'Cd74')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot)+ format
ggsave("MCA_Interstitial_Macro.png", width = 8.5, height = 7) 
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# DC-Naaa high - 15 (except not Naaa high really)
markers.to.plot <- c('Cst3',
                     'Naaa',
                     'Plbd1',
                     'Irf8',
                     'Cd74',
                     'H2-DMa',
                     'Ppt1',
                     'H2-Eb1',
                     'Ckb',
                     'H2-Ab1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_DC_Naaa_High.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Stromal cell Dcn high - 14
markers.to.plot <- c('Dcn',
                     'Col3a1',
                     'Col1a2',
                     'Col1a1',
                     'Gsn',
                     'Serping1',
                     'Clec3b',
                     'Igfbp4',
                     'Mgp',
                     'Serpinf1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Stromal_Cell_Dcn_High.png", width = 8.5, height = 7)

# Stromal cell inmt high - 14 also (20 ish)
markers.to.plot <- c('Inmt',
                     'Mgp',
                     'Cxcl14',
                     'Sparcl1',
                     'Gpx3',
                     'Sod3',
                     'Pcolce2',
                     'Mfap4',
                     'Gsn',
                     'Bgn')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Stromal_Cell_Inmt_High.png", width = 8.5, height = 7)

# Clara cell - 16 likely
markers.to.plot <- c('Scgb1a1',
                     'Scgb3a2',
                     'Cyp2f2',
                     'Scgb3a1',
                     'Hp',
                     'Aldh1a1',
                     'Retnla',
                     'Cp',
                     'Cldn10',
                     'Cbr2')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Clara_Cell.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# cDC Gngt2 high - 6 maybe
markers.to.plot <- c('Gngt2',
                     'Lst1',
                     'Plac8',
                     'Adgre4',
                     'Clec4a1',
                     'Ifitm3',
                     'Ace',
                     'Ifitm6',
                     'Csf1r',
                     'Clec4a3')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_cDC_Gngt2_high.png", width = 8.5, height = 7)

# AT1 cell - N/A
markers.to.plot <- c('Ager',
                     'Igfbp2',
                     'Hopx',
                     'Gprc5a',
                     'Vegfa',
                     'Cyp2b10',
                     'Spock2',
                     'Krt7',
                     'Emp2',
                     'Aqp5')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_AT1.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# cDC Mgl2 high - 15 (I would call this CD209+ DC)
markers.to.plot <- c('Ccl17',
                     'Cd209a',
                     'Mgl2',
                     'Cd74',
                     'H2-Eb1',
                     'Ccl22',
                     'Cd83',
                     'H2-Ab1',
                     'H2-Aa',
                     'Il1b')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_cDC_Mg12_high.png", width = 8.5, height = 7)
DotPlot(combined, assay = "SCT", features = markers.to.plot)

# Nuocyte - 12 
markers.to.plot <- c('Cxcr6',
                     'Icos',
                     'Thy1',
                     'S100a4',
                     'Cd3g',
                     'Trbc2',
                     'Trac',
                     'Tcrg-C1',
                     'Rora',
                     'Crem')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Nuocyte.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Endothelial cell Tmem100 high - 2 (similar to 13ish, no Plvap in 13)
markers.to.plot <- c('Plvap',
                     'Tmem100',
                     'Tspan7',
                     'Ramp2',
                     'Egfl7',
                     'Cldn5',
                     'Gpihbp1',
                     'Hpgd',
                     'Aqp1',
                     'Ptprb')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Endothelial_Cell_Tmem100_High.png", width = 8.5, height = 7)

# Cilliated cell - N/A
markers.to.plot <- c('Ccdc153',
                     'Tmem212',
                     'Dynlrb2',
                     'Fam183b',
                     'Cfap126',
                     'Tppp3',
                     'Sec14l3',
                     '1110017D15Rik',
                     'Cyp2s1',
                     'Scgb1a1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Cilliated_Cell.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Endothelial cell Kdr high - 13
markers.to.plot <- c('Cldn5',
                     'Kdr',
                     'Car4',
                     'Cyp4b1',
                     'Cdh5',
                     'Ramp2',
                     'Calcrl',
                     'Thbd',
                     'Cyp1a1',
                     'Ly6a')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Endothelial_Cell_Kdr_High.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Endothelial cells Vwf high - 2, 8-ish
markers.to.plot <- c('Vwf',
                     'Lyve1',
                     'Plvap',
                     'Aqp1',
                     'Prss23',
                     'Cytl1',
                     'Tm4sf1',
                     'Ptprb',
                     'Vcam1',
                     'Cav1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Endothelial_Cell_Vwf_High.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Neutrophil - 10 probably, competing with eosinophil
markers.to.plot <- c('Ngp',
                     'Camp',
                     'Ltf',
                     'S100a9',
                     'S100a8',
                     'Retnlg',
                     'Wfdc21',
                     'Mmp8',
                     'Lcn2',
                     'G0s2')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Neutrophil.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Dividing DC - N/A
markers.to.plot <- c('Cst3',
                     'Naaa',
                     'Plbd1',
                     'Irf8',
                     'Tubb5',
                     'Hmgb2',
                     'H2afz',
                     'Pclaf',
                     'Birc5',
                     'Tuba1b')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Dividing_DC.png", width = 8.5, height = 7)

# Stromal cell Acta2 high - 14
markers.to.plot <- c('Acta2',
                     'Myl9',
                     'Hhip',
                     'Mustn1',
                     'Actc1',
                     'Tpm2',
                     'Enpp2',
                     'Igfbp5',
                     'Tagln',
                     'Des')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Sromal_Cell_Acta2_High.png", width = 8.5, height = 7)

# Dividing T cells - 12 maybe
markers.to.plot <- c('Hmgb2',
                     'Pclaf',
                     'Trac',
                     'Ube2c',
                     'Trbc2',
                     'Birc5',
                     'Ms4a4b',
                     'Lgals1',
                     'Top2a',
                     'Nusap1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Dividing_T_Cells.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Dividing cells - N/A (some similar to 12)
markers.to.plot <- c('Cdc20',
                     'Ube2c',
                     'Stmn1',
                     'Pclaf',
                     'Tubb5',
                     'Tuba1b',
                     'Lockd',
                     'Birc5',
                     'H2afx',
                     'Spc24')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Dividing_Cells.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Alveolar macrophage  Pclaf high - N/A (some similar to 3)
markers.to.plot <- c('Pclaf',
                     'Chil3',
                     'Tuba1b',
                     'Krt79',
                     'Birc5',
                     'Stmn1',
                     'Tubb5',
                     'Rrm2',
                     'Atp6v0d2',
                     'Tlr2')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Alveolar_Macrophage_Pclaf_High.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Alveolar bipotent progenitor - 16 maybe
markers.to.plot <- c('Clu',
                     'Ctsh',
                     'Krt8',
                     'Krt18',
                     'Phlda1',
                     'Sftpd',
                     'Chia1',
                     'Cldn18',
                     'Sftpa1',
                     'Ndnf')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot)  + format
ggsave("MCA_Alveolar_Bipotent_Progenitor.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Ig producing B cell - No J chain but otherwise 5
markers.to.plot <- c('Jchain',
                     'Igha',
                     'Igkc',
                     'Ighm',
                     'Igkv2-109',
                     'Iglv1',
                     'Iglc2',
                     'Iglc1',
                     'Mzb1',
                     'Txndc5')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Ig_Producing_B_Cell.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# cDC H2-M2 high - N/A (15 some similar)
markers.to.plot <- c('Fscn1',
                     'Ccl22',
                     'Nudt17',
                     'Il4i1',
                     'Zmynd15',
                     'Tmem123',
                     'Apol7c',
                     'Ccr7',
                     'Tspan3',
                     'Ccl17')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_cDC_H2-M2_High.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Basophill - 17
markers.to.plot <- c('Ccl4',
                     'Ccl3',
                     'Il6',
                     'Osm',
                     'Ptgs2',
                     'Ifitm1',
                     'Ier3',
                     'Cks2',
                     'Il4',
                     'Nfkbia')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Basophil.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# Monocyte progenitor - N/A
markers.to.plot <- c('Elane',
                     'Mpo',
                     'Ctsg',
                     'Prtn3',
                     'Ms4a3',
                     'Pclaf',
                     'Ly6c2',
                     'Plac8',
                     'Cmtm7',
                     'I830127L07Rik')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_Monocyte_Progenitor.png", width = 8.5, height = 7)
DotPlot(combined, assay = "RNA", features = markers.to.plot)

# cDC Tubb5 high - N/A (some shared with 15)
markers.to.plot <- c('Ccl17',
                     'Cd209a',
                     'Tubb5',
                     'Pclaf',
                     'Tuba1b',
                     'Ube2c',
                     'Cd74',
                     'Top2a',
                     'H2-Eb1',
                     'H2-Ab1')
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot) + format
ggsave("MCA_cDC_Tubb5_High.png", width = 8.5, height = 7)

# Markers to define metaclusters: ---------
markers.to.plot <- c('Ptprc', 'Pecam1', 'Prox1', 'Hopx', 'Sftpc',
                     'Trbc2', 'Ccl5', 'Retnlg', 'Mcpt8', 'Mcpt4', 'F13a1', 
                     'Epcam', 'Cdh5', 'Col1a2', 'Scgb3a2', 'Foxj1', 'Enpp2', 
                     'Mfap4', 'Ly6c2', "Akap5", "Lamp3", "Foxp3")
DefaultAssay(combined) <- "SCT"
for (marker in markers.to.plot){
  FeaturePlot(combined, features = marker, reduction = 'wnn.umap')
  filename <- paste("FeaturePlot_", marker, ".png", sep = "")
  ggsave(filename = filename, width = 5, height = 5)
}
DefaultAssay(combined) <- "integratedSCT_"

# PARAMITA MARKERS ---------
#avaleor macrophage - 
markers.to.plot <- c("Chi3l3", "Ccl6", "Mrc1", "Fabp1", "Ctsd", "Ear1", "Ear2", "Cd9", "Il18", "Krt77", "Cd44")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot,cols = c("blue", "red"), dot.scale = 8) +  RotatedAxis()
ggsave('paramita_alveolar_macro.png', width = 8.5, height = 7)

# Fn+ macrophage - ?                           
markers.to.plot <- c("Saa3", "Prg4", "Fcna", "Alox15", "Cd5l", "Selp", "Fabp7", "Vsig4", "Gbp2b", "Padi4", "Msr1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +  RotatedAxis() 
ggsave('paramita_Fn+_macro.png', width = 8.5, height = 7)

# Bcells -                                            
markers.to.plot <- c("Ms4a1", "Ly6d", "Cd79b", "Iglc2", "Bank1", "Fcmr", "Cd79a", "Cd37", "Cd19", "Siglecg", "Cd74")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +  RotatedAxis()
ggsave('paramita_Bcells.png', width = 8.5, height = 7)
#VlnPlot(combined, features = c("Ms4a1"), assay = "RNA")
#VlnPlot(combined, features = c("CD19-TotalA"), assay = "integratedADT_")

# Capilary endothelial - 21                                             
markers.to.plot <- c("Ednrb", "Car4", "Fibin", "Kdr", "Kitl", "Tbx3", "Icam2", "Prx", "Tmeff2", "Tmcc2", "Ly6c1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +  RotatedAxis()
ggsave('paramita_capilary_endo.png', width = 8.5, height = 7)

# Dendritic cells - N/A
markers.to.plot <- c("Fscn1", "Ccl22", "Il4i1", "Mreg", "Cacnb3", "Nr4a3", "Pmaip1", "Nudt17", "Relb", "Samsn1", "Ccr7")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('paramita_dendritic_cells.png', width = 8.5, height = 7)

# CD103+ dendritic cells - 18 transcriptionally, 17 and 46 ab
markers.to.plot <- c("Plbd1", "Ifi205", "Irf8", "Xcr1", "Qpct", "Sept3", "Clec9a", "Itgae", "Gcsam", "Klrb1b", "Ckb")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +  RotatedAxis()
ggsave('paramita_CD103+_DC.png', width = 8.5, height = 7)

FeaturePlot(combined, features = c("CD103-TotalA"))
VlnPlot(combined, features = c("CD103-TotalA"))

# CD209+dendritic cells - 18 looks more obvious, 39 also
markers.to.plot <- c("Ccl17", "Mgl2", "Cd209a", "Clec4b1", "H2-Eb1", "Cfp", "S100a4", "Il1rl1", "Slamf7", "Fgfr1", "Wfdc17")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +  RotatedAxis()
ggsave('paramita_CD209+_DC.png', width = 8.5, height = 7)

VlnPlot(combined, features = c("Cd209a"), assay = "integratedSCT_")

# CD4 Tcells - 5, 12 (icos+), 28(icos+), 35, 39-ish, 42, 46 (icos+), 50, 51
markers.to.plot <- c('Cd4', "Trbc2", "Il7r", "Cd3g", "Icos", "Cxcr6", "Cd3d", "Trac", "Thy1", "Cd3e", "Skap1", "Trbc1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, dot.scale = 8) + RotatedAxis()
ggsave('paramita_CD4+_T.png', width = 8.5, height = 7)

# CD8 Tcells - 17, 28-ish? (maybe 28 is mixed)
markers.to.plot <- c("Ms4a4b", "Cd8b1", "Ccl5", "Il7r", "Ly6c2", "Ms4a6b", "Cd3d", "Trbc2", "Nkg7", "Dapl1", "Cd8a")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_CD8+_T.png', width = 8.5, height = 7)

VlnPlot(combined, features = c("CD8a-TotalA"))
ggsave("CD8a_res1_ASINH_GEOM_vlnPlot.png", width = 10, height = 5)
VlnPlot(combined, features = c("CD8b-TotalA"))
ggsave("CD8b_res1_ASINH_GEOM_vlnPlot.png", width = 10, height = 5)

# Ciliated markers - 52 definitely
markers.to.plot <- c("Ccdc153", "Dynlrb2", "Sec14l3", "Tmem212", "Fam183b", "Tppp3", "Rsph1", "Ccdc39", "Riiad1", "Pcp4l1", "Cyp2s1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_cilliated_.png', width = 8.5, height = 7)

# Classical monocytes - 13 most obviously
markers.to.plot <- c("Ccr2", "Plac8", "F13a1", "Ms4a6c", "Gm9733", "Ifitm6", "Vcan", "S100a4", "Fn1", "Ccl9", "Gpr141")
DotPlot(combined, assay = "integratedSCT_",features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_classical_mono.png', width = 8.5, height = 7)

# Club cells - 52
markers.to.plot <- c("Scgb1a1", "Cyp2f2", "Scgb3a2", "Hp", "Aldh1a1", "Cbr2", "Cp", "Fmo3", "Cldn10", "Gsta3", "Prdx6")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_club.png', width = 8.5, height = 7)

# Eosinophil- a lot... 0, 1, 11 (probably neutrophil), 14 (low expressing gene markers), 24, 33, 37, 50, 53
markers.to.plot <- c("S100a9", "S100a8", "Retnlg", "Il1b", "Cxcr2", "Csf3r", "Stfa2l1", "Clec4d", "Il1f9", "H2-Q10", "Il1r2", "Cxcl2", "Acod1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('eosinophil_paramita.png', width = 8.5, height = 7)

FeaturePlot(combined, features = c("Ly6c2"), reduction = "wnn.umap")

#Fibroblasts - 40 definitely, maybe 20? idk
markers.to.plot <- c("Dcn", "Col1a2", "Col3a1", "Serping1", "Clec3b", "Dpt", "Mfap5", "Mmp3", "Col1a1", "Serpinf1", "Cygb")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
#VlnPlot(combined, assay = "integratedSCT_", features = c("Dcn"))
ggsave('paramita_fibroblasts.png', width = 8.5, height = 7)

# Interstitial macrophage - 19 def., maybe 18
markers.to.plot <- c("C1qb", "C1qc", "C1qa", "Ms4a7", "Pf4", "Aif1", "Ccl12", "Cd163", "Mgl2", "Mafb", "C3ar1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_interstitial_macro.png', width = 8.5, height = 7)

# lipo fribroblast - 16, 29, 38 (most likely AT1 cells), 40, 43
markers.to.plot <- c("Inmt", "Gpx3", "Ogn", "Mfap4", "Pcolce2", "Tcf21", "Sod3", "Cdo1", "Itga8", "Pdgfra", "Angpt1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_lipo_fibro.png', width = 8.5, height = 7)

# Lymphatic fibroblast - 49 definitely
markers.to.plot <- c("Mmrn1", "Prox1", "Gm525", "Reln", "Flt4", "Klhl4", "Gng11", "F8", "Fxyd6", "Sema3d", "Thrsp")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_lymphatic_fibro.png', width = 8.5, height = 7)

# Megakaryocytes - N/A
markers.to.plot <- c("Ppbp", "Nrgn", "Clec1b", "Alox12", "Gp9", "Itga2b", "Tmem40", "Tubb1", "Treml1", "Gp1ba", "Gp5")
DotPlot(combined, assay = "RNA", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_megakarocyte.png', width = 8.5, height = 7)

# NK cells - 9, 47, 28ish? (28 is hard to figure out)
markers.to.plot <- c("Gzma", "Ccl5", "Klra8", "Nkg7", "Klrb1c", "Klre1", "Klra4", "Prf1", "Gzmb", "Klra9", "Klrg1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
#FeaturePlot(combined, features = c("Ncr1"), reduction = "wnn.umap")
ggsave('paramita_NK.png', width = 8.5, height = 7)

# Neutrophils - 11 based on RNAseq
markers.to.plot <- c("Lcn2", "G0s2", "Camp", "Ngp", "Ltf", "Cd177", "Wfdc21", "Mmp8", "Itgb2l", "Ms4a3", "Cebpe", "Ly6g", "Mpo", "Retnlg")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
#FeaturePlot(combined, features = c("Ly6g.Ly6c-TotalA"))
#VlnPlot(combined, assay = "integratedSCT_", features = c("Ly6g"))
ggsave('paramita_neutro.png', width = 8.5, height = 7)

# Nonclassical Monocyte - 15, 30
markers.to.plot <- c("Nr4a1", "Gngt2", "Lst1", "Apoc2", "Ifitm6", "Emr4", "Cd300e", "Clec4a1", "Eno3", "Plac8", "Clec4a3")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_nonclassical_mono.png', width = 8.5, height = 7)

# Plasma cells - N/A
markers.to.plot <- c("Igj", "Mzb1", "Iglc2", "Prg2", "Eaf2", "Derl3", "Tnfrsf17", "Igkj2", "Txndc5", "Iglv3", "Oosp1")
DotPlot(combined, assay = "RNA", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_plasma_cell.png', width = 8.5, height = 7)

# Goblet cells - N/A (52 kinda)
markers.to.plot <- c("Scgb3a1", "Reg3g", "Bpifb1", "Tff2", "Lypd2", "Sult1d1", "Ltf", "Muc5b", "Chad", "Pigr", "Bpifa1")
DotPlot(combined, assay = "RNA", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_goblet.png', width = 8.5, height = 7)

# AT2 cells - 8, 25, 27 ish but not really
markers.to.plot <- c("Sftpc", "Sftpa1", "Sftpb", "Lyz2", "Cxcl15", "Slc34a2", "Lamp3", "Chil1", "Sftpd", "Napsa", "Scd1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_AT2.png', width = 8.5, height = 7)

# AT1 cells - 31, 38 (obvious)
markers.to.plot <- c("Ager", "Igfbp2", "Hopx", "Gprc5a", "Vegfa", "Cyp2b10", "Spock2", "Krt7", "Emp2", "Aqp5", "Pdpn")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_AT1.png', width = 8.5, height = 7)

# Dendritic cells - 18
markers.to.plot <- c("Ccl17", "Cd209a", "Mgl2", "Cd74", "H2-Eb1", "Ccl22", "Cd83", "H2-Ab1", "H2-Aa", "Il1b", "Bcl2a1d")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_dendritic_cells.png', width = 8.5, height = 7)

# Endothelial cells_Kdr high - 21
markers.to.plot <- c("Kdr", "Car4", "Cyp4b1", "Cdh5", "Ramp2", "Calcrl", "Thbd", "Cyp1a1", "Ly6a", "Icam2", "Tspan7")
DotPlot(combined, assay = "integratedSCT_",features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_endothelial_kdr.png', width = 8.5, height = 7)

# Endothelial cells_Vwf high - 26
markers.to.plot <- c("Lyve1", "Vwf", "Plvap", "Aqp1", "Prss23", "Cytl1", "Tm4sf1", "Ptprb", "Vcam1", "Cav1", "Cpe")
DotPlot(combined, assay = "integratedSCT_",features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_endothelial_vwf.png', width = 8.5, height = 7)

# Alveolar bipotent progenitor - 8? idk
markers.to.plot <- c("Clu", "Ctsh", "Krt8", "Krt18", "Phlda1", "Chia1", "Cldn18", "Ndnf", "Ccdc107", "Anxa3", "Areg")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_alv_bipotent.png', width = 8.5, height = 7)

# Myeloid cells - 34
markers.to.plot <- c("Ccl4", "Ccl3", "Il6", "Osm", "Ptgs2", "Ifitm1", "Ier3", "Cks2", "Il4", "Nfkbia", "Lilr4b")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave('paramita_myeloid.png', width = 8.5, height = 7)

# ???
markers.to.plot <- c("Tspan7", "Acta2", "Cd8b1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
VlnPlot(combined, assay = "integratedSCT_", features = c("Acta2"))

## from fibroblast paper Liu et al. bioRxiv Fig. 1

# Lipofibroblast - N/A
markers.to.plot <- c("Tcf21", "Plin2", "Fgf10", "G0s2", "Gyg", "Macf1", "Wnt2", "Col13a1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('liu_lipofibroblast.png', width = 8.5, height = 7)

# Myofibroblast - 20 probably, from fibroblast paper
markers.to.plot <- c("Acta2", "Myh11", "Tagln", "Pdgfra", "Tgfbi", "Hhip", "Enpp2", "Wnt5a")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('liu_myofibroblast.png', width = 8.5, height = 7)

# Ebf1+ fibroblasts - 22 clearly
markers.to.plot <- c("Pdgfrb", "Higd1b", "Cox4i2", "Notch3", "Ebf1", "Gucy1a3", "Pdzd2", "Postn")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('liu_Ebf14+_fibro.png', width = 8.5, height = 7)

# Intermediate fibroblasts - N/A
markers.to.plot <- c("Agtr2", "Prss35", "Igfbp7", "Fbln5", "Ptn", "Heyl", "Fstl1", "Tm4sf1")
DotPlot(combined, assay = "RNA", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('liu_intermediate_fibro.png', width = 8.5, height = 7)

# Mesothelial cells - N/A
markers.to.plot <- c("Upk3b", "Wt1", "Msln", "Upk1b", "Aldh1a2", "Cpe", "Gm12840", "Aqp1")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('liu_mesothelial.png', width = 8.5, height = 7)

# Proliferative fibroblasts - N/A
markers.to.plot <- c("Hmmr", "Mki67", "Pcna", "Top2a", "Spc25", "Cdca3", "Ccnb2", "Hist1h2ap")
DotPlot(combined, assay = "integratedSCT_", features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave('liu_proliferative_fibro.png', width = 8.5, height = 7)

## from endothelial paper Gillich et al. Extended Data Fig. 1

# aCap -  (Same as Endothelial kdr high)
markers.to.plot <- c("Car4", "Ednrb", "Fibin", "Tbx2", "Cdkn2b", "Rprml", "Chst1", "Apln")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)
ggsave('Gillich_aCap.png', width = 8.5, height = 7)

# gCap - 2 and 3 most strongly, 10 perhaps (debatable)
markers.to.plot <- c("Cd93", "Ptprb", "Plvap", "Gpihbp1", "H2-Ab1", "Tek", "Kit", "Aplnr")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)
ggsave('Gillich_gCap.png', width = 8.5, height = 7)

# Artery - N/A
markers.to.plot <- c("Mgp", "Cdh13", "Htra1", "Bmx", "Gja5", "Fbln2", "Sulf1")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)
ggsave('Gillich_artery.png', width = 8.5, height = 7)

# Vein - Same as endothelial cells Vwf hi, already labeled
markers.to.plot <- c("Vwf", "Slc6a2", "Bst1", "Car8", "Amigo2", "Mustin1", "Vegfc", "Csrp2", "Nr2f2")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)
ggsave('Gillich_vein_endo.png', width = 8.5, height = 7)

# Lymphatics - Same as lymphatic fibroblasts, already labeled
markers.to.plot <- c("Mmrn1", "Fxyd6", "Reln", "Pdpn", "Thy1", "Nrp2", "Tbx1", "Gja1")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)
ggsave('Gillich_lymphatics.png', width = 8.5, height = 7)

## Mouse Cell Atlas: Table S4

# Clara cell - N/A
markers.to.plot <- c("Scgb1a1", "Aldh1a1", "Cyp2f2", "Scgb3a1", "Hp")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)

# Nuocyte - Probably 23. debatable
# see Barlow et al. JLB 2011, Figure 2:
markers.to.plot <- c("Cxcr6", "Icos", "Thy1", "S100a4", "Il7r")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)

# Basophil - N/A
markers.to.plot <- c("Ccl4", "Ccl3", "Il6", "Cd69", "Cd200r30")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)

# Plasmacytoid DC - N/A
markers.to.plot <- c("Ms4a6c", "Plac8", "Bst2", "Irf7", "Irf5")
DotPlot(combined, assay = "integratedSCT_", features  = markers.to.plot)

# Neutrophil Markers ---------
# Markers based on Zillionis paper
markers.to.plot <- c('Mmp8', 'Ifit1', 'Cxcl3', 'Pald1', 'Ccl3', 'Ctsc')
DefaultAssay(neutro) <- "SCT"
for (marker in markers.to.plot){
  FeaturePlot(neutro, features = marker, reduction = 'sct.umap')
  filename <- paste("neutro_FeaturePlot_", marker, ".png", sep = "")
  ggsave(filename = filename, width = 5, height = 5)
}
DefaultAssay(neutro) <- "integratedSCT_"

# Zillionis dot plots
# mN1
markers.to.plot <- c('Retnlb', 'Mmp8', 'Retnlg', 'Rsad2', 'Ifit3', 'Ifit1', 
                     'Cxcl10', 'Stat1', 'Irf7', 'Cxcl3', 'Tgm2', 'Cd14', 'Cass4', 
                     'Xbp1', 'Pald1', 'Exoc4', 'Gpnmb', 'Ccl3', 'Cstb', 'Ctsb', 'Cd63',
                     'Irak2', 'Ngp', 'Adamdec1', 'Chil3', 'Ctsc', 'Fcnb')
DotPlot(neutro, assay = "integratedSCT_", features  = markers.to.plot)
ggsave('Zillionis_dot.png', width = 8.5, height = 7)
