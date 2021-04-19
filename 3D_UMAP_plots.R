# 3D UMAP Projections
# Based on the following script:
# https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0/blob/master/3D%20UMAP%20Plotting%20v1.3.R
library(plotly)

# load Combined
combined <- readRDS('combined_04192021.rds')

combined <- RunUMAP(combined, nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap3", reduction.key = "wnnUMAP3_", 
                    n.components = 3L)
umap_1 <- combined[["wnn.umap3"]]@cell.embeddings[,1]
umap_2 <- combined[["wnn.umap3"]]@cell.embeddings[,2]
umap_3 <- combined[["wnn.umap3"]]@cell.embeddings[,3]

head(Embeddings(object = combined, reduction = "wnn.umap3"))

plot.data <- FetchData(object = combined, vars = c("wnnUMAP3_1", "wnnUMAP3_2", "wnnUMAP3_3", "celltype"))
#plot.data$label <- paste(rownames(plot.data))
plot.data$label <- paste(Idents(combined))

rownames(color_df) <- color_df$Celltype
color_df <- color_df[order(color_df$Freq, decreasing = T), ]

fig <- plot_ly(data = plot.data, 
               x = ~wnnUMAP3_1, y = ~wnnUMAP3_2, z = ~wnnUMAP3_3, 
               color = ~celltype,
               colors = color_df$Color, # from main script,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
fig

# looking at Mmp8

plot.data <- FetchData(object = combined, vars = c("wnnUMAP3_1", "wnnUMAP3_2", "wnnUMAP3_3", "Mmp8"), slot = 'data')
plot.data$label <- paste(Idents(combined))

plot_ly(data = plot.data, 
        x = ~wnnUMAP3_1, y = ~wnnUMAP3_2, z = ~wnnUMAP3_3, 
        color = ~Mmp8, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text")

# looking at Ccl3

plot.data <- FetchData(object = combined, vars = c("wnnUMAP3_1", "wnnUMAP3_2", "wnnUMAP3_3", "Ccl3"), slot = 'data')
plot.data$label <- paste(Idents(combined))

plot_ly(data = plot.data, 
        x = ~wnnUMAP3_1, y = ~wnnUMAP3_2, z = ~wnnUMAP3_3, 
        color = ~Ccl3, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text")

# load Neutro
# WORKSPACE SAVED 04/19/2021 ~5PM APPROX

neutro <- RunUMAP(neutro, reduction = 'pca', dims = 1:30,
                  assay = 'integratedSCT_', reduction.name = 'sct.3d',
                  reduction.key = 'sct3d_', n.components = 3L)

umap_1 <- neutro[["sct.3d"]]@cell.embeddings[,1]
umap_2 <- neutro[["sct.3d"]]@cell.embeddings[,2]
umap_3 <- neutro[["sct.3d"]]@cell.embeddings[,3]

Embeddings(object = neutro, reduction = "sct.3d")

plot.data <- FetchData(object = neutro, vars = c("sct3d_1", "sct3d_2", "sct3d_3", "seurat_clusters"))
plot.data$label <- paste(rownames(plot.data))

fig <- plot_ly(data = plot.data, 
               x = ~sct3d_1, y = ~sct3d_2, z = ~sct3d_3, 
               color = ~seurat_clusters,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

plot.data <- FetchData(object = neutro, vars = c("sct3d_1", "sct3d_2", "sct3d_3", "celltype"))
plot.data$label <- paste(rownames(plot.data))

plot_ly(data = plot.data, 
        x = ~sct3d_1, y = ~sct3d_2, z = ~sct3d_3, 
        color = ~celltype,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)

fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig
fig_cube

DefaultAssay(neutro) <- "SCT"
plot.data <- FetchData(object = neutro, vars = c("sct3d_1", "sct3d_2", "sct3d_3", "Mmp8"), slot = 'data')

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plot.data$changed <- ifelse(test = plot.data$Mmp8 <1, yes = plot.data$Mmp8, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plot.data$label <- paste(rownames(plot.data)," - ", plot.data$Mmp8, sep="")

plot_ly(data = plot.data, 
        x = ~sct3d_1, y = ~sct3d_2, z = ~sct3d_3, 
        color = ~Mmp8, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text")
        