#------------------- Includes  ------------------- 

rm(list = ls())
source('bin/includes.R')
setwd('figures/fig2/B-E-F')

#------------------- Parameters  ------------------- 

fname <- 'perplexity-20-no-box'
no.box <- TRUE

#------------------- Loadings  ------------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
spots.table <- append.tsne.to.spots.table(spots.table, tsne.2d.path)

#------------------- Plotting  ------------------- 

#Vivid color
spots.table <- append.common.vivid.colors(spots.table)
.plot.tsne(spots.table, file.path = paste('tsne',fname,'vivid-color.pdf',sep='-'), cex = 0.3, no.box = no.box)

#3D t-SNE
df.col <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
spots.table$color <- NULL
spots.table$color <- mapvalues(as.character(spots.table$cluster),
                               from = df.col$cluster.id,
                               to = df.col$clusters.colors)
.plot.tsne(spots.table, file.path = paste('tsne',fname,'3dtsne-color.pdf',sep='-'), cex = 0.3, no.box = no.box)


#ARA PLOT
plot.tsne.all(tsne.2d.path, spots.table, 'ARA', paste('tsne',fname,'ARA-color.pdf',sep='-'),  cex = 0.3, no.box = no.box)
