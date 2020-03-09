#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig5/B-C')

#------------- PARAMETERS -----------------

fname <- 'perplexity-20'
tsne.path <- paste(path.matrices, 'tsne_20_palette.tsv', sep = '/')
no.box <- TRUE

#------------------- Loadings  ------------------- 

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), paste(path.matrices, '181-clusters-266-genes.tsv', sep = '/'), min.cluster.size = 10)
spots.table <- append.tsne.to.spots.table(spots.table, tsne.path)
spots.table.original <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)

#------------------- Getting colors from palette  ------------------- 

df.colors.palette <- unique(spots.table[,c('cluster','clusters.named')])
rownames(df.colors.palette) <- df.colors.palette$clusters.named
df.colors.palette <- df.colors.palette[order(df.colors.palette$cluster),]
df.colors.palette$clusters.named <- NULL
df.colors.palette$clusters.colors <- gray.color

for(cl in rownames(df.colors.palette)){
  cnt <- plyr::count(spots.table.original[rownames(subset(spots.table, clusters.named == cl)),'clusters.named'])
  df.colors.palette[cl, 'clusters.colors'] <- df.colors.vivid[as.character(cnt[which.max(cnt$freq), 'x']),'clusters.colors']
}
colnames(df.colors.palette)[1] <- 'cluster.id'

#------------------- Plotting  ------------------- 

spots.table$color <- mapvalues(as.character(spots.table$clusters.named), 
                               from = rownames(df.colors.palette),
                               to = df.colors.palette$clusters.colors)
.plot.tsne(spots.table, file.path = paste('tsne',fname,'vivid-color.pdf',sep='-'), cex = 0.3, no.box = no.box)

# #3D t-SNE
# df.col <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
# spots.table$color <- NULL
# spots.table$color <- mapvalues(as.character(spots.table$cluster),
#                                from = df.col$cluster.id,
#                                to = df.col$clusters.colors)
# .plot.tsne(spots.table, file.path = paste('tsne',fname,'3dtsne-color.pdf',sep='-'), cex = 0.3, no.box = no.box)

#ARA PLOT
plot.tsne.all(tsne.path, spots.table, 'ARA', paste('tsne',fname,'ARA-color.pdf',sep='-'),  cex = 0.3, no.box = no.box)


