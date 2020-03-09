#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig2/D')

#--------------- LOADING --------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)

#--------------- PLOTTING --------------- 

df.col <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
spots.table$color <- NULL
spots.table$color <- mapvalues(as.character(spots.table$cluster),
                                from = df.col$cluster.id,
                                to = df.col$clusters.colors)


.plot.3d.glassbrain.col(spots.table, dim = c(-1920,0,0,2000), HD = TRUE)
plot.3d.glassbrain.setview('medial', zoom = 0.5)
rgl.snapshot('cluster-3dtsne-medial.png', top = TRUE)
plot.3d.glassbrain.setview('dorsal', zoom = 0.55)
rgl.snapshot('cluster-3dtsne-dorsal.png', top = TRUE)
plot.3d.glassbrain.setview('3d', zoom = 0.55)
rgl.snapshot('cluster-3dtsne-3d.png', top = TRUE)
rgl.close()
