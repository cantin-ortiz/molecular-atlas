#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig1/D')

#--------------- LOADING --------------- 

spots.table <- add.parent.acronym(load.spots.table())

#--------------- PLOTTING --------------- 

spots.table$color <- color.from.acronym(spots.table$acronym)

.plot.3d.glassbrain.col(spots.table, dim = c(-1920,0,0,2000), HD = TRUE)
plot.3d.glassbrain.setview('medial', zoom = 0.5)
rgl.snapshot('allen-3d-medial.png', top = TRUE)
plot.3d.glassbrain.setview('dorsal', zoom = 0.55)
rgl.snapshot('allen-3d-dorsal.png', top = TRUE)
plot.3d.glassbrain.setview('3d', zoom = 0.55)
rgl.snapshot('allen-3d-3d.png', top = TRUE)
rgl.close()

spots.table$color <- 'gray'
.plot.3d.glassbrain.col(spots.table, dim = c(-1920,0,0,2000), HD = TRUE)
plot.3d.glassbrain.setview('medial', zoom = 0.5)
rgl.snapshot('gray-3d-medial.png', top = TRUE)
plot.3d.glassbrain.setview('dorsal', zoom = 0.55)
rgl.snapshot('gray-3d-dorsal.png', top = TRUE)
plot.3d.glassbrain.setview('3d', zoom = 0.55)
rgl.snapshot('gray-3d-3d.png', top = TRUE)
rgl.close()
