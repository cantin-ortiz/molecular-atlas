#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig4/B')

library(Rvcg)

#--------------- PARAMETERS --------------- 

clusters.v1 <- setdiff(clusters.v1,2)

M1 <- rbind(c(0.865400373935699,0,0.501081049442291,0),c(0.211013734340668,0.907006323337555,-0.364434778690338,0),c(-0.454483687877655,0.421116977930069,0.784923613071442,0),c(0,0,0,1))
M2 <- rbind(c(0,0,-1,0),c(-1,0,0,0),c(0,1,0,0),c(0,0,0,1))
zoom = 0.43

#--------------- LOADING --------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)

load(path.all.meshes)

#--------------- Graphical parameters ---------------

mesh3d.new.window()

par3d(userMatrix = M2, zoom = zoom)

for(cl in clusters.v1){
  mesh <- l$data[[which(l$clusters == cl)]]
  wire3d(vcgSmooth(mesh,'HC', iteration = 5), col = 'red', materials = list(shininess = 128, specular = 'black'))
}

for(cl in clusters.alm){
  mesh <- l$data[[which(l$clusters == cl)]]
  wire3d(vcgSmooth(mesh,'HC', iteration = 5), col = '#3A3292', materials = list(shininess = 128, specular = 'gray10'))
}

rgl.snapshot('screenshot-ground-truth.png')
