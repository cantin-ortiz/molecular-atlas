#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig2/J')

#------------- PARAMETERS -----------------

spots.table <- load.spots.table()
smoothing.itt <- 5

param.list <- list(list(coord = 1.75, mode = 'saggital'),
                   list(coord = -2, mode = 'horizontal'))

#------------- LOADING -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

load(path.all.meshes)
mesh.outline <- mesh3d.allen.annot.from.id(get.id.from.acronym('root'))

#------------- Rendering sections -----------------
#-------------*** Straight cuts -----------------

for (i in 1:length(param.list)){
  
  
  coord <- param.list[[i]]$coord
  mode  <- param.list[[i]]$mode
  
  l.output <- get.plane.equation(coord, mode = mode)
  N <- l.output$N
  d <- l.output$d
  rm(l.output)
  
  list.intersect <- get.meshes.list.plane.intersection(l$data, smoothing.itt, N, d)
  list.intersect.allen <- get.meshes.list.plane.intersection(list(mesh.outline), 0, N, d)
  
  polygon.list <- get.polygon.cut(list.intersect, mode = mode)
  polygon.list.allen <- get.polygon.cut(list.intersect.allen, mode = mode)

  fname <- sprintf('polygons_outline_%s_%.3f', mode, coord)
  fname <- gsub('[.]', 'mm', fname)
  fname <- sprintf('%s.pdf', fname)
  
  plot.polygon.cut(c(list(polygon.list.allen[[1]]),polygon.list), 
                   c('lightgray', l$color), 
                   l.lim = get.lim.mode(param.list[[i]]$mode),
                   fname = fname, 
                   main = paste(mode,coord,sep='_'))
}


