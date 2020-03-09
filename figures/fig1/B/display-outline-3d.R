#------------- INCLUDES -----------------

rm(list = ls())
source('bin/includes.R')
setwd('figures/fig1/B')

#------------- PARAMETERS -----------------
M <- rbind(c(-1,0,0,0),c(0,0.296854227781296,0.954922616481781,0),c(0,0.954922795295715,-0.296854168176651,0),c(0,0,0,1))

d <- 0.397303260

AP.list <- c(-2.88)
color.list <- c('red', 'blue', 'green')
list.allen <- c('root', 'TH','HY','MB','HB',
                'CB','STR','PAL','Isocortex','OLF',
                'HIP','RHP', 'CTXsp','fiber tracts','VS')
list.allen <- 'root'

obs <- list(x = 0, y = -0.3, z = 12)
zoom <- 0.55
M <- rbind(c(-1,0,0,0),c(0,0.296854227781296,0.954922616481781,0),c(0,0.954922795295715,-0.296854168176651,0),c(0,0,0,1))

#------------- Getting outlines -----------------

allen.id <- sapply(list.allen, get.id.from.acronym)
allen.color <- color.from.acronym(list.allen)
allen.color[1] <- 'black'

mesh.list <- sapply(allen.id, mesh3d.allen.annot.from.id, simplify = FALSE)

outlines.list <- NULL

for(i in 1:length(AP.list)){
  
  l.plane <- get.plane.equation(AP.list[i])
  N <- l.plane$N
  d <- l.plane$d
  rm(l.plane)
  
  outlines.list[[i]] <- get.meshes.list.plane.intersection(mesh.list, 0, N, d)

}

#----------- V2 ----------- 

mesh3d.new.window(T)
ids.to.keep <- rgl.ids()
rgl.pop(id = setdiff(rgl.ids()$id, ids.to.keep))

for(i in 1:length(outlines.list)){
  for(j in 1:length(outlines.list[[i]])){

    df.cur.outline <- data.frame(x = numeric(),
                                 y = numeric(),
                                 z = numeric())

    o <- outlines.list[[i]][[j]]
    if(is.null(o))
      next
    if(length(o) == 0)
      next

    for(k in 1:length(o)){

      if(is.null(o[[k]]))
        next

      df.cur.outline <- rbind(df.cur.outline,
                              o[[k]])

    }

    segments3d(df.cur.outline, col = 'blue', lwd = 2.5)

  }
}

par3d(userMatrix = M)
par3d(zoom = zoom)
observer3d(0,-0.3,13)
rgl.snapshot('section-3d-outline.png')
