# ------------------- Includes -------------------

source('bin/includes.R')
setwd('figures/fig3/O')

# ------------------- Parameters -------------------

M.list <- list(
  rbind(c(0.404086917638779,0,-0.914720594882965,0),c(-0.0524234883487225,0.998356401920319,0.0231585968285799,0),c(0.913217186927795,-0.0573109313845634,0.403422772884369,0),c(0,0,0,1)),
  rbind(c(1,0,0,0),c(0,0.965925812721252,-0.258819043636322,0),c(0,0.258819043636322,0.965925812721252,0),c(0,0,0,1)),
  rbind(c(0,0,1,0),c(1,0,0,0),c(0,1,0,0),c(0,0,0,1))
  )
zoom.list <- c(0.5,0.4,0.44)

observer.list <- list(c(0,-0.2,13),
                      c(0,-0.5,13),
                      c(0,0,13))
                      
show.ara <- c(FALSE, TRUE, TRUE)
  
col.str <- '#4d63ab'   

# ------------------- Loadings -------------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.cl.id <- load.df.cl.id.name(cl.id.names.path)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
load(path.all.meshes)

all.str.clusters <- levels(spots.table$clusters.named)[startsWith(levels(spots.table$clusters.named), 'Striatum')]

groups <- list(all.str.clusters,
               all.str.clusters,
               all.str.clusters)

# ------------------- Initial computation -------------------

allen.acronyms <- 'STR'
allen.ids <- as.numeric(sapply(allen.acronyms, get.id.from.acronym))

allen.meshes <- NULL
for(i in 1:length(allen.ids)){
  allen.meshes[[i]] <- mesh3d.allen.annot.from.id(allen.ids[i])
}

allen.mesh.one.hemisphere <- list()
for(i in 1:length(allen.meshes)){
  
  cur.mesh <- allen.meshes[[i]]
  vertices.discarded <- which(cur.mesh$vb[1,] >0)
  
  discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
  
  cur.mesh$it <- cur.mesh$it[,!discarded.it]
  
  allen.mesh.one.hemisphere[[i]] <- cur.mesh
}

# ------------------- Plot -------------------

for(cam.id in 1:length(groups)){
  
  cnt <- 0
  
  mesh3d.new.window(T)
  
  to.plot <- as.numeric(mapvalues(groups[[cam.id]],
                                  to = df.cl.id$cluster,
                                  from = df.cl.id$clusters.named,
                                  warn_missing = F))
  for (cl in to.plot){
    
    cnt <- cnt + 1
    id <- which(l$clusters == cl)
    if(length(id) == 0)
      next
    mesh <- l$data[[id]]
    
    mesh$vb[1,] <- -mesh$vb[1,]
    
    smooth <- vcgSmooth(mesh, 'HC', 5)
    
    to.discard <- c(which(smooth$vb[2,] > -2),
                    which(smooth$vb[1,] <  0.4 & smooth$vb[3,] < 0.5 & smooth$vb[2,] < -5))
    
    if(length(to.discard) > 0){
      disc.col <- which(is.element(smooth$it[1,],to.discard) | is.element(smooth$it[2,],to.discard) | is.element(smooth$it[3,],to.discard))
      smooth$it <- smooth$it[,setdiff(1:ncol(smooth$it), disc.col)]
    }

    wire3d(smooth, col = df.colors.vivid[df.cl.id[df.cl.id$cluster == cl, 'clusters.named'],'clusters.colors'], material = list(shininess = 128, specular = 'gray30'))
    
    
  }
  
  if(show.ara[cam.id]){
    shade3d(allen.mesh.one.hemisphere[[i]], color = col.str,
            alpha = 0.4, override = F, material = list(shininess = 50,
                                                       ambiant = col.str,
                                                       emission = col.str,
                                                       lit = T))
  }
  
  par3d(userMatrix = M.list[[cam.id]],
        zoom = zoom.list[cam.id])
  observer3d(observer.list[[cam.id]][1],
             observer.list[[cam.id]][2],
             observer.list[[cam.id]][3])
  rgl.snapshot(paste('striatum-3d-view-', cam.id, '.png', sep = ''))
  
  rgl.close()
}

