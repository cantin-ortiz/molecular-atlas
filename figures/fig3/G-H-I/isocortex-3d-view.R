# ------------------- Includes -------------------

source('bin/includes.R')
setwd('figures/fig3/G-H-I')

# ------------------- Parameters -------------------

M.list <- list(
  rbind(c(0.652781248092651,0,0.757546484470367,0),c(0.279805481433868,0.929287374019623,-0.241109654307365,0),c(-0.703978359699249,0.369357496500015,0.606621384620667,0),c(0,0,0,1)),
  rbind(c(0.999997913837433,0,-0.0020315118599683,0),c(7.15604473953135e-05,0.999379396438599,0.0352251417934895,0),c(0.00203025108203292,-0.0352252162992954,0.999377310276031,0),c(0,0,0,1))
)
zoom.list <- c(0.37, 0.35)
observer.list <- list(c(0,-0.25,13),
                      c(0,0,13))


# ------------------- Loadings -------------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.cl.id <- load.df.cl.id.name(cl.id.names.path)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
load(path.all.meshes)
df.cl.id <- load.df.cl.id.name(cl.id.names.path)

# ------------------- Plot -------------------
# ------------------- *** Somatosentory  -------------------

M.list <- list(
  rbind(c(0.578163266181946,0,-0.815921068191528,0),c(-0.459940314292908,0.825974941253662,-0.325914621353149,0),c(0.67393034696579,0.563706874847412,0.477548360824585,0),c(0,0,0,1)),
  rbind(c(1,0,-0,0),c(0,1,0,0),c(0,-0,1,0),c(0,0,0,1)),
  rbind(c(-0,0,-1,0),c(-1,0,0,0),c(0,1,-0,0),c(0,0,0,1))
)

observer.list <- list(c(0,-0.25,13),
                      c(0,-0.25,13),
                      c(0,0,13))
zoom.list <- c(0.4,0.32,0.44)
acr.list <- list(c('SS','MOp', 'PTLp'),
                 c('SS','MOp', 'PTLp'),
                 c('SS','MOp', 'PTLp'))
hemipshere.list <- list('both',
                        'right',
                        'both')

to.plot.chr <- c('Isocortex-02', 'Isocortex-03','Isocortex-06','Isocortex-15', 'Isocortex-18', 'Isocortex-21', 'Isocortex-25', 'Isocortex-34')
to.plot <- as.numeric(mapvalues(to.plot.chr,
                                to = df.cl.id$cluster,
                                from = df.cl.id$clusters.named,
                                warn_missing = F))

for(cam.id in 1:length(M.list)){
  
  mesh3d.new.window(T)
  cnt <- 0

  for (cl in to.plot){
    
    cnt <- cnt + 1
    id <- which(l$clusters == cl)
    if(length(id) == 0)
      next
    mesh <- l$data[[id]]
    
    #Deleting one ventral artefact
    mesh.to.del <- which(mesh$vb[2,] < -5)
    if (length(mesh.to.del) > 0 ){
      it.to.del <- is.element(mesh$it[1,], mesh.to.del) | is.element(mesh$it[2,], mesh.to.del) | is.element(mesh$it[3,], mesh.to.del)
      mesh$it <- mesh$it[,!it.to.del]
    }
  
    mesh$vb[1,] <- -mesh$vb[1,]
    
    smooth <- vcgSmooth(mesh, 'HC', 5)
    wire3d(smooth, col = df.colors.vivid[df.cl.id[df.cl.id$cluster == cl, 'clusters.named'],'clusters.colors'], material = list(shininess = 128, specular = 'gray30'))
  }
  
  for(acr in acr.list[[cam.id]]){

    col <- color.from.acronym(acr)
    mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
    
    cur.mesh <- mesh
    vertices.discarded <- which(cur.mesh$vb[1,] <0)
    discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
    cur.mesh$it <- cur.mesh$it[,!discarded.it]
    mesh.left <- cur.mesh
    
    cur.mesh <- mesh
    vertices.discarded <- which(cur.mesh$vb[1,] >0)
    discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
    cur.mesh$it <- cur.mesh$it[,!discarded.it]
    mesh.right <- cur.mesh
    
    if(is.element(hemipshere.list[cam.id],c('right','both'))){
      shade3d(mesh.right, color = col,
              alpha = 0.4, override = F, material = list(shininess = 50,
                                                         ambiant = col,
                                                         emission = col,
                                                         lit = FALSE))
    }
    
    if(is.element(hemipshere.list[cam.id],c('left','both'))){
      wire3d(mesh.left, col =col, material = list(shininess = 128, specular = 'grey30', alpha = 0.2))
    }
  }
  
  par3d(userMatrix = M.list[[cam.id]],
        zoom = zoom.list[cam.id])
  observer3d(observer.list[[cam.id]][1],
             observer.list[[cam.id]][2],
             observer.list[[cam.id]][3])
  rgl.snapshot(paste('isocortex-3d-view-SS-camera-',cam.id, '.png', sep = ''))
  
  rgl.close()
}


# ------------------- *** RSP  -------------------

M.list <- list(
  rbind(c(0.7438557,0,-0.6683403,0),c(-0.4914499,0.6777104,-0.5469786,0),c(0.4529411,0.7353289,0.5041187,0),c(0,0,0,1)),
  rbind(c(1,0,-0,0),c(0,1,0,0),c(0,-0,1,0),c(0,0,0,1)),
  rbind(c(-0,0,-1,0),c(-1,0,0,0),c(0,1,-0,0),c(0,0,0,1))
)

observer.list <- list(c(0,-0.25,13),
                      c(0,-0.25,13),
                      c(0,0,13))
zoom.list <- c(0.47,0.32,0.44)
acr.list <- list(c('RSP','POST', 'PRE'),
                 c('RSP','POST', 'PRE'),
                 c('RSP','POST', 'PRE'))
hemipshere.list <- list('both',
                        'right',
                        'both')

to.plot.chr <- c('Isocortex-11', 'Isocortex-23','Isocortex-27','Isocortex-37')
to.plot <- as.numeric(mapvalues(to.plot.chr,
                                to = df.cl.id$cluster,
                                from = df.cl.id$clusters.named,
                                warn_missing = F))

for(cam.id in 1:length(M.list)){
  
  mesh3d.new.window(T)
  cnt <- 0
  
  for (cl in to.plot){
    
    cnt <- cnt + 1
    id <- which(l$clusters == cl)
    if(length(id) == 0)
      next
    mesh <- l$data[[id]]
    
    #Deleting one ventral artefact
    mesh.to.del <- which(mesh$vb[2,] < -5)
    if (length(mesh.to.del) > 0 ){
      it.to.del <- is.element(mesh$it[1,], mesh.to.del) | is.element(mesh$it[2,], mesh.to.del) | is.element(mesh$it[3,], mesh.to.del)
      mesh$it <- mesh$it[,!it.to.del]
    }
    
    mesh$vb[1,] <- -mesh$vb[1,]
    
    smooth <- vcgSmooth(mesh, 'HC', 5)
    wire3d(smooth, col = df.colors.vivid[df.cl.id[df.cl.id$cluster == cl, 'clusters.named'],'clusters.colors'], material = list(shininess = 128, specular = 'gray30'))
  }
  
  for(acr in acr.list[[cam.id]]){
    
    # if(is.null(acr))
    #   next
    
    col <- color.from.acronym(acr)
    mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
    
    cur.mesh <- mesh
    vertices.discarded <- which(cur.mesh$vb[1,] <0)
    discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
    cur.mesh$it <- cur.mesh$it[,!discarded.it]
    mesh.left <- cur.mesh
    
    cur.mesh <- mesh
    vertices.discarded <- which(cur.mesh$vb[1,] >0)
    discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
    cur.mesh$it <- cur.mesh$it[,!discarded.it]
    mesh.right <- cur.mesh
    
    if(is.element(hemipshere.list[cam.id],c('right','both'))){
      shade3d(mesh.right, color = col,
              alpha = 0.4, override = F, material = list(shininess = 50,
                                                         ambiant = col,
                                                         emission = col,
                                                         lit = FALSE))
    }
    
    if(is.element(hemipshere.list[cam.id],c('left','both'))){
      wire3d(mesh.left, col =col, material = list(shininess = 128, specular = 'grey30', alpha = 0.2))
    }
  }
  
  par3d(userMatrix = M.list[[cam.id]],
        zoom = zoom.list[cam.id])
  observer3d(observer.list[[cam.id]][1],
             observer.list[[cam.id]][2],
             observer.list[[cam.id]][3])
  rgl.snapshot(paste('isocortex-3d-view-RSP-camera-',cam.id, '.png', sep = ''))
  
  rgl.close()
}


# ------------------- *** Anterior  -------------------

M.list <- list(
  rbind(c(-0.866702675819397,0,-0.498825073242188,0),c(-0.347282767295837,0.717846393585205,0.603399693965912,0),c(0.358079791069031,0.696201503276825,-0.622159361839294,0),c(0,0,0,1)),
  rbind(c(-1,0,0,0),c(-0,1,-0,0),c(-0,-0,-1,0),c(0,0,0,1)),
  rbind(c(-0,0,-1,0),c(-1,0,0,0),c(0,1,-0,0),c(0,0,0,1))
)

observer.list <- list(c(0,-0.75,13),
                      c(0,-0.35,13),
                      c(0,0,13))
zoom.list <- c(0.49,0.47,0.44)
acr.list <- list(c('ACA', 'PL', 'ILA', 'ORB','MOs','MOp'),
                 c('ACA', 'PL', 'ILA', 'ORB','MOs','MOp'),
                 c('ACA', 'PL', 'ILA', 'ORB','MOs','MOp'))
hemipshere.list <- list('right',
                        'right',
                        'both')

to.plot.chr <- c('Isocortex-04', 'Isocortex-05', 'Isocortex-10','Isocortex-14','Isocortex-29','Isocortex-30','Isocortex-33')
# to.plot.chr <- c('Isocortex-29', 'Isocortex-33', 'Isocortex-14')
to.plot <- as.numeric(mapvalues(to.plot.chr,
                                to = df.cl.id$cluster,
                                from = df.cl.id$clusters.named,
                                warn_missing = F))

for(cam.id in 1:length(M.list)){
  
  mesh3d.new.window(T)
  cnt <- 0
  
  for (cl in to.plot){
    
    cnt <- cnt + 1
    id <- which(l$clusters == cl)
    if(length(id) == 0)
      next
    mesh <- l$data[[id]]
    
    #Deleting one ventral artefact
    mesh.to.del <- which(mesh$vb[2,] < -5)
    if (length(mesh.to.del) > 0 ){
      it.to.del <- is.element(mesh$it[1,], mesh.to.del) | is.element(mesh$it[2,], mesh.to.del) | is.element(mesh$it[3,], mesh.to.del)
      mesh$it <- mesh$it[,!it.to.del]
    }
    
    mesh$vb[1,] <- -mesh$vb[1,]
    
    smooth <- vcgSmooth(mesh, 'HC', 5)
    wire3d(smooth, col = df.colors.vivid[df.cl.id[df.cl.id$cluster == cl, 'clusters.named'],'clusters.colors'], material = list(shininess = 128, specular = 'gray30'))
  }
  
  for(acr in acr.list[[cam.id]]){

    col <- color.from.acronym(acr)
    mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
    
    cur.mesh <- mesh
    vertices.discarded <- which(cur.mesh$vb[1,] <0)
    discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
    cur.mesh$it <- cur.mesh$it[,!discarded.it]
    mesh.left <- cur.mesh
    
    cur.mesh <- mesh
    vertices.discarded <- which(cur.mesh$vb[1,] >0)
    discarded.it <- apply(cur.mesh$it, 2, function(x){return(any(is.element(x, vertices.discarded)))})
    cur.mesh$it <- cur.mesh$it[,!discarded.it]
    mesh.right <- cur.mesh
    
    if(is.element(hemipshere.list[cam.id],c('right','both'))){
      shade3d(mesh.right, color = col,
              alpha = 0.4, override = F, material = list(shininess = 50,
                                                         ambiant = col,
                                                         emission = col,
                                                         lit = FALSE))
    }
    
    if(is.element(hemipshere.list[cam.id],c('left','both'))){
      wire3d(mesh.left, col =col, material = list(shininess = 128, specular = 'grey30', alpha = 0.2))
    }
  }
  
  par3d(userMatrix = M.list[[cam.id]],
        zoom = zoom.list[cam.id])
  observer3d(observer.list[[cam.id]][1],
             observer.list[[cam.id]][2],
             observer.list[[cam.id]][3])
  rgl.snapshot(paste('isocortex-3d-view-anterior-camera-',cam.id, '.png', sep = ''))
  
  rgl.close()
}
