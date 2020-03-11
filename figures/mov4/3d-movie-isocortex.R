# ------------ INCLUDES --------------

source('bin/includes.R')
setwd("figures/mov4")

# ------------ LOADING --------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()),cl.file, min.cluster.size = 10)
df.cl.id <- load.df.cl.id.name(cl.id.names.path)

# ------------ FUNCTION --------------

#Function to make all the rotations
function.rotation <- function(t){
  
  is.same.ts <- function(t,ts){
    return(abs(t-ts) < 1e-4)
  }
  
  userMatrix <- NULL
  
  #Rotations
  for (i in 1:length(list.actions)){
    
    action <- names(list.actions)[i]
    
    if(action == 'interpolate' || action == 'spin'){
    
      l <- list.actions[[i]]
        
      #Start and end time
      st <- l$st
      et <- l$et
      
      if (t >= st & t <= et){
        
        userMatrix <- l$f(t)$userMatrix

        if (is.same.ts(t,et))
          M <<- userMatrix
      
      }
    }
  }
  
  #Adding points
  for (i in 1:length(adding.points.ts)){
    if (is.same.ts(t,adding.points.ts[i])){
      
      mesh <- adding.points.data[[i]]
      smooth <- vcgSmooth(mesh, 'HC', 5)
   
      tmp <- wire3d(smooth, col = adding.points.color[[i]], material = list(shininess = 128, specular = 'gray30', alpha = adding.points.alpha[[i]]))
      to.delete <<- c(tmp, to.delete)
    }
  }
  
  #Deleting points
  for (i in 1:length(delete.points.ts)){
    if (is.same.ts(t,delete.points.ts[i])){
      rgl.pop(id = to.delete)
      to.delete <<- NULL
    }
  }
  
  #User matrix in case no modification
  if (is.null(userMatrix))
    userMatrix <- M
  
  
  return(list(userMatrix = userMatrix))
  
}

# ------------ CODE --------------

load(path.all.meshes)

M1 <- rbind(c(0.0164027772843838,0,-0.999865472316742,0),c(-0.707439124584198,0.706678986549377,-0.0116055281832814,0),c(0.706583917140961,0.707534313201904,0.0115914978086948,0),c(0,0,0,1))
M2 <- rbind(c(-1,0,-0,0),c(-0,0,1,0),c(0,1,-0,0),c(0,0,0,1))
M3 <- rbind(c(-0.595849752426147,0,-0.803095936775208,0),c(-0.322527050971985,0.915813148021698,0.239296019077301,0),c(0.735485792160034,0.401604622602463,-0.545687019824982,0),c(0,0,0,1))
M4 <- rbind(c(1,0,-0,0),c(-0,1,-0,0),c(0,0,1,0),c(0,0,0,1))
M5 <- rbind(c(0.843839049339294,0,-0.536596417427063,0),c(-0.243507966399193,0.891103088855743,-0.382934957742691,0),c(0.478162735700607,0.453800946474075,0.751947581768036,0),c(0,0,0,1))
M6 <- rbind(c(-1,0,0,0),c(0,0.837429583072662,0.546500027179718,0),c(-0.0107711674645543,0.546545267105103,-0.837360322475433,0),c(0,0,0,1))

to.plot.ss <- c('Isocortex-25', 'Isocortex-02','Isocortex-03','Isocortex-34', 'Isocortex-21', 'Isocortex-06', 'Isocortex-18', 'Isocortex-15')
to.plot.rsp <- c('Isocortex-27', 'Isocortex-23','Isocortex-37','Isocortex-11')
to.plot.ant <- c('Isocortex-33', 'Isocortex-29', 'Isocortex-10', 'Isocortex-04', 'Isocortex-30', 'Isocortex-14', 'Isocortex-05')

acr.list.ss <- c('MOp', 'SS', 'PTLp')
acr.list.rsp <- rev(c('RSP','POST', 'PRE'))
acr.list.ant <- c('ILA', 'ORB', 'PL','ACA', 'MOs','MOp')

acr.list <- list(acr.list.ss, acr.list.rsp, acr.list.ant)
to.plot.list <- list(to.plot.ss, to.plot.rsp, to.plot.ant)

adding.points.data <- NULL
adding.points.color <- NULL
adding.points.alpha <- NULL

cnt <- 0

for(i in 1:length(acr.list)){
  
  for(j in 1:length(acr.list[[i]])){
    
    cnt <- cnt + 1
    
    adding.points.data[[cnt]] <- mesh3d.allen.annot.from.id(get.id.from.acronym(acr.list[[i]][j]))
    adding.points.color[cnt] <- color.from.acronym(acr.list[[i]][j])
    adding.points.alpha[cnt] <- 0.2
    
  }
  
  cl.list <- as.numeric(mapvalues(to.plot.list[[i]],
                                  from = df.cl.id$clusters.named,
                                  to = df.cl.id$cluster,
                                  warn_missing = F))
  
  for(j in 1:length(cl.list)){
    
    cnt <- cnt + 1
    
    cl <- cl.list[j]
    mesh <- l$data[[cl]]
    
    
    #Deleting one ventral artefact
    mesh.to.del <- which(mesh$vb[2,] < -5)
    if (length(mesh.to.del) > 0 ){
      it.to.del <- is.element(mesh$it[1,], mesh.to.del) | is.element(mesh$it[2,], mesh.to.del) | is.element(mesh$it[3,], mesh.to.del)
      mesh$it <- mesh$it[,!it.to.del]
    }
    
    
    mesh$vb[1,] <- -mesh$vb[1,]
    
    adding.points.data[[cnt]] <- mesh
    adding.points.color[cnt] <- df.colors.vivid[df.colors.vivid$cluster.id == cl,'clusters.colors']
    adding.points.alpha[cnt] <- 1
  }
}

list.actions <- list(
  interpolate = list(st = 11, et = 14, M1 = M1, M2 = M2),
  interpolate = list(st = 15, et = 18, M1 = M2, M2 = M3),
  interpolate = list(st = 26, et = 29, M1 = M3, M2 = M4),
  interpolate = list(st = 30, et = 33, M1 = M4, M2 = M5),
  interpolate = list(st = 45.5, et = 48.5, M1 = M5, M2 = M6)
)

for(i in 1:length(list.actions)){
  
  l <- list.actions[[i]]
  
  if (names(list.actions)[i] == 'interpolate')
    list.actions[[i]]$f <- par3dinterp(times = c(l$st,l$et),userMatrix = c(l$M1,l$M2), method = 'spline', extrapolate = 'constant')
  
  if (names(list.actions)[i] == 'spin')
    list.actions[[i]]$f <- spin3d(axis = l$axis, rpm = l$rpm)
}

adding.points.ts <- c(1,1.5,2,
                      3,4,5,6,7,8,9,10,
                      20,20.5,21,
                      22,23,24,25,
                      35,35.5,36,36.5,37,37.5,
                      38.5,39.5,40.5,41.5,42.5,43.5,44.5)

delete.points.ts <- c(19,34)

#---------------------- Code ----------------------

mesh3d.new.window()
rgl.viewpoint(userMatrix = M1)
par3d(zoom = 0.75)

M <- M1
to.delete <<- NULL

dir.create('movie')

movie3d( function.rotation, 
         duration = 50,
         dir = 'movie',
         fps = 60,
         convert = FALSE,
         type = 'png',
         top = FALSE)
