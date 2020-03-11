# ------------ INCLUDES --------------

source('bin/includes.R')
setwd("figures/mov5")

# ------------ LOADING --------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()),cl.file, min.cluster.size = 10)

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
      
      mesh$vb[1,] <- -mesh$vb[1,]
      
      smooth <- vcgSmooth(mesh, 'HC', 5)
      
      to.discard <- c(which(smooth$vb[2,] > -2),
                      which(smooth$vb[1,] <  0.4 & smooth$vb[3,] < 0.5 & smooth$vb[2,] < -5))
      
      if(length(to.discard) > 0){
        disc.col <- which(is.element(smooth$it[1,],to.discard) | is.element(smooth$it[2,],to.discard) | is.element(smooth$it[3,],to.discard))
        smooth$it <- smooth$it[,setdiff(1:ncol(smooth$it), disc.col)]
      }
  
      wire3d(smooth, col =adding.points.color[[i]], material = list(shininess = 128, specular = 'gray30'))

    }
  }
  
  #User matrix in case no modification
  if (is.null(userMatrix))
    userMatrix <- M
  
  
  return(list(userMatrix = userMatrix))
  
}

# ------------ CODE --------------

M1 <-  rbind(c(0.761785268783569,0,-0.647829592227936,0),c(-0.336158335208893,0.854835331439972,-0.395289868116379,0),c(0.553787648677826,0.518899321556091,0.651200950145721,0),c(0,0,0,1))
M2 <- rbind(c(0.739074349403381,0,0.673623859882355,0),c(0.125216752290726,0.982571482658386,-0.137383028864861,0),c(-0.661883592605591,0.185885265469551,0.726193368434906,0),c(0,0,0,1))
M3 <- rbind(c(-0.999984622001648,0,0.0055404300801456,0),c(0.00501331873238087,0.425707370042801,0.904847025871277,0),c(-0.00235860189422965,0.90486091375351,-0.425700813531876,0),c(0,0,0,1))
M4 <- rbind(c(0.841581583023071,0,-0.540130019187927,0),c(-0.0777203291654587,0.989593386650085,-0.121096760034561,0),c(0.534509122371674,0.143891885876656,0.832823574542999,0),c(0,0,0,1))

zoom = 0.75
df.cl <- unique(spots.table[,c('cluster','clusters.named')])

cl.list.2 <- as.numeric(df.cl[startsWith(as.character(df.cl$clusters.named), 'Striatum'),'cluster'])
cl.list.2 <- sort(cl.list.2)

cl.list <- cl.list.2[c(3,2,1,13,7,15,8,10,12,16,14,11,4,5,9,6)]

load(path.all.meshes)

adding.points.data <- NULL
adding.points.color <- NULL

for(i in 1:length(cl.list)){
  
  cl <- cl.list[i]
  
  adding.points.data[[i]] <-l$data[[cl]]
  adding.points.color[[i]] <- df.colors.vivid[df.colors.vivid$cluster.id == cl,'clusters.colors']
}

list.actions <- list(
  interpolate = list(st = 23, et = 26, M1 = M1, M2 = M2),
  interpolate = list(st = 29, et = 32, M1 = M2, M2 = M3),
  interpolate = list(st = 34, et = 37, M1 = M3, M2 = M4)
)

for(i in 1:length(list.actions)){
  
  l <- list.actions[[i]]
  
  if (names(list.actions)[i] == 'interpolate')
    list.actions[[i]]$f <- par3dinterp(times = c(l$st,l$et),userMatrix = c(l$M1,l$M2), method = 'spline', extrapolate = 'constant')
  
  if (names(list.actions)[i] == 'spin')
    list.actions[[i]]$f <- spin3d(axis = l$axis, rpm = l$rpm)
}


adding.points.ts <- c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5,16,17.5,19,20.5,22,27)

mesh.allen <- mesh3d.allen.annot.from.id(id.from.acronym('STR'))

mesh3d.new.window()
rgl.viewpoint(userMatrix = M1)
par3d(zoom = zoom)
wire3d(mesh.allen, col = color.from.acronym('STR'), alpha = 0.2, material =  list(shininess = 128, specular = 'gray30'))

M <- M1

dir.create('movie')

movie3d( function.rotation, 
         duration = 39,
         dir = 'movie',
         fps = 60,
         convert = FALSE,
         type = 'png',
         top = FALSE)

# ---------------- Isocortex  ----------------
