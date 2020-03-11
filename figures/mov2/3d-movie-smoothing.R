# ------------ INCLUDES --------------

source('bin/includes.R')
setwd("figures/mov2")

# ------------ LOADING --------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()),cl.file, min.cluster.size = 10)
df.cl.id <- load.df.cl.id.name(cl.id.names.path)
load(path.all.meshes)
l.meshes <- l
rm(l)

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
  
  #Adding point
  for (i in 1:length(adding.points.ts)){
    if (is.same.ts(t,adding.points.ts[i])){
      
      mesh <- adding.points.data[[i]]
      mesh$vb[1,] <- -mesh$vb[1,]
      smooth <<- vcgSmooth(mesh, 'HC', 5)
      col <<- adding.points.color[[i]]
      
      alpha.list <<- alpha.base
      cur.ids <- rgl.ids()$id
      
      wire3d(smooth, col = col, material = list(shininess = 128, specular = 'gray30', alpha = alpha.list[2]))
      alpha.list <<- alpha.list[2:length(alpha.list)]
      shape.being.added <<- setdiff(rgl.ids()$id, cur.ids)
      
    }
  }
  
  #Looping through the sphere table
  for (i in 1:nrow(df.spheres.add)){
    
    #Adding sphere instantly if relevant
    if (is.same.ts(t,df.spheres.add[i,'add.time'])){
      
      cur.id <- rgl.ids()$id
      sub.sp <<- subset(spots.table, cluster == df.spheres.add[i, 'cluster'])
      spheres3d(sub.sp$ML, sub.sp$DV, sub.sp$AP, radius = 0.075, col = df.spheres.add[i, 'color'])
      df.spheres.add[i,'id'] <<- setdiff(rgl.ids()$id, cur.id)
    }
    
    #Deleting sphere instantly if relevant
    if (is.same.ts(t,df.spheres.add[i,'delete.time'])){
      rgl.pop(id = df.spheres.add[i,'id'])
    }
  }
  
  if(length(alpha.list) != 0){
    
    par3d(skipRedraw=TRUE) 
    rgl.pop(id = shape.being.added)
    
    cur.ids <- rgl.ids()$id
    alpha <- alpha.list[1]
    
      wire3d(smooth, col = col, material = list(shininess = 128, specular = 'gray30', alpha = alpha))
    par3d(skipRedraw=FALSE) 
    alpha.list <<- setdiff(alpha.list, alpha.list[1])    
    shape.being.added <<- setdiff(rgl.ids()$id, cur.ids)
    
  }
  
  
  #User matrix in case no modification
  if (is.null(userMatrix))
    userMatrix <- M
  
  return(list(userMatrix = userMatrix))
}

# ------------ Parameters --------------

M0 <- rbind(c(0.984246373176575,0,0.176802426576614,0),c(0.0990536957979202,0.828323066234589,-0.55142480134964,0),c(-0.146449521183968,0.560250759124756,0.815274000167847,0),c(0,0,0,1))
M1 <- rbind(c(0.0744546800851822,0,0.997224390506744,0),c(0.519449412822723,0.853620588779449,-0.0387830883264542,0),c(-0.851251244544983,0.52089524269104,0.063556045293808,0),c(0,0,0,1))
M2 <- rbind(c(0.90315043926239,0,-0.429324269294739,0),c(-0.196877881884575,0.888655126094818,-0.414163261651993,0),c(0.381521224975586,0.458576172590256,0.802589237689972,0),c(0,0,0,1))
M3 <- rbind(c(-0.0805743187665939,0,-0.996748626232147,0),c(0.108674563467503,0.994038581848145,-0.00878494139760733,0),c(0.990806579589844,-0.109029054641724,-0.0800939798355103,0),c(0,0,0,1))
M4 <- rbind(c(-1,0,-0,0),c(-0,0,1,0),c(0,1,-0,0),c(0,0,0,1))
M5 <- rbind(c(0.845470309257507,0,0.534022450447083,0),c(0.104178078472614,0.980786979198456,-0.164935901761055,0),c(-0.523762285709381,0.195081830024719,0.82922625541687,0),c(0,0,0,1))

zoom = 0.75

alpha.base <<- exp(seq(from = 0, to = 4, length.out = 61))
alpha.base <<- alpha.base / max(alpha.base)

# alpha.base <- c(0,1)


# ------------ CODE -------------

cl.list.chr <- c('Hippocampal region-02', 'Midbrain-06', 'Striatum-06', 'Hypothalamus-02', 'Thalamus-01', 'Striatum-04','Olfactory areas-03')
cl.list <- as.numeric(mapvalues(cl.list.chr,
                                from = df.cl.id$clusters.named,
                                to = df.cl.id$cluster,
                                warn_missing = F))

adding.points.data <- NULL
adding.points.color <- NULL

for(i in 1:length(cl.list)){
  
  cl <- cl.list[i]
  adding.points.data[[i]] <-l.meshes$data[[cl]]
  adding.points.color[[i]] <- df.colors.vivid[df.colors.vivid$cluster.id == cl,'clusters.colors']
}

list.actions <- list(
  interpolate = list(st = 4, et = 6, M1 = M0, M2 = M1),  
  interpolate = list(st = 16, et = 19, M1 = M1, M2 = M2),
  interpolate = list(st = 26, et = 29, M1 = M2, M2 = M3),
  interpolate = list(st = 33, et = 36, M1 = M3, M2 = M4),
  interpolate = list(st = 38, et = 41, M1 = M4, M2 = M5)
)

for(i in 1:length(list.actions)){
  
  l <- list.actions[[i]]
  
  if (names(list.actions)[i] == 'interpolate')
    list.actions[[i]]$f <- par3dinterp(times = c(l$st,l$et),userMatrix = c(l$M1,l$M2), method = 'spline', extrapolate = 'constant')
  
  if (names(list.actions)[i] == 'spin')
    list.actions[[i]]$f <- spin3d(axis = l$axis, rpm = l$rpm)
}


#------------------------------------------ MOVIE ------------------------------------------ 

dir.create('movie')

adding.points.ts <- c(2,8,11,14,21,24,31)

df.spheres.add <- data.frame(
  add.time = adding.points.ts - 1,
  delete.time = adding.points.ts+1,
  id = rep(NA,length(adding.points.data)),
  cluster = cl.list,
  color = df.colors.vivid[cl.list.chr, 'clusters.colors'],
  stringsAsFactors = F
)

sphere.being.deleted <- NA
shape.being.added <- NA

mesh3d.new.window()
rgl.viewpoint(userMatrix = M1)
par3d(zoom = zoom)

M <- M0

smooth <- NULL
alpha.list <- numeric()
alpha.list.sphere <- numeric()
col <- character()
col.sphere <- character()

movie3d( function.rotation, 
         duration = 43,
         dir = 'movie',
         fps = 60,
         convert = FALSE,
         type = 'png',
         top = FALSE)
