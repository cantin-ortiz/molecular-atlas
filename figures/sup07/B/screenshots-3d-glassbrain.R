# ------------ INCLUDES --------------

library(Rvcg)
source('bin/includes.R')
setwd('figures/sup07/B')

# ------------ PARAMETERS --------------

vector.order <- c(1,2,3)
data.dir <- path.matrices
selected <- '10_0.3_correlation'

n.itt <- 5
M.view1 <- rbind(c(0.821424901485443,0,-0.570316672325134,0),c(-0.174864053726196,0.951835632324219,-0.251856029033661,0),c(0.542847752571106,0.306608706712723,0.781861484050751,0),c(0,0,0,1))
M.saggital <- rbind(c(0,0,-1,0),c(0,1,0,0),c(1,-0,0,0),c(0,0,0,1))
M.horizontal <- rbind(c(0,0,-1,0),c(-1,0,0,0),c(0,1,0,0),c(0,0,0,1))

set.obs <- function(view){
  if(view == 'view1'){
    M <- M.view1
    zoom <- 0.4
    observer3d(0,-0.33,13)
  }
  if(view == 'saggital'){
    M <- M.saggital
    zoom <- 0.36
    observer3d(0.1,-0.15,13)
  }
  if(view == 'horizontal'){
    M <- M.horizontal
    zoom <- 0.43
    observer3d(0.1,0,13)
  }
  
  par3d(zoom = zoom,
        userMatrix = M)
}



# ------------ LOADING --------------

spots.table.raw <- add.parent.acronym(load.spots.table())
spots.table.0 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 0)
spots.table.10 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 10)

# ------------ DF COLORS --------------

path.sel <- paste(data.dir,'/UMAP_3_',selected,'.txt',sep='')

#Reading the tsne
t <-  read.table(path.sel,
                 sep =  '\t',
                 row.names = 1,
                 header = T)
colnames(t) <- c('tsne1', 'tsne2', 'tsne3')
rownames(t) <- rownames(spots.table.0)

#Converting into rgb proporitions
t1 <- t[, vector.order[1]]
t2 <- t[, vector.order[2]]
t3 <- t[, vector.order[3]]

#Normalizing vectors between 0 and 1
vect.1 <- (t1 - min(t1)) / abs(max(t1) - min(t1))
vect.2 <- (t2 - min(t2)) / abs(max(t2) - min(t2))
vect.3 <- (t3 - min(t3)) / abs(max(t3) - min(t3))

#Initializing data frame
df.colors <- unique(spots.table.10[, c('clusters.named', 'cluster')])
df.colors$clusters.named <- as.character(df.colors$clusters.named)
rownames(df.colors) <- df.colors$clusters.named
df.colors$clusters.named <- NULL
colnames(df.colors) <- 'cluster.id'

#Initializing data frame
df.colors <- unique(spots.table.10[, c('clusters.named', 'cluster')])
df.colors$clusters.named <- as.character(df.colors$clusters.named)
rownames(df.colors) <- df.colors$clusters.named
df.colors$clusters.named <- NULL
colnames(df.colors) <- 'cluster.id'

#Filling it with median coordinates of the cluster
for (cl in rownames(df.colors)) {
  sel.spots <- is.element(rownames(t), rownames(spots.table.10[spots.table.10$clusters.named == cl, ]))
  df.colors[cl, 'clusters.colors'] <- rgb(median(vect.1[sel.spots]),
                                          median(vect.2[sel.spots]),
                                          median(vect.3[sel.spots]))
}

df.colors <- df.colors[order(rownames(df.colors)),]

spots.table <- spots.table.10
spots.table$colors <- as.character(mapvalues(spots.table$cluster,
                                   from = df.colors$cluster.id,
                                   to = df.colors$clusters.colors))

# ------------ PLOT --------------

mesh3d.new.window(T)
spheres3d(-spots.table$ML, spots.table$DV, spots.table$AP, col = spots.table$colors,
          radius = 0.04, alpha = 1)

set.obs('view1')
rgl.snapshot('clusters-3d-spots-view1.png')

set.obs('saggital')
rgl.snapshot('clusters-3d-spots-saggital.png')

set.obs('horizontal')
rgl.snapshot('clusters-3d-spots-horizontal.png')

rgl.close()
