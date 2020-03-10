#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/sup07/B')

#------------------- Parameters  ------------------- 

vector.order <- c(1,2,3)
data.dir <- path.matrices
selected <- '10_0.3_correlation'
camera <- rbind(c(-0.584145545959473,0,0.811648964881897,0),c(0.237170696258545,0.95635461807251,0.170692265033722,0),c(-0.776224255561829,0.292208462953568,-0.558650314807892,0),c(0,0,0,1))

#------------------- Loadings  ------------------- 

spots.table.raw <- add.parent.acronym(load.spots.table())
spots.table.0 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 0)
spots.table.10 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 10)

#------------------- Appending tsne coordinates to spots table  ------------------- 

path.sel <- paste(data.dir,'/UMAP_3_',selected,'.txt',sep='')

#Reading the tsne
t <-  read.table(path.sel,
                 sep =  '\t',
                 row.names = 1,
                 header = T)
colnames(t) <- c('tsne1', 'tsne2', 'tsne3')
rownames(t) <- rownames(spots.table.0)

sel.rows <- dplyr::intersect(rownames(t),rownames(spots.table.10))
spots.table.10[sel.rows,'tsne1'] <- t[sel.rows,1]
spots.table.10[sel.rows,'tsne2'] <- t[sel.rows,2]
spots.table.10[sel.rows,'tsne3'] <- t[sel.rows,3]
spots.table.10 <- spots.table.10[sel.rows,]

#------------------- Creating color data frame  ------------------- 

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

spots.table.10$colors.tsne <- mapvalues(as.character(spots.table.10$clusters.named),
                                        from = rownames(df.colors),
                                        to = df.colors$clusters.colors)

spots.table.10$red <- col2rgb(spots.table.10$colors.tsne)[1,]/256
spots.table.10$green <- col2rgb(spots.table.10$colors.tsne)[2,]/256
spots.table.10$blue <- col2rgb(spots.table.10$colors.tsne)[3,]/256

#------------------- Plotting  ------------------- 

spots.table <- spots.table.10

red.segment <- seq(from = min(spots.table$tsne1), to = max(spots.table$tsne1), length.out = 258)
green.segment <- seq(from = min(spots.table$tsne2), to = max(spots.table$tsne2), length.out = 258)
blue.segment <- seq(from = min(spots.table$tsne3), to = max(spots.table$tsne3), length.out = 258)
red.col <- rgb(0:256/256,0,0)
green.col <- rgb(0,0:256/256,0)
blue.col <- rgb(0,0,0:256/256)
           
           
rgl.open()
rgl.bg(color = 'white')
par3d(windowRect = c(0, 0, 1920, 1080))
par3d(zoom = 0.2)
par3d(userMatrix = camera)
spheres3d(spots.table[,'tsne1'],spots.table[,'tsne2'],spots.table[,'tsne3'], radius = 0.05, col = spots.table$colors.tsne)

cur.seg <- red.segment
segments3d(cbind(cur.seg[1:(length(cur.seg)-1)], cur.seg[2:(length(cur.seg))]),min(spots.table$tsne2),min(spots.table$tsne3), col = red.col, lwd = 10)
cur.seg <- green.segment
segments3d(min(spots.table$tsne1), cbind(cur.seg[1:(length(cur.seg)-1)], cur.seg[2:(length(cur.seg))]),min(spots.table$tsne3), col = green.col, lwd = 10)
cur.seg <- blue.segment
segments3d(min(spots.table$tsne1), min(spots.table$tsne2), cbind(cur.seg[1:(length(cur.seg)-1)], cur.seg[2:(length(cur.seg))]), col = blue.col, lwd = 10)


rgl.bbox(color=c('gray40'), emission="gray60",
         specular="gray1", shininess=5, alpha=0.7,
         xlen = 0, ylen = 0, zlen = 0, expand = 1) 
observer3d(0,-1,50)
rgl.snapshot('3d-umap-illustration.png')
rgl.close()

