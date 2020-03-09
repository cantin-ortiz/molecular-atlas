#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/fig2/C')

#------------------- Parameters  ------------------- 

camera <- rbind(c(-0.584145545959473,0,0.811648964881897,0),c(0.237170696258545,0.95635461807251,0.170692265033722,0),c(-0.776224255561829,0.292208462953568,-0.558650314807892,0),c(0,0,0,1))

#------------------- Loadings  ------------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)

t <- read.table(tsne.3d.path, sep='\t', stringsAsFactors = F, header = TRUE, row.names = 1)

#------------------- Plotting  ------------------- 

sel.rows <- dplyr::intersect(rownames(t),rownames(spots.table))
spots.table[sel.rows,'tsne1'] <- t[sel.rows,1]
spots.table[sel.rows,'tsne2'] <- t[sel.rows,2]
spots.table[sel.rows,'tsne3'] <- t[sel.rows,3]
spots.table <- spots.table[sel.rows,]

df.colors.tsne <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
spots.table$colors.tsne <- mapvalues(as.numeric(as.character(spots.table$cluster)),
                                     from = as.numeric(df.colors.tsne$cluster.id),
                                     to = df.colors.tsne$clusters.colors)

spots.table$red <- col2rgb(spots.table$colors.tsne)[1,]/256
spots.table$green <- col2rgb(spots.table$colors.tsne)[2,]/256
spots.table$blue <- col2rgb(spots.table$colors.tsne)[3,]/256

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
spheres3d(spots.table[,'tsne1'],spots.table[,'tsne2'],spots.table[,'tsne3'], radius = 0.2, col = spots.table$colors.tsne)

cur.seg <- red.segment
segments3d(cbind(cur.seg[1:(length(cur.seg)-1)], cur.seg[2:(length(cur.seg))]),min(spots.table$tsne2),min(spots.table$tsne3), col = red.col, lwd = 10)
cur.seg <- green.segment
segments3d(min(spots.table$tsne1), cbind(cur.seg[1:(length(cur.seg)-1)], cur.seg[2:(length(cur.seg))]),min(spots.table$tsne3), col = green.col, lwd = 10)
cur.seg <- blue.segment
segments3d(min(spots.table$tsne1), min(spots.table$tsne2), cbind(cur.seg[1:(length(cur.seg)-1)], cur.seg[2:(length(cur.seg))]), col = blue.col, lwd = 10)


rgl.bbox(color=c('gray40'), emission="gray60",
         specular="gray1", shininess=5, alpha=0.7,
         xlen = 0, ylen = 0, zlen = 0, expand = 1) 
observer3d(0,-4,200)
rgl.snapshot('3d-tsne-illustration.png')
