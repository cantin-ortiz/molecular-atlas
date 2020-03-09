#--------------- INCLUDES ---------------  

rm(list=ls())
source('bin/includes.R')
setwd('figures/fig1/F')

#--------------- Parameters --------------- 

df.parameters <- data.frame(gene = c('Gpr88', 'Rora', 'Rasgrf2'),
                            AP = c(-0.655, -2.355, 1.845),
                            col = c('#3333CC',sapply(c('TH'),color.from.acronym), '#1E9D5A'),
                            stringsAsFactors = F)

use.scran <- T

list.um <- list(rbind(c(0.680512845516205,0,-0.732736170291901,0),c(-0.228889971971512,0.949958205223083,-0.212576612830162,0),c(0.69606876373291,0.312377065420151,0.646458745002747,0),c(0,0,0,1)))
list.zoom <- list(0.45)

#--------------- Loading --------------- 

spots.table <- add.parent.acronym(load.spots.table())
all.AP <- sort(unique(spots.table$AP))
load.expr.mat()
st.data.raw <- st.data
rm(st.data)
load(normalized.matrix.path)

#--------------- 2D plots --------------- 

# cur.pal <- colorRampPalette(brewer.pal(9, 'YlOrRd'))(100)

# cur.pal <- colorRampPalette(c('#ffffff','#fcf4c0', '#ff0000'))(100)
# cur.pal <- colorRampPalette(c('#ffffff','#fcf4c0', '#800026'))(100)

cur.pal <- colorRampPalette(c('#ffffff', '#FDF9E0', '#fcf4c0', '#ff0000', '#800026'))(100)

for(i in 1:nrow(df.parameters)){
  
  marker <- df.parameters[i,'gene']
  cur.AP <- df.parameters[i,'AP']
  
  # for(cur.AP in all.AP[all.AP > 0]){

  id <- which.min(abs(cur.AP - atlas.stereo$AP))
  df.seg <- get.segments.outlines(atlas.stereo$outlines[[id]], 'ML', 'DV')
  cur.sp <- subset(spots.table, AP == cur.AP)
  
  if(use.scran){
    expr <- st.data[,marker]
    names(expr) <- rownames(st.data)
  }else{
    expr <- st.data.raw[,marker]
    names(expr) <- rownames(st.data.raw)
  }
  
  breaks <- seq(from = min(expr),
                to = max(expr),
                length.out = 101)
  breaks[1] <- breaks[1] - 0.01
  breaks[length(breaks)] <- breaks[length(breaks)] + 0.01
  
  cur.sp$expr <- expr[rownames(cur.sp)]
  cur.sp$col <- cur.pal[.bincode(cur.sp$expr, breaks)]
  
  
  pdf.name <- paste(marker, '_', gsub('[.]', 'mm',cur.AP), '.pdf', sep = '')
  pdf(generate.appropriate.file.name(pdf.name), width = 10, height = 7, useDingbats = F)
  plot(1, type = 'n', xlim = c(-6,6), ylim = c(-8,0), asp = 1, main = paste(marker, cur.AP, sep = '/'))
  symbols(cur.sp$ML, cur.sp$DV, circles = rep(0.05, nrow(cur.sp)), inches = F, add = T, bg = cur.sp$col, fg = NA)
  segments(df.seg$x0, df.seg$y0, df.seg$x1, df.seg$y1, lwd = 0.5)
  dev.off()
}

#--------------- Legend --------------- 

pdf('legend1.pdf', useDingbats = F)
plot(1, ylim = c(0,4), type = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n', asp = 1)
points(c(1,1,1), c(1,2,3), pch = 19, col = df.parameters$col)
text(1, 1:3, df.parameters$gene, pos = 4)
dev.off()

pdf('legend2.pdf', useDingbats = F)
plot(1, xlim = c(0,1), ylim = c(0,1), asp = 1, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
rect(0:(length(cur.pal)-1)/length(cur.pal), 0.55, 1:length(cur.pal)/length(cur.pal), 0.45, col = cur.pal, border = NA)
rect(0, 0.45, 1, 0.55)
text(c(0,0.5,1), 0.45, c(0,0.5,1), pos = 1, cex = 2)
dev.off()

#------------ 3D ------------
#------------ *** Pre computation --------------

list.spots <- NULL

for(i in 1:nrow(df.parameters)){

  marker <- df.parameters[i,'gene']

  if(use.scran){
    expr <- st.data[,marker]
    names(expr) <- rownames(st.data)
  }else{
    expr <- st.data.raw[,marker]
    names(expr) <- rownames(st.data.raw)
  }

  cutoff <- quantile(expr, 0.95)
  list.spots[[i]] <- expr[expr > cutoff]

}


#------------ *** Gpr88 --------------

M <- rbind(c(0.845373868942261,0,-0.534175097942352,0),c(-0.358341842889786,0.741609215736389,-0.567103981971741,0),c(0.396149188280106,0.67083215713501,0.626937031745911,0),c(0,0,0,1))
zoom <- 0.48

acr.list <- c('STRd', 'STRv')

mesh3d.new.window(T)
par3d(userMatrix = M, zoom = zoom)
observer3d(0,-0.45,13)
for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  # to.del <- which(mesh$vb[1,] < 0)
  to.del <- NA
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  wire3d(mesh, col = color.from.acronym(acr), material = list(shininess = 128, specular = 'grey30'), alpha = 0.2)
}
sub.spots.table <- spots.table[names(list.spots[[1]]),]
spheres3d(x = -sub.spots.table$ML, sub.spots.table$DV, sub.spots.table$AP, radius = 0.05, col = df.parameters[1,'col'])
rgl.snapshot('3d-gpr88.png')


#------------ *** Rora --------------

M <-rbind(c(0.771081924438477,0,0.636735916137695,0),c(0.461576819419861,0.688842594623566,-0.558965682983398,0),c(-0.438610821962357,0.724910914897919,0.531154096126556,0),c(0,0,0,1))
zoom <- 0.48

acr.list <- c('TH')

mesh3d.new.window(T)
par3d(userMatrix = M, zoom = zoom)
observer3d(0,-0.45,13)
for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  # to.del <- which(mesh$vb[1,] < 0)
  to.del <- NA
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  wire3d(mesh, col = color.from.acronym(acr), material = list(shininess = 128, specular = 'grey30'), alpha = 0.2)
}
sub.spots.table <- spots.table[names(list.spots[[2]]),]
spheres3d(x = sub.spots.table$ML, sub.spots.table$DV, sub.spots.table$AP, radius = 0.05, col = df.parameters[2,'col'])
rgl.snapshot('3d-rora.png')

#------------ *** Rasgrf2 --------------

M <-rbind(c(0.91052520275116,0,-0.413453608751297,0),c(-0.139787063002586,0.94111156463623,-0.307845026254654,0),c(0.389105975627899,0.338096112012863,0.856905817985535,0),c(0,0,0,1))
zoom <- 0.42

acr.list <- list()

mesh3d.new.window(T)
par3d(userMatrix = M, zoom = zoom)
observer3d(0,-0.45,13)
for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  # to.del <- which(mesh$vb[1,] < 0)
  to.del <- NA
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  wire3d(mesh, col = color.from.acronym(acr), material = list(shininess = 128, specular = 'grey30'), alpha = 0.2)
}
sub.spots.table <- spots.table[names(list.spots[[3]]),]
spheres3d(x = -sub.spots.table$ML, sub.spots.table$DV, sub.spots.table$AP, radius = 0.05, col = df.parameters[3,'col'])
rgl.snapshot('3d-rasgrf2.png')

#------------ *** All genes --------------

M1 <- rbind(c(0.460364758968353,0,-0.88772988319397,0),c(-0.207692578434944,0.972246408462524,-0.10770657658577,0),c(0.863092184066772,0.233959212899208,0.447587996721268,0),c(0,0,0,1))
M2 <- rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1))
zoom1 <- 0.35

sub.spots.table <- spots.table[unique(names(unlist(list.spots))),]

for(i in 1:length(list.spots)){
  sub.spots.table[names(list.spots[[i]]), 'color'] <- df.parameters[i, 'col']
}

mesh3d.new.window(T)
par3d(userMatrix = M1,
      zoom = zoom1)
observer3d(0,-0.26,13)
spheres3d(x = -sub.spots.table$ML, sub.spots.table$DV, sub.spots.table$AP, radius = 0.05, col = sub.spots.table$color)
rgl.snapshot('3d-all-1.png')

acr.list <- c('STRd', 'STRv', 'TH')

for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  to.del <- which(mesh$vb[1,] < 0)
  # to.del <- NA
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  wire3d(mesh, col = color.from.acronym(acr), material = list(shininess = 128, diffuse = 'white', ambiant = 'white'), alpha = 0.3)
}
par3d(userMatrix = M2)
rgl.snapshot('3d-all-2.png')
