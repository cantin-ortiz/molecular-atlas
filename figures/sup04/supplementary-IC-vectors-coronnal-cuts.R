#--------------- INCLUDES ---------------  
 
source('bin/includes.R')
setwd('figures/sup04')

#--------------- LOADING --------------- 

spots.table <- add.parent.acronym(load.spots.table())

ic.mat <- as.matrix(read.table(paste(path.matrices, 'ic-matrix-scores.tsv', sep = '/'), sep = '\t', row.names = 1, header = T))

spots.table <- spots.table[intersect(rownames(spots.table), rownames(ic.mat)),]

#--------------- PLOTTING --------------- 

ic.list <- paste('IC',ic.kept, sep='')

for(cur.ic in ic.list){
  
  spots.table$ic.load <- NULL
  
  fname.pdf <- paste(cur.ic, '.pdf', sep = '') 
  fname.png <- paste(cur.ic, '.png', sep = '') 
  
  #Chosing 2D AP coordinates for plotting the IC
  ic.vect <- ic.mat[,cur.ic]
  spots.table[names(ic.vect), 'ic.load'] <- ic.vect
  df.max.ic <- data.frame(AP = NULL, ic.mean = NULL)
  for(AP in unique(spots.table$AP)){
    loads <- spots.table[spots.table$AP == AP, 'ic.load']
    df.max.ic <- rbind(df.max.ic, data.frame(AP = AP, ic.mean =  mean(sort(abs(loads), decreasing = T)[1:100])))
  }
  
  sel.AP <- df.max.ic[which.max(df.max.ic$ic.mean),'AP']
 
  atlas.id <- which.min(abs(atlas.stereo$AP - sel.AP))
  outlines <- atlas.stereo$outlines[[atlas.id]]
  df.seg <- get.segments.outlines(outlines, 'ML','DV')
  
  sp <- subset(spots.table, AP == sel.AP)

  #Obtaining the proper color for each spot
  bin.sep <- seq(from = -max(abs(ic.vect))-0.1, to = max(abs(ic.vect))+0.1, length.out = 202)
  sp$ic.bin <- .bincode(sp$ic.load, breaks = bin.sep)
  
  function.cp <- colorRampPalette(c('#0004ff', '#ffffff', '#ff0000'))
  cp <- function.cp(201)
  sp$color <- cp[sp$ic.bin]
  
  for(i in 1:1){
  
    if(i == 1)
      pdf(fname.pdf, useDingbats = F, fillOddEven = T, width = 10, height = 7)
    if(i == 2)
      png(fname.png, width = 10, height = 7, units = 'in', res = 300)
    
    par(mar= c(0.2, 0.2, 0.2, 0.2),
        xpd=TRUE)
    plot(1, type = 'n', xlim = c(-6,6), ylim = c(-8,0), asp = 1, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    symbols(-sp$ML, sp$DV, circles = rep(0.05, nrow(sp)), inches = F, add = T, bg = sp$col, fg = NA)
    segments(-df.seg$x0, df.seg$y0, -df.seg$x1, df.seg$y1, lwd = lwd.atlas.segments)
    
    #Plotting the inside filling
    for (x in 1:length(outlines)){
      
      #Current outline
      o <- outlines[[x]]
      
      #Plotting polygons for colors
      polygon(o$ML, o$DV, col = atlas.stereo$plate.info[[atlas.id]][x,'col'], border = NA)
      
    }
    segments(df.seg$x0, df.seg$y0, df.seg$x1, df.seg$y1, lwd = lwd.atlas.segments) 
  
  
    dev.off()
    
  }
  
  
  #Screenshot from 3D view
  plot.3d.glassbrain.ic(spots.table, ic.mat, ic.name = cur.ic, cex = 3, top.quantile = 0.95, dim = c(-1920,0,0,2000), HD = TRUE)
  plot.3d.glassbrain.setview('3d', zoom = 0.55)
  plot.3d.glassbrain.setview('medial', zoom = 0.5)
  rgl.snapshot(paste(cur.ic,'-',1,'.png',sep=''), top = TRUE)
  plot.3d.glassbrain.setview('dorsal', zoom = 0.55)
  rgl.snapshot(paste(cur.ic,'-',2,'.png',sep=''), top = TRUE)
  plot.3d.glassbrain.setview('3d', zoom = 0.55)
  rgl.snapshot(paste(cur.ic,'-',3,'.png',sep=''), top = TRUE)
  rgl.close()
 
}








