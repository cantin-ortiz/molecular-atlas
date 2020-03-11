#--------------- INCLUDES ---------------  

rm(list=ls())
source('bin/includes.R')
setwd('figures/mov1')

#--------------- PARAMETERS --------------- 

path.svm <- paste(path.matrices, 'smoothed-atlas/all-genes', sep = '/')
all.files <- list.files(path.svm)
AP.list <- sort(as.numeric(gsub('mm','.',unlist(strsplit(unlist(strsplit(all.files, '[.]'))[seq(from = 1, to = 2*length(all.files), by = 2)],'_'))[seq(from = 2, to = 2*length(all.files), by = 2)])), decreasing = T)

#--------------- LOADING --------------- 

df.cl.id <- load.df.cl.id.name(cl.id.names.path)
spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
AP.list.spots <- sort(unique(spots.table$AP), decreasing = T)


df.colors.vivid$adjusted.color <- df.colors.vivid$clusters.colors

while(length(unique(df.colors.vivid$adjusted.color)) != nrow(df.colors.vivid)){
    
  cnt <- plyr::count(df.colors.vivid$adjusted.color)
  cnt$x <- as.character(cnt$x)
  
  for(col in cnt[cnt$freq > 1, 'x']){
    sel <- rownames(subset(df.colors.vivid, adjusted.color == col))
    sel.ordered <- sel[order(df.colors.vivid[sel, 'cluster.id'])]
    bits <- sapply(0:(length(sel.ordered)-1), intToBits)
    cols.change <- matrix(nrow = 3, ncol = ncol(bits))
    cols.change[2,] <- as.numeric(bits[2,])
    cols.change[3,] <- as.numeric(bits[3,])
    cols.change[1,] <- as.numeric(bits[1,]) + as.numeric(bits[4,])
    
    new.cols <- t(col2rgb((df.colors.vivid[sel.ordered,'adjusted.color'])) + 1*cols.change)
    while(max(new.cols[,1]) > 255)
      new.cols[,1] <- new.cols[,1] - 1
    while(max(new.cols[,2]) > 255)
      new.cols[,2] <- new.cols[,2] - 1
    while(max(new.cols[,3]) > 255)
      new.cols[,3] <- new.cols[,3] - 1
    
    df.colors.vivid[sel.ordered, 'adjusted.color'] <- rgb(new.cols, maxColorValue = 255)
    
  }
}


#--------------- Plotting spots with vivid colors --------------- 

dir.create('movie')
counter <- 0
for(cur.AP in AP.list.spots){
  
  counter <- counter + 1

  png(sprintf('movie/movie%04d.png',counter),
      width = 1920,
      height = 1080,
      res = 300,
      type = 'cairo')
  
  par(mar= c(2.2, 1.2, 1.2, 1.2),
      xpd=TRUE)

  plate.id <- which.min(abs(cur.AP - atlas.stereo$AP))
  sp.section <- subset(spots.table, AP == cur.AP)

  outlines <- atlas.stereo$outlines[[plate.id]]

  plot(1, type = 'n', asp = 1, ylim = c(-8,0), xlim = c(-6,5.3),
       bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='')
 
    for(i in 1:length(outlines)){
      polygon(outlines[[i]]$ML, outlines[[i]]$DV)
      polygon(-outlines[[i]]$ML, outlines[[i]]$DV, col = atlas.stereo$plate.info[[plate.id]][i,'col'])
    }
  
  symbols(sp.section$ML, sp.section$DV, circles = rep(0.05, nrow(sp.section)), inches = F, add = T, bg = df.colors.vivid[as.character(sp.section$clusters.named), 'adjusted.color'], fg = NA)
  
  axis(1, at = -5:5, pos = -8)
  axis(2, at = 0:-8, pos = -5.5, las = 2) 
  mtext('Medial-lateral (mm)', side = 1, at = 0, line = 1.25)
  mtext('Dorsal-ventral (mm)', side = 2, at = -7.75/2, line = -4)
  mtext(sprintf('Bregma = %.2f mm', cur.AP), side = 3, at = 0)

  dev.off()
  
  tmp <- counter
  for(j in 1:19){
    counter <- counter + 1
    file.copy(from = sprintf('movie/movie%04d.png',tmp),
              to = sprintf('movie/movie%04d.png',counter))
  }
}
  
