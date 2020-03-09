#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig3/F')

#--------------- PARAMETERS --------------- 

path.svm <- paste(path.matrices, 'smoothed-atlas/all-genes/', sep = '/')

#--------------- LOADING --------------- 

df.cl.id <- load.df.cl.id.name(cl.id.names.path)

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

all.ap <- unique(spots.table$AP)

#--------------- *** 2/ Plotting spots with vivid colors --------------- 

AP.list.to.use <- c(1.80,0.79,-1.33,-2.34,-3.15,-3.76)

for(cur.AP in AP.list.to.use){
  
  param <- strsplit(as.character(sprintf('%.2f', cur.AP)), '[.]')[[1]]
  file.to.load <- sprintf('%s/predictions_%smm%s.RData', path.svm, param[1], param[2])
  load(file.to.load)
  
  df.colors.isocortex <- df.colors.vivid
  df.colors.isocortex[!startsWith(rownames(df.colors.isocortex), 'Isocortex'), 'clusters.colors'] <- gray.color

  closest.ap <- which.min(abs(all.ap - grid.parameters$AP))
  sp.section <- subset(spots.table, AP == all.ap[closest.ap])
  
  df.spots <- df.colors.isocortex
  df.spots$clusters.colors[df.spots$clusters.colors == gray.color] <- NA
  
  sp.section$color <- mapvalues(as.character(sp.section$clusters.named),
                                from = rownames(df.spots),
                                to = df.spots$clusters.colors,
                                warn_missing = F)
  
  sp.section$ML <- -sp.section$ML
  sp.empty <- subset(sp.section, is.na(color))
  sp.spots <- subset(sp.section, !is.na(color))
  
  outlines <- atlas.stereo$outlines[[which.min(abs(atlas.stereo$AP - cur.AP))]]
  df.segments <- get.segments.outlines(outlines, 'ML', 'DV', 1)  
  
  pdf(generate.appropriate.file.name(sprintf('isocortex-all-in-one-vivid_%s.pdf', paste(param, collapse='mm'))),
      width = 10, height = 7, useDingbats = F)
  
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  
  plot(1, type = 'n', asp = 1, ylim = c(-8,0), xlim = c(-6,6),
       bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='')
  
  plot.raster.atlas(grid.2d, grid.parameters, df.colors.isocortex, show.edges = F, increase.factor = 1, file.path = NULL, add = TRUE, right.hemisphere = T)
  
  x.coeff <- 1
  segments(x.coeff*df.segments$x0, df.segments$y0, x.coeff*df.segments$x1, df.segments$y1, lwd = lwd.atlas.segments)
  x.coeff <- -1
  segments(x.coeff*df.segments$x0, df.segments$y0, x.coeff*df.segments$x1, df.segments$y1, lwd = lwd.atlas.segments)
  
  if(nrow(sp.spots)>1)
    symbols(sp.spots$ML, sp.spots$DV, circles = rep(0.05, nrow(sp.spots)), inches = F, add = T, bg = sp.spots$color, fg = NA)
  
  if(nrow(sp.empty)>1)
    symbols(sp.empty$ML, sp.empty$DV, circles = rep(0.05, nrow(sp.empty)), inches = F, add = T, bg = NA, fg = 'black', lwd = 0.5)
  dev.off()
  
}
