#------------- INCLUDES -----------------

library(alphahull)
source('bin/includes.R')
setwd('figures/fig4/B')


#------------- PARAMETERS -----------------

dir.predictions <-  paste(path.matrices, 'smoothed-atlas/all-genes/', sep = '/')
ground.truth.v1 <- clusters.v1
ground.truth.alm <- clusters.alm

alm.file <- 'predictions_1mm40.RData'
v1.file <- 'predictions_-3mm56.RData'

files <- c(alm.file, v1.file)

#------------- LOADING -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)

#------------- PLOTTING vivid color -----------------

fnames <- paste('ground-truth-2d-vivid-',
               c('alm', 'v1'),
               '.pdf',
               sep = '')

for(i in 1:2){
  
  load(paste(dir.predictions, files[i], sep=''))
  
  pdf(generate.appropriate.file.name(fnames[i]), width = 10, height = 7, useDingbats = F)

  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  
  df.colors <- df.colors.vivid
  if(i == 1)
    df.colors[!is.element(df.colors$cluster.id, ground.truth.alm), 'clusters.colors'] <- NA
  if (i == 2)
    df.colors[!is.element(df.colors$cluster.id, ground.truth.v1), 'clusters.colors'] <- NA

  stereo.plate <- which.min(abs(grid.parameters$AP - atlas.stereo$AP))
  plate.infos <- atlas.stereo$plate.info[[stereo.plate]]
  outlines <- atlas.stereo$outlines[[stereo.plate]]
  
  alm.ara.contour <- which(startsWith(plate.infos$acronym, 'MOs'))
  v1.ara.contour <- grep('^(VISp)([[:digit:]]).*', plate.infos$acronym)
  
  alm.seg <- NULL
  v1.seg <- NULL
  if(length(alm.ara.contour) > 0)
    alm.seg <- get.segments.outlines(outlines[alm.ara.contour], 'ML', 'DV')
  if(length(v1.ara.contour) > 0)
    v1.seg <- get.segments.outlines(outlines[v1.ara.contour], 'ML', 'DV')
  
  df.segments <- get.segments.outlines(outlines, 'ML', 'DV', 1)  
  x.coeff <- 1
  
  df.colors[is.na(df.colors$clusters.colors), 'clusters.colors'] <- gray.color
  
  plot(1, type = 'n', asp = 1, ylim = c(-8,0), xlim = c(-6,6), 
       main = '', bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='') 
  
  plot.raster.atlas(grid.2d, grid.parameters, df.colors, increase.factor = 1, show.edges = F, file.path = NULL, add = TRUE, right.hemisphere = T)
  
  #Plotting the contours
  segments(x.coeff*df.segments$x0, df.segments$y0, x.coeff*df.segments$x1, df.segments$y1, lwd = lwd.atlas.segments)

  if(!is.null(alm.seg)){
    all.points <- unique(data.frame(x = c(alm.seg$x0, alm.seg$x1), y = c(alm.seg$y0, alm.seg$y1)))
    p <- ahull(all.points$x, all.points$y, alpha = 20)
    plot(p, wpoints = F, add = T, lwd = 4, col = '#3A3292')
  }
  
  if(!is.null(v1.seg)){
    all.points <- unique(data.frame(x = c(v1.seg$x0, v1.seg$x1), y = c(v1.seg$y0, v1.seg$y1)))
    p <- ahull(all.points$x, all.points$y, alpha = 5)
    plot(p, wpoints = F, add = T, lwd = 4, col = 'red')
  }

  dev.off()

}
