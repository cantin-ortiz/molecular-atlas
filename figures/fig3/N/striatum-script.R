#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig3/N')

#--------------- PARAMETERS --------------- 

path.svm <- paste(path.matrices, 'smoothed-atlas/all-genes/', sep = '/')
all.files <- list.files(path.svm)

AP.list <- sort(as.numeric(gsub('mm','.',unlist(strsplit(unlist(strsplit(all.files, '[.]'))[seq(from = 1, to = 2*length(all.files), by = 2)],'_'))[seq(from = 2, to = 2*length(all.files), by = 2)])), decreasing = T)

cex.spots <- 0.6

sel.AP <- c(0.09, -0.32, -0.62, 0.79, -0.93, -1.13, 1.40, -1.74, 1.80)

#--------------- LOADING --------------- 

df.cl.id <- load.df.cl.id.name(cl.id.names.path)

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

all.ap <- unique(spots.table$AP)
AP.list.sp <- AP.list[sapply(unique(spots.table$AP), function(x){return(which.min(abs(x-AP.list)))})]

#Alternative 2
df.colors.vivid <- unique(spots.table[,c('cluster','clusters.named')])
id.first <- which(startsWith(as.character(df.colors.vivid$clusters.named), 'Striatum'))
df.colors.vivid <- df.colors.vivid[c(id.first, setdiff(1:nrow(df.colors.vivid), id.first)),]
df.colors.vivid[1:length(id.first),] <- df.colors.vivid[order(df.colors.vivid$clusters.named[1:length(id.first)]),]
rownames(df.colors.vivid) <- as.character(df.colors.vivid$clusters.named)
colnames(df.colors.vivid) <- c('cluster.id', 'clusters.colors')
df.colors.vivid$clusters.colors <- 'lightgray'
df.colors.vivid$clusters.colors[1:length(id.first)] <- c('#E58606','#5D69B1','#99C945','#CC61B0','#24796C','#2F8AC4','#764E9F','#ED645A','#CC3A8E','#A5AA99', brewer.pal(10, 'Set3'))[1:length(id.first)]


#--------------- *** 2/ Plotting spots with vivid colors --------------- 

for(cur.AP in sel.AP){
  
  param <- strsplit(as.character(sprintf('%.2f', cur.AP)), '[.]')[[1]]
  file.to.load <- sprintf('%s/predictions_%smm%s.RData', path.svm, param[1], param[2])
  load(file.to.load)
  
  df.colors.striatum <- df.colors.vivid
  df.colors.striatum[!startsWith(rownames(df.colors.striatum), 'Striatum'), 'clusters.colors'] <- gray.color
  
  closest.ap <- which.min(abs(all.ap - grid.parameters$AP))
  sp.section <- subset(spots.table, AP == all.ap[closest.ap])
  
  df.spots <- df.colors.striatum
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
  
  pdf(generate.appropriate.file.name(sprintf('striatum-all-in-one-vivid_%s.pdf', paste(param, collapse='mm'))),
      width = 10, height = 7, useDingbats = F)
  
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  
  plot(1, type = 'n', asp = 1, ylim = c(-8,0), xlim = c(-6,6),
       bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='')
  
  plot.raster.atlas(grid.2d, grid.parameters, df.colors.striatum, show.edges = F, increase.factor = 1, file.path = NULL, add = TRUE, right.hemisphere = T)
  
  x.coeff <- 1
  segments(x.coeff*df.segments$x0, df.segments$y0, x.coeff*df.segments$x1, df.segments$y1, lwd = lwd.atlas.segments)
  x.coeff <- -1
  segments(x.coeff*df.segments$x0, df.segments$y0, x.coeff*df.segments$x1, df.segments$y1, lwd = lwd.atlas.segments)

  symbols(sp.spots$ML, sp.spots$DV, circles = rep(0.05, nrow(sp.spots)), inches = F, add = T, bg = sp.spots$color, fg = NA)
  symbols(sp.empty$ML, sp.empty$DV, circles = rep(0.05, nrow(sp.empty)), inches = F, add = T, bg = NA, fg = 'black', lwd = 0.5)
  dev.off()
  
}

