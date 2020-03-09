#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/fig2/A')

#------------------- Parameters  ------------------- 

ic.biological <- c(14,17,69)
ic.technical <- c(15,27,75)

df.plots <- data.frame(
  id = c(14,17,69,15,27,75),
  type = c(rep('biological',3),rep('technical',3)),
  AP = c(1,-2,-1.8,-1.65,-2.88,0.14),
  zoom = c(0.5846796, 0.6446092, 0.6446092, 0.5846796, 0.5303217, 0.4810174))
rownames(df.plots) <- paste('IC', df.plots$id, sep = '')

l.userMatrix <- list(
  rbind(c(-0.636359930038452,-0.0348060056567192,-0.770599842071533,0),c(0.371204376220703,-0.889520943164825,-0.266362696886063,0),c(-0.676197290420532,-0.455554902553558,0.578978657722473,0),c(0,0,0,1)),
  rbind(c(-0.625040054321289,-0.205728679895401,-0.75294703245163,0),c(0.546067237854004,-0.80449914932251,-0.233490154147148,0),c(-0.5577312707901,-0.55711841583252,0.615208268165588,0),c(0,0,0,1)),
  rbind(c(0.721286177635193,0.296135991811752,-0.626120746135712,0),c(0.547816872596741,-0.797064006328583,0.254093855619431,0),c(-0.423816949129105,-0.52627968788147,-0.73714804649353,0),c(0,0,0,1)),
  rbind(c(0.999930858612061,-0.00418952852487564,-0.0100989090278745,0),c(0.00526350736618042,-0.625126540660858,0.780494630336761,0),c(-0.00958330743014812,-0.780501306056976,-0.625067114830017,0),c(0,0,0,1)),
  rbind(c(-0.986751139163971,-0.0400493294000626,-0.157138615846634,0),c(0.0697717741131783,-0.979584097862244,-0.188465297222137,0),c(-0.146384567022324,-0.196935087442398,0.969413936138153,0),c(0,0,0,1)),
  rbind(c(0.956236720085144,0.0803724750876427,-0.281284183263779,0),c(0.0379251353442669,-0.98744785785675,-0.153215914964676,0),c(-0.290072441101074,0.135844811797142,-0.947297930717468,0),c(0,0,0,1)))


df.plots <- df.plots[1:4,]

#------------------- Loadings  ------------------- 

spots.table <- load.spots.table()
all.AP <- sort(unique(spots.table$AP))
load(seurat.object.path)
ic.mat <- get.ic.mat(seur.obj, 'fiftypercents')

#------------------- Plot ------------------- 

for(i in 1:dim(df.plots)[1]){

  fname.2d <- sprintf('ic-%s-%02d-2d.pdf', df.plots[i,'type'], df.plots[i,'id'])
  fname.3d <- sprintf('ic-%s-%02d-3d.png', df.plots[i,'type'], df.plots[i,'id'])
  
  sel.AP <- all.AP[which.min(abs(df.plots[i,'AP']-all.AP))]
  atlas.id <- which.min(abs(atlas.stereo$AP - sel.AP))
  outlines <- atlas.stereo$outlines[[atlas.id]]
  df.seg <- get.segments.outlines(outlines, 'ML','DV')
  
  sp <- subset(spots.table, AP == sel.AP)
  
  #Only keeping shared spots between ic.mat and spots.table
  shared.spots <- intersect(rownames(sp), rownames(ic.mat))
  sp <- sp[shared.spots,]
  
  #Vector of ic expression
  ic.vect <- ic.mat[,rownames(df.plots)[i]]
  sp$ic.load <- ic.vect[rownames(sp)]
  
  #Obtaining the proper color for each spot
  bin.sep <- seq(from = -max(abs(ic.vect))-0.1, to = max(abs(ic.vect))+0.1, length.out = 202)
  sp$ic.bin <- .bincode(sp$ic.load, breaks = bin.sep)

  function.cp <- colorRampPalette(c('#0004ff', '#ffffff', '#ff0000'))
  cp <- function.cp(201)
  sp$color <- cp[sp$ic.bin]
  
  pdf(fname.2d, useDingbats = F, fillOddEven = T, width = 10, height = 7)
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
  

  plot.3d.glassbrain.ic(spots.table, ic.mat, ic.name = rownames(df.plots)[i], cex = 3, top.quantile = 0.95, dim = c(-1920,0,0,2000), HD = TRUE)
  par3d(userMatrix = l.userMatrix[[i]],
        zoom = df.plots[i, 'zoom'])
  rgl.snapshot(fname.3d, top = TRUE)
  rgl.close()
}



#----------------- Plotting legend -----------------------

function.cp <- colorRampPalette(c('#0004ff', '#ffffff', '#ff0000'))
cp <- function.cp(201)

pdf('legend-ic.pdf', useDingbats = F)
plot(1, type = 'n', xlim = c(0,200), ylim = c(0,1), bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
rect(0:(length(cp)-1), 0.45, 1:length(cp), 0.55, border = NA, col = cp)
rect(0, 0.45, length(cp), 0.55, border = 'black', col = NA)
text(c(0,(length(cp))/2,length(cp)), 0.45, c(-1,0,1), pos = 1, cex = 1.5)
text(length(cp)/2, 0.55, 'Normalized expression level', pos = 3, cex = 1.5)
dev.off()