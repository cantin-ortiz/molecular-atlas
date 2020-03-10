#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/sup03')

#------------------- Parameters  ------------------- 

cutoff <- 20

#------------------- Loadings  ------------------- 

load(seurat.object.path)
ic.mat <- seur.obj@dr$fiftypercents@gene.loadings
t <- ic.mat[,ic.kept]

#-------------------  ------------------- 

pdf(generate.appropriate.file.name('ic-loads.pdf'), paper = 'a4', width = 9, height = 12, useDingbats = F)
par(mfcol= c(8,6),
    mar = c(1,5,0.75,1),
    mgp = c(1,0.15,0), #3 1 0
    xpd = NA)


for(i in 1:ncol(t)){
  vect <- t[,i]
  vect <- vect[order(abs(vect), decreasing = T)]
  sub.vect <- vect[1:cutoff]
  
  plot(sort(sub.vect),1:cutoff, pch = 16, xaxt = 'n', yaxt= 'n', xlab = '', ylab = '', bty = 'n', 
       xlim = c(-max(abs(sub.vect)), +max(abs(sub.vect))), main = colnames(t)[i], cex.main = 1, cex = 0.8)
  axis(1, at = c(-max(abs(sub.vect)), +max(abs(sub.vect))), tick = T, labels = NA, cex.axis=0.75, lwd.ticks = 0, lwd = 0.5)
  axis(1, cex.axis=0.6, tck=-0.05, lwd = 0.5, lwd.tick = 0.5, lab = NA)
  axis(1, cex.axis=0.6, tck=-0.05, lwd = 0, lwd.tick = 0, pos = 1)
  axis(2, at = 1:cutoff, tck=-0.025, las = 2, labels = names(sub.vect[order(sub.vect)]), cex.axis = 0.6, lwd = 0.5, lwd.tick = 0.5)
  segments(0,0,0,20, lty = 2, lwd = 0.5)
}

dev.off()
