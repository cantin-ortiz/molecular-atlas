#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig2/H')
library(parallel)
library(foreach)
library(itertools)
library(raster)

#------------- PARAMETERS -----------------

cluster.list.1 <- c('Isocortex-04','Isocortex-05','Isocortex-09','Isocortex-10','Isocortex-12','Isocortex-14','Isocortex-29','Isocortex-30','Isocortex-33')
cluster.list.2 <- c('Striatum-04','Striatum-06','Isocortex-02','Isocortex-06','Isocortex-35','Olfactory areas-03','Striatum-08','Hypothalamus-02','Isocortex-08')
cluster.list.3 <- c('Thalamus-01', 'Hippocampal region-02', 'Cortical subplate-02','Hypothalamus-06', 'Thalamus-03', 'Isocortex-07', 'Isocortex-11', 'Thalamus-08', 'Thalamus-06', 'Striatum-03')
cluster.list.4 <- c('Hippocampal region-04', 'Hippocampal region-06', 'Mixed-03 (Pallidum)', 'Midbrain-06', 'Midbrain-08', 'Midbrain-15', 'Hippocampal region-09', 'Retrohippocampal region-04', 'Hippocampal region-07', 'Hippocampal region-02', 'Hippocampal region-10', 'Retrohippocampal region-10','Retrohippocampal region-01')

cl.list <- list(cluster.list.1, cluster.list.2, cluster.list.3, cluster.list.4)

path.list <- c(paste(path.matrices,'smoothed-atlas/all-genes/predictions_1mm80.RData',sep = '/'),
               paste(path.matrices,'smoothed-atlas/all-genes/predictions_0mm39.RData',sep = '/'),
               paste(path.matrices,'smoothed-atlas/all-genes/predictions_-1mm53.RData',sep = '/'),
               paste(path.matrices,'smoothed-atlas/all-genes/predictions_-3mm56.RData',sep = '/'))
  
#------------- LOADING -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

#------------- Plotting ----------------

for(i in 1:length(cl.list)){
  
  load(path.list[i])
  
  df.colors.isocortex <- df.colors.vivid
  df.colors.isocortex[setdiff(rownames(df.colors.isocortex), cl.list[[i]]),'clusters.colors'] <- gray.color
  
  all.ap <- unique(spots.table$AP)
  closest.ap <- which.min(abs(all.ap - grid.parameters$AP))
  sp.section <- subset(spots.table, AP == all.ap[closest.ap])
  
  df.spots <- df.colors.isocortex
  sp.section$color <- mapvalues(as.character(sp.section$clusters.named),
                                from = rownames(df.spots),
                                to = df.spots$clusters.colors,
                                warn_missing = F)
  
  # sp.section$ML <- sp.section$ML
  sp.empty <- subset(sp.section, color == gray.color)
  sp.spots <- subset(sp.section, color != gray.color)
  
  pdf(generate.appropriate.file.name(paste('smoothing-illustration-',i,'.pdf',sep='')),
      width = 10, height = 7, useDingbats = F)
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  plot.2d.mixed(array(grid.2d, dim = c(nrow(grid.2d), ncol(grid.2d),1)),
                1, df.colors.isocortex, grid.parameters, right.hemisphere = T, increase.factor = 1,
                show.box = FALSE, show.title = FALSE, pdf.name = NULL, show.saggital = F, fill.atlas = F)
  
  symbols(-sp.spots$ML, sp.spots$DV, circles = rep(0.05, nrow(sp.spots)), inches = F, add = T, bg = sp.spots$color, fg = NA)
  symbols(-sp.empty$ML,  sp.empty$DV, circles = rep(0.05, nrow(sp.empty)), inches = F, add = T, bg = NA, fg = 'black', lwd = 0.5)
  dev.off()
  
}
