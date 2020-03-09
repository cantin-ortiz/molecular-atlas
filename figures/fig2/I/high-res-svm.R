#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig2/I/')

library(raster)

#------------- LOADING -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

#------------- *** Plotting ----------------

lf <- c(paste(path.matrices,'smoothed-atlas/all-genes/predictions_1mm80.RData',sep = '/'),
        paste(path.matrices,'smoothed-atlas/all-genes/predictions_0mm39.RData',sep = '/'),
        paste(path.matrices,'smoothed-atlas/all-genes/predictions_-0mm62.RData',sep = '/'),
        paste(path.matrices,'smoothed-atlas/all-genes/predictions_-3mm56.RData',sep = '/')) 

cnt <- 0

for(f in lf){
  cnt <- cnt + 1
  load(f)
  pdf(generate.appropriate.file.name(paste('coronal',cnt,'.pdf',sep='')),
      width = 10, height = 7, useDingbats = F)
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  plot.2d.mixed(array(grid.2d, dim = c(nrow(grid.2d), ncol(grid.2d),1)),
                1, df.colors, grid.parameters, right.hemisphere = T, increase.factor = 1,
                show.box = FALSE, show.title = FALSE, pdf.name = NULL, show.saggital = F)
  dev.off()
}
