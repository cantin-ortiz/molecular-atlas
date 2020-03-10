#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/sup6')

#------------- PARAMETERS -----------------

data.path <- paste(path.matrices, 'smoothed-atlas/all-genes/', sep = '/')

#------------- LOADING -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

#------------- PREDICTIONS -----------------

lf <- list.files(data.path, pattern = ('^(prediction).*(.RData)$'))
lf.ap <- as.numeric(gsub('mm', '.', unlist(strsplit(unlist(strsplit(lf, '_'))[seq(from = 2, to = 2*length(lf), by = 2)], '[.]'))[seq(from = 1, to = 2*length(lf), by = 2)]))

spots.ap <- sort(unique(spots.table$AP), decreasing = T)
lf.selected <- lf[sapply(spots.ap, function(x){return(which.min(abs(lf.ap-x)))})]

cnt <- 0

for(f in lf.selected){
  cnt <- cnt + 1
  load(paste(data.path, f, sep = ''))
  pdf(generate.appropriate.file.name(sprintf('coronal%02d.pdf',cnt)),
      width = 10, height = 7, useDingbats = F)
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  plot.2d.mixed(array(grid.2d, dim = c(nrow(grid.2d), ncol(grid.2d),1)),
                1, df.colors, grid.parameters, right.hemisphere = T, increase.factor = 2,
                show.box = FALSE, show.title = FALSE, pdf.name = NULL, show.saggital = F)
  dev.off()
}
