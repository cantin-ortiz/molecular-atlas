#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/mov3')

#------------- LOADING -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

#------------- Plotting ----------------
  
lf <- list.files(paste(path.matrices, 'smoothed-atlas/all-genes', sep = '/'), pattern = ('^(prediction).*(.RData)$'), full.names = T)
ap.lf <- as.numeric(gsub('mm', '.', unlist(strsplit(unlist(strsplit(unlist(strsplit(lf, '/'))[(1:length(lf))*7],'_'))[(1:length(lf))*2],'[.]'))[(1:length(lf))*2-1]))
ord <- order(ap.lf, decreasing = T)

ap.lf <- ap.lf[ord]
lf <- lf[ord]
cnt <- 0
dir.create('movie')
for(i in 1:length(lf)){

  cnt <- cnt + 1
  f <- lf[i]
  
  png(sprintf('movie/movie%04d.png',cnt),
    width = 1920,
    height = 1080,
    res = 300,
    type = 'cairo')

  load(f)

  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  plot.2d.mixed(array(grid.2d, dim = c(nrow(grid.2d), ncol(grid.2d),1)),
                1, df.colors, grid.parameters, right.hemisphere = T, increase.factor = 1,
                show.box = FALSE, show.title = FALSE, pdf.name = NULL, show.saggital = F)
  
  plot.saggital.view(-5.25, 0.2, -1, ap.lf[i])
  dev.off()
  
  tmp <- cnt
  for(j in 1:19){
    cnt <- cnt + 1
    file.copy(from = sprintf('movie/movie%04d.png',tmp),
              to = sprintf('movie/movie%04d.png',cnt))
  }
}
