#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig5/F')


# ------------- Parameters ----------------

path.cluster.266 <- paste(path.matrices, '181-clusters-266-genes.tsv', sep = '/')
ap.list <- c(-5.28,-4.47,-4.06,-3.56,+2.92,-2.24,+2.21,-1.53,+1.50,-1.03,-0.52,+0.29)

# ------------- Loading ----------------

spots.table <- add.parent.acronym(load.spots.table())
spots.table.original <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
spots.table.266 <- append.cluster.to.spots.table(spots.table, path.cluster.266, min.cluster.size = 10)

df.colors.original <- get.color.from.3d.tsne(tsne.3d.path, spots.table.original)
df.colors.266 <- get.color.from.3d.tsne(tsne.3d.path, spots.table.266)

# ------------- Plotting ----------------
  
lf.all <- sort(list.files(paste(path.matrices, 'smoothed-atlas/all-genes', sep = '/'), full.names = T))
lf.266 <- sort(list.files(paste(path.matrices, 'smoothed-atlas/266-genes', sep = '/'), full.names = T))

for(i in 1:(length(lf.all))){
  
  f.original <- lf.all[i]
  f.palette <- lf.266[i]
  
  coord <- strsplit(f.original, '_')[[1]][2]
  coord <- strsplit(coord, '[.]')[[1]][1]
  
  pdf(generate.appropriate.file.name(sprintf('palette-svm-%s.pdf', coord)),
      width = 10, height = 7, useDingbats = F)
  
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  
  plot(1, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
       asp = 1, xlim = c(-6,6), ylim = c(-8,0))
  
  load(f.original)
  plot.raster.atlas(grid.2d, 
                    grid.parameters, 
                    df.colors.original, 
                    increase.factor = 2,
                    file.path = NULL, 
                    add = TRUE, 
                    right.hemisphere = TRUE)
  load(f.palette)
  plot.raster.atlas(grid.2d, 
                    grid.parameters, 
                    df.colors.266, 
                    increase.factor = 2,
                    file.path = NULL, 
                    add = TRUE, 
                    right.hemisphere = FALSE)
  
  dev.off()

}

