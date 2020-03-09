#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig5/H')

# ------------- Parameters ----------------

path.cluster.266 <- paste(path.matrices, '181-clusters-266-genes.tsv', sep = '/')

ap.min <- -3.88
ap.max <- -1.055

# ------------- Loading ----------------

spots.table <- add.parent.acronym(load.spots.table())
spots.table.original <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
spots.table.266 <- append.cluster.to.spots.table(spots.table, path.cluster.266, min.cluster.size = 10)

# ------------- Computing ----------------
  
clusters.266 <- as.character(unique(spots.table.266$clusters.named))
hipp.266 <- clusters.266[startsWith(clusters.266, 'Hippocampal')]
palette.clusters <- hipp.266

clusters.all <- character(length(hipp.266))

a <- 0
for(cl in hipp.266){
  a <- a +1 
  cnt <- plyr::count(spots.table.original[rownames(subset(spots.table.266, clusters.named == cl)),'clusters.named'])
  clusters.all[a] <- as.character(cnt$x)[which.max(cnt$freq)]
}


lf.all <- sort(list.files(paste(path.matrices, 'smoothed-atlas/all-genes', sep = '/'), full.names = T))
lf.266 <- sort(list.files(paste(path.matrices, 'smoothed-atlas/266-genes', sep = '/'), full.names = T))

ap.files <- as.numeric(gsub('mm', '.', unlist(strsplit(unlist(strsplit(lf.all, '_'))[seq(from = 2, to = 2*length(lf.all), by = 2)], '[.]'))[seq(from = 1, to = 2*length(lf.all), by = 2)]))

df.colors.original <- df.colors.vivid
df.colors.original[!is.element(rownames(df.colors.original), clusters.all), 'clusters.colors'] <- gray.color

original.clusters <- clusters.all

df.colors.palette <- unique(spots.table.266[,c('cluster','clusters.named')])
rownames(df.colors.palette) <- df.colors.palette$clusters.named
df.colors.palette <- df.colors.palette[order(df.colors.palette$cluster),]
df.colors.palette$clusters.named <- NULL
df.colors.palette$clusters.colors <- gray.color

df.colors.original[!is.element(rownames(df.colors.original),original.clusters), 'clusters.colors'] <- gray.color

for(cl in palette.clusters){
  cnt <- plyr::count(spots.table.original[rownames(subset(spots.table.266, clusters.named == cl)),'clusters.named'])
  df.colors.palette[cl, 'clusters.colors'] <- df.colors.original[as.character(cnt[which.max(cnt$freq), 'x']),'clusters.colors']
}
colnames(df.colors.palette)[1] <- 'cluster.id'


# ------------- Plotting ----------------

loop.param <- which(ap.files > ap.min & ap.files < ap.max)

selected <- c(21, 30, 37)

for(i in selected){
  
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
                    df.colors.palette,
                    increase.factor = 2,
                    file.path = NULL,
                    add = TRUE,
                    right.hemisphere = FALSE)

  dev.off()

}

