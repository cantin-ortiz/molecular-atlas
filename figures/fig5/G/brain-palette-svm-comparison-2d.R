#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig5/G')

# ------------- Parameters ----------------

path.cluster.266 <- paste(path.matrices, '181-clusters-266-genes.tsv', sep = '/')
palette.clusters <- c('Isocortex-31','Isocortex-08','Retrohippocampal region-07','Isocortex-03','Isocortex-26','Isocortex-19','Retrohippocampal region-06','Isocortex-24','Retrohippocampal region-01','Retrohippocampal region-03','Isocortex-12','Isocortex-16','Isocortex-27','Isocortex-36','Isocortex-06','Isocortex-32','Isocortex-20','Mixed-05 (Isocortex)','Isocortex-13','Isocortex-04','Retrohippocampal region-04','Isocortex-11','Isocortex-14','Isocortex-01','Isocortex-15','Mixed-22 (Isocortex)','Mixed-15 (Isocortex)','Olfactory areas-04','Olfactory areas-01','Olfactory areas-06','Olfactory areas-02','Olfactory areas-07')
orignal.clusters <- c('Olfactory areas-10','Olfactory areas-02','Olfactory areas-06','Olfactory areas-09','Retrohippocampal region-09','Olfactory areas-03','Olfactory areas-05','Olfactory areas-01','Isocortex-10','Isocortex-02','Isocortex-25','Isocortex-17','Isocortex-36','Isocortex-27','Retrohippocampal region-03','Isocortex-08','Isocortex-26','Isocortex-12','Mixed-05 (Cortical subplate)','Isocortex-07','Isocortex-41','Isocortex-06','Isocortex-21','Isocortex-40','Isocortex-28','Isocortex-16','Isocortex-34','Retrohippocampal region-02','Retrohippocampal region-07','Retrohippocampal region-05','Isocortex-11','Retrohippocampal region-06','Isocortex-23','Isocortex-37','Isocortex-14','Isocortex-30','Retrohippocampal region-08','Isocortex-04','Isocortex-29')

# ------------- Loading ----------------

spots.table <- add.parent.acronym(load.spots.table())
spots.table.original <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
spots.table.266 <- append.cluster.to.spots.table(spots.table, path.cluster.266, min.cluster.size = 10)

# ------------- Computing ----------------
  
lf.all <- sort(list.files(paste(path.matrices, 'smoothed-atlas/all-genes', sep = '/'), full.names = T))
lf.266 <- sort(list.files(paste(path.matrices, 'smoothed-atlas/266-genes', sep = '/'), full.names = T))

df.colors.original <- df.colors.vivid
df.colors.original[is.na(df.colors.original$clusters.colors) & is.element(rownames(df.colors.original),orignal.clusters),'clusters.colors'] <- 'black'

df.colors.palette <- unique(spots.table.266[,c('cluster','clusters.named')])
rownames(df.colors.palette) <- df.colors.palette$clusters.named
df.colors.palette <- df.colors.palette[order(df.colors.palette$cluster),]
df.colors.palette$clusters.named <- NULL
df.colors.palette$clusters.colors <- gray.color

df.colors.original[!is.element(rownames(df.colors.original),orignal.clusters), 'clusters.colors'] <- gray.color

for(cl in palette.clusters){
  cnt <- plyr::count(spots.table.original[rownames(subset(spots.table.266, clusters.named == cl)),'clusters.named'])
  df.colors.palette[cl, 'clusters.colors'] <- df.colors.original[as.character(cnt[which.max(cnt$freq), 'x']),'clusters.colors']
}
colnames(df.colors.palette)[1] <- 'cluster.id'


# ------------- Plotting ----------------

loop.param <- c(18,23,36,45,63,73)

for(i in loop.param){
  
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

