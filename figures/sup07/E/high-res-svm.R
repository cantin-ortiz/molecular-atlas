#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/sup07/E')
library(parallel)
library(foreach)
library(itertools)
library(raster)

#------------- PARAMETERS -----------------

vector.order <- c(1,2,3)
# data.dir <- 'C:/Users/MatLab/Desktop/transcripBrainAtlas/figures/Vector atlas SVM/pixels/all-genes-fullHD/'
data.dir <- paste(path.matrices, 'smoothed-atlas/all-genes-HD/pixels', sep = '/')
selected <- '10_0.3_correlation'

sel.files <- c('2mm25','1mm84', '0mm74', '0mm14', '-0mm66', '-1mm05', '-1mm55', '-1mm96', '-2mm25', '-2mm88', '-3mm58', '-3mm78', '-4mm46', '-5mm25')

#------------- LOADING -----------------

spots.table.raw <- add.parent.acronym(load.spots.table())
spots.table.0 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 0)
spots.table.10 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 10)

df.colors.tsne <- get.color.from.3d.tsne(tsne.3d.path, spots.table.10)

#------------- DF COLORS -----------------

path.sel <- paste(data.dir,'/UMAP_3_',selected,'.txt',sep='')

#Reading the tsne
t <-  read.table(path.sel,
                 sep =  '\t',
                 row.names = 1,
                 header = T)
colnames(t) <- c('tsne1', 'tsne2', 'tsne3')
rownames(t) <- rownames(spots.table.0)

#Converting into rgb proporitions
t1 <- t[, vector.order[1]]
t2 <- t[, vector.order[2]]
t3 <- t[, vector.order[3]]

#Normalizing vectors between 0 and 1
vect.1 <- (t1 - min(t1)) / abs(max(t1) - min(t1))
vect.2 <- (t2 - min(t2)) / abs(max(t2) - min(t2))
vect.3 <- (t3 - min(t3)) / abs(max(t3) - min(t3))

#Initializing data frame
df.colors <- unique(spots.table.10[, c('clusters.named', 'cluster')])
df.colors$clusters.named <- as.character(df.colors$clusters.named)
rownames(df.colors) <- df.colors$clusters.named
df.colors$clusters.named <- NULL
colnames(df.colors) <- 'cluster.id'

#Initializing data frame
df.colors <- unique(spots.table.10[, c('clusters.named', 'cluster')])
df.colors$clusters.named <- as.character(df.colors$clusters.named)
rownames(df.colors) <- df.colors$clusters.named
df.colors$clusters.named <- NULL
colnames(df.colors) <- 'cluster.id'

#Filling it with median coordinates of the cluster
for (cl in rownames(df.colors)) {
  sel.spots <- is.element(rownames(t), rownames(spots.table.10[spots.table.10$clusters.named == cl, ]))
  df.colors[cl, 'clusters.colors'] <- rgb(median(vect.1[sel.spots]),
                                          median(vect.2[sel.spots]),
                                          median(vect.3[sel.spots]))
}

df.colors.umap <- df.colors[order(rownames(df.colors)),]


#------------- *** Plotting ----------------
  
spots.table <- spots.table.10

for(f in sel.files){
  
  fname <- paste(data.dir, '/predictions_', f, '.RData', sep = '')
  load(fname)
  
  pdf(generate.appropriate.file.name(sprintf('coronal-umap-tsne-%s.pdf', f)),
      width = 10, height = 7, useDingbats = F)

  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)

  plot(1, type = 'n', asp = 1, ylim = c(-8,0), xlim = c(-6,6),
       bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='')

  plot.raster.atlas(grid.2d, grid.parameters, df.colors.tsne, show.edges = T, increase.factor = 1, file.path = NULL, add = TRUE, right.hemisphere = T)
  plot.raster.atlas(grid.2d, grid.parameters, df.colors.umap, show.edges = T, increase.factor = 1, file.path = NULL, add = TRUE, right.hemisphere = F)
 
  dev.off()
}
