#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/sup07/A-C')

#------------------- Parameters  ------------------- 

data.dir <- path.matrices
selected <- NULL
selected <- '10_0.3_correlation'

vector.order <- c(1,2,3)
no.box <- TRUE
all.files <- list.files(data.dir, pattern = '^(UMAP_2).*')

umap.path <- paste('UMAP_2_',selected, '.txt', sep = '')


#------------------- Loadings  ------------------- 

spots.table.raw <- add.parent.acronym(load.spots.table())
spots.table.0 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 0)
spots.table.10 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 10)

#------------------- Plotting 2D umap with ABA/vivid color code ------------------- 

for(path in umap.path){

  split <- strsplit(path, '[.]')[[1]]
  fname <- paste(split[1],split[2],sep='-')

  t <- read.table(paste(data.dir, path, sep = '/'), sep='\t', stringsAsFactors = F, header = TRUE, row.names = 1)
  spots.table.0$tsne1 <- t$X0
  spots.table.0$tsne2 <- t$X1

  spots.table.10[,c('tsne1','tsne2')] <- spots.table.0[rownames(spots.table.10),c('tsne1','tsne2')]

  spots.table.10 <- append.common.vivid.colors(spots.table.10)
  .plot.tsne(spots.table.10, file.path = paste(fname,'-vivid-color.pdf',sep=''), cex = 0.3, no.box = no.box)

  #ARA PLOT
  spots.table.10$color <- color.from.acronym(spots.table.10$acronym)
  .plot.tsne(spots.table.10, file.path = paste(fname,'-aba-color.pdf',sep=''), cex = 0.3, no.box = no.box)

}

#------------------- Plotting 2D umap with 3D colors codes------------------- 

path.sel <- paste(data.dir,'/UMAP_3_',selected,'.txt',sep='')
path.sel2D <- paste(data.dir,'/UMAP_2_',selected,'.txt',sep='')

t <- read.table(path.sel2D, sep='\t', stringsAsFactors = F, header = TRUE, row.names = 1)
spots.table.0$tsne1 <- t$X0
spots.table.0$tsne2 <- t$X1

spots.table.10[,c('tsne1','tsne2')] <- spots.table.0[rownames(spots.table.10),c('tsne1','tsne2')]

df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table.10)
df.colors.tsne <- df.colors
spots.table.10$color <- mapvalues(as.character(spots.table.10$clusters.named),
                                  from = rownames(df.colors),
                                  to = df.colors$clusters.colors)
.plot.tsne(spots.table.10, file.path = paste('UMAP_',gsub('[.]', '-', selected),'-3dtsne-color.pdf',sep=''), cex = 0.3, no.box = no.box)

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

#Filling it with median coordinates of the cluster
for (cl in rownames(df.colors)) {
  sel.spots <- is.element(rownames(t), rownames(spots.table.10[spots.table.10$clusters.named == cl, ]))
  df.colors[cl, 'clusters.colors'] <- rgb(median(vect.1[sel.spots]),
                                          median(vect.2[sel.spots]),
                                          median(vect.3[sel.spots]))
}

df.colors <- df.colors[order(rownames(df.colors)),]
df.colors.umap <- df.colors

spots.table.10$color <- mapvalues(as.character(spots.table.10$clusters.named),
                                  from = rownames(df.colors),
                                  to = df.colors$clusters.colors)

.plot.tsne(spots.table.10, file.path = paste('UMAP_',gsub('[.]', '-', selected),'-3dumap-color.pdf',sep=''), cex = 0.3, no.box = no.box)

#------------------- Plotting 2D t-SNE with 3D colors codes------------------- 

spots.table.tsne <- append.tsne.to.spots.table(spots.table.10, tsne.2d.path)
df.colors <- df.colors.tsne
spots.table.tsne$color <- mapvalues(as.character(spots.table.tsne$clusters.named),
                                    from = rownames(df.colors),
                                    to = df.colors$clusters.colors)
.plot.tsne(spots.table.tsne, file.path = paste('TSNE-3dtsne-color.pdf',sep=''), cex = 0.3, no.box = no.box)

df.colors <- df.colors.umap
spots.table.tsne$color <- mapvalues(as.character(spots.table.tsne$clusters.named),
                                    from = rownames(df.colors),
                                    to = df.colors$clusters.colors)
.plot.tsne(spots.table.tsne, file.path = paste('TSNE-3dumap-color.pdf',sep=''), cex = 0.3, no.box = no.box)

