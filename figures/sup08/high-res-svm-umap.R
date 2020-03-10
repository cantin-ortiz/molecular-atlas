#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/sup08')

#------------- PARAMETERS -----------------

data.path <- paste(path.matrices, 'smoothed-atlas/all-genes', sep = '/')
data.dir <- path.matrices
selected <- '10_0.3_correlation'
path.sel <- paste(data.dir,'/UMAP_3_',selected,'.txt',sep='')
vector.order <- c(1,2,3)

#------------- LOADING -----------------

spots.table.raw <- add.parent.acronym(load.spots.table())
spots.table.0 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 0)
spots.table.10 <- append.cluster.to.spots.table(spots.table.raw, cl.file, min.cluster.size = 10)
spots.table <- spots.table.10

#Reading the umap
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

#------------- PREDICTIONS -----------------

lf <- list.files(data.path, pattern = ('^(prediction).*(.RData)$'))
lf.ap <- as.numeric(gsub('mm', '.', unlist(strsplit(unlist(strsplit(lf, '_'))[seq(from = 2, to = 2*length(lf), by = 2)], '[.]'))[seq(from = 1, to = 2*length(lf), by = 2)]))

spots.ap <- sort(unique(spots.table$AP), decreasing = T)
lf.selected <- lf[sapply(spots.ap, function(x){return(which.min(abs(lf.ap-x)))})]

cnt <- 0

for(f in lf.selected){
  cnt <- cnt + 1
  load(paste(data.path, f, sep = '/'))
  pdf(generate.appropriate.file.name(sprintf('coronal%02d.pdf',cnt)),
      width = 10, height = 7, useDingbats = F)
  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  plot.2d.mixed(array(grid.2d, dim = c(nrow(grid.2d), ncol(grid.2d),1)),
                1, df.colors, grid.parameters, right.hemisphere = T, increase.factor = 2,
                show.box = FALSE, show.title = FALSE, pdf.name = NULL, show.saggital = F)
  dev.off()
}
