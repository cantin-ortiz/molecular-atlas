#------------- INCLUDES -----------------

library(alphahull)
source('bin/includes.R')
setwd('figures/fig4/D')

#------------- PARAMETERS -----------------

ground.truth.v1 <- clusters.v1
ground.truth.alm <- clusters.alm

n.points <- 100000
x.range <- c(0, 6)
y.range <- c(-7,0)

dir.data <- 'C:/Users/MatLab/Desktop/transcripBrainAtlas/figures/main-figures/figure4'
fname.pred <- paste(path.matrices, 'predicted_classes_nn.tsv', sep = '/')
fname.cells <- paste(path.matrices, 'cells-to-classify.tsv', sep = '/')

v1.colors <- df.colors.vivid[c("Isocortex-01","Isocortex-19","Isocortex-39","Isocortex-28","Isocortex-20","Isocortex-17","Isocortex-26"),'clusters.colors']
alm.colors <- df.colors.vivid[c("Isocortex-33","Isocortex-09","Isocortex-05","Isocortex-14","Isocortex-04","Isocortex-10","Isocortex-08"),'clusters.colors']

alm.xlim <- c(0,2.1)
v1.xlim <- c(1.7,3.8)

alm.ylim <- c(-3.25,-1)
v1.ylim <- c(-2.25,0)

#------------- LOADINGS -----------------

load(paste(path.matrices, 'svm-remap-v1.RData', sep = '/'))
load(paste(path.matrices, 'svm-remap-alm.RData', sep = '/'))

df.cl.name <- load.df.cl.id.name(cl.id.names.path)

meta <- sc.load.meta.file(fname.cells)
l <- sc.load.prediction.file(fname.pred, meta, df.cl.name)
meta.pred <- l$meta
rm(l)

t <- read.table(paste(path.matrices, 'ordered-cells-tassic-paper.txt', sep = '/'), sep='\n', row.names = NULL, stringsAsFactors = F)[,1]
ordered.labels <- sapply(t, function(x){if(substr(x, nchar(x), nchar(x)) == ' ')return(substr(x, 1, nchar(x)-1))else return(x)})
cluster.list <- ordered.labels[1:55]

#---------- Contour to color ---------------

all.points <- unique(data.frame(x = c(l.svm.alm$seg$x0, l.svm.alm$seg$x1), y = c(l.svm.alm$seg$y0, l.svm.alm$seg$y1)))
p.alm <- ahull(all.points$x, all.points$y, alpha = 20)

all.points <- unique(data.frame(x = c(l.svm.v1$seg$x0, l.svm.v1$seg$x1), y = c(l.svm.v1$seg$y0, l.svm.v1$seg$y1)))
p.v1 <- ahull(all.points$x, all.points$y, alpha = 5)

#------------- *** Randomized points for the single cell clusters  -----------------
#------------- ****** Random points in ALM -----------------

for (cluster.names in c("L2/3 IT VISp Adamts2","L4 IT VISp Rspo1","L5 IT ALM Gkn1 Pcdh19","L6b ALM Olfr111 Spon1")){
# for(cluster.names in cluster.list){
  
  l.svm <- l.svm.alm
  df.points <- data.frame(ML = runif(n.points, min = x.range[1], max = x.range[2]),
                          DV = runif(n.points, min = y.range[1], max = y.range[2]))
  
  for(i in 1:length(l.svm$outline)){
    df.points[as.logical(point.in.polygon(df.points$ML, df.points$DV, l.svm$outline[[i]]$x, l.svm$outline[[i]]$y)), 'cluster.id'] <- l.svm$id[i]
  }
  
  df.points <- df.points[!is.na(df.points$cluster.id),]
  df.points <- df.points[df.points$cluster.id != 0,]
  df.points.alm <- df.points
  rm(df.points)
  
  #------------- ****** Random points in V1 -----------------
  
  l.svm <- l.svm.v1
  df.points <- data.frame(ML = runif(n.points, min = x.range[1], max = x.range[2]),
                          DV = runif(n.points, min = y.range[1], max = y.range[2]))
  
  for(i in 1:length(l.svm$outline)){
    df.points[as.logical(point.in.polygon(df.points$ML, df.points$DV, l.svm$outline[[i]]$x, l.svm$outline[[i]]$y)), 'cluster.id'] <- l.svm$id[i]
  }
  
  df.points <- df.points[!is.na(df.points$cluster.id),]
  df.points <- df.points[df.points$cluster.id != 0,]
  df.points.v1 <- df.points
  rm(df.points)
  
  #------------- *** Predicting clusters -----------------
  
  cur.cluster <- subset(meta.pred, cell_cluster == cluster.names)
  
  cur.cluster$is.alm <- is.element(cur.cluster$predicted, ground.truth.alm)
  cur.cluster$is.v1 <- is.element(cur.cluster$predicted,  ground.truth.v1)

  cur.cluster.alm <- subset(cur.cluster, is.alm)
  df.points <- df.points.alm
  for(cl in unique(cur.cluster.alm$predicted)){
    n.add <- dim(cur.cluster.alm[cur.cluster.alm$predicted == cl,])[1]
    cur.cluster.alm[cur.cluster.alm$predicted == cl,'ML'] <- df.points[df.points$cluster.id == cl, 'ML'][1:n.add]
    cur.cluster.alm[cur.cluster.alm$predicted == cl,'DV'] <- df.points[df.points$cluster.id == cl, 'DV'][1:n.add]
  }
  
  cur.cluster.v1 <- subset(cur.cluster, is.v1)
  df.points <- df.points.v1
  for(cl in unique(cur.cluster.v1$predicted)){
    n.add <- dim(cur.cluster.v1[cur.cluster.v1$predicted == cl,])[1]
    cur.cluster.v1[cur.cluster.v1$predicted == cl,'ML'] <- df.points.v1[df.points.v1$cluster.id == cl, 'ML'][1:n.add]
    cur.cluster.v1[cur.cluster.v1$predicted == cl,'DV'] <- df.points.v1[df.points.v1$cluster.id == cl, 'DV'][1:n.add]
  }
  
  rm(df.points)
  
  #------------- PLOTTING -----------------
  
  fname <- paste(cluster.names, '.pdf', sep ='')
  fname <- gsub('/', '-', fname)
  pdf(generate.appropriate.file.name(fname), useDingbats = F,
      paper = 'a4', width = 8, height = 11)
  par(mfrow = c(3,1),
      mai = c(0.01, 0.01, 0.01, 0.01))
  
  #ALM
  l.svm <- l.svm.alm
  plot(1, type = 'n', xlim = alm.xlim, ylim = alm.ylim,
       asp = 1, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
  for(i in 1:length(l.svm$outline)){
    polygon(l.svm$outline[[i]]$x, l.svm$outline[[i]]$y, col = alm.colors[[i]], border = NA)
  }
  segments(l.svm$seg$x0, l.svm$seg$y0, l.svm$seg$x1, l.svm$seg$y1, lwd = 1)
  points(cur.cluster.alm$ML, cur.cluster.alm$DV, cex = 0.5, pch = 16)
  text(0.75,-1,cluster.names, cex = 1.5)
  plot(p.alm, wpoints = F, add = T, lwd = 6, col = '#3A3292')
  
  #V1
  l.svm <- l.svm.v1
  plot(1, type = 'n', xlim = v1.xlim, ylim = v1.ylim, 
       asp = 1, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
  for(i in 1:length(l.svm$outline)){
    polygon(l.svm$outline[[i]]$x, l.svm$outline[[i]]$y, col = v1.colors[i], border = NA)
  }
  segments(l.svm$seg$x0, l.svm$seg$y0, l.svm$seg$x1, l.svm$seg$y1, lwd = 1)
  points(cur.cluster.v1$ML, cur.cluster.v1$DV, cex = 0.5, pch = 16)
  plot(p.v1, wpoints = F, add = T, lwd = 6, col = 'red')
  
  #Proportion
  n.ALM <- sum(cur.cluster$is.alm)
  n.v1 <- sum(cur.cluster$is.v1)
  cut.1 <- n.ALM/dim(cur.cluster)[1]
  cut.2 <- (n.ALM+n.v1)/dim(cur.cluster)[1]
  
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
  
  rect.param <- list(x0 = 0.2,
                     y0 = 0.4,
                     width = 0.6,
                     height = 0.2)
  
  rect(rect.param$x0, rect.param$y0, rect.param$x0 + rect.param$width, rect.param$y0 + rect.param$height)
  rect(rect.param$x0, rect.param$y0, rect.param$x0 + cut.1*rect.param$width, rect.param$y0 + rect.param$height, col =  '#3A3292')
  rect(rect.param$x0 + cut.1*rect.param$width, rect.param$y0, rect.param$x0 + cut.2*rect.param$width, rect.param$y0 + rect.param$height, col = 'red')
  rect(rect.param$x0 + cut.2*rect.param$width, rect.param$y0, rect.param$x0 + rect.param$width, rect.param$y0 + rect.param$height, col = 'gray')
  
  text(mean(c(rect.param$x0, rect.param$x0 + cut.1*rect.param$width)),
       rect.param$y0 + 1.5*rect.param$height,
       sprintf('ALM = %.1f%%', cut.1*100), srt = 0)
  
  text(mean(c(rect.param$x0 + cut.1*rect.param$width,  rect.param$x0 + cut.2*rect.param$width)),
       rect.param$y0 + 1.5*rect.param$height,
       sprintf('V1 = %.1f%%', 100*n.v1/nrow(cur.cluster)), srt = 0)
  
  text(mean(c(rect.param$x0 + cut.2*rect.param$width,  rect.param$x0 + rect.param$width)),
       rect.param$y0 + 1.5*rect.param$height,
       sprintf('Other = %.1f%%', (nrow(cur.cluster) - n.ALM - n.v1)*100/nrow(cur.cluster)), srt = 0)
  
  dev.off()
}

