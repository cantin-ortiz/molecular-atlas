#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/sup10/C-D')

#------------- PARAMETERS -----------------

fname.pred <-  paste(path.matrices, 'predicted_classes_nn.tsv', sep = '/')
fname.cells <- paste(path.matrices, 'cells-to-classify.tsv', sep = '/')
cutoff <- 0.3

col.lines <- 'black'
lwd.segments <- 0.5

#------------- LOADINGS -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.cl.name <- load.df.cl.id.name(cl.id.names.path)
meta <- sc.load.meta.file(fname.cells)
l <- sc.load.prediction.file(fname.pred, meta, df.cl.name)
meta.pred <- l$meta
rm(l)

v1.chr <- mapvalues(clusters.v1, 
                    from = df.cl.name$cluster,
                    to = df.cl.name$clusters.named,
                    warn_missing = F)

alm.chr <- mapvalues(clusters.alm, 
                    from = df.cl.name$cluster,
                    to = df.cl.name$clusters.named,
                    warn_missing = F)

#------------- Initial computation -----------------

meta.gluta <- subset(meta.pred, cell_class == 'Glutamatergic')
df.truth <- data.frame(sc.clusters =  setdiff(sort(unique(meta.gluta$cell_cluster)), 'CR Lhx5'),
                       stringsAsFactors = F)
df.truth$layer <- unlist(lapply(strsplit(df.truth$sc.clusters, ' '), function(x){return(x[[1]])}))
df.truth$region <- NA
df.truth[grep('ALM', df.truth$sc.clusters),'region'] <- 'ALM'
df.truth[grep('VISp', df.truth$sc.clusters),'region'] <- 'VISp'

df.truth <- df.truth[!is.na(df.truth$region),]

meta.gluta <- meta.gluta[is.element(meta.gluta$cell_cluster, df.truth$sc.clusters),]

layer.count <- get.layer.count(spots.table[,c('clusters.named','acronym')])
layer.count <- layer.count[c(alm.chr,v1.chr),]
layer.count[,3] <- layer.count[,2] + layer.count[,3]
layer.count[,2] <- NULL
layer.count[,ncol(layer.count)] <- NULL

colnames(layer.count) <- paste('L',colnames(layer.count),sep='')


df.st <- data.frame(st.cluster = rownames(layer.count),
                    region = NA,
                    stringsAsFactors = F)
df.st[is.element(df.st$st.cluster, v1.chr), 'region'] <- 'VISp'
df.st[is.element(df.st$st.cluster, alm.chr), 'region'] <- 'ALM'
df.st$cluster.id <- as.numeric(mapvalues(df.st$st.cluster,
                              to = df.cl.name$cluster,
                              from = df.cl.name$clusters.named,
                              warn_missing = F))
df.st$layer.main <- ''
df.st$layer.secondary <- ''
df.st$layer.tertiary <- ''

for(i in 1:nrow(df.st)){
  
  sort.layer <- sort(layer.count[df.st[i, 'st.cluster'],],
                     decreasing = T)
  
  if(sort.layer[1] > cutoff)
    df.st[i, 'layer.main'] <- names(sort.layer)[1]
  if(sort.layer[2] > cutoff)
    df.st[i, 'layer.secondary'] <- names(sort.layer)[2]  
  if(sort.layer[3] > cutoff)
    df.st[i, 'layer.tertiary'] <- names(sort.layer)[3]  
}

#------------- Computation -----------------

cnt <- 0

for(i in 1:nrow(df.truth)){
  sc.cl <- df.truth[i, 'sc.clusters']
  sub <- subset(meta.gluta, cell_cluster == sc.cl)

  valid.cluster <- df.st[(df.truth[i,'region'] == df.st$region) & 
                         (df.truth[i, 'layer'] == df.st$layer.main | df.truth[i, 'layer'] == df.st$layer.secondary | df.truth[i, 'layer'] == df.st$layer.tertiary),
                         'cluster.id']
  
  cnt <- cnt + sum(is.element(sub$predicted,valid.cluster))
}

acc <- cnt / nrow(meta.gluta)

#------------- Plotting accuracy -----------------

pdf(generate.appropriate.file.name('accuracy-summary-suppl.pdf'),
    useDingbats = F)
plot(1, type = 'n', xlim = c(-0.5,1), ylim = c(0,1), xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
rect(1,0.51,0,0.61, col = 'gray90', border = 'black')
rect(0,0.51,acc,0.61, col = 'black', border = NA)

rect(1,0.49,0,0.39, col = 'gray90', border = 'black')
rect(0,0.49,0.65,0.39, col = 'black', border = NA)

axis(2, pos = 0, at = c(0.56,0.44), lab = c(sprintf('Glutamatergic, layers (%.1f%%)', 100*acc),'Glutamatergic, no layers (65.0%)'), las = 2, lwd = 0)
dev.off()

#------------- Plotting spot view -----------------

layer.count.same.order <- get.layer.count(spots.table[,c('clusters.named','acronym')],F)
layer.count.same.order <- layer.count.same.order[c(alm.chr, v1.chr),]
layer.count.same.order[,3] <- layer.count.same.order[,3] + layer.count.same.order[,2]
layer.count.same.order[,2] <- NULL

layer.count.same.order <- layer.count.same.order[rev(c(6,7,2,5,1,4,3,8,10,14,15,13,11,9,12)),]

df.grid <- expand.grid(0.5:(nrow(layer.count.same.order)-0.5), 0.5:(ncol(layer.count.same.order)-0.5))
radius <-  0.45*sqrt(as.vector(as.matrix(layer.count.same.order)))/sqrt(max(layer.count.same.order))

col.labels <- paste('L',colnames(layer.count.same.order),sep='')
col.labels[7] <-  'Non cortical'

radius.scale <- 0.45*sqrt(c(25,100,200))/sqrt(max(layer.count.same.order))

sel <- which(radius != 0)

n.row <- length(unique(df.grid[,1]))
n.col <- length(unique(df.grid[,2]))

back.color <- df.st
rownames(back.color) <- back.color$st.cluster
back.color <- back.color[rownames(layer.count.same.order),]

df.approved <- data.frame(y = 1:nrow(layer.count.same.order),
                          x = sapply(back.color$layer.main, function(x){return(which(x == paste('L', colnames(layer.count.same.order), sep='')))})) 
df.approved <- rbind(df.approved,
                     data.frame(y = (1:nrow(layer.count.same.order))[back.color$layer.secondary != ''],
                                x = sapply(back.color$layer.secondary[back.color$layer.secondary != ''], function(x){return(which(x == paste('L', colnames(layer.count.same.order), sep='')))}))) 
df.approved <- rbind(df.approved,
                     data.frame(y = (1:nrow(layer.count.same.order))[back.color$layer.tertiary != ''],
                                x = sapply(back.color$layer.tertiary[back.color$layer.tertiary != ''], function(x){return(which(x == paste('L', colnames(layer.count.same.order), sep='')))}))) 

df.approved$sel <- (df.approved$x - 1) * 15 + df.approved$y

fname <- 'dotplot-isocortex-layer.pdf'

pdf(generate.appropriate.file.name(fname),
    useDingbats = F,
    paper = 'a4',
    width = 7,
    height = 11)

plot(1, type = 'n', xlim = c(-10, n.col+5), ylim = c(-4.5,n.row), asp = 1, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
rect(df.approved$x, df.approved$y, df.approved$x-1, df.approved$y-1, border = NA, col = 'lightgreen')
rect(-7, 0, 0, 8, col = 'red', border = NA)
rect(-7, 8, 0, 15, col = '#3A3292', border = NA)
# col = 'red'

segments(-7, 0:n.row, n.col, 0:n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0:n.col, -7, 0:n.col, n.row, lty = 1, col = col.lines, lwd = lwd.segments)

sel.grey <- setdiff(sel, df.approved$sel)
symbols(x = df.grid[sel.grey,2], y = df.grid[sel.grey,1], circles = radius[sel.grey], inches = F, add = T, lwd = 0.25, bg = 'darkgray')
symbols(x = df.grid[df.approved$sel,2], y = df.grid[df.approved$sel,1], circles = radius[df.approved$sel], inches = F, add = T, lwd = 0.25, bg = 'darkgreen')

symbols(rep(n.col + 2,3), c(3,4,5), circles = radius.scale, inches = F, add = T, lwd = 0.25, bg = 'darkgray')
text(n.col+2.5, c(3,4,5), c(25,100,200), adj = c(0,0.5), cex =  1)
text(n.col+2.5,10,'Spots',adj = c(0,0.5), cex= 1)
# points(df.approved$x, df.approved$y, pch = 16, col = 'red')

axis(2, pos = 0, at = 0.5:n.row, labels= rownames(layer.count.same.order), las = 2, lwd= 0, cex.axis = 1.5)
# axis(4, pos = n.col-1, at = 0.5:n.row, labels= paste('(',cnt[rownames(layer.count.same.order),'freq'],')',sep=''), las = 2, lwd= 0, cex.axis = 0.8)
axis(1, pos = 0, at = 0.5:(n.col-0.5), labels = col.labels, las = 2, lwd = 0, cex.axis = 1.5)


dev.off()

