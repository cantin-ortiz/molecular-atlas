#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig3/D-E')


#------------- PARAMETERS -----------------

max.radius <- 373
acronym.depth <- get.acronym.child('Isocortex')
acronym.depth.detailed <- c(setdiff(acronym.depth, c('SS', 'MO')), 'SSp', 'SSs', 'MOs', 'MOp')
list.markers <- c('Slc1a2', 'Cux2', 'Rorb', 'Bcl6', 'Foxp2', 'Ctgf', NA)
zscore.gene.expression <- T
col.lines <- 'black'
lwd.segments <- 0.5
layers.color <- c('#F48D2D', 'deeppink1', '#B79AC8', 'slateblue3','cyan3', 'darkgreen', 'lightgray')

showed.3d <- c(c('Isocortex-02', 'Isocortex-03','Isocortex-06','Isocortex-15', 'Isocortex-18', 'Isocortex-21', 'Isocortex-25', 'Isocortex-34'),
              to.plot.chr <- c('Isocortex-11', 'Isocortex-23','Isocortex-27','Isocortex-37'),
              to.plot.chr <- c('Isocortex-04', 'Isocortex-05', 'Isocortex-10','Isocortex-14','Isocortex-29','Isocortex-30','Isocortex-33'))

#------------- LOADINGS -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
load(normalized.matrix.path)

#------------- Initial computation -----------------

cl.isocortex <- levels(spots.table$clusters.named)[startsWith(levels(spots.table$clusters.named), 'Isocortex')]
sp.iso <- subset(spots.table, is.element(clusters.named, cl.isocortex))
sp.iso <- add.all.acronyms(sp.iso)

#------------- Layers -----------------

layer.count <- get.layer.count(spots.table[,c('clusters.named', 'acronym')], FALSE)
layer.count$`2/3` <- layer.count$`2/3`+layer.count$`2`
layer.count$`2` <- NULL

#------------- Acronym fine ------------- 
#------------- *** Execution -----------------

df.bool <- as.data.frame(matrix(FALSE, nrow = nrow(sp.iso), ncol = length(acronym.depth.detailed)))

colnames(df.bool) <- acronym.depth.detailed
rownames(df.bool) <- rownames(sp.iso)

for(acr in colnames(df.bool))
  df.bool[,acr] <- sapply(sp.iso$all.acronyms, function(x){return(is.element(acr, x))})

count.cluster <- as.data.frame(matrix(0, nrow = length(cl.isocortex), ncol = ncol(df.bool)))
colnames(count.cluster) <- colnames(df.bool)
rownames(count.cluster) <- cl.isocortex

for(cl in cl.isocortex){
  sub.sp <- subset(sp.iso, clusters.named == cl)
  count.cluster[cl,] <- colSums(df.bool[rownames(sub.sp),])
}

cnt <- plyr::count(spots.table$clusters.named)
rownames(cnt) <- as.character(cnt$x)
non.cortical <- cnt[rownames(count.cluster), 'freq'] - rowSums(count.cluster)
count.cluster <- cbind(count.cluster,non.cortical)

prop.cluster <- count.cluster
prop.cluster <- t(apply(count.cluster, 1, function(x){return(x/sum(x))}))

ordered.cols <- c('MOp','MOs','SSp','SSs','GU','VISC','AUD','VIS','ACA','PL','ILA','ORB','AI','RSP','PTLp','TEa','PERI','ECT', 'non.cortical')


prop.cluster <- prop.cluster[,ordered.cols]

#Group by peak value
prop.cluster <- prop.cluster[order(apply(prop.cluster, 1, which.max)),]

max.val <- unique(apply(prop.cluster, 1, which.max))
for(val in max.val){
  sel <- apply(prop.cluster, 1, which.max) == val
  reorder <- order(prop.cluster[sel, val], decreasing = T)
  prop.cluster[which(sel),] <- prop.cluster[rownames(prop.cluster)[sel][reorder],]
  rownames(prop.cluster)[which(sel)] <- rownames(prop.cluster)[which(sel)][reorder]
}

prop.cluster <- prop.cluster[nrow(prop.cluster):1,]

breaks <- seq(from = -0.01, to = 1.01, by = 0.01)
col <- colorRampPalette(brewer.pal(9,'YlOrRd'))(length(breaks)-1)

count.cluster <- count.cluster[rownames(prop.cluster),]
count.cluster <- count.cluster[,colnames(prop.cluster)]

#------------- *** Plotting -----------------

fname <- 'dotplot-isocortex-region.pdf'

df.grid <- expand.grid(0.5:(nrow(prop.cluster)-0.5), 0.5:(ncol(prop.cluster)-0.5))
radius <-  0.45*sqrt(as.vector(as.matrix(count.cluster)))/sqrt(max.radius)

col.labels <- colnames(prop.cluster)
col.labels[length(col.labels)] <- 'Non cortical'

sel <- which(radius != 0)

n.row <- length(unique(df.grid[,1]))
n.col <- length(unique(df.grid[,2]))

rect.y <- which(is.element(rownames(count.cluster), showed.3d))
col.rect <- df.colors.vivid[rownames(count.cluster)[rect.y], 'clusters.colors']

cols.x.rect <- color.from.acronym(col.labels)
cols.x.rect[length(cols.x.rect)] <- 'lightgray'


pdf(generate.appropriate.file.name(fname),
    useDingbats = F,
    paper = 'a4',
    width = 7,
    height = 11)

plot(1, type = 'n', xlim = c(0, n.col), ylim = c(-4.5,n.row), asp = 1, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
rect(0:(n.col-1), 0, 1:n.col,-6.75, col = cols.x.rect, border = NA)
rect(0, rect.y, -6.5, rect.y-1, col = col.rect, border = NA)


segments(-6.5, 0:n.row, n.col, 0:n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0:n.col, -6.75, 0:n.col, n.row, lty = 1, col = col.lines, lwd = lwd.segments)

segments(-6.5, 0, -6.5, n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0, -6.75, n.col, -6.75, lty = 1, col = col.lines, lwd = lwd.segments)

# symbols(x = df.grid[sel,2], y = df.grid[sel,1], circles = radius[sel], inches = F, add = T, lwd = 0.25, bg = col[.bincode(as.vector(as.matrix(prop.cluster))[sel], breaks)])
symbols(x = df.grid[sel,2], y = df.grid[sel,1], circles = radius[sel], inches = F, add = T, lwd = 0.25, bg = 'darkgray')

axis(2, pos = 1, at = 0.5:n.row, labels= rownames(prop.cluster), las = 2, lwd= 0, cex.axis = 0.8)
axis(1, pos = 1, at = 0.5:(n.col-0.5), labels = col.labels, las = 2, lwd = 0, cex.axis = 0.8)
dev.off()
  
#-------------- Layers -------------- 
#-------------- *** Execution -------------- 

layer.count.same.order <- layer.count[rownames(prop.cluster),]
prop.layer.count <- t(apply(layer.count.same.order, 1, function(x){return(x/sum(x))}))

#Group by peak value, not considering not cortical
prop.layer.count <- prop.layer.count[order(apply(prop.layer.count[,1:ncol(prop.layer.count)], 1, which.max)),]
max.val <- unique(apply(prop.layer.count, 1, which.max))

for(val in max.val){
  sel <- apply(prop.layer.count[,1:ncol(prop.layer.count)], 1, which.max) == val
  reorder <- order(prop.layer.count[sel, val], decreasing = T)
  prop.layer.count[which(sel),] <- prop.layer.count[rownames(prop.layer.count)[sel][reorder],]
  rownames(prop.layer.count)[which(sel)] <- rownames(prop.layer.count)[which(sel)][reorder]
}

prop.layer.count <- prop.layer.count[nrow(prop.cluster):1,]
layer.count.same.order <- layer.count.same.order[rownames(prop.layer.count),]

#-------------- *** Plotting -------------- 

df.grid <- expand.grid(0.5:(nrow(layer.count.same.order)-0.5), 0.5:(ncol(layer.count.same.order)-0.5))
radius <-  0.45*sqrt(as.vector(as.matrix(layer.count.same.order)))/sqrt(max.radius)
radius.non.cortical <- 0.45*sqrt(non.cortical)/sqrt(max.radius)

col.labels <- colnames(layer.count.same.order)
col.labels[7] <-  'Non cortical'

sel <- which(radius != 0)

n.row <- length(unique(df.grid[,1]))
n.col <- length(unique(df.grid[,2]))

df.markers <- matrix(NA, nrow = nrow(layer.count.same.order), ncol = ncol(layer.count.same.order))
rownames(df.markers) <- rownames(layer.count.same.order)
colnames(df.markers) <- colnames(layer.count.same.order)

for(i in 1:length(list.markers)){
  mark <- list.markers[i]
  if(is.na(mark))
    next
  
  expr <- sapply(rownames(df.markers), function(x){return(mean(st.data[rownames(spots.table[spots.table$clusters.named == x,]),mark]))})
  breaks <- seq(from = min(expr)-0.01, to = max(expr)+0.01, length.out = 101)

  if(zscore.gene.expression){
    all.expr <- st.data[rownames(subset(spots.table, is.element(clusters.named, rownames(layer.count.same.order)))),mark]
    expr <- (expr - mean(all.expr))/sd(all.expr)
    breaks <- seq(from = -max(-min(expr), max(expr)),
                  to = max(-min(expr), max(expr)),
                  length.out = 101)
    breaks <- c(-1000,seq(from = -1, to = 1, length.out = 101),1000)
  }
  
  col <- colorRampPalette(rev(brewer.pal(9,'RdBu')))(length(breaks)-1)
  df.markers[,i] <- col[.bincode(expr, breaks)]
}


fname <- 'dotplot-isocortex-layer.pdf'

pdf(generate.appropriate.file.name(fname),
    useDingbats = F,
    paper = 'a4',
    width = 7,
    height = 11)

plot(1, type = 'n', xlim = c(0, n.col), ylim = c(-4.5,n.row), asp = 1, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
rect(0:6, 0, 1:7,-6.75, col = layers.color, border = NA)

for(i in 1:ncol(df.markers)){
  rect(i-1, 0:(n.row-1), i, 1:n.row, col = df.markers[,i], border = NA)
}
segments(-6.5, 0:n.row, n.col, 0:n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0:n.col, -6.75, 0:n.col, n.row, lty = 1, col = col.lines, lwd = lwd.segments)

segments(-6.5, 0, -6.5, n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0, -6.75, n.col, -6.75, lty = 1, col = col.lines, lwd = lwd.segments)

symbols(x = df.grid[sel,2], y = df.grid[sel,1], circles = radius[sel], inches = F, add = T, lwd = 0.25, bg = 'darkgray')

print.col <- paste(list.markers,col.labels, sep = ' - L')
print.col[length(print.col)] <- 'Non cortical'

axis(2, pos = 1, at = 0.5:n.row, labels= rownames(layer.count.same.order), las = 2, lwd= 0, cex.axis = 0.8)
axis(1, pos = 1, at = 0.5:(n.col-0.5), labels = print.col, las = 2, lwd = 0, cex.axis = 0.8)
dev.off()




