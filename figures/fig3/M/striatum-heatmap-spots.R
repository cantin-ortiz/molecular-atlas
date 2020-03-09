#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig3/M')

#------------- PARAMETERS -----------------

max.radius <- 373
acronym.depth <- c('STRd',
                   get.acronym.child('STRv'),
                   'LSX',
                   get.acronym.child('sAMY'))

acronym.depth <- setdiff(acronym.depth, c('FS', 'IA', 'BA', 'AAA'))

list.markers <- c('Gpr88', 'Cartpt', 'Ndrg4', 'Zic1', 'Nr2f2', 'Ahi1', NA)
zscore.gene.expression <- T
col.lines <- 'black'
lwd.segments <- 0.5
cex.axis <- 0.8
str.color <- c("darkorange","dodgerblue3","firebrick3" , "springgreen3", 'lightgray')
n.legends <- c(25,100,300)

#------------- LOADINGS -----------------

load(paste(path.matrices, 'vivid-colors.RData', sep='/'))
spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
df.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
load(normalized.matrix.path)

#------------- Initial computation -----------------

cl.striatum <- levels(spots.table$clusters.named)[startsWith(levels(spots.table$clusters.named), 'Striatum')]
sp.str <- subset(spots.table, is.element(clusters.named, cl.striatum))
sp.str <- add.all.acronyms(sp.str)

#------------- Execution -----------------

df.bool <- as.data.frame(matrix(FALSE, nrow = nrow(sp.str), ncol = length(acronym.depth)))

colnames(df.bool) <- acronym.depth
rownames(df.bool) <- rownames(sp.str)

for(acr in colnames(df.bool))
  df.bool[,acr] <- sapply(sp.str$all.acronyms, function(x){return(is.element(acr, x))})
  
count.cluster <- as.data.frame(matrix(0, nrow = length(cl.striatum), ncol = ncol(df.bool)))
colnames(count.cluster) <- colnames(df.bool)
rownames(count.cluster) <- cl.striatum

prop.cluster <- count.cluster

for(cl in cl.striatum){
  sub.sp <- subset(sp.str, clusters.named == cl)
  count.cluster[cl,] <- colSums(df.bool[rownames(sub.sp),])
}

prop.cluster <- t(apply(count.cluster, 1, function(x){return(x/sum(x))}))

all.max.val <- apply(prop.cluster, 1, which.max)
cnt <- plyr::count(all.max.val)
ordered.cols <- cnt$x[order(cnt$freq, decreasing = T)]
not.max.cols <- setdiff(1:ncol(prop.cluster), ordered.cols)
if(length(not.max.cols) > 1)
  not.max.cols <- not.max.cols[order(colSums(prop.cluster[,not.max.cols]), decreasing = T)]
ordered.cols <- c(ordered.cols, not.max.cols)
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

cnt <- plyr::count(spots.table$clusters.named)
rownames(cnt) <- as.character(cnt$x)

non.cortical <- cnt[rownames(count.cluster), 'freq'] - rowSums(count.cluster)

#------------- Pre-plotting -----------------

show.non.striatal <- T
  
if(show.non.striatal){
  fname <- 'dotplot-striatum-region-non-striatal.pdf'
}else{
  fname <- 'dotplot-striatm-region.pdf'
}

df.grid <- expand.grid(0.5:(nrow(prop.cluster)-0.5), 0.5:(ncol(prop.cluster)-0.5))
# radius <-  0.45*sqrt(as.vector(as.matrix(count.cluster)))/sqrt(max(count.cluster))
# radius.non.cortical <- 0.45*sqrt(non.cortical)/sqrt(max(count.cluster))
# radius.legend <- 0.45*sqrt(n.legends)/sqrt(max(count.cluster))

radius <-  0.45*sqrt(as.vector(as.matrix(count.cluster)))/sqrt(max.radius)
radius.non.cortical <- 0.45*sqrt(non.cortical)/sqrt(max.radius)
radius.legend <- 0.45*sqrt(n.legends)/sqrt(max.radius)

col.labels <- colnames(prop.cluster)

if(show.non.striatal){
  radius <- c(radius,radius.non.cortical)
  df.grid <-  rbind(df.grid, data.frame(Var1 = sort(unique(df.grid[,1])), Var2 = rep(max(df.grid[,2]+1)))) 
  col.labels <- c(col.labels, 'Other')
}
sel <- which(radius != 0)

n.row <- length(unique(df.grid[,1]))
n.col <- length(unique(df.grid[,2]))

df.markers <- matrix(NA, nrow = n.row, ncol = n.col)
rownames(df.markers) <- rownames(prop.cluster)
colnames(df.markers) <- col.labels

for(i in 1:length(list.markers)){
  mark <- list.markers[i]
  if(is.na(mark))
    next
  expr <- sapply(rownames(df.markers), function(x){return(mean(st.data[rownames(spots.table[spots.table$clusters.named == x,]),mark]))})
  breaks <- seq(from = min(expr)-0.01, to = max(expr)+0.01, length.out = 101)
  
  if(zscore.gene.expression){
    all.expr <- st.data[rownames(subset(spots.table, is.element(clusters.named, rownames(df.markers)))),mark]
    expr <- (expr - mean(all.expr))/sd(all.expr)
    breaks <- seq(from = -max(-min(expr), max(expr)),
                  to = max(-min(expr), max(expr)),
                  length.out = 101)
    breaks <- c(-1000,seq(from = -1, to = 1, length.out = 101),1000)
  }
  
  col <- colorRampPalette(rev(brewer.pal(9,'RdBu')))(length(breaks)-1)
  df.markers[,i] <- col[.bincode(expr, breaks)]
}

#------------- Plotting without number of spots -----------------

fname <- 'dotplot-striatum-region-non-striatal-V2.pdf'
pdf(generate.appropriate.file.name(fname),
    useDingbats = F,
    paper = 'a4',
    width = 7,
    height = 11)

plot(1, type = 'n', xlim = c(0, 7), ylim = c(-4.5,55), asp = 1, bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
rect(0,0:(n.row-1),-6.5,1:n.row, border = NA, col = mapvalues(rownames(prop.cluster),
                                                                from = rownames(df.colors.vivid),
                                                                to = df.colors.vivid$clusters.colors,
                                                                warn_missing = F))
rect(c(0,1,3,4,6), -6.75, c(1,3,4,6,7), 0, col = str.color, border = NA)
                                                                                                                                
for(i in 1:ncol(df.markers)){
  rect(i-1, 0:(n.row-1), i, 1:n.row, col = df.markers[,i], border = NA)
}
segments(-6.5, 0:n.row, n.col, 0:n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0:n.col, -6.75, 0:n.col, n.row, lty = 1, col = col.lines, lwd = lwd.segments)

segments(-6.5, 0, -6.5, n.row, lty = 1, col = col.lines, lwd = lwd.segments)
segments(0, -6.75, n.col, -6.75, lty = 1, col = col.lines, lwd = lwd.segments)

symbols(x = df.grid[sel,2], y = df.grid[sel,1], circles = radius[sel], inches = F, add = T, lwd = 0.25, bg = 'darkgray')


#Legend

symbols(x = c(8,8,8), y = c(9,8,7), circles =  radius.legend, inches = F, add = T, lwd = 0.25, bg = 'darkgray')
text(8.75,c(9,8,7),n.legends, adj = c(0,0.5), cex =  cex.axis)
text(8.75,10,'Spots',adj = c(0,0.5), cex= cex.axis)

print.col <- paste(list.markers,col.labels, sep = ' - ')
print.col[length(print.col)] <- 'Other'

axis(2, pos = 1, at = 0.5:n.row, labels= rownames(prop.cluster), las = 2, lwd= 0, cex.axis = cex.axis)
axis(1, pos = 1, at = 0.5:(n.col-0.5), labels = print.col, las = 2, lwd = 0, cex.axis = cex.axis)
dev.off()

