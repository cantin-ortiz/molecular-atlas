#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/sup12')

#------------- PARAMETERS -----------------

gene.path  <- paste(path.matrices, 'genes-list-266.tsv', sep = '/')

method = 'ward.D2'
breaks <- seq(from = -0.01, to = 1.01, by = 0.01)
col <- colorRampPalette(brewer.pal(9,'YlOrRd'))(length(breaks)-1)
length.colors <- 4

col.hipp <- '#ffa6f8'
col.block <- '#f3ffad'
x.offset <- 38.5

#------------- LOADINGS -----------------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
list.genes <- read.table(gene.path, stringsAsFactors = F)[,1]
load(paste(path.matrices, 'normalized_matrix.RData', sep = '/'))

load(seurat.object.path)
mat.ic.all <- get.ic.mat(seur.obj, 'fiftypercents')
mat.ic.all.2 <- mat.ic.all[,ic.kept]
mat.ic.cl.avg.all <- get.ic.cluster.average(spots.table, mat.ic.all.2)
mat.ic.cl.dist <- get.ic.cluster.dist(mat.ic.cl.avg.all)

#------------- Initial computation -----------------

df.mean <- as.data.frame(matrix(nrow = length(unique(spots.table$clusters.named)),
                                ncol = length(list.genes)))

rownames(df.mean) <- unique(spots.table$clusters.named)
colnames(df.mean) <- list.genes

st.data.zscore <- st.data[,colnames(df.mean)]
st.data.zscore <- apply(st.data.zscore, 2, function(x){return((x-mean(x))/sd(x))})

for(cl in rownames(df.mean)){
  r.sel <- rownames(subset(spots.table, clusters.named == cl))
  df.mean[cl,] <- colMeans(st.data.zscore[r.sel,])
}

#Getting clustering tree
cl.tree <- get.ic.hiearchical.tree(mat.ic.cl.dist, spots.table, method)
ordered.labels <- cl.tree$labels[cl.tree$order]

df.mean <- df.mean[ordered.labels,]

#Group by peak value
df.mean <- df.mean[,order(apply(df.mean, 2, which.max), decreasing = F)]
max.val <- unique(apply(df.mean, 2, which.max))

for(val in max.val){
  sel <- apply(df.mean, 2, which.max) == val
  reorder <- order(df.mean[val,sel], decreasing = F)
  df.mean[,which(sel)] <- df.mean[,colnames(df.mean)[sel][reorder]]
  colnames(df.mean)[which(sel)] <- colnames(df.mean)[which(sel)][reorder]
}

breaks <- c(-1000,seq(from = -1, to = 1, length.out = 101),1000)
col <- colorRampPalette(rev(brewer.pal(9,'RdBu')))(length(breaks)-1)

pdf(generate.appropriate.file.name('gene-clusters-heatmap.pdf'),
    paper = 'a4',
    width = 7.5,
    height = 11.7)
image(x = 1:nrow(df.mean),
      y = 1:ncol(df.mean),
      z = as.matrix(df.mean),
      asp = 1,
      breaks = breaks,
      col = col,
      xaxt = 'n',
      yaxt = 'n',
      xlab = '',
      ylab = '',
      bty = 'n')
axis(1, 1:nrow(df.mean), rownames(df.mean), las = 2, cex.axis = 0.25, lwd = 0, pos = 5.5)
axis(2, 1:ncol(df.mean), colnames(df.mean), las = 1, cex.axis = 0.25, lwd = 0, pos = 5.5)
dev.off()

#------------- Computation -----------------
# 
# df.count <- as.data.frame(matrix(0, nrow = length(unique(sp$cluster.all)), ncol = length(unique(sp$cluster.266))))
# colnames(df.count) <- sprintf('palette-%03d',sort(unique(as.numeric(as.character(sp$cluster.266)))))
# rownames(df.count) <- sprintf('original-%03d',sort(unique(as.numeric(as.character(sp$cluster.all)))))
# 
# for(cl in unique(as.numeric(as.character(sp$cluster.all)))){
#   
#   tmp <- subset(sp, cluster.all == cl)
#   tmp$cluster.266 <- as.numeric(as.character(tmp$cluster.266))
#   cnt <- plyr::count(tmp$cluster.266)
#   df.count[sprintf('original-%03d',cl),
#            sprintf('palette-%03d',cnt$x)] <- cnt$freq
# }
# 
# tmp <- unlist(strsplit(rownames(df.count), '-'))
# tmp <- as.numeric(tmp[seq(from = 2, to = length(tmp), by = 2)])
# df.cl.id <- unique(spots.table[,c('cluster','clusters.named')])
# rownames(df.count) <- paste('original-', 
#                             mapvalues(tmp, from = as.numeric(as.character(df.cl.id$cluster)), to = as.character(df.cl.id$clusters.named)),
#                             sep = '')
# 
# tmp <- unlist(strsplit(colnames(df.count), '-'))
# tmp <- as.numeric(tmp[seq(from = 2, to = length(tmp), by = 2)])
# df.cl.id <- unique(spots.table.266[,c('cluster','clusters.named')])
# colnames(df.count) <- paste('palette-', 
#                             mapvalues(tmp, from = as.numeric(as.character(df.cl.id$cluster)), to = as.character(df.cl.id$clusters.named)),
#                             sep = '')
# 
# norm.original <- apply(df.count, 1, function(x){return(x/sum(x))})
# norm.palette <- apply(df.count, 2, function(x){return(x/sum(x))})
# 
# 
# #------------- Preparing plot ----------------
# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# 
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# labels.palette <- unlist(strsplit(rownames(cur.norm), 'palette-'))[seq(from = 2, to = 2*nrow(cur.norm), by = 2)]
# labels.original <- unlist(strsplit(colnames(cur.norm), 'original-'))[seq(from = 2, to = 2*ncol(cur.norm), by = 2)]
# 
# cut.palette.1 <- which(rownames(cur.norm) == 'palette-Isocortex-31')
# cut.palette.2 <- which(rownames(cur.norm) == 'palette-Olfactory areas-07')
# 
# cut.original.1 <- which(colnames(cur.norm) == 'original-Olfactory areas-10')
# cut.original.2 <- which(colnames(cur.norm) == 'original-Isocortex-29')
# 
# clusters.original <- c(labels.original[cut.original.1:cut.original.2], labels.original[startsWith(labels.original, 'Hippocampal')]) 
# clusters.palette <- c(labels.palette[cut.palette.1:cut.palette.2], labels.palette[startsWith(labels.palette, 'Hippocampal')]) 
# 
# all.col <- unique(df.colors.vivid$clusters.colors)
# hash.col <- 1:length(all.col)-1000.5
# col.used <- c(NA, all.col, col)
# breaks.used <- c(-3000,-1500,hash.col,breaks[2:length(breaks)])
# 
# #Adding original colors
# row.append <- rep(-2000, ncol(cur.norm))
# names(row.append) <- colnames(cur.norm)
# row.append[paste('original-', clusters.original, sep = '')] <- mapvalues(df.colors.vivid[clusters.original, 'clusters.colors'],
#                                                                          from = all.col,
#                                                                          to = hash.col,
#                                                                          warn_missing = F)
# row.append <- as.numeric(row.append)   
# plot.norm <- rbind(t(replicate(length.colors, row.append)),cur.norm)
# 
# col.append <- rep(-2000, nrow(plot.norm))
# names(col.append) <- rownames(plot.norm)
# for(cl in clusters.palette){
#   cnt <- plyr::count(spots.table[rownames(spots.table.266[spots.table.266$clusters.named == cl,]),'clusters.named'])
#   col.append[paste('palette-',cl,sep='')] <- mapvalues(df.colors.vivid[as.character(cnt$x)[which.max(cnt$freq)],'clusters.colors'],
#                                                        from = all.col,
#                                                        to = hash.col,
#                                                        warn_missing = F)
# }
# 
# col.append <- as.numeric(col.append)
# plot.norm <- cbind(replicate(length.colors, col.append),plot.norm)
# 
# labels.palette <- c(rep('',nrow(plot.norm) - length(labels.palette)),labels.palette)
# labels.original <-c(rep('',ncol(plot.norm) - length(labels.original)),labels.original)
# 
# x.rect.palette <- which(is.element(labels.palette, clusters.palette)) + x.offset - 0.5
# cols.rect.palette <- ifelse(startsWith(labels.palette[x.rect.palette+0.5-x.offset], 'Hippocampal'),
#                             col.hipp,
#                             col.block)
# 
# y.rect.original <- which(is.element(labels.original, clusters.original)) - 0.5
# cols.rect.original <- ifelse(startsWith(labels.original[y.rect.original+0.5], 'Hippocampal'),
#                             col.hipp,
#                             col.block)
#---------- Plot ---------- 
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram.pdf'),
#     useDingbats = F,
#     paper = 'a4',
#     width = 8.27,
#     height = 11.69)
# par(mar= c(0.2, 0.2, 0.2, 0.2),
#     xpd=TRUE)
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(plot.norm)+0.5+x.offset), ylim = c(-0.5, ncol(plot.norm) + 0.5), asp = 1)
# 
# rect(x.rect.palette,0.5, x.rect.palette+1, -25, border = NA, col = cols.rect.palette)
# rect(0.5+x.offset, y.rect.original, -25+x.offset, y.rect.original + 1, border = NA, col = cols.rect.original)
# 
# image(x = 1:nrow(plot.norm)+x.offset,
#       y = 1:ncol(plot.norm),
#       z = plot.norm, 
#       breaks = breaks.used, 
#       col = col.used,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(plot.norm)+x.offset, labels = labels.palette, las = 2, lwd = 0, cex.axis = 0.25, pos = 6)
# axis(2,at = 1:ncol(plot.norm), labels = labels.original, las = 2, lwd = 0, cex.axis = 0.25, pos = x.offset+6)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.5)
# 
# rect(cut.palette.1+x.offset-0.5+length.colors,
#      cut.original.1-0.5 + length.colors,
#      cut.palette.2+x.offset+0.5+length.colors, 
#      cut.original.2+0.5 + length.colors)
# 
# # mtext('Clusters from reduced gene palette', side = 3, line = -5)
# # mtext('Clusters from all genes', side = 4, line = -1)
# text(x.offset+nrow(plot.norm)+2, (ncol(plot.norm)+1)/2, 'Clusters from all genes', srt = -90, pos= 3)
# text(x.offset+(nrow(plot.norm)+1)/2, ncol(plot.norm), 'Clusters from reduced gene palette', pos= 3)
# dev.off()
# 

# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,sort(colnames(cur.norm), decreasing = T)]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# # for(val in max.val){
# #   sel <- apply(cur.norm, 1, which.max) == val
# #   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][order(cur.norm[sel, val], decreasing = T)],]
# # }
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-name.pdf'),
#     useDingbats = F)
# image(x = 1:nrow(cur.norm),
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '')
# 
# axis(1,at = 1:nrow(cur.norm), labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2)
# dev.off()

# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# #Coloring tree if required
# clus.col <- cutree(cl.tree.all, 13)
# names(clus.col) <- paste('original', names(clus.col), sep = '-')
# clus.col <- clus.col[ordered.labels.all]
# x.splits <- which(as.logical(diff(clus.col) != 0)) + 0.5
# 

# 
# x.offset <- 38.5
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram-lines.pdf'),
#     useDingbats = F)
# 
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(cur.norm)+0.5+x.offset), ylim = c(-0.5, ncol(cur.norm) + 0.5), asp = 1)
# 
# image(x = 1:nrow(cur.norm)+x.offset,
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(cur.norm)+x.offset, labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = 4)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = x.offset+6)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.5)
# 
# segments(x.offset+0.5, x.splits, x.offset+nrow(cur.norm)+0.5, x.splits, lwd = 0.5, col = 'black', lend='butt')
# segments(x.offset+0.5, x.splits, -50, x.splits, lwd = 0.5, lty = 2, col = 'black', lend='butt')
# 
# dev.off()
# 

# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# #Coloring tree if required
# clus.col <- cutree(cl.tree.all, 13)
# names(clus.col) <- paste('original', names(clus.col), sep = '-')
# clus.col <- clus.col[ordered.labels.all]
# x.splits <- which(as.logical(diff(clus.col) != 0)) + 0.5
# m.all.val <- apply(cur.norm, 1, which.max)
# y.splits <- sapply(x.splits, function(x){return(which(x > m.all.val)[1])}) - 0.5
# 

# 
# x.offset <- 38.5
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram-lines-vertical.pdf'),
#     useDingbats = F)
# 
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(cur.norm)+0.5+x.offset), ylim = c(-0.5, ncol(cur.norm) + 0.5), asp = 1)
# 
# image(x = 1:nrow(cur.norm)+x.offset,
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(cur.norm)+x.offset, labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = 8)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = x.offset+6)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.5)
# 
# segments(x.offset+0.5, x.splits, x.offset+nrow(cur.norm)+0.5, x.splits, lwd = 0.25, col = 'black', lend='butt')
# segments(x.offset+0.5, x.splits, -50, x.splits, lwd = 0.25, lty = 2, col = 'black', lend='butt')
# 
# segments(y.splits+x.offset,0.5, y.splits+x.offset, ncol(cur.norm)+0.5, lwd = 0.25, lty = 1, col = 'black', lend='butt')
# segments(y.splits+x.offset,0.5, y.splits+x.offset, -60, lwd = 0.25, lty = 2, col = 'black', lend='butt')
# 
# dev.off()
# 

# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# cur.norm <- cur.norm[ordered.labels.266,]
# 

# 
# x.offset <- 38.5
# y.offset <- -39
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram-both.pdf'),
#     useDingbats = F)
# 
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(cur.norm)+0.5+x.offset), ylim = c(-34, ncol(cur.norm) + 0.5), asp = 1)
# 
# image(x = 1:nrow(cur.norm)+x.offset,
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(cur.norm)+x.offset, labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = 8)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = x.offset+7)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.25,
#                     xx.offset = 3)
# plot.phylo.modified(as.phylo(cl.tree.266), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.25, direction = 'upwards',
#                     xx.offset = x.offset, yy.offset = y.offset)
# dev.off()
# 
# ##------------ LEGEND ------------ 
# 
# breaks.legend <- (breaks + 0.01) / 1.02
# 
# pdf('legend.pdf', useDingbats = F)
# plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
# rect(breaks.legend[1:(length(breaks.legend)-1)],0.4,breaks.legend[2:length(breaks.legend)],0.6, col = col, border = NA)
# rect(0,0.4,1,0.6)
# axis(1, pos = 0.3)
# dev.off()
# 
# 
