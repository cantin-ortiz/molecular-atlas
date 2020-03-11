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
