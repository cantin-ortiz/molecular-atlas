#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig3/A-B/')

#--------------- PARAMETERS --------------- 

sort.by.depth <- FALSE
col.non.selected <- 'lightgray'
layers.list <- c("-1", "1", "2/3", "4", "5", "6a", "6b")
layers.color <- c(col.non.selected, '#F48D2D', 'deeppink1', '#B79AC8', 'slateblue3','cyan3', 'darkgreen')
x.left.prop <- 1
x.right.prop <- 4

x.spatial <- seq(from = 4.5,
                 to = 6.5,
                 by = 0.5)

x.3d.tsne <- c(7,7.5)

x.count <- 8

#--------------- LOADING --------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
spots.table <- append.tsne.to.spots.table(spots.table, tsne.2d.path)
# load(grid.svm.path)

load(paste(path.matrices, 'layerstable.RData', sep = '/'))
layers.table[layers.table$layer == '2', 'layer'] <- '2/3' 

#--------------- Computation --------------- 

cl.list <- unique(as.character(spots.table$clusters.named))
iso.list <- cl.list[grep('^(Isocortex).*', cl.list)]

spots.table$is.iso <- is.element(as.character(spots.table$clusters.named),iso.list)
df.3dtsne.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

missing.acr <- setdiff(unique(spots.table$acronym),layers.table$acronym)
layers.table <- rbind(layers.table, data.frame(acronym = missing.acr, 
                                               name = '',
                                               layer = -1))       
spots.table$layer <- mapvalues(spots.table$acronym,
                               from = layers.table$acronym,
                               to = layers.table$layer)

#----------------*** t-SNE with vivid colors  --------------- 

spots.table <- append.common.vivid.colors(spots.table)
spots.table[!spots.table$is.iso, 'color'] <- col.non.selected
.plot.tsne(spots.table, file.path = 'tsne-isocortex-vivid.pdf', cex = 0.3, no.box = T)

#----------------*** t-SNE with layers  --------------- 

pal <- brewer.pal(9,'YlOrRd')


n.pal <- length(unique(layers.table$layer))
pal <- pal[(length(pal)-n.pal+1):length(pal)]

spots.table$color <- col.non.selected

spots.table$color <- mapvalues(spots.table$layer,
                               from = layers.list,
                               to = layers.color)

sp <- spots.table[order(as.numeric(mapvalues(spots.table$layer, from = sort(unique(spots.table$layer)), to = 8-c(7,2,6,3,5,4,1)))),]

.plot.tsne(sp, file.path = 'tsne-isocortex-layers.pdf', cex = 0.3, no.box = T)

#----------------*** Legend  --------------- 

pdf('legend isocortex layers.pdf', useDingbats = F)
plot(1, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
legend('center',c('Not cortex', paste('L', c("1", "2/3", "4", "5", "6a", "6b"), sep ='')), 
       pch = 19, col = layers.color,
       cex = 0.75)
dev.off()