#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/fig3/J-K')

#--------------- PARAMETERS --------------- 

col.non.selected <- 'lightgray'
str.color <- c("darkorange","dodgerblue3","firebrick3" , "springgreen3",col.non.selected)
main.region.striatum <- c('STRd', 'STRv', 'LSX', 'sAMY', 'not striatum')

#--------------- LOADING --------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
spots.table <- append.tsne.to.spots.table(spots.table, tsne.2d.path)
spots.table <- add.all.acronyms(spots.table)

acr <- sapply(spots.table$all.acronyms, function(x){return(unlist(sapply(x, function(y){return(which(y == main.region.striatum))})))})
spots.table$region.interest <- unlist(lapply(acr, function(x){return(ifelse(length(x)==0, 'not striatum', main.region.striatum[x]))}))

#--------------- Computation --------------- 

cl.list <- unique(as.character(spots.table$clusters.named))
str.list <- cl.list[grep('^(Striatum).*', cl.list)]

spots.table$is.str <- is.element(as.character(spots.table$clusters.named),str.list)
df.3dtsne.colors <- get.color.from.3d.tsne(tsne.3d.path, spots.table)

#--------------- Plots: t-SNE --------------- 
#----------------*** t-SNE with vivid colors  --------------- 

spots.table <- append.common.vivid.colors(spots.table)
spots.table[!spots.table$is.str, 'color'] <- col.non.selected
.plot.tsne(spots.table, file.path = 'tsne-striatum-vivid.pdf', cex = 0.3, no.box = T)

#----------------*** t-SNE with layers  --------------- 

spots.table$color <- col.non.selected
spots.table$color <- mapvalues(spots.table$region.interest, from = main.region.striatum, to = str.color)

.plot.tsne(spots.table, file.path = 'tsne-striatum-regions.pdf', cex = 0.3, no.box = T)

#----------------*** Legend  --------------- 

pdf('legend striatum.pdf', useDingbats = F)
plot(1, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
legend('center',main.region.striatum, 
       pch = 19, col = str.color,
       cex = 0.75)
dev.off()
