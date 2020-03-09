#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/fig4/C')

#------------------- Parameters  ------------------- 

fname.tassic.order <- paste(path.matrices, 'ordered-cells-tassic-paper.txt', sep = '/')
fname.pred <- paste(path.matrices, 'predicted_classes_nn.tsv', sep = '/')
fname.cells <- paste(path.matrices, 'cells-to-classify.tsv', sep = '/')
fname <- 'accuracy-summary.pdf'

gluta.range <- 1:55
gaba.range <- 56:117
nn.range <- 118:133
range.list <- list(gluta.range, gaba.range, nn.range)

color.correct <- 'black'
color.wrong <- 'gray90'

#------------------- Loadings ------------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
tassic.order <- read.table(fname.tassic.order, sep='\n', row.names = NULL, stringsAsFactors = F)[,1]
tassic.order <- sapply(tassic.order, function(x){if(substr(x, nchar(x), nchar(x)) == ' ')return(substr(x, 1, nchar(x)-1))else return(x)})

#------------------- Initial computation ------------------- 
    
df.cl.name <- load.df.cl.id.name(cl.id.names.path)
meta <- sc.load.meta.file(fname.cells)
l <- sc.load.prediction.file(fname.pred, meta, df.cl.name)
meta.pred <- l$meta
rm(l)

meta.pred[meta.pred$dissected_region == 'VISp', 'correct.mapping'] <- is.element(meta.pred[meta.pred$dissected_region == 'VISp', 'predicted'], clusters.v1)
meta.pred[meta.pred$dissected_region == 'ALM', 'correct.mapping'] <- is.element(meta.pred[meta.pred$dissected_region == 'ALM', 'predicted'], clusters.alm)


#------------------- Computation  ------------------- 

y.list <- list(c(0.7, 0.8), c(0.55,0.65), c(0.4,0.5))
prop <- list()
text.list <- c('Glutamatergic', 'GABAergic', 'Non-neuronal')
x1 <- 0.3
x2 <- 1

for(i in 1:3){
  sub.meta <- subset(meta.pred, is.element(cell_cluster, tassic.order[range.list[[i]]]))
  prop[[i]] <- sum(sub.meta$correct.mapping)/nrow(sub.meta)
}

pdf(generate.appropriate.file.name(fname), useDingbats = F)

plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab ='')

for(i in 1:3){
  rect(x1, y.list[[i]][1], (x2-x1)*prop[[i]] + x1, y.list[[i]][2], col = color.correct)
  rect(x2, y.list[[i]][1], (x2-x1)*prop[[i]] + x1, y.list[[i]][2], col = color.wrong)
  text(x1-0.2, mean(y.list[[i]]), sprintf('%s\n%0.1f%%', text.list[[i]], 100*prop[[i]]))
}

dev.off()



