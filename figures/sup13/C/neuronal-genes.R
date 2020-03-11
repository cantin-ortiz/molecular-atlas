#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/sup13/C')

#------------------- Parameters  ------------------- 

path.genes <- paste(path.matrices, 'genes-list-266.tsv', sep = '/')

#------------------- Loadings  ------------------- 

t <- read.table(paste(path.matrices, 'enrichR/neuronal-genes.txt', sep = '/'), sep = '\t')
colnames(t) <- c('ID', 'name')
t <- t[order(t$name),]

genes <- sort(as.character(read.table(path.genes)[,1]))

#------------------- Written stats ------------------- 

is.in.list <- sapply(genes, function(x){return(is.element(x, t$name))})

genes.neuron <- genes[is.in.list]
genes.not.neuron <- genes[!is.in.list]

paste(genes.not.neuron, collapse = ', ')
paste(genes.neuron, collapse = ', ')

#------------------- Plotting table ------------------- 

nrows <- 40
cex <- 0.5
ncols.neur <- ceiling(length(genes.neuron)/nrows)
ncols.not.neur <- ceiling(length(genes.not.neuron)/nrows)
ncols <- ncols.neur + ncols.not.neur

pdf('table-genes-neuronal.pdf', useDingbats = F, paper = 'a4', width = 9, height = 7)

plot(1, type = 'n', bty = 'n', ylim = c(nrows,-2), xlim = c(0, ncols + 1), axes = F, xlab = '', ylab = '')

for(i in 1:length(genes.neuron)){
  x <- (i-1) %% ncols.neur + 1
  y <- floor((i-1) / ncols.neur)
  text(x,y,genes.neuron[i], cex = cex)
}
y.max <- y

for(i in 1:length(genes.not.neuron)){
  x <- (i - ncols.not.neur - 1) %% ncols.not.neur + ncols.neur + 1
  y <- floor((i-1) / ncols.not.neur)
  text(x,y,genes.not.neuron[i], cex = cex)
}

y.max <- max(y, y.max)
x.seg <- setdiff(1:(ncols-1), ncols.neur) + 0.5
y.seg <- 0.5:(y.max-0.5)

segments(x.seg, -0.5, x.seg, y.max + 0.5, lend = 'butt', lwd = 0.5)
segments(0.5, y.seg, ncols+0.5, y.seg,  lend = 'butt', lwd = 0.5)

segments(c(0, ncols.neur, ncols) + 0.5, y.max+0.5, c(0, ncols.neur, ncols) + 0.5, -2, lwd = 2,  lend = 'butt')
segments(0.5, c(-2, -0.5, nrows-0.5), ncols+0.5, c(-2, -0.5, nrows-0.5), lwd = 2,  lend = 'butt')
# segments(0.5, -2, ncols+0.5, -2,  lend = 'butt')

text(ncols.neur / 2 + 0.5, -1.25, 'Neuronal genes', cex = 0.75)
text(ncols.not.neur / 2 + 0.5 + ncols.neur, -1.25, 'Non neuronal genes', cex = 0.75)
graphics.off()
