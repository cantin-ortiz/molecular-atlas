#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/fig4/E')

#------------------- Initial computation  ------------------- 

spots.table <- add.parent.acronym(load.spots.table())
spots.table <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)

df.colors <- unique(spots.table[,c('cluster','clusters.named')])
rownames(df.colors) <- df.colors$clusters.named
df.colors$clusters.named <- NULL
colnames(df.colors) <- 'cluster.id'
df.colors$clusters.colors <- sapply(rownames(df.colors), function(x){
  cnt <- count(spots.table[spots.table$clusters.named == x,'acronym.parent'])
  cnt$x <- as.character(cnt$x)
  return(color.from.acronym(cnt[which.max(cnt$freq),'x']))})

#------------------- Plotting  ------------------- 

bool.dendogram <- TRUE
i <- 2

fname <- 'detailed-accuracy.pdf'

fname.pred <- paste(path.matrices, 'predicted_classes_nn.tsv', sep = '/')
fname.cells <- paste(path.matrices, 'cells-to-classify.tsv', sep = '/')

rect.x.acc.left <- 0.18
rect.x.acc.right <- 0.34

rect.x.left <- 0.34
rect.x.right <- 0.5

text.x.count <- 0.55
text.x.cluster.lab <- 0.53

segment.x.right <- 0.75

ranges.group <- c(0.5,55.5,117.5,133.5)
ranges.titles <- c('Glutamatergic', 'GABAergic', 'Non neuronal')

#------------------- *** Loadings  ------------------- 

df.cl.name <- load.df.cl.id.name(cl.id.names.path)
meta <- sc.load.meta.file(fname.cells)
l <- sc.load.prediction.file(fname.pred, meta, df.cl.name)
meta.pred <- l$meta
rm(l)

meta.pred[meta.pred$dissected_region == 'VISp', 'correct.mapping'] <- is.element(meta.pred[meta.pred$dissected_region == 'VISp', 'predicted'], clusters.v1)
meta.pred[meta.pred$dissected_region == 'ALM', 'correct.mapping'] <- is.element(meta.pred[meta.pred$dissected_region == 'ALM', 'predicted'], clusters.alm)

#------------------- *** Computation  ------------------- 

spots.table$color <- mapvalues(as.numeric(spots.table$cluster),
                               from = df.colors$cluster.id,
                               to = df.colors$clusters.colors)

sc.cl.list <- sort(unique(meta.pred$cell_cluster))

list.prop <- NULL
prop.correct <- numeric(length(sc.cl.list))

#Looping through clusters
for(i in 1:length(sc.cl.list)){
  
  sc.cl <- sc.cl.list[i]
  
  #Counting the parent acronym for every cluster
  meta.subset <- subset(meta.pred, cell_cluster == sc.cl)
  cnt <- count(meta.subset$predicted.name)
  cnt$x <- as.character(cnt$x)
  cnt$prop <- cnt$freq / sum(cnt$freq)
  cnt <- cnt[order(cnt$freq, decreasing = T),]
  cnt$color <- df.colors[cnt$x,'clusters.colors']
  list.prop[[i]] <- cnt
  names(list.prop)[i] <- sc.cl
  prop.correct[i] <- sum(meta.subset$correct.mapping)/dim(meta.subset)[1] 
  names(prop.correct)[i] <- sc.cl
  
}

t <- read.table(paste(path.matrices, 'ordered-cells-tassic-paper.txt', sep = '/'), sep='\n', row.names = NULL, stringsAsFactors = F)[,1]
ordered.labels <- sapply(t, function(x){if(substr(x, nchar(x), nchar(x)) == ' ')return(substr(x, 1, nchar(x)-1))else return(x)})

list.prop <- list.prop[ordered.labels]
prop.correct <- prop.correct[ordered.labels]
#------------------- *** Plotting  ------------------- 

pdf.name <- generate.appropriate.file.name(fname)
pdf(pdf.name, paper = 'a4', width = 8.25, height = 11.5)
par(mar = c(1.1,4.1,1.1,2.1))
plot(1, type = 'n', xlim = c(0,1), ylim = c(length(list.prop)+1,0), bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')

rect(rect.x.acc.left,
     0.5,
     rect.x.acc.right,
     length(list.prop) + 0.5,
     border = NA,
     col = 'gray90')

x.lines <- seq(from = rect.x.acc.left, to = rect.x.acc.right, length.out = 6)
segments(x.lines[2:5], 0.5, x.lines[2:5], length(list.prop)+0.5, col = 'white')

#General plot
for(i in 1:length(list.prop)){
  
  to.plot <- list.prop[[i]]
  cl.name <- names(list.prop)[i]
  
  x.list.end <- sapply(1:dim(to.plot)[1],function(x){return(sum(to.plot$prop[1:x]))})*(rect.x.right-rect.x.left) + rect.x.left
  
  if(length(x.list.end) == 1)
    x.list.start <- rect.x.left
  else
    x.list.start <- c(rect.x.left,x.list.end[1:(length(x.list.end)-1)])
  
  #Proportion filling rectangle
  rect(x.list.start, i-0.5, x.list.end, i+0.5, col = to.plot$color, lwd = 0.25, border = NA)
  
  #Accuracy outline rectangle
  # rect(rect.x.acc.left, i-0.5, rect.x.acc.right, i+0.5)
  
  #Accuracy filling rectangle
  x.sep <- rect.x.acc.right - prop.correct[i]*(rect.x.acc.right-rect.x.acc.left)
  rect(c(rect.x.acc.right,x.sep),
       rep(i-0.5,2),
       c(x.sep,rect.x.acc.left),
       rep(i+0.5,2),
       col = c('black',	NA),
       lwd = 0.25,
       border = NA)

  #Count of cells
  text(text.x.count, i, sum(to.plot$freq), pos = 2, cex = 0.3)  
  
  #Text label
  text(text.x.cluster.lab, i, cl.name, pos = 4, cex = 0.3)
}

text(text.x.count, -1, 'n cells', cex = 0.4, pos = 2)
text((rect.x.left + rect.x.right) / 2, -1, 'ST clusters\n(% SC cluster)', cex = 0.4)
text((rect.x.acc.left + rect.x.acc.right) / 2, -1, 'Accuracy', cex = 0.4)
text(x.lines, 0, c(1,0.8,0.6,0.4,0.2,0), cex = 0.4)

if(bool.dendogram){
  segments(rect.x.acc.left, ranges.group, segment.x.right, ranges.group, lty = 2, col = 'darkgray')
  text(segment.x.right, (ranges.group[1:(length(ranges.group)-1)] + ranges.group[2:length(ranges.group)])/2, ranges.titles, srt = -90)
}
dev.off()


