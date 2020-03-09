#-------------- Includes ----------

source('bin/includes.R')
setwd('figures/fig4/C')

#-------------- Files paths  and constant ---------

fname.plot.list <- 'probability-threshold.pdf'
fname.pred <- paste(path.matrices, 'predicted_classes_nn.tsv', sep = '/')
fname.cells <- paste(path.matrices, 'cells-to-classify.tsv', sep = '/')

#Returns accuracy depending on cutoff
get.acc <- function(cutoff){
  return(sum(corre.val > cutoff) / (sum(corre.val > cutoff ) + sum(wrong.val > cutoff)))
}
#Returns proporiton of discarded depending on cutoff
get.discard <- function(cutoff){
  return((sum(corre.val <= cutoff) + sum(wrong.val <= cutoff))/(length(corre.val) + length(wrong.val)))
}

#-------------- Loadings  ---------

spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file ,min.cluster.size = 10)
df.cl.name <- load.df.cl.id.name(cl.id.names.path)

meta <- sc.load.meta.file(fname.cells)
l <- sc.load.prediction.file(fname.pred, meta, df.cl.name)

pred.meta <- l$meta
pred.coeff <- l$prob
rm(l)

pred.prob <- as.data.frame(t(apply(pred.coeff, 1, function(x){return(exp(x)/sum(exp(x)))})))

#-------------- Computations  ---------

#Checking if classification is indeed argmax
max.clust <- colnames(pred.coeff)[apply(pred.coeff,1,which.max)]
identical(max.clust, pred.meta$predicted.name)

#Finding max coeff
max.val.coeff <- apply(pred.coeff,1,max)
max.val.prob <- apply(pred.prob,1,max)

pred.meta[names(max.val.coeff), 'pred.coeff'] <- max.val.coeff
pred.meta[names(max.val.prob), 'pred.prob'] <- max.val.prob

#Appending the coeff with coarse accuracy
pred.gluta <- subset(pred.meta, cell_class == 'Glutamatergic')

pred.gluta[pred.gluta$dissected_region == 'VISp', 'correct.mapping'] <- is.element(pred.gluta[pred.gluta$dissected_region == 'VISp', 'predicted'], clusters.v1)
pred.gluta[pred.gluta$dissected_region == 'ALM', 'correct.mapping'] <- is.element(pred.gluta[pred.gluta$dissected_region == 'ALM', 'predicted'], clusters.alm)

#-------- ****** Accuracy threshold, more point in the end  -------- 

pdf(generate.appropriate.file.name('probability-threshold.pdf'), useDingbats = F)

corre.val <- unlist(subset(pred.gluta, correct.mapping, pred.prob))
wrong.val <- unlist(subset(pred.gluta, !correct.mapping, pred.prob))

median.prob <- median(pred.gluta$pred.prob)

proba.list<- c(seq(from = 0, to = 0.99, by = 0.01),
               seq(from = 0.99, to = median.prob, by = 0.00001))
# proba.list <- proba.list[1:(length(proba.list)-1)]
y.disc <- sapply(proba.list, get.discard)
y.acc <- sapply(proba.list, get.acc)

par(mar = c(5.1,5.1,4.1,4.1))
plot(proba.list,y.acc, col = 'black', type = 'l', lty = 1, 
     xlab = 'Threshold value', ylab = 'Accuracy', yaxt = 'n',
     ylim = c(0,1), main = 'Accuracy function of prediction probability threshold',
     cex.axis = 1.5, cex.lab = 1.5)
lines(proba.list, y.disc, col = 'red')
lines(proba.list, y.acc, col = 'black', type = 'l', lty = 1)
axis(4, col = 'red', col.ticks = 'red', col.axis = 'red', cex.axis = 1.5, las = 2)
axis(2, col = 'black', las = 2, cex.axis = 1.5)
mtext('Proportion of discarded cells', 4, col = 'red', line = 3, cex = 1.5)


dev.off()

print(sprintf('Prob discarded: %.2f, Accuracy: %.4f, Threshold: %.4f', max(y.disc), y.acc[length(y.acc)], median.prob))

      