rm(list = ls())
source('bin/includes.R')
setwd('figures/sup2/A-B-D-E')

load.expr.mat()

genes.expressed <- rowSums(st.data != 0)
reads.count <- rowSums(st.data)

spot.expressed <- colSums(st.data != 0)
reads.gene <- colSums(st.data)

title <- ''

#--------- Gene per spot ---------- 

pdf('genes-per-spot.pdf', width = 5, height = 5)

hist(genes.expressed,
     breaks = 50,
     xlab = 'Genes',
     ylab = 'Spots',
     main = title,
     col = 'lightgray')

segments(1000, -100000, 1000, +10000, lwd = 1, lty = 2, col = 'black')

dev.off()

#--------- Read per spot ---------- 

pdf('reads-per-spot.pdf', width = 5, height = 5)

hist(reads.count,
     breaks = 50,
     xlab = 'Reads',
     ylab = 'Spots',
     main = title,
     col = 'lightgray')

dev.off()

#--------- READ COUNT ---------- 

pdf('spots-per-gene.pdf', width = 5, height = 5)

hist(log10(spot.expressed),
     breaks = 50,
     xlab = 'Spots',
     ylab = 'Genes',
     main = title,
     xaxt = 'n',
     col = 'lightgray',
     xlim = c(0,log10(100000)))
axis(1, 0:5, label = c(1,10,100,1000,10000,100000), cex.axis = 0.75)
segments(log10(100), -100000, log10(100), +10000, lwd = 1, lty = 2, col = 'black')
dev.off()

#--------- READ COUNT ---------- 

pdf('reads-per-gene.pdf', width = 5, height = 5)

hist(log10(reads.gene) ,
     breaks = 50,
     xlab = 'Reads',
     ylab = 'Genes',
     main = title,
     xaxt = 'n',
     xlim = c(0,log10(1000000)),
     col = 'lightgray')
axis(1, 0:6, label = c(1,10,100,1000,10000,100000,1000000), cex.axis = 0.75)

dev.off()
