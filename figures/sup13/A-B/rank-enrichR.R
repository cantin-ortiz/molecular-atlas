source('bin/includes.R')
setwd('figures/sup13/A-B')

t.bio <- read.table(paste(path.matrices, 'enrichR/enriched_genes_GO_biological_.tsv', sep = '/'),
                    stringsAsFactors = F)
t.cel <- read.table(paste(path.matrices, 'enrichR/enriched_genes_GO_cellular_.tsv', sep = '/'),
                    stringsAsFactors = F)
t.mol <- read.table(paste(path.matrices, 'enrichR/enriched_genes_GO_molecular_.tsv', sep = '/'),
                    stringsAsFactors = F)
t.keg <- read.table(paste(path.matrices, 'enrichR/enriched_genes_KEGG_.tsv', sep = '/'),
                    stringsAsFactors = F) 

t.bio$group <- 'bio'
t.cel$group <- 'cel'
t.mol$group <- 'mol'
t.keg$group <- 'keg'

t.merged <- rbind(t.bio, t.cel, t.mol, t.keg)
t.merged <- subset(t.merged, group != 'keg')
t.merged <- t.merged[order(t.merged$Combined.Score, decreasing = T),]

selected <- t.merged[1:6,]
all.genes <- unique(unlist(strsplit(selected$Genes, ';')))


df.plot <- data.frame(matrix(nrow = length(all.genes),
                             ncol = 6))
colnames(df.plot) <- paste('path', 1:6, sep = '-')
rownames(df.plot) <- all.genes

for(r in rownames(df.plot)){
  for(i in 1:ncol(df.plot)){
    df.plot[r,i] <- is.element(r, unlist(strsplit(selected[i,'Genes'], ';')))
  }
}

df.plot <- as.matrix(df.plot)
breaks <- c(-0.5,0.5,1.5)
col <- c('#FDF9E0', '#800026')

pdf('top-combined-score.pdf', paper = 'a4r', width = 11.7, heigh = 7.4)
par("mar" = c(5.1,14.1,4.1,2.1))
plot(1, type = 'n', 
     xlim = c(0,80), 
     ylim = c(0,60),       
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     bty = 'n')

segments(c(0,20,40,60),-10,c(0,20,40,60),70, col = 'gray')

rect(0,
     seq(from = 0.5, to = 50.5, by = 10),
     rev(selected$Combined.Score),
     seq(from = 9.5, to = 59.5, by = 10),
     col = 'lightblue')



axis(1, at = c(0,20, 40, 60), las = 1, cex.axis = 0.7, lwd = 0)
axis(2, at =  seq(from = 5, to = 55, by = 10), gsub('\\(', '\n\\(', rev(paste(selected$Term, selected$Overlap))), las = 2, cex.axis = 0.7, lwd = 0)
dev.off()




df.plot.2 <- df.plot[,1:(ncol(df.plot)-1)]
df.plot.2 <- df.plot.2[rowSums(df.plot.2) != 0,]


pdf('gene-per-term.pdf', paper = 'a4r', width = 11.7, heigh = 7.4)
par("mar" = c(5.1,14.1,4.1,2.1))
plot(1, type = 'n', 
     xlim = c(0.5,(nrow(df.plot.2)+0.5)), 
     ylim = rev(c(0.5,(ncol(df.plot.2)+0.5))),       
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     bty = 'n')
image(x = 1:nrow(df.plot.2),
      y = 1:ncol(df.plot.2),
      as.matrix(df.plot.2),
      add = T,
      breaks = breaks,
      col = col)

axis(1, at = 1:nrow(df.plot.2), rownames(df.plot.2), las = 2, cex.axis = 0.7, col = NA, lwd = 0, pos = ncol(df.plot.2)+0.35)
axis(2, at = 1:ncol(df.plot.2), gsub('\\(', '\n\\(', selected$Term[1:ncol(df.plot.2)]), las = 2, cex.axis = 0.7, lwd= 0, pos= 0.65)
dev.off()



