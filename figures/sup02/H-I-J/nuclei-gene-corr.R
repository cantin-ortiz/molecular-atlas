#--------------- INCLUDES ---------------  

source('bin/includes.R')
setwd('figures/sup02/H-I-J')

#--------------- PARAMETERS --------------- 

x0 <- -10
x1 <- 100

#--------------- LOADING --------------- 

load.expr.mat()
spots.table <- load.spots.table()

#--------------- Computation --------------- 

n.counts <- unique(spots.table$nuclei)

df <- data.frame(reads = apply(st.data, 1, sum),
                 genes = apply(st.data > 0, 1, sum))

df[rownames(spots.table),'nuclei'] <- spots.table$nuclei

cor.genes <- cor.test(df$nuclei, df$genes, method = 'pearson')
cor.reads <- cor.test(df$nuclei, df$reads, method = 'pearson')

lm.reads <- lm(reads ~ nuclei , data=df)
lm.genes <- lm(genes ~ nuclei, data=df)

#--------------- Plots --------------- 
#----------------*** Correlation nuclei genes --------------- 

pdf('correlation-nuclei-genes.pdf',
    useDingbats = F)

boxplot(formula = genes ~ nuclei, data = df, ylab = 'Genes/spot', xlab = 'Nuclei/spot', xaxt = 'n')
# legend('topright', sprintf("Pearson corr coeff: %.3f", cor.genes$estimate))
y0 <- lm.genes$coefficients['(Intercept)'] + x0 * lm.genes$coefficients['nuclei']
y1 <- lm.genes$coefficients['(Intercept)'] + x1 * lm.genes$coefficients['nuclei']
segments(x0,y0,x1,y1, col = 'red', lwd = 2)
axis(1, at = seq(from = 1, to = 21, by = 2), labels = 2*(0:10))
axis(1, at = seq(from = 1, to = 22, by = 1), labels = rep('',22))
dev.off()

#----------------*** Correlation nuclei reads --------------- 

pdf('correlation-nuclei-reads.pdf',
    useDingbats = F)

boxplot(formula = reads ~ nuclei, data = df, ylab = 'Reads/spot', xlab = 'Nuclei/spot', xaxt = 'n')
# legend('topright', sprintf("Pearson corr coeff: %.3f", cor.reads$estimate))
y0 <- lm.reads$coefficients['(Intercept)'] + x0 * lm.reads$coefficients['nuclei']
y1 <- lm.reads$coefficients['(Intercept)'] + x1 * lm.reads$coefficients['nuclei']
segments(x0,y0,x1,y1, col = 'red', lwd = 2)
axis(1, at = seq(from = 1, to = 21, by = 2), labels = 2*(0:10))
axis(1, at = seq(from = 1, to = 22, by = 1), labels = rep('',22))
dev.off()

#----------------*** Nuclei per spots --------------- 

pdf('nuclei-per-spots.pdf',
    useDingbats = F)

cnt <- plyr::count(spots.table$nuclei)

plot(1, type = 'n', xlim = c(-0.5,21.5), ylim = c(0,4000), ylab = 'Spots', xlab = 'Nuclei', xaxt = 'n')
rect(-0.45:20.55, 0, 0.45:21.45, cnt$freq, col = 'lightgray')
axis(1, at = seq(from = 0, to = 21, by = 2), labels = 2*(0:10))
axis(1, at = seq(from = 0, to = 21, by = 1), labels = rep('',22))

# plot(cnt$x, cnt$freq, t = 'o', xlab = 'Number of nuclei', ylab = 'Number of spots')

dev.off()

