# NOTE: Unlike other figure scripts, this one was not updated to run instantly.

#!/usr/bin/env Rscript
# Script created by Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
# Run it like this:
# ./st_sc_correlation.R st_counts.tsv sc_counts.tsv st_clusters.tsv sc_clusters.tsv
#
library(ggrepel)
library(data.table)
library(ggplot2)
library(repr)
library(ggpubr)

NORM = TRUE
LOG = FALSE

args = commandArgs(trailingOnly=TRUE)

# Load the counts
st_counts = read.delim(args[1], sep="\t", header=T, row.names=1)
sc_counts = read.delim(args[2], sep="\t", header=T, row.names=1)

# Load spot meta info
st_clusters = read.delim(args[3], sep="\t", header=F, row.names=1)
sc_clusters = read.delim(args[4], sep="\t", header=F, row.names=1)

unique_sc_cluster = unique(c(sc_clusters$V2, sc_clusters$V3, sc_clusters$V4, 
                           sc_clusters$V5, sc_clusters$V6, sc_clusters$V7))
unique_sc_cluster = unique_sc_cluster[1:length(unique_sc_cluster)-1]

st_clusters$V3 = rownames(st_clusters)
st_clusters_index = st_clusters[st_clusters$V2 %in% unique_sc_cluster,]$V3

st_counts_filtered = st_counts[st_clusters_index,]
st_counts_filtered = na.omit(st_counts_filtered)

genes_intersect = intersect(colnames(st_counts_filtered), colnames(sc_counts))

sc_counts = sc_counts[,genes_intersect]
st_counts_filtered = st_counts_filtered[,genes_intersect]

print(dim(st_counts_filtered))
print(dim(sc_counts))

st_gene_counts = data.frame(ST=colMeans(st_counts_filtered), row.names=colnames(st_counts_filtered))
sc_gene_counts = data.frame(SC=colMeans(sc_counts), row.names=colnames(sc_counts))


both = data.frame(merge(st_gene_counts, sc_gene_counts, by = "row.names"), row.names = 1)
both = both[rowSums(both) > 0,]
both[is.na(both)] = 0

if (NORM) {
  both = t(t(both)/colSums(both)*mean(colSums(both)))
}

if (LOG) {
  both = log2(both + 1)
}

both = as.data.frame(both)

pdf("st_sc_correlation1.pdf")
plot(both$ST, both$SC)
dev.off()

fm <- lm(SC ~ ST, both)
both$resid = resid(fm)

# Subset to exclude genes with close to 0 expression
# Higher threhsold => labels more highly expressed genes
both_subset <- subset(both, ST > 1 & SC > 1)

# Select most/least extreme values
# Change to any number of top + bottom genes to label
high <- order(both_subset$resid, decreasing = T)[1:10]
low <- order(both_subset$resid, decreasing = F)[1:10]
both_subset$label <- ""
both_subset$label[high] <- rownames(both_subset)[high]
both_subset$label[low] <- rownames(both_subset)[low]

both$label <- both_subset[rownames(both), ]$label

options(repr.plot.width = 10, repr.plot.height = 8)
cols <- rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 9))
pdf("st_sc_correlation2.pdf")
ggplot(both, aes(ST, SC, label = label)) +
  geom_point(size = 1) +
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = F, n = 200, alpha = 0.6) +
  scale_fill_gradientn(colours = c(NA, cols)) +
  guides(fill = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_text_repel(data = subset(both, label != "")) +
  labs(title = "Gene-gene scatter plot: SC vs ST",
       x = "log [ST]",
       y = "log [SC]") +
  stat_cor(method = "pearson", label.x = 0, label.y = 2)
dev.off()
