library(scran)

sc = read.delim("~/Projects/st/st_brain_atlas/single_cell_mapping/sc_counts.tsv", sep="\t", header=T, row.names=1)
counts = as.matrix(t(sc))
sce = SingleCellExperiment(assays=list(counts=counts))

# remove all genes with no/low expression, here set to expression in more than 100 spots with > 0 count
sce = sce[rowSums(counts(sce) > 0) > 100,]

# remove all spots with low number of genes, 1000 genes with > 0 count
sce = sce[,colSums(counts(sce) > 0) > 1000]

# Apply Scran normalization 
clusters = quickCluster(counts(sce), method="igraph", min.mean=0.1, min.size=100)
sce = scran::computeSumFactors(sce, clusters=clusters, sizes=seq(50, 100, 10))
sce = scater::normalize(sce)
log_counts = logcounts(sce)
write.table(t(log_counts), file="~/Projects/st/st_brain_atlas/single_cell_mapping/sc_counts_norm.tsv", sep="\t")