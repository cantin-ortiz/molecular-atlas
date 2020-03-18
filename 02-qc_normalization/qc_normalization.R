#!/usr/bin/env Rscript
# Script created by Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
# Run it like this:
# ./qc_normalization.R counts.tsv spots_meta.tsv slices_meta.tsv FALSE 0
#
suppressMessages(library(BiocParallel))
register(MulticoreParam(8))
suppressMessages(library(scater))
suppressMessages(library(ggfortify))
suppressMessages(library(scran))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

qc_analysis <- function(sce, slice_meta, name, log_total_pca=TRUE) {
  
  pdf(name)
  
  sce = calculateQCMetrics(sce)
  
  # Plot highest expressed genes.
  print(plotHighestExprs(sce, colour_cells_by="animal"))
  
  # Plot frequency of expression (number of spots with detection) vs mean normalised expression.
  print(plotExprsFreqVsMean(sce))
  
  # Plot log10 total count vs number of spots.
  print(plotRowData(sce, x="n_cells_by_counts", y="log10_total_counts"))
  
  # Plot the percentage of expression accounted for by feature controls against total_features_by_X.
  print(plotColData(sce, x="total_features_by_counts", y="pct_counts_in_top_500_features", colour_by="acronym"))
  print(plotColData(sce, x="total_features_by_counts", y="pct_counts_in_top_500_features", colour_by="animal"))
  print(plotColData(sce, x="total_features_by_counts", y="pct_counts_in_top_500_features", colour_by="slice_index"))
  print(plotColData(sce, x="total_features_by_counts", y="pct_counts_in_top_500_features", colour_by="seq_dep"))
  
  # PCA - with different coloring
  print(plotPCA(sce, ncomponents=2, colour_by="animal") + geom_point(aes_string(fill = "colour_by"), size = 2, shape = 21, colour = "gray70", alpha = 0.25))
  print(plotPCA(sce, ncomponents=2, colour_by="seq_dep"))
  print(plotPCA(sce, ncomponents=2, colour_by="acronym"))
  print(plotPCA(sce, ncomponents=2, colour_by="slice_index"))
  print(plotPCA(sce, ncomponents=2, colour_by="total_features_by_counts"))
  print(plotPCA(sce, ncomponents=2, colour_by="total_counts"))
  
  # Shows how much of the data variation is explained by a single variable.
  print(plotExplanatoryPCs(sce, variables=c("total_features_by_counts", "total_counts", "animal", "seq_dep", "acronym", "slice_index")))
  print(plotExplanatoryPCs(sce, variable="total_features_by_counts", plot_type="pairs-pcs"))
  print(plotExplanatoryPCs(sce, variable="total_counts", plot_type="pairs-pcs"))
  print(plotExplanatoryPCs(sce, variable="animal", plot_type="pairs-pcs"))
  print(plotExplanatoryPCs(sce, variable="seq_dep", plot_type="pairs-pcs"))
  print(plotExplanatoryPCs(sce, variable="acronym", plot_type="pairs-pcs"))
  print(plotExplanatoryPCs(sce, variable="slice_index", plot_type="pairs-pcs"))
  
  # PCA total counts 
  log_counts = normcounts(sce)
  indexes = 1:length(rownames(slice_meta))
  sum_log_counts = matrix(nrow=length(rownames(slice_meta)), 
                          ncol=length(rownames(log_counts)))
  rownames(sum_log_counts) = rownames(slice_meta)
  colnames(sum_log_counts) = rownames(log_counts)
  for(slice in indexes) {
    sum_log_counts[slice,] = rowSums(log_counts[,grepl(slice, colnames(log_counts))])
  }
  sce_total = SingleCellExperiment(assays=list(counts=t(sum_log_counts)),
                                   colData=slice_meta[,c("animal", "AP", "seq_dep")])
  if (log_total_pca) {
    logcounts(sce_total) = log2(counts(sce_total) + 1)
  } else {
    logcounts(sce_total) = counts(sce_total)
  }
  
  sce_total = calculateQCMetrics(sce_total)
  print(plotPCA(sce_total, colour_by="AP", shape_by="animal", run_args=list(method="prcomp")))
  print(plotPCA(sce_total, colour_by="seq_dep", shape_by="animal", run_args=list(method="prcomp")))
  print(plotPCA(sce_total, colour_by="animal", run_args=list(method="prcomp")))
  print(plotPCA(sce_total, colour_by="seq_dep", run_args=list(method="prcomp")))
  print(plotPCA(sce_total, colour_by="AP", run_args=list(method="prcomp")))
  print(plotPCA(sce_total, colour_by="total_features_by_counts", run_args=list(method="prcomp")))
  print(plotPCA(sce_total, colour_by="total_counts", run_args=list(method="prcomp")))
  print(autoplot(prcomp(t(logcounts(sce_total))), data=slice_meta, colour="animal"))
  print(autoplot(prcomp(t(logcounts(sce_total))), data=slice_meta, colour="seq_dep"))
  print(autoplot(prcomp(t(logcounts(sce_total))), data=slice_meta, colour="AP"))
  
  dev.off()
}


MIN_GENES_SPOT = 1000
MIN_SPOTS_GENE = 100

contaminated_genes = c("Ttr", "Pmch", "Lars2", "Prkcd", "Enpp2", "Nrgn")
contaminated_genes_long = c("Ttr", "Pmch", "Lars2", "Prkcd", "Enpp2", 
                            "Ptgds", "Uba52", "Rpl9.ps6", "Supt7l",
                            "Gm10126", "Rpl37", "Gpx4", "Klf9", "Zbtb20", 
                            "Mt3", "Eno1", "Penk", "Pcsk1n", "Hspa8", "Rpsa", "Eef1a1",
                            "Arhgap5", "Celf4", "Hba.a2", "Rplp2", "Rps16", 
                            "Rps27rt", "Rps23", "Rpl28", "Tpt1", "Lrrc58", "Rps8",
                            "Gm10073", "Rpl35a", "Rps27", "Hbb.bt", "Elavl3", 
                            "Hba.a1", "Rgs7bp", "Calm1", "Sphkap", "Gm11808", "Sv2a",
                            "Cbx6", "Brsk1", "Pura")

args = commandArgs(trailingOnly=TRUE)

# Load the counts
counts = read.delim(args[1], sep="\t", header=T, row.names=1)

# Load spot meta info
meta = read.delim(args[2], sep="\t", header=T, row.names=1)

# Load the slice meta
slice_meta = read.delim(args[3], sep="\t", header=T, row.names=1)

LONG_GENE_LIST = as.logical(args[4])

MIN_COUNT = as.numeric(args[5])

# Remove contaminated genes
if (LONG_GENE_LIST) {
  counts = counts[,!colnames(counts) %in% contaminated_genes_long]
} else {
  counts = counts[,!colnames(counts) %in% contaminated_genes]
}

# Genes as rows
counts = t(counts)

# Perrm a small filtering
# Remove all genes with low number of spots detected
counts = counts[rowSums(counts > MIN_COUNT) > MIN_SPOTS_GENE,]
# Remove all spots with low number of genes detected
counts = counts[,colSums(counts > MIN_COUNT) > MIN_GENES_SPOT]

# Update the meta 
meta = meta[colnames(counts),]
meta = as.data.frame(meta)

# Create SCE object
sce = SingleCellExperiment(assays=list(counts=counts), colData=meta[,c("slice_index", 
                                                                       "acronym", 
                                                                       "animal", 
                                                                       "seq_dep")])

# RAW counts
normcounts(sce) = counts
logcounts(sce) = log2(counts + 1)
qc_analysis(sce, slice_meta, "qc_metrics.pdf", log_total_pca=TRUE)

# Apply Scran normalization
clusters = quickCluster(counts(sce), method="igraph", min.mean=0.1, min.size=100)
sce = scran::computeSumFactors(sce, clusters=clusters, sizes=seq(50, 100, 10))
sce = scater::normalize(sce)
write.table(t(logcounts(sce)), file="log_counts_normalized.tsv", sep="\t")
qc_analysis(sce, slice_meta, "qc_metrics_norm.pdf", log_total_pca=FALSE)

# Apply Scran Batch correction (normalize each batch separately)
counts = counts(sce)
a1 = counts[,colData(sce)$animal == "A1"]
sce_a1 = SingleCellExperiment(assays=list(counts=a1))
clusters_a1 = quickCluster(counts(sce_a1), method="igraph", min.mean=0.1, min.size=100)
sce_a1 = scran::computeSumFactors(sce_a1, clusters=clusters_a1, sizes=seq(50, 100, 10))
sce_a1 = scater::normalize(sce_a1)
a1_norm = logcounts(sce_a1)

a2 = counts[,colData(sce)$animal == "A2"]
sce_a2 = SingleCellExperiment(assays=list(counts=a2))
clusters_a2 = quickCluster(counts(sce_a2), method="igraph", min.mean=0.1, min.size=100)
sce_a2 = scran::computeSumFactors(sce_a2, clusters=clusters_a2, sizes=seq(50, 100, 10))
sce_a2 = scater::normalize(sce_a2)
a2_norm = logcounts(sce_a2)

a3 = counts[,colData(sce)$animal == "A3"]
sce_a3 = SingleCellExperiment(assays=list(counts=a3))
clusters_a3 = quickCluster(counts(sce_a3), method="igraph", min.mean=0.1, min.size=100)
sce_a3 = scran::computeSumFactors(sce_a3, clusters=clusters_a3, sizes=seq(50, 100, 10))
sce_a3 = scater::normalize(sce_a3)
a3_norm = logcounts(sce_a3)

batch_corrected = mnnCorrect(a1_norm, a2_norm, a3_norm, k=20, cos.norm.in=TRUE, cos.norm.out=FALSE)
a1_corrected = as.data.frame(t(batch_corrected$corrected[[1]]))
a2_corrected = as.data.frame(t(batch_corrected$corrected[[2]]))
a3_corrected = as.data.frame(t(batch_corrected$corrected[[3]]))
merged = join(a1_corrected, a2_corrected, by=NULL, type="full", match="all")
log_counts_corrected = join(merged, a3_corrected, by=NULL, type="full", match="all")
rownames(log_counts_corrected) = c(rownames(a1_corrected), 
                                   rownames(a2_corrected), 
                                   rownames(a3_corrected))
write.table(log_counts_corrected, file="log_counts_corrected.tsv", sep="\t")
counts = t(log_counts_corrected)

# Update the meta 
meta = meta[colnames(counts),]
meta = as.data.frame(meta)

sce = SingleCellExperiment(assays=list(counts=counts), colData=meta[,c("slice_index", 
                                                                       "acronym", 
                                                                       "animal",
                                                                       "seq_dep")])
normcounts(sce) = counts
logcounts(sce) = counts
qc_analysis(sce, slice_meta, "qc_metrics_norm_bc.pdf", log_total_pca=FALSE)

