#!/usr/bin/env Rscript
# Script created by Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
# Run it like this:
# ./dea_analysis_edgeR_FT.R counts.tsv clusters.tsv spots_meta.tsv comparisons.tsv TRUE
#
#!/usr/bin/env Rscript
suppressMessages(library(edgeR))
suppressMessages(library(scran))
suppressMessages(library(BiocParallel))
register(MulticoreParam(4))

REDUCE = TRUE

args = commandArgs(trailingOnly=TRUE)

# Load the counts
counts = read.delim(args[1], sep="\t", header=T, row.names=1)

# Remove contaminated genes
contaminated_genes = c("Ttr", "Pmch", "Lars2", "Prkcd", "Enpp2", "Nrgn")
counts = counts[,!colnames(counts) %in% contaminated_genes]

# Load the meta clusters
meta = read.delim(args[2], sep="\t", header=T, row.names=1)
colnames(meta) = c("spot", "cluster")
rownames(meta) = meta$spot

# Load spot meta info
meta_info = read.delim(args[3], sep="\t", header=T, row.names=1)

# Add animal to meta file
meta$animal = meta_info[rownames(meta),]$animal

# Load the DE tests
de_tests = read.delim(args[4], header=FALSE)

# Number of clusters to reduce (random sampling)
REDUCE_NUMB = as.numeric(args[5])

# Parse the de tests
de_tests_list = list()
for (i in 1:nrow(de_tests)) {
  c1 = unlist(strsplit(as.character(de_tests$V1[i]), ","))
  c2 = unlist(strsplit(as.character(de_tests$V2[i]), ","))
  if (all(c1 %in% meta$cluster) && all(c2 %in% meta$cluster)) {
    de_tests_list[[i]] = list(c1,c2)
  }
}
unique_clusters = unique(unlist(de_tests_list))

# Keep only clusters in the file
meta = meta[meta$cluster %in% unique_clusters,]

if (REDUCE) {
  # Sample spots per cluster to reduce computational time
  to_keep = list()
  for(c in unique(meta$cluster)) {
    rows = rownames(meta[meta$cluster == c,])
    to_keep = c(to_keep, sample(rows, min(REDUCE_NUMB, length(rows))))
  }
  meta = meta[unlist(to_keep),]
}

# keep only spots in the meta file
counts = counts[rownames(meta),]

# create the Scran object
print("Creating Scran object")
sce = SingleCellExperiment(assays=list(counts=as.matrix(t(counts))))
# remove all genes with no/low expression, here set to expression in more than 100 spots with > 0 count
sce = sce[rowSums(counts(sce) > 0) > 100,]
# remove all spots with low number of genes, 1000 genes with > 0 count
sce = sce[,colSums(counts(sce) > 0) > 1000]

# apply Scran normalization 
print("Normalizing")
clusters = quickCluster(counts(sce), method="igraph", min.size=100)
sce = scran::computeSumFactors(sce, clusters=clusters, sizes=seq(50, 100, 10))
sce = scater::normalize(sce)
dge = scran::convertTo(sce, type="edgeR")

# Converting to edgeR seems to remove some genes/spots
meta = meta[rownames(dge$samples),]

# Iterate tests
for (i in 1:nrow(de_tests)) {
  c1 = unlist(de_tests_list[[i]][1])
  c2 = unlist(de_tests_list[[i]][2])
  print(paste(paste(unlist(c1),sep="-",collapse="-"), "_vs_", 
              paste(unlist(c2),sep="-",collapse="-"), sep=""))
  
  # Slice and create the design matrix
  meta_slice = meta[meta$cluster %in% unique(c(c1,c2)),]
  animal = meta_slice$animal
  cluster = as.character(meta_slice$cluster)
  cluster[cluster %in% c1] = "A"
  cluster[cluster %in% c2] = "B"

  # Slice EdgeR object
  print("Estimating dispersions and fitting")
  dge_slice = dge[,rownames(meta_slice)]
  dge_slice$samples$group = as.factor(cluster)
  dge_slice = estimateCommonDisp(dge_slice)
  dge_slice = estimateTagwiseDisp(dge_slice)
	  
  # Do test
  print("Testing")
  res = exactTest(dge_slice, pair=c("B", "A"))

  if (length(rownames(res$table)) > 0) {
    res = topTags(res, n=length(rownames(res$table)))
    res = res[!is.na(res$table$FDR),]
    res = res[order(res$table$FDR),]
    write.table(res$table, 
                file=paste(i, ".tsv", sep="_"), 
                sep='\t', quote=FALSE, col.names=NA)
  }
}
