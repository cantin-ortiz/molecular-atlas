#!/usr/bin/env Rscript
# Script created by Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
# Run it like this:
# ./st_sc_correlation.R st_counts.tsv sc_counts.tsv st_clusters.tsv sc_clusters.tsv
#
library(enrichR)
library(stringr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GOplot)
library(pheatmap)
library(ggplot2)

heatmap_pathways <- function(pathways, name, max=25) {
  all_pathways = pathways$Term[1:min(nrow(pathways),max)]
  all_genes = getPathwaysGenes(pathways, all_pathways)
  data = data.frame(matrix(ncol=length(all_pathways), nrow=length(all_genes)))
  rownames(data) = all_genes
  colnames(data) = all_pathways
  
  for(pathway in all_pathways) {
    genes = getPathwaysGenes(pathways, pathway)
    data[genes,pathway] = 1
  }
  
  ind = apply(data, 1, function(x) all(is.na(x)))
  data = data[!ind,] 
  data[is.na(data)] = 0

  pdf(file=paste0(name, ".pdf"), height=12)
  print(pheatmap(data, cluster_rows=FALSE, cellwidth=5, cellheight=5, show_rownames=TRUE, legend=FALSE,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize_row=4, fontsize_col=4, angle_col=45))
  dev.off()
}

barplot_pathways <- function(pathways, name, max=25) {
  all_genes = getPathwaysGenes(pathways, pathways$Term)
  pathways = pathways[1:min(nrow(pathways),max),]
  pathways$Count = str_count(pathways$Genes, ";") + 1
  pathways$Ratio = pathways$Count / length(all_genes)
  pathways$Term = factor(pathways$Term, levels=pathways$Term[order(pathways$Combined.Score)])
  p = ggplot(pathways[1:min(nrow(pathways),max),], aes(x=Term, y=Combined.Score, fill=Ratio), position=position_stack(reverse=TRUE)) + geom_bar(stat="identity") + theme_minimal() + coord_flip() + theme(legend.position="top")
  pdf(file=paste0(name, ".pdf"))
  print(p)
  dev.off()
}

getPathwaysGenes <- function(enrichment_results, pathways) {
  pathway_genes = c()
  for(pathway in pathways) {
    genes = c()
    if (length(strsplit(enrichment_results[pathway,]$Genes, ";")) > 0) {
      genes = c(genes, strsplit(enrichment_results[pathway,]$Genes, ";")[[1]])
    }
    pathway_genes = c(pathway_genes, unique(genes))
  }
  pathway_genes = na.omit(pathway_genes)
  pathway_genes = unique(sapply(pathway_genes, tolower))
  firstup <- function(x) {
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    return(x)
  }
  pathway_genes = sapply(pathway_genes, firstup)
  return(pathway_genes)
}

extract_go_id <- function(go_term) {
  s = unlist(strsplit(go_term, " "))
  go = gsub("\\(", "", s[length(s)])
  go = gsub("\\)", "", go)
  return(go)
}

extract_go_term <- function(go_term) {
  s = unlist(strsplit(go_term, " "))
  term = do.call(paste, c(as.list(s[1:length(s)-1]), sep = " "))
  return(term)
}

convert_genes <- function(genes) {
  return(gsub(";", ", ", genes))
}


# The EnrichR databases
dbs = c("GO_Biological_Process_2018", "KEGG_2019_Mouse", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")

args = commandArgs(trailingOnly=TRUE)

# Load the genes
genes = read.delim(args[1], sep="\t", header=T, row.names=1)
genes = rownames(genes)
substring = "enriched_genes"
min_fdr = 0.05
min_genes = 3

# Compute pathways
enriched = enrichr(genes, dbs)
go_biological = enriched[[1]]
go_biological = go_biological[which(go_biological$Adjusted.P.value < min_fdr & str_count(go_biological$Genes, ";") >= min_genes + 1), c(1,2,3,4,7,8,9)]
write.table(go_biological, file=paste(substring, "GO_biological", ".tsv", sep="_"), sep="\t", row.names=T, col.names=T)

kegg = enriched[[2]]
kegg = kegg[which(kegg$Adjusted.P.value < min_fdr & str_count(kegg$Genes, ";") >= min_genes - 1), c(1,2,3,4,7,8,9)]
write.table(kegg, file=paste(substring, "KEGG", ".tsv", sep="_"), sep="\t", row.names=T, col.names=T)

go_molecular = enriched[[3]]
go_molecular = go_molecular[which(go_molecular$Adjusted.P.value < min_fdr & str_count(go_molecular$Genes, ";") >= min_genes + 1), c(1,2,3,4,7,8,9)]
write.table(go_molecular, file=paste(substring, "GO_molecular", ".tsv", sep="_"), sep="\t", row.names=T, col.names=T)

go_cellular = enriched[[4]]
go_cellular = go_cellular[which(go_cellular$Adjusted.P.value < min_fdr & str_count(go_cellular$Genes, ";") >= min_genes + 1), c(1,2,3,4,7,8,9)]
write.table(go_cellular, file=paste(substring, "GO_cellular", ".tsv", sep="_"), sep="\t", row.names=T, col.names=T)

# Create figures for GO biological 
new_ids = unlist(lapply(go_biological$Term, extract_go_id))
new_terms = unlist(lapply(go_biological$Term, extract_go_term))
go_biological$Term = new_terms
go_biological$ID = new_ids
rownames(go_biological) = go_biological$Term

heatmap_pathways(go_biological, "go_biological_heatmap_top10", max=10)
barplot_pathways(go_biological, "go_biological_barplot_top10", max=10)
heatmap_pathways(go_biological, "go_biological_heatmap_top25", max=25)
barplot_pathways(go_biological, "go_biological_barplot_top25", max=25)
heatmap_pathways(go_biological, "go_biological_heatmap_top50", max=50)
barplot_pathways(go_biological, "go_biological_barplot_top50", max=50)

# Create figures for GO molecular 
new_ids = unlist(lapply(go_molecular$Term, extract_go_id))
new_terms = unlist(lapply(go_molecular$Term, extract_go_term))
go_molecular$Term = new_terms
go_molecular$ID = new_ids
rownames(go_molecular) = go_molecular$Term

heatmap_pathways(go_molecular, "go_molecular_heatmap_top10", max=10)
barplot_pathways(go_molecular, "go_molecular_barplot_top10", max=10)
heatmap_pathways(go_molecular, "go_molecular_heatmap_top25", max=25)
barplot_pathways(go_molecular, "go_molecular_barplot_top25", max=25)
heatmap_pathways(go_molecular, "go_molecular_heatmap_top50", max=50)
barplot_pathways(go_molecular, "go_molecular_barplot_top50", max=50)

# Create figures for GO celullar  
new_ids = unlist(lapply(go_cellular$Term, extract_go_id))
new_terms = unlist(lapply(go_cellular$Term, extract_go_term))
go_cellular$Term = new_terms
go_cellular$ID = new_ids
rownames(go_cellular) = go_cellular$Term

heatmap_pathways(go_cellular, "go_cellular_heatmap_top10", max=10)
barplot_pathways(go_cellular, "go_cellular_barplot_top10", max=10)
heatmap_pathways(go_cellular, "go_cellular_heatmap_top25", max=25)
barplot_pathways(go_cellular, "go_cellular_barplot_top25", max=25)
heatmap_pathways(go_cellular, "go_cellular_heatmap_top50", max=50)
barplot_pathways(go_cellular, "go_cellular_barplot_top50", max=50)

# Create figures for KEGG
kegg$ID = kegg$Term 
rownames(kegg) = kegg$Term

heatmap_pathways(kegg, "kegg_heatmap_top10", max=10)
barplot_pathways(kegg, "kegg_barplot_top10", max=10)
heatmap_pathways(kegg, "kegg_heatmap_top25", max=25)
barplot_pathways(kegg, "kegg_barplot_top25", max=25)
heatmap_pathways(kegg, "kegg_heatmap_top50", max=50)
barplot_pathways(kegg, "kegg_barplot_top50", max=50)

# Create clusters
entrez = mapIds(org.Mm.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
data = data.frame(Entrez=entrez, group="Palette1")
res1 = compareCluster(Entrez~group, data=data, fun='enrichKEGG', organism="mmu", pvalueCutoff=min_fdr)
res2 = compareCluster(Entrez~group, data=data, fun='enrichGO', OrgDb='org.Mm.eg.db', pvalueCutoff=min_fdr)
res3 = compareCluster(Entrez~group, data=data, fun='groupGO', OrgDb='org.Mm.eg.db')

pdf("cluster_enrichKEGG.pdf")
dotplot(res1, title="KEGG", font.size=10)
dev.off()
pdf("cluster_enrichGO.pdf",  width = 10)
dotplot(res2, title="GO", font.size=10)
dev.off()
pdf("cluster_groupGO.pdf",  width = 10)
dotplot(res3, title="GROUPS GO", font.size=10)
dev.off()
