#------------------- INCLUDES AND PATH ------------------- 

source('bin/includes.R')
setwd('03-clustering')

#------------------- INITIAL TSV READS ------------------- 

#Path to the normalized expression matrix
path.to.normalized.expression.matrix <- 'expr_normalized_table.tsv'

#Based on manual classification of IC vectors
ic.kept <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,26,28,29,30,32,33,34,36,37,38,39,42,43,44,47,50,54,61,67,69,70,71)

#Clustering parameter
resolution <- 29

#Parsing table
st.data <- read.table(path.to.normalized.expression.matrix, sep = '\t', stringsAsFactors = F, header = T, row.names = 1)

#------------------- CREATING SEURAT OBJECT -------------------

library(Seurat)
library(data.table)

st.data.transp <- transpose(st.data)
colnames(st.data.transp) <- rownames(st.data)
rownames(st.data.transp) <- colnames(st.data)

seur.obj <- CreateSeuratObject(st.data.transp , project = "Few genes no strict")

#------------------- ICA -------------------

seur.obj@scale.data <- seur.obj@raw.data
seur.obj <- RunICA(seur.obj, ic.genes = rownames(st.data.transp), ics.print = 1:80,
                   ics.compute = 80, genes.print = 30, ica.function = 'icajade')


seur.obj <- FindVariableGenes(seur.obj, selection.method = 'dispersion', top.genes= length(rownames(seur.obj@data))/2)


seur.obj <- RunICA(seur.obj, ic.genes = seur.obj@var.genes, ics.print = 1:80,
                   ics.compute = 80, genes.print = 30, ica.function = 'icajade', reduction.name = 'fiftypercents')

save(seur.obj, file = 'seurobj-few-genes-no-strict.RData')

#------------------- RUN THE CLUSTERING ---------------------

#Clustering
seur.obj@dr$ica <- seur.obj@dr$fiftypercents
seur.obj <- FindClusters(object = seur.obj, reduction.type = "ica", dims.use = ic.kept, resolution = resolution,
                         save.SNN = FALSE, n.start = 100, nn.eps = 0, print.output = FALSE, force.recal = TRUE)  

#Experting clusters
table.export <- data.frame(cluster = as.numeric(seur.obj@ident))
rownames(table.export) <- names(seur.obj@ident)
write.table(table.export, file=path, sep='\t', quote=FALSE, col.names = FALSE)
