#------------------- INCLUDES AND PATH ------------------- 
rm(list=ls())
source('C:/Users/MatLab/Desktop/transcripBrainAtlas/bin/includes.R')

#------------------- INITIAL TSV READS ------------------- 

st.data <- read.table('E:/20181203-finalMatrices/log_counts_corrected_few_genes_no_strict.tsv')

#------------------- CREATING SEURAT OBJECT -------------------

library(Seurat)
library(data.table)

st.data.transp <- transpose(st.data)
colnames(st.data.transp) <- rownames(st.data)
rownames(st.data.transp) <- colnames(st.data)

seur.obj <- CreateSeuratObject(st.data.transp , project = "Few genes no strict seuratV2")

#------------------- ICA -------------------

seur.obj@scale.data <- seur.obj@raw.data
seur.obj <- RunICA(seur.obj, ic.genes = rownames(st.data.transp), ics.print = 1:80,
                   ics.compute = 80, genes.print = 30, ica.function = 'icajade')


seur.obj <- FindVariableGenes(seur.obj, selection.method = 'dispersion', top.genes= length(rownames(seur.obj@data))/2)


seur.obj <- RunICA(seur.obj, ic.genes = seur.obj@var.genes, ics.print = 1:80,
                   ics.compute = 80, genes.print = 30, ica.function = 'icajade', reduction.name = 'fiftypercents')
save(seur.obj, file = 'seurobj-few-genes-no-strict.RData')

#------------------- CREATING 3D IC PLOTS ------------------- 

spots.table <- add.parent.acronym(load.spots.table())
lf <- names(seur.obj@dr)
lf <- sort(lf, decreasing = TRUE)

for(f in lf){

  seur.obj@dr[[names(seur.obj@dr)[2]]]@cell.embeddings

  dir.name <- f
  dir.name <- paste(dir.name, 'ica3dplots', sep = '-')
  new.dir <- paste(dir,dir.name,sep='/')
  dir.create(new.dir)

  ica.matrix <-  seur.obj@dr[[f]]@cell.embeddings

  for (i in 1:80){
    p <- plot.3d.ica(spots.table, ica.matrix, i)
    htmlwidgets::saveWidget(p, file = paste(new.dir,'/vect-', i, '.html', sep=''), selfcontained = FALSE, title = i)
  }
}

#------------------- RUN THE CLUSTERING ---------------------

seur.obj.2@dr$ica <- seur.obj@dr$fiftypercents
generate.cluster.range(spots.table, seur.obj.2, 29,
                       dir.path = 'C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/190109-seuratV2FewGenesNoStrict/fifty-percents',
                       n.ica = ic.kept.50.percent)

#------------------- PLOTTING ATLAS ---------------------
# 
# rm(list=ls())
# source('C:/Users/MatLab/Desktop/transcripBrainAtlas/bin/includes.R')
# 
# spots.table <- add.parent.acronym(load.spots.table())


# generate.smoothed.atlas('C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/190109-seuratV2FewGenesNoStrict/fifty-percents-new-wo-31/res-29.tsv',
#                         cost.list = 100,
#                         gamma.list = 1,
#                         resolution = 0.015)
# 
# plot.smooth.atlas.tsne.colors(file.path = 'atlas.pdf',
#                               grid.3d = 'E:/20181102-continuousClusters/20190204-080528-11/grid-3d.RData', 
#                               spots.table,
#                               'C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/190109-seuratV2FewGenesNoStrict/fifty-percents-new-wo-31/res-29.tsv',
#                               path.to.tsne = 'C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/181114-tsneFromJose/tsne_50_3joint_000.txt',
#                               colors.order = c(3,2,1))











#------------------- EXPORT CLUSTERS ---------------------
# rm(list=ls())
# source('C:/Users/MatLab/Desktop/transcripBrainAtlas/bin/includes.R')
# setwd('C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/181217-fewGenesNoStrict')
# 
# spots.table <- add.parent.acronym(load.spots.table())
# spots.table <- append.cluster.to.spots.table(spots.table,
#                                              'C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/181217-fewGenesNoStrict/fewgenes-no-strict/res-32.tsv',
#                                              min.cluster.size = 10)
# 
# table.export <- data.frame(cluster = spots.table[,'cluster'])
# rownames(table.export) <- rownames(spots.table)
# write.table(table.export, file='few-genes-no-strict-32.tsv', sep='\t', quote=FALSE, col.names = FALSE)
# 


plot.smooth.atlas.tsne.colors()
#------------------- CREATE SMOOTH OUTLINE ---------------------

# 
# rm(list=ls())
# gamma <- 0.2
# cost.list <- c(1,10,100,1000,10000)
# source('C:/Users/MatLab/Desktop/transcripBrainAtlas/bin/includes.R')
# setwd('C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/181217-fewGenesNoStrict')
# spots.table <- add.parent.acronym(load.spots.table())
# sp <- append.cluster.to.spots.table(spots.table, 'fewgenes-no-strict/res-32.tsv')
# source('svm-script.R')
