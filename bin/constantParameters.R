#Clustering file
cl.file <- paste(path.matrices ,'res-29.tsv', sep='/')

#File with correspondance between clusters ID and paths
cl.id.names.path <- paste(path.matrices ,'cluster-correspondance-id-name.csv', sep='/')

#Path to the 3d t-SNE
tsne.3d.path <-  paste(path.matrices ,'tsne_20_3joint_000.txt',sep = '/')

#Path to the 2d t-SNE
tsne.2d.path <-  paste(path.matrices ,'tsne_20_2joint_000.txt',sep = '/')

#Path to the seurat object
seurat.object.path <- paste(path.matrices , 'seurobj-few-genes-no-strict.RData',sep = '/')

#Ic kept as biological
ic.kept <- setdiff(1:80,
                   c(15,25,27,31,35,40,41,45,46,48,49,51:53,55:60,62:66,68,72:80))

grid.svm.path.mesh.outline <- paste(path.matrices , 'grid-3d-mesh-outline.RData', sep = '/')
grid.svm.path.palette.mesh.outline <- paste(path.matrices , 'grid-3d-mesh-outline-palette.RData' ,sep = '/')

#Allen annotation path 
allen.annot.path <- paste(path.matrices , 'allen-mesh', sep = '/')

#Contains a list with all meshes already computed
path.all.meshes <- paste(path.matrices, 'all-meshes-outline-cut.RData', sep = '/')

#Color for non-white background
gray.color <- 'gray93'

lwd.atlas.segments <- 0.75

#Ground truth
clusters.v1 <- c(2,44,49,51,76,82,121,168)
clusters.alm <- c(7,11,19,30,37,21,95)

#Normalized matrix
normalized.matrix.path <- paste(path.matrices, 'normalized_matrix.RData', sep = '/')
