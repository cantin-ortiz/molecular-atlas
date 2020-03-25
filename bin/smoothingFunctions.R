#------------- Generate smooth atlas ------------- 

#Generate a smoothen atlas using SVM
generate.smoothed.atlas <- function(cluster.path, #(str): path to a clustering file
                                    cost.list, #(numeric): list of cost value to use for the svm
                                    gamma.list, #(numeric): list of gamma value to use for the svm
                                    resolution = 0.05, #(numeric): size in mm of every voxel
                                    ncores = 2, #(int): number of cores to use for parallel computing
                                    main.dir = 'E:/20181102-continuousClusters', #(str): main directory in which atlas will be saved
                                    min.cluster.size = 10,  #(int): clusters with less spots than this will be discarded
                                    use.aligned.atlas.coordinates = TRUE, #Should the matrix be computed in the whole atlas, or only at the sections coordinates?
                                    model.to.use = NULL, #If null, building a new model. Otherwise loading another one
                                    AP.meshes = NULL,
                                    kernel.type = 'radial') #Only relevant if use.aligned.atlas.coordinates = TRUE. If not null, using meshes at the specificied AP coordinates to compute outline
  {

  library(e1071)
  library(pracma)  
  library(doParallel)
  
  #Loading and appending clusters
  spots.table <- add.parent.acronym(load.spots.table())
  spots.table <- append.cluster.to.spots.table(spots.table, cluster.path, min.cluster.size = min.cluster.size)
  sp <- spots.table
  
  #Getting unique AP coordinates and initializing volume
  AP.list.model <- sort(unique(spots.table$AP), decreasing = TRUE)
  volume.3d.model <- sp[sapply(sp$AP, function(x){return(is.element(x,AP.list.model))}),]
  volume.3d.model[,c('somaX','somaY')] <- atlas.spots$spots[rownames(volume.3d.model),c('somaX','somaY')]
  volume.3d.model <- volume.3d.model[c('somaY','somaX','ML','DV','AP','cluster')]
  
  #All coordinates used are the ones from the model
  if (!use.aligned.atlas.coordinates){
    AP.list <- AP.list.model
    volume.3d <- volume.3d.model
    
  #Case were using more AP coordinate
  }else{
    
    #Just atlas
    if(is.null(AP.meshes)){
    
      AP.list <- aligned.atlas$AP
      volume.3d <- NULL
      for(i in 1:length(aligned.atlas$soma))
        volume.3d <- rbind(volume.3d,aligned.atlas$soma[[i]])
      
    #From meshes
    }else{
      AP.list <- AP.meshes
      volume.3d <- NULL
    }
  }
  
  #List of all pairs of parameters values
  list.all.parameters.combination <- expand.grid(cost.list, gamma.list)
  colnames(list.all.parameters.combination) <- c('cost','gamma')
  
  #Directory name
  list.all.parameters.combination$dir.name <- paste(format(Sys.time(), "%Y%m%d-%H%M%S"),1:(dim(list.all.parameters.combination)[1]), sep='-')
  
  #Printing statuts
  print(paste('Starting parallel processing using ', ncores, ' cores ...', sep=''))
  print(paste(dim(list.all.parameters.combination)[1], ' loops to do.', sep=''))
  
  cl.par <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl.par)
  
  #For each loop
  l <- foreach(i=1:(dim(list.all.parameters.combination)[1]),
               .export = c('sub.function.smooth.atlas','aligned.atlas','ML.to.pixels', 'DV.to.pixels',
                           'get.ratio.stereo.to.pixels', 'point.in.polygon', 'predict', 'svm', 'atlas.spots',
                           'get.coronal.outline.from.allen.meshes', 'get.plane.equation',
                           'vector.cross.product.3d', 'get.id.from.acronym',  'mesh3d.allen.annot.from.id',
                           'get.meshes.list.plane.intersection', 'dist.from.plane', 'allen.annot.path', 'path.bin',
                           'get.triangle.plane.intersection', 'get.segment.plane.intersection', 'get.polygon.cut'),
               .packages = c('igraph', 'Rvcg', 'plyr', 'rjson', 'rgl')) %dopar% {

    sub.function.smooth.atlas(volume.3d, volume.3d.model,
                              list.all.parameters.combination[i,'cost'],
                              list.all.parameters.combination[i,'gamma'],
                              list.all.parameters.combination[i,'dir.name'],
                              spots.table, main.dir, cluster.path, resolution, model.to.use,
                              AP.list, AP.list.model, use.aligned.atlas.coordinates, AP.meshes,kernel.type)
  }
  
  # l <- lapply(1:(dim(list.all.parameters.combination)[1]),function(i){return(
  #             sub.function.smooth.atlas(volume.3d, volume.3d.model, 
  #                                       list.all.parameters.combination[i,'cost'],
  #                                       list.all.parameters.combination[i,'gamma'],
  #                                       list.all.parameters.combination[i,'dir.name'],
  #                                       spots.table, main.dir, cluster.path, resolution, model.to.use,
  #                                       AP.list, AP.list.model, use.aligned.atlas.coordinates, AP.meshes)
  # )})
  
  stopCluster(cl.par)
  return(l)
}

#Sub-function to generate a smooth atlas using SVM. Should always be called through "generate.smoothed.atlas"
sub.function.smooth.atlas <- function(volume.3d, #the volume of voxels to study
                                      volume.3d.model, #the volume with training data
                                      cost, #cost for svm
                                      gamma, #gamma for svm
                                      dir.name, #directory name
                                      sp,  #spot table
                                      main.dir, #path to main directory
                                      cluster.path,  #path to clustering file
                                      resolution = 0.05, #voxel size (mm)
                                      model.to.use = NULL,  #path to an already computed model
                                      AP.list,  #List of AP to compute the voxel
                                      AP.list.model, #List of AP used for training
                                      use.aligned.atlas.coordinates, #Boolean regarding the use of all atlas coordinates or not)
                                      AP.meshes = NULL,
                                      kernel.type = 'radial'){ #Only relevant if use.aligned.atlas.coordinates = TRUE. If not null, using meshes at the specificied AP coordinates to compute outline{
  
  #Base directory
  setwd(main.dir)

  #Parameters
  xlim.stereo <- c(-1,10)
  ylim.stereo <- c(-8,0)

  dir.create(dir.name)
  setwd(dir.name)
  
  #Managing the case with model having / not having to be built
  if(is.null(model.to.use)){
    print('Building model')
    model.to.use <- paste(getwd(),'/svm-model.RData',sep='')
    build.model <- TRUE
  }else{
    print(paste('Loading model: ', model.to.use, sep=''))
    build.model <- FALSE
    load(model.to.use)
  }
  
  #Generating a constant colorpalette
  colorpalette <- RColorBrewer::brewer.pal(12,'Set3')
  cp <- colorRampPalette(colorpalette)
  cp <- cp(length(unique(sp$cluster)))
  set.seed(42)
  cp <- sample(cp)
  
  #File with infos
  file.create('config.txt')
  fileConn<-file('config.txt')
  writeLines(c(dir.name,
               paste('gamma:', gamma),
               paste('cost:', cost),
               paste('AP.list:', paste(AP.list.model, collapse=';')),
               paste('Kernel:', kernel.type),
               paste('Model used:', model.to.use),
               paste('Cluster file:', cluster.path),
               paste('Resolution:', resolution),
               paste('Use atlas coordinates:', use.aligned.atlas.coordinates),
               paste('AP meshes:', paste(AP.meshes, collapse=';'))),
             fileConn)
  close(fileConn)
  
  #If required to build the model
  if(build.model){
    
    #Fitting model and saving it
    svmfit = svm(cluster ~ ML+DV+AP , data = volume.3d.model, kernel = kernel.type, cost = cost, gamma = gamma, scale = FALSE)
    save(svmfit, file = 'svm-model.RData')
  }
  
  #Data frame for remapping clusters from ID to name
  df.map.clusters <- unique(sp[,c('cluster','clusters.named')])
  df.map.clusters <- df.map.clusters[order(df.map.clusters$cluster),]
  df.map.clusters$cluster <- as.numeric(df.map.clusters$cluster)
  df.map.clusters$clusters.named <- as.character(df.map.clusters$clusters.named)
  rownames(df.map.clusters) <- df.map.clusters$cluster
  df.map.clusters$col <- cp[1:(dim(df.map.clusters)[1])]
  
  #Grid ticks
  grid.ML <- seq(from = xlim.stereo[1], to = xlim.stereo[2], by = resolution) #0.02
  grid.DV <- seq(from = ylim.stereo[1], to = ylim.stereo[2], by = resolution)
  
  #Creating the 3D grid
  grid.3d <- array(0, dim= c(length(grid.ML),length(grid.DV),length(AP.list)))
  storage.mode(grid.3d) <- "integer"
  
  #Plotting for every sections the clusters
  for (cur.AP in AP.list){
    
    #Displaying time to see progress
    print(Sys.time())

    #Outline from current section
    if(!use.aligned.atlas.coordinates){
      
      outlines.cur.vol.info <- atlas.spots$outlines$AP[atlas.spots$outlines$AP$AP == cur.AP,]
      outlines.cur.vol <- atlas.spots$outlines$outlines[[outlines.cur.vol.info$i]]
      
    }else{
      if(is.null(AP.meshes)){
        outlines.cur.vol.info <- which(aligned.atlas$AP == cur.AP)
        outlines.cur.vol <- aligned.atlas$outlines[[outlines.cur.vol.info]]
      }else{
        outlines.cur.vol <- get.coronal.outline.from.allen.meshes(cur.AP)
        
        if(length(outlines.cur.vol) == 0){
          next
        }
        for(i in 1:length(outlines.cur.vol)){
          outlines.cur.vol[[i]] <- as.data.frame(outlines.cur.vol[[i]])
          colnames(outlines.cur.vol[[i]]) <- c('x', 'y')
        }
      }
    }

    #Spots from current sections
    cur.vol <- volume.3d[volume.3d$AP == cur.AP,]
    
    #Creating a grid (change "by" for acccuracy)
    grid <- expand.grid(ML = grid.ML, 
                        DV = grid.DV)
    grid$AP <- cur.AP
    
    #Converting from a stereotactic gird to pixel grid
    grid.plot <- grid[,c('ML','DV','AP')]
    grid.plot <- subset(grid.plot, ML >= 0)
    
    if(!is.null(AP.meshes)){
      ML.to.pixels <- function(x,y){return(x)}
      DV.to.pixels <- function(x,y){return(x)}
    }
    
    grid.plot$somaX <- ML.to.pixels(grid.plot$ML, cur.vol)
    grid.plot$somaY <- DV.to.pixels(grid.plot$DV, cur.vol)
    xlim.pixels <- ML.to.pixels(xlim.stereo, cur.vol)
    ylim.pixels <- DV.to.pixels(ylim.stereo, cur.vol)
    
    #Check for every point if they are inside a contour
    matrix.inside <- matrix(nrow = dim(grid.plot)[1], ncol = length(outlines.cur.vol))
    
    if(is.null(AP.meshes)){
      for (i in 1:length(outlines.cur.vol)){
        matrix.inside[,i] <- point.in.polygon(round(grid.plot$somaX), round(grid.plot$somaY), 
                                              pol.x = round(as.numeric(outlines.cur.vol[[i]]$x)), 
                                              pol.y = round(as.numeric(outlines.cur.vol[[i]]$y)))
      }
    }else{
      for (i in 1:length(outlines.cur.vol)){
        matrix.inside[,i] <- point.in.polygon(grid.plot$somaX, grid.plot$somaY, 
                                              pol.x = as.numeric(outlines.cur.vol[[i]]$x), 
                                              pol.y = as.numeric(outlines.cur.vol[[i]]$y))
      }
    }
    
    l.inside <- apply(matrix.inside, 1, sum) > 0
    
    #Keeping only points inside contour
    grid.plot <- grid.plot[l.inside,]
    
    #Predicting cluster for each point
    p.plot <- predict(svmfit, grid.plot[,c('ML','DV','AP')])
    
    #Saving prediction in the grid
    idx <- data.frame(ML = sapply(grid.plot$ML, function(x){return(which(x == grid.ML))}),
                      DV = sapply(grid.plot$DV, function(x){return(which(x == grid.DV))}),
                      AP = which(cur.AP == AP.list))
    
    for (i in (1:(dim(idx)[1]))){
      grid.3d[idx[i,'ML'],idx[i,'DV'],idx[i,'AP']] <- p.plot[i]
    }
    
    #Plot name
    sign.AP <- '-'
    if (cur.AP > 0)
      sign.AP <- '+'
    plot.name <- sprintf('section_%02d_%s%.3f.pdf',which(AP.list == cur.AP), sign.AP, abs(cur.AP))
    
    #Creating plot
    pdf(paste('no-text_',plot.name,sep=''), fillOddEven = TRUE, useDingbats = F)
    plot(grid.plot[,c('somaX','somaY')], pch = 15, col = df.map.clusters[p.plot, 'col'], asp = 1, ylim = ylim.pixels, xlim = xlim.pixels, 
         main = paste('Bregma:', cur.AP), bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='', cex = 0.35)
    

    for (i in 1:length(outlines.cur.vol)){
      lines(as.numeric(outlines.cur.vol[[i]]$x), as.numeric(outlines.cur.vol[[i]]$y))
    }
    
    #Selecting clusters that are displayed
    cl <- unique(as.character(p.plot))
    
    #Df containing the legend
    df.legend <- data.frame(cluster = sort(as.numeric(cl)))
    df.legend$clusters.named <- df.map.clusters[df.legend$cluster, 'clusters.named']
    df.legend$col <- df.map.clusters[df.legend$cluster, 'col']
    df.legend$size <- sapply(df.legend$cluster, function(x){return(sum(p.plot == x))})
    #Displaying the legend
    points(rep(ML.to.pixels(7,cur.vol),length(cl)),seq(from = ylim.pixels[2], to = ylim.pixels[1], length.out = length(cl)),
           cex = 0.6, bg = df.legend$col, pch = 21)
    text(rep(ML.to.pixels(7,cur.vol),length(cl)),seq(from = ylim.pixels[2], to = ylim.pixels[1], length.out = length(cl)),
         paste(df.legend$cluster, '-', df.legend$clusters.named, ' (', df.legend$size, ')',sep=''), cex = 0.35, pos = 4)
    
    
    #Displaying the axis
    axis(3, pos = DV.to.pixels(0, cur.vol), at = ML.to.pixels(-1:6, cur.vol),labels = -1:6, tck  = 0.02, mgp=c(3, .5, 0))
    axis(2, pos = ML.to.pixels(-1, cur.vol), at = DV.to.pixels(-8:0, cur.vol),labels = -8:0, tck  = 0.02, mgp=c(3, .5, 0), las= 1)
    axis(1, pos = DV.to.pixels(-8, cur.vol), at = ML.to.pixels(-1:6, cur.vol),labels = -1:6, tck  = 0.02, mgp=c(3, .5, 0))
    axis(4, pos = ML.to.pixels(6, cur.vol), at = DV.to.pixels(-8:0, cur.vol),labels = -8:0, tck  = 0.02, mgp=c(3, .5, 0), las= 1)
    
    dev.off()
    
    #Creating plot
    pdf(plot.name, useDingbats = F)
    plot(1, type = 'n', asp = 1, ylim = ylim.pixels, xlim = xlim.pixels, 
         main = paste('Bregma:', cur.AP), bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='')
    
    points(grid.plot[,c('somaX','somaY')],col = df.map.clusters[p.plot, 'col'], pch = 15, cex = 0.35)
    
    for (i in 1:length(outlines.cur.vol)){
      lines(as.numeric(outlines.cur.vol[[i]]$x), as.numeric(outlines.cur.vol[[i]]$y))
    }
    
    #Selecting clusters that are displayed
    cl <- unique(as.character(p.plot))
    
    #Plotting the cluster ID in the median coordinates
    text(sapply(cl, function(x){return(median(grid.plot[p.plot == x, 'somaX']))}),
         sapply(cl, function(x){return(median(grid.plot[p.plot == x, 'somaY']))}),
         cl, font = 2)
    
    #Df containing the legend
    df.legend <- data.frame(cluster = sort(as.numeric(cl)))
    df.legend$clusters.named <- df.map.clusters[df.legend$cluster, 'clusters.named']
    df.legend$col <- df.map.clusters[df.legend$cluster, 'col']
    df.legend$size <- sapply(df.legend$cluster, function(x){return(sum(p.plot == x))})
    #Displaying the legend
    points(rep(ML.to.pixels(7,cur.vol),length(cl)),seq(from = ylim.pixels[2], to = ylim.pixels[1], length.out = length(cl)),
           cex = 0.6, bg = df.legend$col, pch = 21)
    text(rep(ML.to.pixels(7,cur.vol),length(cl)),seq(from = ylim.pixels[2], to = ylim.pixels[1], length.out = length(cl)),
         paste(df.legend$cluster, '-', df.legend$clusters.named, ' (', df.legend$size, ')',sep=''), cex = 0.35, pos = 4)
    
    
    #Displaying the axis
    axis(3, pos = DV.to.pixels(0, cur.vol), at = ML.to.pixels(-1:6, cur.vol),labels = -1:6, tck  = 0.02, mgp=c(3, .5, 0))
    axis(2, pos = ML.to.pixels(-1, cur.vol), at = DV.to.pixels(-8:0, cur.vol),labels = -8:0, tck  = 0.02, mgp=c(3, .5, 0), las= 1)
    axis(1, pos = DV.to.pixels(-8, cur.vol), at = ML.to.pixels(-1:6, cur.vol),labels = -1:6, tck  = 0.02, mgp=c(3, .5, 0))
    axis(4, pos = ML.to.pixels(6, cur.vol), at = DV.to.pixels(-8:0, cur.vol),labels = -8:0, tck  = 0.02, mgp=c(3, .5, 0), las= 1)
    
    dev.off()
  }
  
  grid.parameters <- list(ML = grid.ML, DV = grid.DV, AP = AP.list)
  save(grid.3d, grid.parameters, file = 'grid-3d.RData')
}

# ----------- Analysis of smoothing function ----------- 

get.accuracy.of.smoothing <- function(spots.table, path.to.grid.file = NULL){
  
  #Loading if needed
  if(!is.null(path.to.grid.file))
    load(path.to.grid.file)
  
  #Data frame containing the accuracy of smoothing for every cluster
  df.cluster.accuracy <- unique(spots.table[,c('cluster', 'clusters.named')])
  df.cluster.accuracy$cluster <- as.numeric(df.cluster.accuracy$cluster)
  df.cluster.accuracy$clusters.named <- as.character(df.cluster.accuracy$clusters.named)
  rownames(df.cluster.accuracy) <- df.cluster.accuracy$clusters.named
  df.cluster.accuracy <- df.cluster.accuracy[order(df.cluster.accuracy$clusters.named),]
  
  #Initializing relevent fields
  df.cluster.accuracy$matching.spots.and.shape <- 0
  df.cluster.accuracy$spots.in.the.wrong.shape <- 0
  df.cluster.accuracy$shape.arround.wrong.sp <- 0
  
  #Looping through all spots
  for (i in (1:dim(spots.table)[1])){
    
    #Selecting row
    s <- spots.table[i,]
    
    #Finding closest voxel
    idx.AP <- which.min(abs(s$AP - grid.parameters$AP))
    idx.ML <- which.min(abs(s$ML - grid.parameters$ML))
    idx.DV <- which.min(abs(s$DV - grid.parameters$DV))
    
    #Appending to appropriate list depending if inside/outside smoothed cluster
    if (grid.3d[idx.ML, idx.DV, idx.AP] != s$cluster)
      df.cluster.accuracy[df.cluster.accuracy$cluster == grid.3d[idx.ML, idx.DV, idx.AP] , 'shape.arround.wrong.sp'] <- df.cluster.accuracy[df.cluster.accuracy$cluster == grid.3d[idx.ML, idx.DV, idx.AP] , 'shape.arround.wrong.sp']  + 1
  }
  
  #Looking at every cluster
  for(cl.idx in 1:(dim(df.cluster.accuracy)[1])){
    
    #Corresponding cluster and spots
    cl <- df.cluster.accuracy[cl.idx, 'cluster']
    sp <- spots.table[spots.table$cluster == cl,]
    
    #Contains spots that are inside/outside the matching cluster after smoothing
    list.spots.in <- NULL
    list.spots.out <- NULL
    
    #For all spots from the cluster
    for (i in (1:dim(sp)[1])){
      
      #Selecting row
      s <- sp[i,]
      
      #Finding closest voxel
      idx.AP <- which.min(abs(s$AP - grid.parameters$AP))
      idx.ML <- which.min(abs(s$ML - grid.parameters$ML))
      idx.DV <- which.min(abs(s$DV - grid.parameters$DV))
      
      #Appending to appropriate list depending if inside/outside smoothed cluster
      if (grid.3d[idx.ML, idx.DV, idx.AP] == cl)
        list.spots.in <- c(list.spots.in, rownames(s))
      else
        list.spots.out <- c(list.spots.out, rownames(s))
      
    }
    
    #Saving the count
    df.cluster.accuracy[cl.idx,'matching.spots.and.shape'] <- length(list.spots.in) 
    df.cluster.accuracy[cl.idx,'spots.in.the.wrong.shape'] <- length(list.spots.out)
    
  }
  
  #Percentage of spots put in the proper shape
  df.cluster.accuracy$sensitivty <- df.cluster.accuracy$matching.spots.and.shape / (df.cluster.accuracy$matching.spots.and.shape + df.cluster.accuracy$spots.in.the.wrong.shape)
  
  #Percentage of spots in the shape that should not be there
  df.cluster.accuracy$false.classification <- df.cluster.accuracy$matching.spots.and.shape / (df.cluster.accuracy$matching.spots.and.shape + df.cluster.accuracy$shape.arround.wrong.sp)
  
  #Replacing NaN values
  df.cluster.accuracy[is.nan(df.cluster.accuracy$sensitivty), 'sensitivty'] <- 0
  df.cluster.accuracy[is.nan(df.cluster.accuracy$false.classification), 'false.classification'] <- 0
  
  #Plotting averaged information
  print(sprintf('Sensitivty: %.2f', mean(df.cluster.accuracy$sensitivty)))
  print(sprintf('False clas: %.2f', mean(df.cluster.accuracy$false.classification)))
  
  #Output
  return(df.cluster.accuracy)
}

# ----------- Plotting after generating ----------- 

#Generate the plot 
plot.smooth.atlas.tsne.colors <- function(file.path = 'atlas.pdf',
                                          grid.3d, 
                                          spots.table, 
                                          path.to.cluster,
                                          path.to.tsne = tsne.3d.path,
                                          colors.order = c(3,2,1),
                                          min.cluster.size = 10){
  
  #file.path(string): path or file name, with extension. If null, plots are not saved in a file.
  #grid.3d(string): path to a grid.3d object
  #spots.table
  #path.to.cluster(string): path to the clustering file
  #path.to.tsne(string): path to a 3D tsne object
  #colors.order(numerical): triplet that gives order of converison r,g,b/tsne1,2,3 for the get.colors.from.3d.tsne function
  #min.cluster.size(numerical): smaller clusters will be discarded
  
  #If saving the plot in a file
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path) 
  }
  
  #Loading the grid
  load(grid.3d)
  
  #Appending clusters
  spots.table <- append.cluster.to.spots.table(spots.table, path.to.cluster, min.cluster.size = min.cluster.size)
  
  #Getting proper colors
  df.colors <- get.color.from.3d.tsne(path.to.tsne, spots.table, use.hsv = FALSE, colors.order)
  
  #Plotting every section
  for (id.grid in 1:(dim(grid.3d)[3])){
    plot.2d.mixed(grid.3d, id.grid, df.colors, grid.parameters = grid.parameters)
  }
  
  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
}

#Plot the atlas as a raster 
plot.raster.atlas <- function(grid, 
                              grid.parameters, 
                              df.colors, 
                              title = '', 
                              hide.box = TRUE, 
                              show.edges = TRUE, 
                              increase.factor = 3, 
                              file.path = 'raster-atlas.pdf', 
                              xlim = c(-0.5,6), 
                              ylim = c(-8,1),
                              add = FALSE,
                              right.hemisphere = TRUE){
  
  library(raster)
  
  df.colors$cluster.id <- as.numeric(df.colors$cluster.id)
  df.colors <- df.colors[order(df.colors$cluster.id),]
  
  increased_matrix <- grid[rep(1:nrow(grid), each = increase.factor), 
                           rep(1:ncol(grid), each = increase.factor)]
  
  if(!right.hemisphere){
    grid.parameters$ML <- -grid.parameters$ML
    xlim <- c(-6,6)
    ylim <- c(-8,8)
  }

  r <- raster(t(increased_matrix),
              xmn = min(grid.parameters$ML),
              xmx = max(grid.parameters$ML),
              ymn = min(grid.parameters$DV),
              ymx = max(grid.parameters$DV))
  
  r <- raster::flip(r, 'y')
  
  if(!right.hemisphere)
    r <- raster::flip(r, 'x')
  
  #Detecting edges if required
  if(show.edges){
    
    #Detecting edges
    grid.edges <- matrix(FALSE, nrow = dim(r)[1], ncol = dim(r)[2])
    for (i in 1:(dim(grid.edges)[1]-1)){
      grid.edges[i,] <- (r[i+1,] - r[i,]) != 0 | grid.edges[i,]
    }
    for (i in 1:(dim(grid.edges)[2]-1)){
      grid.edges[,i] <- ((r[,i+1] - r[,i]) != 0) | grid.edges[,i]
    } 
    
    pix <- which(grid.edges, arr.ind = T)
    r[pix] <- -1
  }
  
  #Creating pdf if required
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path, useDingbats = F) 
  }
  
  plot(r,
       breaks = seq(from = -1.5, to = dim(df.colors)[1] + 0.5, by = 1),
       col = c('black', NA, df.colors$clusters.colors),
       maxpixels=5000000, interpolate = FALSE,
       legend = FALSE, main = title, xlim = xlim, ylim = ylim, 
       axes = FALSE, box = !hide.box, add = add)    

  
  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
  
}

#Read the config file
read.config.file <- function(config.path){
  
  list.param <- list()
  
  con <- file(config.path,open="r")
  
  idx <- 0
  while (TRUE){
    line <- readLines(con, n = 1)
    idx <- idx + 1
    if (length(line) == 0)
      break
    
    #Dir
    if (idx == 1){
      dir <- line
      list.param <-  c(list.param, dir = dir)
    }
    
    #Gamma
    if (idx == 2){
      gamma <- as.numeric(strsplit(line, '[:]')[[1]][2])
      list.param <- c(list.param, gamma = gamma) 
    }
    
    #Costs
    if (idx == 3){
      cost <- as.numeric(strsplit(line, '[:]')[[1]][2])
      list.param <- c(list.param, cost = cost) 
    }
    
    if (idx == 4){
      kernel <- strsplit(line, '[:]')[[1]][2]
      list.param <- c(list.param, kernel = kernel) 
    }
  }
  close(con)
  return (list.param)
}