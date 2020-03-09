## ----------- Extracting layer information ----------- 

#Returns a data frame with the number of spots/voxels falling in each cortical layer for each cluster
get.layer.count <- function(df.cl.acronym, bool.return.proportion = T) {
  #df.cl.acronym(data.frame): rows are spots/voxels,
  #                           1st column is 'cluster' with cluster name / ID
  #                           2nd column is 'acronym' with brain region acronym
  #bool.return.proportion(logical): if TRUE, returning the proportion in each layer of spots/voxels.
  #                                 if FALSE, returning absolute count
  

  #Loading the layer table
  load(paste(path.bin, 'layerstable.RData', sep = '/'))
  
  #Mapping the layer in the spots table
  df.cl.acronym$layer <- mapvalues(df.cl.acronym$acronym,
                                   from = layers.table$acronym,
                                   to = layers.table$layer,
                                   warn_missing = FALSE)
  
  #Non cortical region get the "NA" acronym
  df.cl.acronym[!is.element(df.cl.acronym$acronym, layers.table$acronym), 'layer'] <-
    NA
  
  #List of cluster
  list.cl <- unique(as.character(df.cl.acronym$cluster))
  
  #List of layers, replacing NA by clearer name
  list.layers <- unique(df.cl.acronym$layer)
  list.layers <- sort(list.layers[!is.na(list.layers)])
  list.layers <- c(list.layers, 'non.cortical')
  
  #Generating data frame that will contain all information
  df.cl.layers <- as.data.frame(matrix(nrow = length(list.cl), ncol = length(list.layers)))
  
  #Rows: clusters, cols: layers
  rownames(df.cl.layers) <- list.cl
  colnames(df.cl.layers) <- list.layers
  
  #Making a simple loop to fill the data frame
  for (cl in rownames(df.cl.layers)) {
    for (lay in colnames(df.cl.layers)) {
      df.cl.layers[cl, lay] <-
        sum(df.cl.acronym[df.cl.acronym$cluster == cl, 'layer'] == lay, na.rm = T)
    }
    
    #Handling NA case
    df.cl.layers[cl, 'non.cortical'] <-
      sum(is.na(df.cl.acronym[df.cl.acronym$cluster == cl, 'layer']))
  }
  
  #If required, normalizing each cluster by number of spots / voxels
  if (bool.return.proportion) {
    norm.vect <- rowSums(df.cl.layers)
    df.cl.layers <- df.cl.layers / norm.vect
  }
  
  #Returning the generated data frame
  return(df.cl.layers)
  
}

#Returns the layer count for cluster based on the voxels volume (3D) as a list for every AP coordinates
get.layer.from.voxels <- function(grid.3d, grid.parameters, parallel.cores = 0){
  
  #grid.3d(num 3D array): voxels
  #grid.paramaters(num 3 lists): ML/DV/AP coordinates of points in the grid
  #parallel.cores(integer): nunber of cores to use for parallel processing. If 0, no parallel processing.
  
  if (parallel.cores == 0){
    
    #Initializing list that will contain all information
    list.layers.df <- NULL
    
    #Looping through AP coordinates
    for (i in 1:dim(grid.3d)[3]){
      
      #Converting 3d grid into 2d grid
      grid.2d <- grid.3d[,,i]
      
      #Converting parameters to 2D parameters (only one AP to call)
      gp <- grid.parameters
      gp$AP <- grid.parameters$AP[i]
      
      #Using the 2D layer function
      list.layers.df[[i]] <- get.layer.from.pixels(grid.2d, gp)
      
      #Naming the list with current ap
      names(list.layers.df)[i] <- grid.parameters$AP[i]
    }
  }else{
    
    library(foreach)
    library(doParallel)
    
    #Registering cores
    registerDoParallel(cores=parallel.cores)
    
    #For each loop
    list.layers.df <- foreach(i=1:(dim(grid.3d)[3]),
                              .packages=c('wholebrain'),
                              .export=c('get.layer.from.pixels',
                                        'path.bin',
                                        'get.acronym.depth',
                                        'get.layer.count')) %dopar% {
                                          
                                          gp <- grid.parameters
                                          gp$AP <- grid.parameters$AP[i]
                                          get.layer.from.pixels(grid.3d[,,i], gp)
                                        }
  }
  
  #Returning the generated list
  return(list.layers.df)
}

#Returns the layer for count for cluster based on the pixels surface (2D)
get.layer.from.pixels <- function(grid.2d, grid.parameters){
  
  #Loading the stereotaxic atlac
  load(paste(path.bin, 'atlasstereo.RData', sep='/'))
  
  #Expanding grid and converting array into a data frame with 3 cols: ML/DV/cluster
  coord <- expand.grid(grid.parameters$ML, grid.parameters$DV)
  df.pixels <- data.frame(cluster = as.numeric(grid.2d),
                          ML = coord[,1],
                          DV = coord[,2])
  
  #Discarding outside of brain outline
  df.pixels <- df.pixels[df.pixels$cluster != 0,]
  
  #Selecting plate and outline of interest
  plate <- which(atlas.stereo$AP == grid.parameters$AP)
  o <- atlas.stereo$outlines[[plate]]
  
  #Getting depth of the different outlines
  depth <- get.acronym.depth(as.character(atlas.stereo$plate.info[[plate]]$acronym))
  
  #2d points X acronyms. Used to chose the deepest acronym for each point
  mat.acronym <- matrix(FALSE, nrow = dim(df.pixels)[1], ncol = length(o))
  
  #Looping through all outlines
  for (i in 1:length(o)){
    
    #Detecting point inside an outline
    bool.is.inside.outline <- point.in.polygon(df.pixels[,'ML'], df.pixels[,'DV'],o[[i]]$ML, o[[i]]$DV)
    
    #Updating the acronym matrix
    mat.acronym[which(bool.is.inside.outline == TRUE),i] <- TRUE
  }
  
  #For some reason, few points get outside of outlines. Discarding them.
  kept.points <- rowSums(mat.acronym) > 0
  print(sprintf('Points discarded: %d/%d (outside of outline while not supposed to be)',
                sum(rowSums(mat.acronym) <= 0), dim(df.pixels)[1]))
  mat.acronym <- mat.acronym[kept.points,]
  
  
  #For each 2d point, find the deepest acronym possible
  acronym.list <- apply(mat.acronym,1,function(x){names(which.max(depth[x]))})
  
  #counting per layer
  df.layers <- (get.layer.count(data.frame(cluster = df.pixels[kept.points,'cluster'], 
                                           acronym =  acronym.list,
                                           stringsAsFactors = FALSE),
                                FALSE))
  
  #Returning the data frame with info about layers
  return(df.layers)
  
}

#Merged the list of layers information from the voxels layering into a simple data frame (losing AP information)
get.merged.layer.from.voxels <- function(list.layers, min.cortical.proportion = 0.5, bool.normalize.to.proportion = TRUE){
  
  #list.layers(data frame list): list of layer from pixels for each AP coordinates of voxels
  #min.cortical.proportion(num): clusters with more than this proportion of non-cortical layer will be discarded
  #bool.normalize.to.proportion(logical): if TRUE, counts are converted into proportions
  
  #Will contain all clusters and layers part of the list, with repetition
  vect.cl <- numeric()
  layers.vect <- character()
  
  #Looping through every AP coordinates
  for(i in 1:length(list.layers)){
    
    #Appending plate clusters and layers
    vect.cl <- c(vect.cl,as.numeric(rownames(list.layers[[i]])))
    layers.vect <- c(layers.vect, colnames(list.layers[[i]]))
  }
  
  #Sorting and selecting unique clusters/layers
  clusters <- sort(unique(vect.cl))
  layers <- sort(unique(layers.vect))
  
  #Defining data frame that will contain counts
  #Rows: clusters, Cols: layers
  df.layers.count <- as.data.frame(matrix(0, nrow = length(clusters), ncol = length(layers)))
  rownames(df.layers.count) <- clusters
  colnames(df.layers.count) <- layers
  
  #Looping through every AP coordinates
  for(i in 1:length(list.layers)){
    
    #Looping through layers
    for(j in colnames(list.layers[[i]])){
      
      #Adding the count
      df.layers.count[rownames(list.layers[[i]]),j] <- df.layers.count[rownames(list.layers[[i]]),j] + list.layers[[i]][,j]
    }
  }
  
  #Normalizing to get repartition
  df.layers.count.norm <- df.layers.count / rowSums(df.layers.count)
  
  #Normalizing if required
  if(bool.normalize.to.proportion)
    df.to.return <- df.layers.count.norm
  else
    df.to.return <- df.layers.count
  
  #Keeping only relevant clusters
  df.to.return <- df.to.return[df.layers.count.norm$non.cortical <= (1-min.cortical.proportion),]
  
  return(df.to.return)
  
}

#Return the main layer for each cluster based on the spots table
get.cluster.layer.from.spots.table <- function(spots.table, min.cortical.proportion = 0.5){
  
  #spots.table
  #min.cortical.proportion(num): minimum proportion of spots that should be cortical in order for the cluster to receive a layer
  
  #Getting the layer count
  l.count <- get.layer.count(data.frame(cluster = spots.table$clusters.named,
                                        acronym = spots.table$acronym,
                                        stringsAsFactors = F),
                             bool.return.proportion = T)
  
  #Discarding non cortical clusters and the non-cortical column
  l.count <- l.count[l.count$non.cortical < (1-min.cortical.proportion),1:(dim(l.count)[2]- 1)]
  
  #Finding the dominant layer
  l.max <- colnames(l.count)[apply(l.count, 1, which.max)]
  
  #Adding cluster names
  names(l.max) <- rownames(l.count)
  
  #Returning the list
  return(l.max)
  
}
