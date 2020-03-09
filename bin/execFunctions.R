## ----------- Spots table manipualation ----------- 

#Add the clusters and the clusters named columns to a spot table from the seurat object or a .tsv table
append.cluster.to.spots.table <- function(spots.table, cluster.list, name.clusters = TRUE, min.cluster.size = 0, new.names = TRUE){
  
  #spots.table
  #clusters(str or Seurat): clusters are in @ident if Seurat, path to a .tsv file containing clusters otherwise
  #name.clusters(bool): if true, clusters get a name assigned. Otherwise, the ID is used as a name
  #min.cluster.size(int): clusters with less elements than this value will be discarded
  #new.naming(bool): use the new way for naming clusters
  
  #Case where a path to .tsv was sent
  if (is.character(cluster.list)){
    
    t <- read.table(cluster.list,stringsAsFactors = FALSE)
    spots.table <- spots.table[t[,1],]
    spots.table$cluster <- NULL
    spots.table$clusters.named <- NULL
    t[is.na(t[,2]),2] <- -1
    spots.table[t[,1],'cluster'] <- factor(as.numeric(t[,2]))
    
    #Case where a seurat object was sent      
  }else{
    
    seur.obj <- cluster.list
    clusters <- seur.obj@ident
    clusters <- clusters[!is.na(clusters)]
    
    spots.table <- spots.table[names(clusters),]
    
    if (is.null(spots.table$acronym.parent))
      spots.table <- add.parent.acronym(spots.table)
    
    spots.table[names(clusters),'cluster'] <- clusters
  }
  
  #Discarding clusters than are too small
  if(min.cluster.size > 0){
    cnt <- plyr::count(spots.table$cluster)
    cnt <- cnt[cnt$freq >= min.cluster.size,]
    spots.table <- spots.table[sapply(spots.table$cluster,function(x){return(is.element(x,cnt$x))}),]
  }
  
  #Eventually naming clusters
  if (name.clusters){
    
    #old naming
    if(!new.names){
      
      #Adding parents acronyms if not already there
      if (is.null(spots.table$acronym.parent))
        spots.table <- add.parent.acronym(spots.table, c('TH','HY','MB','HB',
                                                         'CB','STR','PAL','Isocortex','OLF',
                                                         'HPF','CTXsp','root','fiber tracts','VS'))
      spots.table <- name.clusters(spots.table, cutoff.majority = 50, new.names = FALSE)
      
    #new naming
    }else{
      if (is.null(spots.table$acronym.parent))
        spots.table <- add.parent.acronym(spots.table, c('TH','HY','MB','HB',
                                                        'CB','STR','PAL','Isocortex','OLF',
                                                        'HIP','RHP', 'CTXsp','root','fiber tracts','VS'))
      spots.table <- name.clusters(spots.table, cutoff.majority = 40, new.names = TRUE)
    }
    
  }else{
    spots.table$clusters.named <- factor(spots.table$cluster)
  }
  
  return(spots.table)
}

#Export a matrix or a vector as a string to be loaded into R 
export.matrix.string <- function(M, type = 'numeric'){
  
  #M: matrix or vector
  #type: if numeric, exported as numerical. If anything else, exported as a string
  
  #Collapsing a numeric matrix / vector
  if (type == 'numeric'){
    
    #Matrix
    if (!is.null(dim(M)[2])){
      str <- paste(apply(M,1,function(x){paste(x,sep='',collapse=',')}), collapse = '),c(')
      str <- paste('rbind(c(',str,'))',sep='')
      
      #Vector
    }else{
      str <- paste(M,collapse = ',')
      str <- paste('c(', str, ')', sep='')
    }
    
    #Collapse chars
  }else
    
    #Matrix
    if (!is.null(dim(M)[2])){
      str <- paste(apply(M,1,function(x){paste(x,sep='',collapse="','")}), collapse = "'),c('")
      str <- paste("rbind(c('",str,"'))",sep='')
      
      #Vector
    }else{
      str <- paste(M,collapse = "','")
      str <- paste("c('", str, "')", sep='')
    }
  
  
  return(str)
}

#Write a table in memory in standardize .tsv format
write.mat <- function(mat, file){
  
  #mat(df): dataframe to save
  #file(char): file path + name to save
  
  #Error if file exists already
  if(file.exists(file))
    stop(paste(file,'File already exists.', sep=' - '))
  
  write.table(mat, file=file, sep='\t', quote=FALSE, col.names = NA)
}


## ----------- Clusters ----------- 

#Returns a default df.colors with the inputed clusters
get.default.df.colors <- function(cl.list){
  
  #cl.list(char): list of clusters that must be included in the data frame. Should be as  the clusters.named field
  
  #If more clsuters than the palette size, interpolating
  if (length(cl.list) > 9){
    f <- colorRampPalette(brewer.pal(9,'Set1'))
    cp <- f(length(cl.list))
    
    #Otherwise directly selecting from the palette
  }else{
    cp <- brewer.pal(length(cl.list),'Set1')
  }
  
  #Creating the data frame
  df.col <- data.frame(cluster = cl.list, color = cp, stringsAsFactors = FALSE)
  rownames(df.col) <- df.col$cluster
  
  #Returning
  return(df.col)
  
}

#Name the clusters
name.clusters <- function(spots.table, cutoff.majority = 50, new.names = TRUE){
  
  cl.list <- unique(as.numeric(spots.table$cluster))
  
  df.cl.name <- data.frame(clusters.named = character(length(cl.list)),
                           cluster = cl.list,
                           count= numeric(length(cl.list)))
  
  df.cl.name$clusters.named <- character(dim(df.cl.name)[1])
  rownames(df.cl.name) <- df.cl.name$cluster
  df.cl.name <- df.cl.name[order(df.cl.name$cluster),]
  df.cl.name$count <- plyr::count(spots.table$cluster)[rownames(df.cl.name),'freq']
 
  #Contains current counts to append id
  cur.count <- numeric()
  
  #Looping through all rows
  for(r in rownames(df.cl.name)){
    
    sp.sub <- subset(spots.table, cluster == r)
    cnt <- plyr::count(sp.sub$acronym.parent)
    cnt$prop <- 100*cnt$freq / sum(cnt$freq)
    cnt <- cnt[order(cnt$freq, decreasing = T),]
    cnt$x <- as.character(cnt$x)
    
    if(cnt[1,'prop'] < cutoff.majority)
      cl.name.prepend <- 'Mixed'
    else
      cl.name.prepend <- as.character(name.from.acronym(cnt[1,'x']))
    
    if(is.element(cl.name.prepend, names(cur.count))){
      id <- cur.count[cl.name.prepend] + 1
      cur.count[cl.name.prepend] <- id
    }else{
      id <- 1
      former.names <- names(cur.count)
      cur.count <- c(cur.count, id)
      names(cur.count) <- c(former.names, cl.name.prepend)
    }
    
    df.cl.name[r,'clusters.named'] <- sprintf('%s-%02d',cl.name.prepend,id)
    
    #If mixed and new naming method
    if(cl.name.prepend == 'Mixed' & new.names){
      
      df.cl.name[r, 'clusters.named'] <- sprintf('%s (%s)', 
                                                 df.cl.name[r, 'clusters.named'], 
                                                 as.character(name.from.acronym(cnt[1,'x'])))
    }
  }
  
  spots.table$clusters.named <- mapvalues(spots.table$cluster,
                                          from = df.cl.name$cluster,
                                          to = df.cl.name$clusters.named)
  spots.table$clusters.named <- factor(spots.table$clusters.named)
  spots.table$clusters.named <- factor(spots.table$clusters.named, levels(spots.table$clusters.named)[order(tolower(levels(spots.table$clusters.named)))])
  
  return(spots.table)
  
}

find.matching.clusters <- function(clust.file.1, clust.file.2, min.cluster.size = 0, spots.table = NULL){
  
  if(is.null(spots.table))
    spots.table <- add.parent.acronym(load.spots.table())
  
  spots.table.1 <- append.cluster.to.spots.table(spots.table,clust.file.1, min.cluster.size = min.cluster.size)
  spots.table.2 <- append.cluster.to.spots.table(spots.table,clust.file.2, min.cluster.size = min.cluster.size)
  
  #Overlap is defined as 1 - number of spots in 1 and not in 2 + number in 2 and not in 1 devided by number of unique spots 
  df <- as.data.frame(matrix(nrow=length(unique(spots.table.1$clusters.named)),ncol=length(unique(spots.table.2$clusters.named))))
  rownames(df) <- unique(spots.table.1$clusters.named)
  colnames(df) <- unique(spots.table.2$clusters.named)
  
  for(cl.1 in unique(spots.table.1$clusters.named)){
    spots.1 <- rownames(spots.table.1[spots.table.1$clusters.named == cl.1,])
    for (cl.2 in unique(spots.table.2$clusters.named)){
      spots.2 <- rownames(spots.table.2[spots.table.2$clusters.named == cl.2,])
      spots.not.in.common <- length(union(setdiff(spots.1,spots.2),setdiff(spots.2,spots.1))) / (length(spots.1) + length(spots.2))
      df[cl.1,cl.2] <- 1-spots.not.in.common
    }
  }
  
  t <- data.frame(apply(df,1,which.max))
  colnames(t) <- 'id'
  t$best.match <- unique(spots.table.2$clusters.named)[t$id]
  for(i in 1:length(rownames(t))){
    t$ratio[i] <- df[rownames(t)[i],t$id[i]]
  }
  
  t$percent <- sprintf('%.1f',100*t$ratio)
  t <- t[order(rownames(t)),]
  return(t)
}



generate.cluster.range <- function(spots.table, seur.obj, resolution.range, dir.path, parent.list = c('VS','TH','STR','RHP','P','PAL','OLF','MY','MB','HY','HIP','fiber tracts','CTX','CB'), n.ica = 80){
  
  if(dir.exists(dir.path))
    stop('Directory already exists')
  
  dir.create(dir.path)
  
  spots.table <- add.parent.acronym(spots.table, list.acronym.parents = c(parent.list,'root'))
  
  clusters.path <- NULL
  
  if (length(n.ica) == 1){
    ica.range <- 1:n.ica
  }else{
    ica.range <- n.ica
  }
  
  for (res in resolution.range){
    
    print(paste('Computation starting for resolution: ',res,sep=''))
    
    seur.obj <- FindClusters(object = seur.obj, reduction.type = "ica", dims.use = ica.range, resolution = res,
                             save.SNN = FALSE, n.start = 100, nn.eps = 0, print.output = FALSE, force.recal = TRUE)  
    
    print(paste('Computation finished for resolution: ',res,'. Exporting results.',sep=''))
    
    p <- paste(dir.path,'/res-',res,'.tsv',sep='')
    export.cluster.table(seur.obj, p)
    clusters.path <- c(clusters.path,p)
    
    print(paste('Results exported for resolution: ',res,sep=''))
    
  }
  return(get.S.value.per.clusters(spots.table,clusters.path,parent.list))
}

#Write the clusters from a Seurat object into a .tsv table
export.cluster.table <- function(seur.obj, path){
  #seur.obj(Seurat): object that contains the clusters in @ident
  #path(str): path + name (including extension) of the table to save.
  
  table.export <- data.frame(cluster = as.numeric(seur.obj@ident))
  rownames(table.export) <- names(seur.obj@ident)
  write.table(table.export, file=path, sep='\t', quote=FALSE, col.names = FALSE)
}

get.expression.cluster <- function(spots.table, genes, cluster.names){
  
  # Loading the expression matrix if not loaded yet
  if (!exists('st.data'))
    load.expr.mat()
  
  n.c <- length(cluster.names)
  n.g <- length(genes)
  n.r <- n.c * n.g
  
  df <- data.frame(gene = character(n.r), cluster = character(n.r), mean = numeric(n.r), sd = numeric(n.r), prop.gene.expr = numeric(n.r), stringsAsFactors = FALSE) 
  
  idx <- 0
  for(j in 1:length(genes)){
    for (i in 1:length(cluster.names)){
      
      idx <- idx + 1
      selected.spots <- spots.table[spots.table$clusters.named == cluster.names[i],]
      data <- st.data[rownames(selected.spots),genes[j]]
      
      df[idx, 'gene'] <- genes[j]
      df[idx, 'cluster'] <- cluster.names[i]
      df[idx, 'mean'] <- mean(data)
      df[idx, 'sd'] <- sd(data)
      df[idx, 'prop.gene.expr'] <- 100*sum(data>0)/(dim(selected.spots)[1])
    }
  }
  
  df$gene <- factor(df$gene)
  df$cluster <- factor(df$cluster)
  
  return(df)
}

## ----------- Expression ----------- 

#Returning spots with coordinates, expression level and color scale
get.expression.slice <- function(gene, spots.table, st.data, min.expr, log.scale, normalized.scale, slice_index, rbPal){
  
  #gene(str): gene of intereset
  #spots.table
  #st.data
  #min.expr(float): between 0 and 1, indicates percentile of expression. Above 1: integer cutoff value
  #log.scale(bool): using a logarithmic scale
  #normalized.scale(bool): normalized expression by peak expression across the brain
  #slice_index(char): vector of slices indexes
  #rbPal(func): palette function from RColorBrewer
  
  dataset <- NULL
  
  for (i in 1:length(slice_index)){
    
    #Selecting appropraite spots
    dataset[[i]] <- spots.table[spots.table$slice_index == slice_index[i],c('ML','DV','AP','acronym','x','y')]
    
    #Extracting the expression of gene of interest
    dataset[[i]]$expr <- st.data[rownames(dataset[[i]]),gene]
    
    #This adds a column of color values
    # based on the y values
    if (min.expr >= 1)
      dataset[[i]] <- dataset[[i]][dataset[[i]]$expr >= min.expr,]
    else if(min.expr >0){
      dataset[[i]] <- dataset[[i]][(dataset[[i]]$expr > as.numeric(quantile(dataset[[i]]$expr,min.expr))) == TRUE,]
    }
    
    if (normalized.scale)
      dataset[[i]]$expr <- (dataset[[i]]$expr / max(st.data[,gene]))*100
    
    m.val <- max(dataset[[i]]$expr)
    if (log.scale){
      dataset[[i]]$expr[dataset[[i]]$expr==0] <- 0.98
      dataset[[i]]$expr <- dataset[[i]]$expr + 0.01
      dataset[[i]]$expr <- ceiling(5*log(dataset[[i]]$expr))
    }
    
    if (normalized.scale){
      if (log.scale)
        dataset[[i]]$col <- rbPal(ceiling(5*log(101))+1)[(as.numeric(dataset[[i]]$expr)+1)]
      else
        dataset[[i]]$col <- rbPal(101)[(as.numeric(dataset[[i]]$expr)+1)]
      m.val <- 100
    }else
      dataset[[i]]$col <- rbPal(max(dataset[[i]]$expr)+1)[(as.numeric(dataset[[i]]$expr)+1)]
  }
  dataset[[(length(dataset)+1)]] <- m.val
  return(dataset)
}

## ----------- Selection ----------- 

#Returns the spots table with only spots that are in the regions from acronym list
select.specific.sub.areas <- function(spots.table, acronyms.list){
  #spots.table
  #acronym.slist(str): list of acronyms to select
  
  spots.table <- add.all.acronyms(spots.table)
  spots.table$selected <- sapply(spots.table$all.acronyms, function(x){return(length(intersect(acronyms.list,x)) > 0)})
  
  return(spots.table[spots.table$selected,])
}

#Returns the closest AP value which is element of the spots.table
select.closest.AP <- function(AP, spots.table){
  #AP(float): desired AP
  #spots.table
  
  AP.list <- unique(spots.table$AP)
  AP.diff <- abs(AP.list - AP)
  return(AP.list[which.min(AP.diff)])
}

## ----------- Generation  ----------- 

#Old version of create.2d.atlas, not relevant anymore
create.web.2d.projection <- function(spots.table, spots.table.wb, seur.obj, st.data, directory, name, discard.posterior.slices = TRUE, n.markers = 5){
  
  library(RColorBrewer)
  
  dir.path <- paste(directory,name,sep='/')
  
  if (dir.exists(dir.path))
    stop('Directory already exists')
  
  dir.create(dir.path)
  
  list.coord <- sort(unique(spots.table$AP),decreasing = TRUE)
  
  #Very posterior slices uses 2 or 3 different registration, the remapping might not work
  if (discard.posterior.slices)
    list.coord <- list.coord[list.coord > -4.9]
  
  col.tab.paired <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff00','#b15928')
  col.vect <- rep(col.tab.paired,40)
  col.vect <- col.vect[1:length(levels(spots.table.wb$clusters.named))]
  p <- plot.3d.all(spots.table.wb, col.vect, title = 'Clustering across the whole brain')
  p <- htmlwidgets::saveWidget(p, sprintf('%s/3d-plot-wholebrain.html',dir.path), selfcontained = FALSE)
  
  # col.vect <- colorRampPalette(brewer.pal(9, "Set1"))(length(levels(spots.table$clusters.named)))
  # col.vect <- col.vect[1:length(unique(spots.table$clusters.named))]
  p <- plot.3d.all(spots.table, col.vect, title = 'Clustering only in specific area')
  p <- htmlwidgets::saveWidget(p, sprintf('%s/3d-plot.html',dir.path), selfcontained = FALSE)
  
  idx <- 0
  
  for (coord in list.coord){
    idx <- idx + 1
    print(sprintf('Computing slice %d/%d',idx,length(list.coord)))
    p <- plot.2d.atlas.clust(spots.table,coord)
    p <- htmlwidgets::saveWidget(p, sprintf('%s/cut-%d.html',dir.path,idx), selfcontained = FALSE)
  }
  
  str.table <- '<tr><th>Link to section</th><th>AP</th><th>Slices</th><th>Spots number</th></tr>'
  idx <- 0
  for (coord in list.coord){
    idx <- idx + 1
    cut <- sprintf('Section %d',idx)
    spots.count <- dim(spots.table[spots.table$AP == coord,])[1]
    slices <- paste(unique(spots.table[spots.table$AP == coord,"slice_index"]),collapse=' & ')
    path.to.plate <- sprintf('cut-%d.html',idx)
    str.table <- c(str.table, sprintf('<tr><th><a target="_blank" href="%s">%s</a></th><th>%s</th><th>%s</th><th>%s</th></tr>',
                                      path.to.plate,cut,coord,slices,spots.count))
  }
  
  # clusters.matching <- find.matching.clusters(spots.table,spots.table.wb)
  # str.matching <- '<tr><th>Cluster from specific area</th><th>Best matching cluster from whole brain</th><th>Overlap percentage</th></tr>'
  # for (i in 1:dim(clusters.matching)[1]){
  #   str.matching <- c(str.matching, sprintf('<tr><th>%s</th><th>%s</th><th>%s</th></tr>',
  #                                           rownames(clusters.matching)[i],clusters.matching$best.match[i],clusters.matching$percent[i]))  
  # }
  
  
  str.gene.expr <- ''
  # str.gene.expr <- '<tr><th>Cluster</th><th>Gene 1</th><th>Gene 2</th><th>Gene 3</th><th>Gene 4</th><th>Gene 5</th></tr>'  
  # 
  # ct <- 0
  # for(i in levels(spots.table$clusters.named)){
  #   
  #   ct <- ct + 1
  #   print(sprintf('Computing markers: %d/%d',ct,length(levels(spots.table$clusters.named))))
  #   
  #   clust.id <- unique(spots.table$cluster[spots.table$clusters.named == i])
  #   clust.id <- as.numeric(as.character(clust.id[!is.na(clust.id)]))
  #   markers <- FindMarkers(seur.obj, clust.id, test.use = 'tobit')
  #   markers <- markers[1:n.markers,]
  #   
  #   p <- plot.3d.diff.expr(rownames(markers), spots.table, st.data, i, percentile = 0.85)
  #   path.to.mark <- sprintf('markers-%d.html',clust.id)
  #   htmlwidgets::saveWidget(p, sprintf('%s/%s',dir.path,path.to.mark), selfcontained = FALSE)
  #   str.gene.expr <- c(str.gene.expr,'<tr>')
  #   str.gene.expr <- c(str.gene.expr,sprintf('<th><a target="_blank" href="%s">%s</a></th>',path.to.mark,i))
  #   for (j in 1:n.markers){
  #     str.gene.expr <- c(str.gene.expr, sprintf('<th>%s, p-val: %.1e</th>',rownames(markers)[j],markers$p_val_adj[j]))
  #   }
  #   
  #   str.gene.expr <- c(str.gene.expr,'</tr>')
  # }
  
  
  
  
  
  fileConn<-file(sprintf('%s/index.html',dir.path))
  writeLines(c('<!DOCTYPE html>',
               '<html>',
               '<style>',
               'table, th, td {',
               '    border: 1px solid black;',
               '    border-collapse: collapse;',
               '}',
               'th, td {',
               '    padding: 5px;',
               '}',
               '</style>',
               '<center>',
               sprintf('<h1>%s</h1>',name),
               '<h2>3d-plot</h2>',
               '<a target="_blank" href="3d-plot.html">Clustering only in specific area</a><br />',
               '<a target="_blank" href="3d-plot-wholebrain.html">Clustering accross the wole brain</a>',
               '<h2>Coronal atlas</h2>',
               '<table>',
               str.table,
               '</table>',
               '<h2>Clusters overlap</h2>',
               '<table>',
               str.matching,
               '</table>',
               '<h2>Top markers</h2>',
               '<table>',
               str.gene.expr,
               '</table>',
               '</center>',
               '</html>')
             , fileConn)
  close(fileConn)
  
  
}

#Generates a 2D projection of the clusters in the atlas, with index file and 3d plot
create.2d.atlas <- function(dir.path, cluster.list = NULL, spots.table = NULL, col.vect = NULL, min.cluster.size = 0, order.by.cluster.size = TRUE, minimum.spots.slice.display = 0, name.clusters = TRUE){
  
  library(plyr)
  library(RColorBrewer)
  
  #Avoiding erasing previous atlas by mistake
  if (dir.exists(dir.path))
    stop('Directory already exists')
  
  #Creating directory
  dir.create(dir.path)
  
  #Loading spots table if not sent
  if (is.null(spots.table))
    spots.table <- add.parent.acronym(load.spots.table())
  
  #Appending the clusters from the cluster file if they are not in the spots table
  if (is.null(spots.table$cluster))
    spots.table <- append.cluster.to.spots.table(spots.table, cluster.list, name.clusters, min.cluster.size)
  
  #Counting cluster size
  cnt <- plyr::count(spots.table$clusters.named)
  cnt <- cnt[order(cnt$freq, decreasing = TRUE),]
  
  #Setting up colors if required
  if (is.null(col.vect)){
    
    #We want similar color for the 3d plot and the 2D plot
    col.vect <- data.frame(clusters.named = as.character(unique(spots.table$clusters.named)), stringsAsFactors = FALSE)
    col.base <- colorRampPalette(brewer.pal(9, "Set1"))(18)
    col.base <- rep(col.base,40)
    col.vect <- data.frame(clusters.named = cnt$x,
                           color = col.base[1:dim(cnt)[1]],
                           stringsAsFactors = FALSE)
  }
  
  #Getting color in proper order
  col <- mapvalues(levels(spots.table$clusters.named),from = col.vect$clusters.named, to = col.vect$color)
  
  #3d plot
  p <- plot.3d.all(spots.table, col.vect = col, title = '3d plot')
  path.3d <- sprintf('%s/%s/3d-plot.html',getwd(),dir.path)
  htmlwidgets::saveWidget(p, path.3d, selfcontained = FALSE, title = '3d-plot', libdir = 'lib')
  
  #List of all AP coordinates ordered
  AP.list <- sort(unique(spots.table$AP), decreasing = TRUE)
  
  #Going through all sections
  for(i in 1:length(AP.list)){
    
    #Selecting appropriate AP and spots
    AP <- AP.list[i]
    sub.spots.table <- spots.table[spots.table$AP == AP,]
    
    #Plotting on the atlas
    p <- plot.2d.atlas.clust(sub.spots.table, AP, col.vect = col.vect, order.by.cluster.size = order.by.cluster.size, minimum.spots = minimum.spots.slice.display)
    
    #Saving the section
    section.path <- sprintf('%s/%s/section-%d.html',getwd(),dir.path,i)
    htmlwidgets::saveWidget(p, section.path, selfcontained = FALSE, title = AP, libdir = 'lib')
    
    #Path to previous and next section
    prev.section.path <- sprintf('section-%d.html',(i-1))
    next.section.path <- sprintf('section-%d.html',(i+1))
    
    #Appending the section links for easier navigation (specific case for first and last section)
    if(i == 1)
      add.navigation.2d.atlas(section.path, next.section = next.section.path)
    else if (i == length(AP.list))
      add.navigation.2d.atlas(section.path, previous.section = prev.section.path)
    else
      add.navigation.2d.atlas(section.path, prev.section.path, next.section.path)
  }
  
  create.index.file(dir.path, spots.table)
  
}

#Generates an index file for the 2D atlas in the dir.path directory
create.index.file <- function(dir.path, spots.table){
  
  #dir.path(str): path to the directory
  #spots.table
  
  #Parsing the name from the directory path
  str <- strsplit(dir.path,'/')[[1]]
  name <- str[length(str)]
  
  #Counting the number of cluster in each region
  list.regions <- unique(as.character(spots.table$full.name.parent))
  list.regions <- c(list.regions, 'Mixed')
  list.clusters <- unique(as.character(spots.table$clusters.named))
  
  #df that gonna contains the value
  df.cl.count <- data.frame(matrix(nrow =length(list.regions), ncol = 2))
  rownames(df.cl.count) <- list.regions
  colnames(df.cl.count) <- c('n.clusters','n.spots')
  
  #Looping through the cluster names and counting
  for(r in rownames(df.cl.count)){
    cl <- list.clusters[startsWith(list.clusters,r)]
    df.cl.count[r,'n.clusters'] <- length(cl)
    df.cl.count[r,'n.spots'] <- sum(sapply(spots.table$clusters.named,function(x){
      return(is.element(x,cl))}))
  }
  
  #Sorting df
  df.cl.count <- df.cl.count[order(df.cl.count$n.clusters, decreasing = TRUE),]
  print(df.cl.count)
  
  str.cluster.number <- '<tr><th>Region</th><th>N clusters</th><th>N spots</th></tr>'
  for (i in 1:(dim(df.cl.count)[1])){
    str.cluster.number <- c(str.cluster.number,paste(
      '<tr><th>',
      rownames(df.cl.count)[i],
      '</th><th>',
      df.cl.count[i, 'n.clusters'],
      '</th><th>',
      df.cl.count[i, 'n.spots'],
      '</th></tr>'
    ))
  }
  
  
  
  #List of all AP coordinates ordered
  AP.list <- sort(unique(spots.table$AP), decreasing = TRUE)
  
  #Initialization
  idx <- 0
  str.table <- '<tr><th>See in atlas</th><th>AP</th><th>Sections</th><th>Number of spots</th></tr>'
  
  #Going through AP and appending text
  for (AP in AP.list){
    idx <- idx + 1
    section <- sprintf('Section %d',idx)
    spots.count <- dim(spots.table[spots.table$AP == AP,])[1]
    slices <- paste(unique(spots.table[spots.table$AP == AP,"slice_index"]),collapse=' & ')
    path.to.plate <- sprintf('section-%d.html',idx)
    str.table <- c(str.table, sprintf('<tr><th><a target="_blank" href="%s">%s</a></th><th>%s</th><th>%s</th><th>%s</th></tr>',
                                      path.to.plate,section,AP,slices,spots.count))
  }
  
  #Openning file
  fileConn<-file(sprintf('%s/index.html',dir.path))
  
  #Writitng lines
  writeLines(c('<!DOCTYPE html>',
               '<html>',
               '<style>',
               'table, th, td {',
               '    border: 1px solid black;',
               '    border-collapse: collapse;',
               '}',
               'th, td {',
               '    padding: 5px;',
               '}',
               '</style>',
               '<center>',
               sprintf('<h1>%s: %d clusters</h1>',name, length(unique(spots.table$clusters.named))),
               '<h2>3d-plot</h2>',
               '<a target="_blank" href="3d-plot.html">3d plot</a><br />',
               '<h2>Repartition of clusters</h2>',
               '<table>',
               str.cluster.number,
               '</table>',
               '<h2>Coronal atlas</h2>',
               '<table>',
               str.table,
               '</table>',
               '</center>',
               '</html>')
             , fileConn)
  
  #Closing file
  close(fileConn)
}

#Adds link to open the previous and the next section from the 2d plotly projection
add.navigation.2d.atlas <- function(html.path, previous.section = NULL, next.section = NULL){
  
  #html.path(str): path to the file to edit
  #previous.section(str): path to the previous section (if null, nothing written)
  #next.section(str): path to the next section (if null, nothing written)
  
  #Loading the html file
  str <- readChar(html.path, file.info(html.path)$size)
  
  #List of breaks and strings to interpolate
  interp.1 <- '<body style="background-color:white;">'
  
  interp.2 <- '<div style = "width:90%;">\n<p style="text-align:left;">'
  
  if(!is.null(previous.section))
    interp.2 <- paste(interp.2,
                      '<a href="',
                      previous.section,
                      '">Previous section</a>',
                      sep='')
  
  if(!is.null(next.section))
    interp.2 <- paste(interp.2,
                      '<span style="float: right"><a href="',
                      next.section,
                      '">Next section</a></span>',
                      sep='')
  
  interp.2 <- paste(interp.2,
                    '</div>\n<div style = "position:absolute; top: 10%; width: 100%; height: 90%">',
                    sep='')
  
  interp.3 <- '</div>'
  interp.4 <- '<script type="application/json"'
  
  
  #Splitting according to the breaks
  split <- strsplit(str,interp.1)[[1]]
  split.1 <- split[[1]]
  split.tmp <- split[[2]]
  
  split <- strsplit(split.tmp,interp.4)[[1]]
  split.2 <- split[[1]]
  split.3 <- split[[2]]
  
  #Merging all strings with linked appended
  merged.str <- paste(split.1,interp.1,interp.2,split.2,interp.3,interp.4,split.3,sep='')
  
  #Writting the file
  write(merged.str, file=html.path)
  
}

#Returns the file name/path with an ID appended behind the file name if it already exists.
generate.appropriate.file.name <- function(file.path){
  
  #file.path(string): path or file name if working in current directory. Should include the file extension.
  
  #Splitting arround slash
  directories.splits <- strsplit(file.path, '/')[[1]]
  
  #Contains only the file name part
  file.name.ext <- directories.splits[length(directories.splits)]
  file.name <- strsplit(file.name.ext, '[.]')[[1]][1]
  file.ext <- strsplit(file.name.ext, '[.]')[[1]][2]
  
  #Case with directory to prefix
  if (length(directories.splits) > 1)
    prefix.dir <- paste(directories.splits[1:(length(directories.splits)-1)], collapse = '/')
  
  #No prefix
  else
    prefix.dir <- NULL
  
  #Appending prefix if required
  file.attempt <- paste(file.name,file.ext,sep='.')
  if (!is.null(prefix.dir))
    file.attempt <- paste(prefix.dir,file.attempt,sep='/')
  
  #Index to append
  i <- 1
  
  #While the file name is not usable, increasing i value
  while(file.exists(file.attempt)){
    
    file.attempt <- paste(file.name,'-',i,'.',file.ext,sep='')
    if (!is.null(prefix.dir))
      file.attempt <- paste(prefix.dir,file.attempt,sep='/')    
    
    #Incrementing counter
    i <- i+1
  }
  
  return(file.attempt)
}

## ----------- Loading matrices ----------- 

#Loading the spots table
load.spots.table <- function(add.default.parents = FALSE){
  
  #add.default.parents(bool): appending the parent acronyms to the spots table
  
  spots.table <- read.table(paste(path.matrices,'spotstable.tsv',sep='/'),
                            sep='\t', header = TRUE, row.names = 1, quote="", stringsAsFactors = FALSE)
  if(add.default.parents)
    spots.table <- add.parent.acronym(spots.table)
  
  return(spots.table)
}

#Loading the expression matrix
load.expr.mat <- function(use.RData = TRUE){
  
  #use.RData(bool): use the RData format instead of parsing the .tsv matrix
  
  if(use.RData){
    load(paste(path.matrices,'exprmat.RData',sep='/'),.GlobalEnv)
    return()
  }else{
    st.data <- read.table(paste(path.matrices,'exprmat.tsv',sep='/'),
                          sep='\t', header = TRUE, row.names = 1, quote="", stringsAsFactors = FALSE)
    return(st.data)
  }
}

#Loading the slice table
load.slice.table <- function(){
  return(read.table(paste(path.matrices,'slicestable.tsv',sep='/'),
                    sep='\t', header = TRUE, row.names = 1, quote="", stringsAsFactors = FALSE))
}

#Loading the data frame with correpsondance cluster ID - names
load.df.cl.id.name <- function(cl.id.names.path){
  
  df.cl.name <- read.table(cl.id.names.path, 
                           header = T,
                           sep = ';',
                           stringsAsFactors = F)
  
  return(df.cl.name)
}
