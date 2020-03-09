get.SI <- function(spots.table,
                   partition.1.field,
                   partition.2.field){
  
  sp.i <- subset(spots.table,TRUE,partition.1.field)
  sp.j <- subset(spots.table,TRUE,partition.2.field)
  
  sp.i[,1] <- as.character(sp.i[,1])
  sp.j[,1] <- as.character(sp.j[,1])
  colnames(sp.i) <- 'cl'
  colnames(sp.j) <- 'cl'
  
  r.i <- sapply(sort(unique(sp.i[,1])),function(x){return(rownames(subset(sp.i, cl == x)))})
  r.j <- sapply(sort(unique(sp.j[,1])),function(x){return(rownames(subset(sp.j, cl == x)))})
  
  P.ij <- matrix(nrow = length(r.i),
                 ncol = length(r.j))
  rownames(P.ij) <- names(r.i)
  colnames(P.ij) <- names(r.j)

  for(i in rownames(P.ij)){
    for(j in colnames(P.ij)){
      P.ij[i,j] <- length(intersect(r.i[[i]],r.j[[j]])) / length(r.j[[j]])
    }
  }
  
  P.ji <- matrix(nrow = length(r.j),
                 ncol = length(r.i))
  rownames(P.ji) <- names(r.j)
  colnames(P.ji) <- names(r.i)
  
  for(i in colnames(P.ji)){
    for(j in rownames(P.ji)){
      P.ji[j,i] <- length(intersect(r.j[[j]],r.i[[i]])) / length(r.i[[i]])
    }
  }
  
  X.ij <- matrix(nrow = length(r.i),
              ncol = length(r.j))
  rownames(X.ij) <- names(r.i)
  colnames(X.ij) <- names(r.j)
  
  for(i in rownames(X.ij)){
    for(j in colnames(X.ij)){
      X.ij[i,j] <- max(P.ij[i,j],P.ji[j,i])
    }
  }
  
  U.ij <- matrix(nrow = length(r.i),
                 ncol = length(r.j))
  rownames(U.ij) <- names(r.i)
  colnames(U.ij) <- names(r.j)
  for(i in rownames(U.ij)){
    for(j in colnames(U.ij)){
      if(X.ij[i,j] > 0)
        U.ij[i,j] <- min(length(r.i[[i]]), length(r.j[[j]]))
      else
        U.ij[i,j] <- 0
    }
  }  
  
  W.ij <- U.ij / sum(U.ij)
  S <- 1 - 4*sum(W.ij*X.ij*(1-X.ij))
  return(S)
}

get.S.value.per.clusters <- function(spots.table, clusters.path = NULL, parent.list = c('VS','TH','STR','RHP','P','PAL','OLF','MY','MB','HY','HIP','fiber tracts','CTX','CB'),
                                     discard.non.clustered = TRUE, min.cluster.size = 0, deep.classification = FALSE, force.acronym = TRUE){
  
  if (deep.classification){
    load(paste(path.bin,'deepClassification.RData',sep='/'))
    parent.list <- brain
  }
  
  if (is.null(clusters.path))
    clusters.path <- '0'
  
  S.vect <- as.data.frame(matrix(0,length(clusters.path),2))
  colnames(S.vect) <- c('N.clusters','S')
  
  if (force.acronym){
    spots.table$full.name.parent <- NULL
    spots.table$acronym.parent <- NULL
    spots.table <- add.parent.acronym(spots.table,c(parent.list,'root'))
    spots.table$full.name.parent <- name.from.acronym(spots.table$acronym.parent)
  }
  
  for (idx in 1:length(clusters.path)){
    if(clusters.path[[1]] != '0'){
      spots.table.cur <- append.cluster.to.spots.table(spots.table,clusters.path[[idx]], min.cluster.size = min.cluster.size)
    }
    
    spots.table.cur <- spots.table.cur[spots.table.cur$acronym.parent != 'root',]
    
    if(discard.non.clustered){
      spots.table.cur <- spots.table.cur[!(spots.table.cur$cluster == -1),]
      spots.table.cur$cluster <- factor(spots.table.cur$cluster)
    }
    
    browser()
    
    r.anat <- lapply(parent.list,function(x){return(rownames(spots.table.cur[spots.table.cur$acronym.parent == x,]))})
    names(r.anat) <- parent.list
    
    r.clust <- lapply(as.character(levels(spots.table.cur$cluster)),function(x){return(rownames(spots.table.cur[spots.table.cur$cluster == x,]))})
    names(r.clust) <- as.character(levels(spots.table.cur$cluster))
    
    P.region.in.cluster <- matrix(0,length(r.anat),length(r.clust))
    rownames(P.region.in.cluster) <- names(r.anat)
    colnames(P.region.in.cluster) <- names(r.clust)
    
    P.cluster.in.region <- P.region.in.cluster
    
    for (i in 1:length(r.anat)){
      for (j in 1:length(r.clust)){
        inter <- intersect(r.anat[[i]],r.clust[[j]])
        P.cluster.in.region[i,j] <- length(inter) / length(r.clust[[j]])
        P.region.in.cluster[i,j] <- length(inter)/ length(r.anat[[i]])
      }
    }
    
    X <- matrix(0,length(r.anat),length(r.clust))
    colnames(X) <- colnames(P.cluster.in.region)
    rownames(X) <- rownames(P.cluster.in.region)
    for (i in 1:dim(X)[1]){
      for (j in 1:dim(X)[2]){
        X[i,j] <- max(P.cluster.in.region[i,j],P.region.in.cluster[i,j],na.rm = TRUE)
      }
    }
    
    U <- matrix(0,length(r.anat),length(r.clust))
    colnames(U) <- colnames(P.cluster.in.region)
    rownames(U) <- rownames(P.cluster.in.region)
    for (i in 1:dim(U)[1]){
      for (j in 1:dim(U)[2]){
        if(X[i,j] > 0)
          U[i,j] <- min(length(r.anat[[i]]),length(r.clust[[j]]))
        else
          U[i,j] <- 0
      }
    }
    
    W <- matrix(0,length(r.anat),length(r.clust))
    colnames(W) <- colnames(P.cluster.in.region)
    rownames(W) <- rownames(P.cluster.in.region)
    for (i in 1:dim(W)[1]){
      for (j in 1:dim(W)[2]){
        W[i,j] <- U[i,j]/sum(U)
      }
    }
    
    S.vect[idx,'S'] <- 1 - 4*sum(W*X*(1 - X))
    S.vect[idx,'N.clusters'] <- length(levels(spots.table.cur$cluster))
  }
  
  return(S.vect)
  
}



plot.overlap.matrices <- function(cluster.path, spots.table, add.names = TRUE, reorder = TRUE){
  
  library(lattice)
  
  t <- read.table(cluster.path, stringsAsFactors = FALSE)
  spots.table.cur <- spots.table[t[,1],]
  spots.table.cur$cluster <- NULL
  spots.table.cur$clusters.named <- NULL
  spots.table.cur[t[,1],'cluster'] <- as.factor(t[,2])
  
  spots.table.cur <- spots.table.cur[spots.table.cur$acronym.parent != 'root',]
  spots.table.cur <- spots.table.cur[!is.na(spots.table.cur$cluster),]
  
  r.anat <- lapply(parent.list,function(x){return(rownames(spots.table.cur[spots.table.cur$acronym.parent == x,]))})
  names(r.anat) <- parent.list
  
  if (add.names){
    sp <- name.clusters(spots.table.cur)
    r.clust <- lapply(as.character(levels(sp$clusters.named)),function(x){return(rownames(sp[sp$clusters.named == x,]))})
    names(r.clust) <- as.character(levels(sp$clusters.named))    
  }else{
    r.clust <- lapply(as.character(levels(spots.table.cur$cluster)),function(x){return(rownames(spots.table.cur[spots.table.cur$cluster == x,]))})
    names(r.clust) <- as.character(levels(spots.table.cur$cluster))
  }
  P.region.in.cluster <- matrix(0,length(r.anat),length(r.clust))
  rownames(P.region.in.cluster) <- names(r.anat)
  colnames(P.region.in.cluster) <- names(r.clust)
  
  P.cluster.in.region <- P.region.in.cluster
  
  print(length(r.clust))
  
  for (i in 1:length(r.anat)){
    for (j in 1:length(r.clust)){
      inter <- intersect(r.anat[[i]],r.clust[[j]])
      P.cluster.in.region[i,j] <- length(inter) / length(r.clust[[j]])
      P.region.in.cluster[i,j] <- length(inter)/ length(r.anat[[i]])
    }
  }
  rgb.palette = brewer.pal(9, "YlOrRd")
  rgb.palette <- rev(rgb.palette)
  rgb.palette <- c('#000000',rgb.palette,'#ffffff')
  
  if (reorder){
    ord.vect <- rep('',dim(P.cluster.in.region)[2])
    main.clust <- apply(P.cluster.in.region, 2, which.max)
    curs <- 1
    for(idx in sort(unique(main.clust),decreasing = TRUE)){
      sub.sel <- P.cluster.in.region[idx,main.clust == idx]
      if (length(sub.sel) == 1)
        ord.vect[curs] <- names(which(main.clust == idx))
      else if(length(sub.sel) > 0)
        ord.vect[curs:(curs+length(sub.sel)-1)] <- names(sort(sub.sel))
      curs <- curs + length(sub.sel)
    }
    P.cluster.in.region <- P.cluster.in.region[,ord.vect]
    
    ord.vect <- rep('',dim(P.region.in.cluster)[1])
    main.clust <- apply(P.region.in.cluster, 2, which.max)
    curs <- 1
    for(idx in sort(unique(main.clust),decreasing = TRUE)){
      sub.sel <- P.region.in.cluster[idx,main.clust == idx]
      if (length(sub.sel) == 1)
        ord.vect[curs] <- names(which(main.clust == idx))
      else if(length(sub.sel) > 0)
        ord.vect[curs:(curs+length(sub.sel)-1)] <- names(sort(sub.sel))
      curs <- curs + length(sub.sel)
    }
    P.region.in.cluster <- P.region.in.cluster[,ord.vect]
  }
  win.graph()
  plot.new()
  p1 <- levelplot(P.region.in.cluster,col.regions=colorRampPalette(rgb.palette)(500),cuts=200, at=seq(0,1,0.005),xlab='Anatomical region',ylab = 'Cluster', main = 'P(cluster | region)', scales=list(x=list(cex=.6, rot = 90),y=list(cex=.6)))
  p2 <- levelplot(P.cluster.in.region,col.regions=colorRampPalette(rgb.palette)(500),cuts=200, at=seq(0,1,0.005),xlab='Anatomical region',ylab = 'Cluster', main = 'P(region | cluster)', scales=list(x=list(cex=.6, rot = 90),y=list(cex=.6)))
  library(gridExtra)
  grid.arrange(p1,p2,ncol=2)
}