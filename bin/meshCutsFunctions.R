#---------------3D GEOMETRY INTERSECTIONS ---------------

library(igraph)

#Computes the cross product
vector.cross.product.3d <- function(a,b){
  x <- a[2] * b[3] - a[3] * b[2]
  y <- a[3] * b[1] - a[1] * b[3]
  z <- a[1] * b[2] - a[2] * b[1]
  
  return(c(x,y,z))
}

#Returns the plane equation parameters N and d in a list
get.plane.equation <- function(coord, mode = 'coronal'){
  
  #Mode : saggital, coronal, horizontal: coord is a numerical
  #Mode : points, named list with P, Q, R triplets of points
  
  if(mode == 'saggital'){
    P <- c(-abs(coord), -3, -1)
    Q <- c(-abs(coord), -1, -2)
    R <- c(-abs(coord), -5, -3)
  }
  
  if(mode == 'coronal'){
    P <- c(-3, -1, coord)
    Q <- c(-1, -2, coord)
    R <- c(-5, -3, coord)
  }
  
  if(mode == 'horizontal'){
    P <- c(-3, -abs(coord), -1)
    Q <- c(-1, -abs(coord), -2)
    R <- c(-5, -abs(coord), -3)
  }
  
  if(mode == 'points'){
    P <- coord$P
    Q <- coord$Q
    R <- coord$R
  }
  
  #Get vectors
  PQ <- Q-P
  PR <- R-P
  
  #Get normal
  N <- vector.cross.product.3d(PQ,PR)
  
  #Normalize
  N <- N/sqrt(N[1]*N[1] + N[2]*N[2] + N[3]*N[3])
  
  #d
  d <- -sum(N*P)
  
  return(list(N = N, d = d))
}

#Signed distance from plane
dist.from.plane <- function(P, N, d){
  return (as.numeric(N %*% P + d))
}

#Get the intersection of a segment with a plane
get.segment.plane.intersection <- function(P1, P2, N, d, d1 = NULL, d2 = NULL, eps = 1e-100 ){
  
  tips <- NULL
  
  if(is.null(d1))
    d1 <- dist.from.plane(P1, N, d)
  
  if(is.null(d2))
    d2 <- dist.from.plane(P2, N, d)
  
  bool.P1.on.plane <- abs(d1) < eps
  bool.P2.on.plane <- abs(d2) < eps
  
  if (bool.P1.on.plane)
    tips <- rbind(tips, P1)
  if (bool.P2.on.plane)
    tips <- rbind(tips, P2)
  
  if (bool.P1.on.plane & bool.P2.on.plane)
    return(tips)
  
  #Same side of a planes
  if (d1 * d2 > eps)
    return(tips)
  
  t <- d1 / (d1 - d2)
  
  tips <- rbind(tips, (P1 + t * (P2 - P1)))
  return(tips)
  
}

#Get the intersection of a triangle with a plane
get.triangle.plane.intersection <- function(P1, P2, P3, N, d, d1 = NULL, d2 = NULL, d3 = NULL){
  
  tips <- rbind(get.segment.plane.intersection(P1,P2,N,d,d1,d2),
                get.segment.plane.intersection(P2,P3,N,d,d2,d3),
                get.segment.plane.intersection(P3,P1,N,d,d3,d1))
  
  return(unique(tips))
  
}

#Get the intersection between a list of meshes and a plane
get.meshes.list.plane.intersection <- function(mesh.list, smoothing.itt = 5, N, d){
  
  time.start <- Sys.time()
  
  list.intersect <- NULL
  
  for(k in 1:length(mesh.list)){
    
    list.intersect[[k]] <- list()
    
    time.diff <- as.numeric(round(difftime(Sys.time(),time.start,units='secs')))
    cat('\r                                                                                   ')
    cat(sprintf('\r%d/%d clusters treated. Ellapsed time: %d s', k, length(mesh.list), time.diff))
    
    mesh <- mesh.list[[k]]
    
    if(is.null(mesh))
      next
    
    mesh <- vcgSmooth(mesh, 'HC', smoothing.itt)
    
    all.dist <- dist.from.plane(mesh$vb[1:3,], N, d)
    
    for(i in 1:dim(mesh$it)[2]){
      
      P1 <- mesh$vb[1:3,mesh$it[1,i]]
      P2 <- mesh$vb[1:3,mesh$it[2,i]]
      P3 <- mesh$vb[1:3,mesh$it[3,i]]
      
      d1 <- all.dist[mesh$it[1,i]]
      d2 <- all.dist[mesh$it[2,i]]
      d3 <- all.dist[mesh$it[3,i]]
      
      inter.coord <- get.triangle.plane.intersection(P1,P2,P3,N,d,d1,d2,d3)
      list.intersect[[k]][[i]] <- inter.coord
    }
  }
  return(list.intersect)
}

#Plot the intersection lines as segments
plot.segments.intersection <- function(list.intersect, list.colors, mode = 'coronal', fname = 'segments-plots.pdf'){
  
  l.lim <- get.lim.mode(mode)
  
  pdf(generate.appropriate.file.name(fname), useDingbats = F)
  plot(1, type = 'n', xlim = l.lim$xlim, ylim = l.lim$ylim, asp = 1)
  
  for(k in 1:length(list.intersect)){
    
    col <- list.colors[k]
    
    if(length(list.intersect[[k]]) == 0)
      next
    
    for(i in 1:length(list.intersect[[k]])){
      
      inter.coord <- list.intersect[[k]][[i]]
      
      if(!is.null(inter.coord)){
        if(dim(inter.coord)[1] == 2){
          if(mode == 'coronal'){
            x0 <- inter.coord[1,1]
            y0 <- inter.coord[1,2]
            x1 <- inter.coord[2,1]
            y1 <- inter.coord[2,2]
          }else if(mode == 'saggital'){
            x0 <- inter.coord[1,3]
            y0 <- inter.coord[1,2]
            x1 <- inter.coord[2,3]
            y1 <- inter.coord[2,2]          
          }else if(mode == 'horizontal'){
            x0 <- inter.coord[1,3]
            y0 <- inter.coord[1,1]
            x1 <- inter.coord[2,3]
            y1 <- inter.coord[2,1]          
          }
          segments(x0,y0,x1,y1,
                   lwd = 0.05,
                   col = col)
        }
      }
    }
  }
  dev.off()
}
  
#Returns appropriate limits for certain plotting mode
get.lim.mode <- function(mode){
  if(mode == 'coronal')
    return(list(xlim = c(-6,6), ylim = c(-8,1)))    
  if(mode == 'saggital')
    return(list(xlim = c(-8,5), ylim = c(-8,1)))
  if(mode == 'horizontal')
    return(list(xlim = c(-8,5), ylim = c(-6,6)))
}

#--------------- CONVERTING INTERSECTION SEGMENTS INTO A POLYGON ---------------

#Plot the polygonal cuts
get.polygon.cut <- function(list.intersect, mode = 'coronal', eps.equality = 1e-10){
  
  #Columns selected
  if (mode == 'coronal'){
    cols.to.select <- c(1,2)
  }else if (mode == 'saggital'){
    cols.to.select <- c(3,2)
  }else if (mode == 'horizontal'){
    cols.to.select <- c(3,1)
  }
  
  polygon.list <- list()
  
  #Looping through clusters
  for (j in 1:length(list.intersect)){
    
    polygon.list[[j]] <- list()
    
    cur.intersect <- list.intersect[[j]]
    
    #If no intersections with that cluster: next
    if (length(cur.intersect) == 0)
      next
    
    #Intializing list of segments
    seg.list <- NULL
    
    #Looping through each mesh intersection
    for(i in 1:length(cur.intersect)){
      
      if(!is.null(cur.intersect[[i]])){
        
        #Kepping only the proper segments, not point
        if(dim(cur.intersect[[i]])[1] == 2)
          seg.list <- rbind(seg.list, cbind(cur.intersect[[i]], i))
      }
    }
    
    #No segment case
    if(is.null(seg.list))
      next

    #Getting unique points 
    points <- unique(seg.list[,cols.to.select])
    
    #Checking for unicity of points
    while(TRUE){
      
      #Pairwise distance between points
      d <- as.matrix(dist(points))
      
      #Finding points that are unique 
      is.unique <- apply(d, 1, function(x){return(sum(x < eps.equality))})
      
      #If they are all unique: continue
      if(sum(is.unique != 1) == 0){
        break
        
        #Deleting a non-unique line 
      }else{
        index.to.keep <- setdiff(1:dim(points)[1],which(is.unique != 1)[1])
        points <- points[index.to.keep,]
      }
    }
    
    #Unique ID based on the points linked
    id.seg <- apply(seg.list, 1, function(x){return(which(abs(points[,1] - x[cols.to.select[1]]) < eps.equality & abs(points[,2] - x[cols.to.select[2]]) < eps.equality))})
    
    #Getting the edge
    edges.mat <- matrix(id.seg, ncol = 2, byrow = T)
    
    #Creating the graph
    graph <- graph.edgelist(edges.mat, directed = FALSE)
    
    #Finding individual groups
    graph.groups <- igraph::groups(components(graph))
    
    #Looking at chains in each group
    for(index in 1:length(graph.groups)){
      
      #Subgraph with current group
      g <- graph.groups[[index]]
      sub.graph <- induced_subgraph(graph, g)
      
      #Finding the order for the circle
      points.order <- g[girth(sub.graph)$circle]
      
      #Case with no loop
      if(length(points.order) == 0){
        
        #Cannot have an interesting graph with less than 3 points anyway
        if(length(graph.groups[[index]]) < 3)
          next
        
        degree.equal.1 <- which(degree(sub.graph) == 1)
        
        if(length(degree.equal.1) < 2)
          next
      
        #Adding edges if required
        sub.graph <- add_edges(sub.graph, V(sub.graph)[degree.equal.1[1:2]])
        
        #Reselecting the points
        points.order <- g[girth(sub.graph)$circle]
        
        #Case where it did not solve everything
        if(length(points.order) == 0)
          warning(paste('No closed circle after adding edge. (cluster = ', j, ', index = ', index, ', degree = ' ,sum(degree(sub.graph) != 2), ')', sep = ''))
        
      }
      polygon.list[[j]][[index]] <- points[points.order,] 
    }
    
    
    #Find contour in another contour
    # updated <- TRUE
    # while(updated){
    # 
    #   updated <- FALSE
    # 
    #   if(length(polygon.list[[j]]) < 2)
    #     break
    # 
    #   for(i1 in 1:(length(polygon.list[[j]])-1)){
    # 
    #     if(i1 > length(polygon.list[[j]]))
    #       next
    # 
    #     if(is.null(polygon.list[[j]][[i1]]))
    #       next
    # 
    #     for(i2 in setdiff(1:length(polygon.list[[j]]), i1)){
    # 
    #       if(i2 > length(polygon.list[[j]]))
    #         next
    # 
    #       if(is.null(polygon.list[[j]][[i2]]))
    #         next
    # 
    #       is.inside <- point.in.polygon(polygon.list[[j]][[i2]][,1],
    #                                     polygon.list[[j]][[i2]][,2],
    #                                     polygon.list[[j]][[i1]][,1],
    #                                     polygon.list[[j]][[i1]][,2])
    # 
    #       if(sum(is.inside)/dim(polygon.list[[j]][[i2]])[1] > 0.3){
    #         polygon.list[[j]][[i1]] <- rbind(polygon.list[[j]][[i1]] ,
    #                                          NA,
    #                                          polygon.list[[j]][[i2]])
    #         polygon.list[[j]][[i2]] <- NULL
    #         print('updated')
    #         updated <- TRUE
    #         next
    #       }
    #     }
    #   }
    # }
  }
  return(polygon.list)
}

#Plot the polygonal cut
plot.polygon.cut <- function(polygon.list, list.colors, l.lim = get.lim.mode('coronal'), fname = 'mesh-cut.pdf', main = ''){
 
  #Opening pdf 
  if(!is.null(fname))
    pdf(generate.appropriate.file.name(fname), useDingbats = F)
  
  #Initaliazing plot
  plot(1, type = 'n', xlim = l.lim$xlim, ylim = l.lim$ylim, asp = 1, main = main)
  
  #Looping through polygons
  for(j in 1:length(polygon.list)){
    
    #Handling colors
    col <- list.colors[j]
    
    #If no polygon here, next
    if(length(polygon.list[[j]]) == 0)
      next
    
    #Showing every contour
    for (id in 1:length(polygon.list[[j]])){
      
      polygon(polygon.list[[j]][[id]][, 1],
              polygon.list[[j]][[id]][, 2],
              col = col,
              lwd = 0.1,
              fillOddEven = T)
      
    }
  }
  
  if(!is.null(fname))
    dev.off()
  
}

#Returns the coronal outline from the brain based on meshes
get.coronal.outline.from.allen.meshes <- function(AP){
  
  l.plane <- get.plane.equation(AP, 'coronal')
  N <- l.plane$N
  d <- l.plane$d
  
  allen.outline.mesh <- mesh3d.allen.annot.from.id(get.id.from.acronym('root'))
  
  list.intersect <- get.meshes.list.plane.intersection(list(allen.outline.mesh), 0, N, d)
  
  outline <- get.polygon.cut(list.intersect, 'coronal')
  
  return(outline[[1]])
  
}
