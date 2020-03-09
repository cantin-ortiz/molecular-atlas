library(ape)

#------------------- ACCESSING DATA -------------------

#Return the spots/IC matrix from a seurat object
get.ic.mat <- function(seur.obj, field.name = 'ica'){
  
  #seur.obj: seurat object with IC run
  #field.name(char): name of the dimensionality reduction field 
  return(seur.obj@dr[[field.name]]@cell.embeddings)
}

#------------------- BASIC PLOTS -------------------

#Plot the top genes of an IC
plot.ic.top.genes <- function(seur.obj, ic.name = 'IC1', field.name = 'ica', n.top = 30, show.ylabs = TRUE){
  
  #seur.obj
  #ic.name(char): name of the IC to plot
  #field.name: name of the dr
  #n.top: number of top genes to plot
  
  #Selecting proper data
  data <- seur.obj@dr[[field.name]]@gene.loadings[,ic.name]
  
  #Keeping only top genes (and sorting them)
  top.genes <- data[order(abs(data), decreasing = TRUE)]
  top.genes <- sort(top.genes[1:n.top])
  
  #Defining xlim
  diff.lim <- max(top.genes) - min(top.genes)
  xlim = c(min(top.genes) - 0.4*diff.lim, max(top.genes) + 0.05*diff.lim)
  left.axis <- min(top.genes) - 0.05*diff.lim

  #Potential step roundings.
  desired.step.roundings <- c(100,250,500,1000,2000,2500,5000,10000)
    
  #Deciding x.ticks: case with different signs
  if((left.axis * xlim[2]) < 0){
    
    #Chosing reference point
    ref <- max(abs(left.axis), xlim[2])
    
    #Desired number of ticks on that side (~10 in total)
    n.desired.ticks <- (round(ref/diff.lim*10))
    
    #Finding the closest desired step
    perfect.step <- ref/n.desired.ticks
    step <- desired.step.roundings[which.min(abs(desired.step.roundings-perfect.step))]
    # print(paste(perfect.step,step,sep='/'))

    #First  the negative ticks if relevant
    if (-step > left.axis)
      left.ticks <- rev(seq(from = 0, to = left.axis, by = -step))
    else
      left.ticks <- NULL
    
    #Positive ticks
    if (step < xlim[2])
      right.ticks <- seq(from = 0, to = xlim[2], by = step)
    else
      right.ticks <- NULL
    
    #Merging together and discarding the 0 if present twice
    if(is.element(0,left.ticks) & is.element(0,right.ticks))
      right.ticks <- sort(setdiff(right.ticks,0))

    ticks <- c(left.ticks,right.ticks)
  
  #Case with only one sign for values
  }else{
    #Finding the closest desired step
    perfect.step <- diff.lim/10
    step <- desired.step.roundings[which.min(abs(desired.step.roundings-perfect.step))]
    start.val <- ceiling(left.axis/step)*step
    ticks <- seq(from = start.val, to = xlim[2], by = step)
  }
  
  #Plotting the top genes
  plot(top.genes, 1:n.top, xaxt = 'n', yaxt = 'n', ylim = c(-5,n.top+1), xlim = xlim, xlab = '', ylab = '', bty = 'n', pch = 19)
  par(xpd = TRUE)
  axis(1, pos = 0, at = c(left.axis,xlim[2]), labels = NA, lwd.ticks = 0)
  axis(2, pos = left.axis, at = c(0,n.top+1), labels = NA, lwd.ticks = 0)
  axis(3, pos = n.top+1, at = c(left.axis,xlim[2]), labels = NA, lwd.ticks = 0)
  axis(4, pos = xlim[2], at = c(0,n.top+1), labels = NA, lwd.ticks = 0)

  axis(1, pos = 0, at = ticks)
  
  if(show.ylabs){
    axis(2, pos = left.axis, at = 1:n.top, labels = names(top.genes), las = 2, cex.axis = 0.65)
  }else{
    text(left.axis, n.top/2, 'Genes', srt = 0, pos = 2)
  }

  #0 line
  segments(0,0,0,n.top+1,lty=2)
  
  #xlab
  text(mean(c(left.axis,xlim[2])),-1000,'Gene loading')
}

#Plot in the atlas, with a "double plate" layout, the expression of an IC
plot.2d.atlas.ic <- function(spots.table, ic.mat, AP = 0, ic.name = 'IC1', cex = 1, lwd = 1){
  
  #spots.table
  #ic.mat
  #AP(num): closest AP desired for plotting
  #ic.name(char): name of the IC to plot
  #cex: cex for spots plot
  #lwd: lwd for spots plot
  
  #Finding closest matching AP
  AP.list <- unique(spots.table$AP)
  AP.sel <- AP.list[which.min(abs(AP.list-AP))]
  
  #Selecting the proper spots
  sp <- spots.table[spots.table$AP == AP.sel,]
  
  #Only keeping shared spots between ic.mat and spots.table
  shared.spots <- intersect(rownames(sp), rownames(ic.mat))
  sp <- sp[shared.spots,]

  #Vector of ic expression
  ic.vect <- ic.mat[,ic.name]
  sp$ic.load <- ic.vect[rownames(sp)]
    
  #Obtaining the proper color for each spot
  bin.sep <- seq(from = -max(abs(ic.vect))-0.1, to = max(abs(ic.vect))+0.1, length.out = 202)
  sp$ic.bin <- .bincode(sp$ic.load, breaks = bin.sep)
  # function.cp <- colorRampPalette(brewer.pal(11, 'RdBu'))
  # cp <- rev(function.cp(201))
  function.cp <- colorRampPalette(c('#0004ff', '#ffffff', '#ff0000'))
  cp <- function.cp(201)
  sp$color <- cp[sp$ic.bin]
  
  #Reordering so strongest expression is plot on top
  sp <- sp[order(abs(sp$ic.bin - 101)),]

  #Plotting  
  .plot.2d.atlas.double.plate.col(sp, cex, lwd)
  
}

#Plot in 3d with glassbrain the expression of an IC
plot.3d.glassbrain.ic <- function(spots.table, ic.mat, ic.name = 'IC1', cex = 1.5, top.quantile = 0.9, dim = c(0, 0, 720,1080), HD = FALSE){
  
  #spots.table
  #ic.mat
  #ic.name(char): name of the IC to plot
  #cex: cex for spots plot
  #top.quantile: only spots in that quantile of the IC will be plotted
  #dim(numeric): dimension of the window
  #HD(logical): use HD contour
  
  #Only keeping shared spots between mat.ic and spots.table
  shared.spots <- intersect(rownames(spots.table), rownames(ic.mat))
  spots.table <- spots.table[shared.spots,]
  ic.mat <- ic.mat[shared.spots,]
  
  #Selecting corresponding component
  vect.ic <- ic.mat[,ic.name]
  
  #Adding the expression
  spots.table[rownames(ic.mat), 'ic.load'] <- vect.ic
  
  #Binning the expression level
  bin.sep <- seq(from = -max(abs(vect.ic)), to = max(abs(vect.ic)), length.out = 202)
  bin.sep[1] <- bin.sep[1] - 0.1
  bin.sep[length(bin.sep)] <- bin.sep[length(bin.sep)] + 0.1
  spots.table$ic.bin <- .bincode(spots.table$ic.load, breaks = bin.sep)
  
  #Getting the color
  function.cp <- colorRampPalette(c('#0004ff', '#bfbfbf', '#ff0000'))
  cp <- function.cp(201)
  spots.table$color <- cp[spots.table$ic.bin]
  
  #Keeping only the proper quantile
  quant <- quantile(abs(vect.ic), top.quantile)
  spots.table <- spots.table[(spots.table$ic.load >= quant) | (spots.table$ic.load <= -quant),]

  #Plotting in 3D
  .plot.3d.glassbrain.col(spots.table, cex = cex, HD = HD, dim = dim)
}

# #Plot a heatmap of the expression level of clusters on every IC
# ploc.ic.clusters.heatmap <- function(spots.table, mat.ic.cl.avg){
#   
# }

#------------------- MATCHING IC BETWEEN DIFFERENT DIMENSIONALITY REDUCTION -------------------

#Returns the correlation matrix of expression between two IC matrix
get.ic.correlation <- function(mat.ic.1, mat.ic.2, names = c('m1.ic.','m2.ic.')){

  #mat.ic.1: first spots/IC matrix (as rows in output)
  #mat.ic.2: second spots/IC matrix (as cols in output)
  #names: prefix of each IC vector name, first one corresponds to mat.ic.1 (rows) and second to mat.ic.2 (cols)
  
  #Taking only the shared spots between the two matrices
  m1 <- mat.ic.1[intersect(rownames(mat.ic.1), rownames(mat.ic.2)),]
  m2 <- mat.ic.2[intersect(rownames(mat.ic.1), rownames(mat.ic.2)),]
  
  #Will contain the results
  df.result <- data.frame(matrix(nrow = dim(m1)[2], ncol = dim(m1)[2]))
  rownames(df.result) <- paste(names[1], 1:80, sep = '.')
  colnames(df.result) <- paste(names[2], 1:80, sep = '.')
  
  #Computing CC for all pairs of IC
  for (i in 1:(dim(m1)[2])){

    for (j in 1:dim(m2)[2]){
      cc <- ccf(m1[,i], m2[,j], plot = FALSE)
      lag.0 <- which(cc$lag == 0)
      #Selecting only the 0 lag (aligned spots)
      df.result[i,j] <- cc$acf[lag.0]
    }
  }
  
  return(df.result)
}

#For each row, returns the col that best matches (in absolute value) from a correlation matrix, i.e closest to 1/-1
get.cor.mat.best.match <- function(correlation.matrix){
  
  #correlation.matrix: correlation.matrix from get.ic.correlation
  df.summary <- data.frame(matching.ic = numeric(dim(correlation.matrix)[1]),
                           cross.corr = numeric(dim(correlation.matrix)[1]))
  rownames(df.summary) <- rownames(correlation.matrix)
  df.summary$matching.ic <- colnames(correlation.matrix)[apply(abs(correlation.matrix), 1, which.max)]
  df.summary$cross.corr <- apply(abs(correlation.matrix), 1, max)
  df.summary <- df.summary[order(df.summary$cross.corr, decreasing = T),]
  return(df.summary)
}

#For each IC vectors in mat.ic.1, finds the best matching vector in mat.ic.2 (closest to 1 in absolute value)
get.ic.mat.best.batch <- function(mat.ic.1, mat.ic.2, names = c('m1.ic.','m2.ic.')){
  
  #mat.ic.1: first spots/IC matrix (as rows in output)
  #mat.ic.2: second spots/IC matrix (as cols in output)
  #names: prefix of each IC vector name
  
  corr.mat <- get.ic.correlation(mat.ic.1, mat.ic.2, names)
  best.match <- get.cor.mat.best.match(corr.mat)
  
  return(best.match)
}

#------------------- AVERAGE EXPRESSION OF IC IN EVERY CLUSTER -------------------

#Returns a matrix with the average expression of every cluster (rows) on every IC (cols)
get.ic.cluster.average <- function(spots.table, mat.ic){
  
  #spots.table
  #mat.ic: IC matrix from get.ic.mat
  
  #Keeping only common spots
  shared.spots <- intersect(rownames(spots.table), rownames(mat.ic))
  spots.table <- spots.table[shared.spots,]
  mat.ic <- mat.ic[shared.spots,]
  
  #Selecting unique clusters
  cl.list <- unique(spots.table$clusters.named)
  
  #Initializing matrix
  mat.cl.ic <- matrix(nrow = length(cl.list), ncol = dim(mat.ic)[2])
  colnames(mat.cl.ic) <- colnames(mat.ic)
  rownames(mat.cl.ic) <- cl.list
  
  #Looping through clusters
  for (cl in cl.list){
    
    #Selecting relevant spots
    sp.cl <- rownames(spots.table[spots.table$clusters.named == cl,])
    
    #Extracting relevant IC matrix and getting mean value
    sub.ica.mat <- mat.ic[sp.cl,]
    m <- colMeans(sub.ica.mat)
    
    #Saving
    mat.cl.ic[cl,names(m)] <- m
  }
  
  #Returning matrix
  return(mat.cl.ic)
}

#Returns the matrix of euclidian distance between clusters based on average IC expression
get.ic.cluster.dist <- function(mat.ic.cl.avg){
  
  #mat.ic.cl.avg: matrix of IC/cluster average
  
  #List of clusters
  cl.list <- rownames(mat.ic.cl.avg)

  #Initializing matrix
  mat.ic.cl.dist <- matrix(nrow = length(cl.list), ncol = length(cl.list))
  colnames(mat.ic.cl.dist) <- cl.list
  rownames(mat.ic.cl.dist) <- cl.list

  #Looping through every pairs (triangular matrix is enough because of symmetry)
  for (i in 1:(length(cl.list)-1)){
    for (j in (i+1):length(cl.list)){
      
      m1 <- mat.ic.cl.avg[i,]
      m2 <- mat.ic.cl.avg[j,]
      # cc <- ccf(m1,m2, plot = FALSE)
      # val <- cc$acf[cc$lag == 0]
      
      val <- sqrt(sum((m1-m2)*(m1-m2)))
      
      mat.ic.cl.dist[i,j] <- val
      mat.ic.cl.dist[j,i] <- val
    }
  }
  
  #Replacing NA values by 0 (diagonal)
  mat.ic.cl.dist[is.na(mat.ic.cl.dist)] <- 0
  
  #Returning  matrix
  return(mat.ic.cl.dist)
}

#Returns a tree object obtained by running hierarchical clustering
get.ic.hiearchical.tree <- function(mat.ic.cl.dist, spots.table, method = 'ward.D2'){
  
  #mat.ic.cl.dist: matrix of IC/cluster distance based on averaged IC expression
  #spots.table
  #method(char): a method valid with hclust
  
  cl.count <- count(spots.table$clusters.named)
  cl.count$x <- as.character(cl.count$x)
  rownames(cl.count) <- cl.count$x
  
  cl.tree <- hclust(as.dist(mat.ic.cl.dist), 
                    method = method, 
                    members = cl.count[rownames(mat.ic.cl.dist), 'freq'])
  
  return(cl.tree)
}

#Plot the contribution of each IC for a specific cluster (based on average value per spot)
plot.ic.cluster.average <- function(mat.ic.cl.avg, file.path = 'main-IC.pdf', cl.list = NULL, spots.table = NULL){
  
  #mat.ic.cl.avg: matrix of IC/cluster average
  #file.path(string): path or file name, with extension. If null, plots are not saved in a file.
  #cl.list(string): list of cluster to plot. If null, all clusters are plotted
  #spots.table: if specified, adds the number of spot in each cluster in the title of the plot
  
  #If saving the plot in a file
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path) 
  }
    
  #Default value for cluster list if it is null
  if (is.null(cl.list))
    cl.list <- sort(rownames(mat.ic.cl.avg))
  
  #Looping through clusters of interest
  for (cl in cl.list){
    
    #Row of interest
    m <- mat.ic.cl.avg.all[cl,]
    
    #Ordering values
    l <- sort(m, decreasing = FALSE)
    
    #5 genes with highest contribution (absolute value)
    top.genes <- l[order(abs(l), decreasing = TRUE)][1:5]
    
    #Definning xlim and ylim
    maxx <- max(l)
    minx <- min(l)
    xdiff <- maxx-minx
    xlim <- c(minx-0.05*xdiff,maxx+0.05*xdiff)
    n.ic <- length(m)
    ylim <- c(0, n.ic+1)
    
    #Sequence of vertical lines
    x.lines <- seq(ceiling(xlim[2]*2)/2, floor(xlim[1]*2)/2, by = -0.5)
    x.lines <- setdiff(x.lines,0)
    
    #Initializing title
    title <- cl
    
    #If spots table is specified, adding the number of spots in the title.
    if (!is.null(spots.table))
      title <- paste(title,'\n(',sum(spots.table$clusters.named == cl), ' spots)',sep='')
    
    #Plotting
    plot(1,
         type = 'n',
         main = title,
         xlim = xlim,
         ylim = ylim,
         xlab = 'IC load',
         ylab = '',
         yaxt = 'n')
    
    segments(x.lines,rep(1000,length(x.lines)),x.lines,rep(-1000,length(x.lines)), lty = 1, 
             col = 'lightgray', lwd = 0.5, lend = 'butt')
    segments(0,1000,0,-1000, lty = 1, col = 'lightgray', lwd = 2.5, lend = 'butt')
    
    points(l,
           1:length(l),
           col = 'darkgray',
           pch = 16)
    points(top.genes,
           sapply(names(top.genes), function(x){return(which(names(l) == x))}),
           pch = 16,
           col = 'red')
    text(l, 1:length(l), names(l), pos = 4, cex = 0.5)
  }
  
  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
}

#Plot the clustering tree of clusters using the euclidian distance between IC as dissimilarity measurement
plot.ic.cluster.tree <- function(mat.ic.cl.dist, 
                                 spots.table, 
                                 method = 'average', 
                                 file.path = 'hierarchical-cluster.pdf', 
                                 ncolors = 50,
                                 plot.type = 'fan',
                                 width = 100, 
                                 height = 100){

  #mat.ic.cl.dist: matrix of IC/cluster distance based on averaged IC expression
  #spots.table: used for normalizing by number of spots per cluster
  #method(string): a method accepted by hclust algorithm
  #file.path(string): path or file name, with extension. If null, plots are not saved in a file.
  #ncolors(numeric): if #0, then the tree is colored with that many colors based on branching
  #plot.type: If string, must be fan, phylogram or unrooted. Coded in respective order by numerical values 1 2 3
  #width(numerical): width of the pdf
  #height(numerical): height of the pdf
  
  library(ape)
    
  #If saving the plot in a file
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path, width = width, height = height) 
  }
  
  #Converting numerical values to usable strings
  if(is.element(plot.type, 1:3))
    plot.type <- c('fan', 'phylogram', 'unrooted')[plot.type]
  
  #Size of clusters
  cl.count <- count(spots.table$clusters.named)
  
  #Comverting into a distance matrix
  # mat.dist <- 1 - abs(mat.ic.cl.dist)
  mat.dist <- mat.ic.cl.dist
  
  #Computing tree
  cl.tree <- hclust(as.dist(mat.dist), method = method, members = cl.count$freq)
  
  #Coloring tree if required
  if(ncolors > 0){
    clus.col <- cutree(cl.tree, ncolors)
    f.col <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Set1'))
    cp <- f.col(ncolors)
  }
  
  #Plotting tree
  plot(as.phylo(cl.tree), type = plot.type, cex = 1, tip.color = cp[clus.col])
  

  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
}

#Plot the distance matrix between clusters ordered by a hierarchical clustering
plot.ic.cluster.dist <- function(mat.ic.cl.dist, spots.table, method = 'average', file.path = 'distance-matrix.pdf', ncolors = 50){

  #mat.ic.cl.dist: matrix of IC/cluster distance based on averaged IC expression
  #spots.table: used for normalizing by number of spots per cluster
  #method(string): a method accepted by hclust algorithm
  #file.path(string): path or file name, with extension. If null, plots are not saved in a file.
  #ncolors(numeric): if #0, then the tree is colored with that many colors based on branching
  
  #If saving the plot in a file
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path, width = 100, height = 100) 
  }else{
    plot.new()
  }

  #Size of clusters
  cl.count <- count(spots.table$clusters.named)
  
  #Comverting into a distance matrix
  mat.dist <- mat.ic.cl.dist
  
  #Computing tree
  cl.tree <- hclust(as.dist(mat.dist), method = method, members = cl.count$freq)
  
  #Reordering the distance matrix based on the hierarchical order
  mat.dist.ordered <- mat.ic.cl.dist[cl.tree$order,]
  mat.dist.ordered <- mat.dist.ordered[,cl.tree$order]
  
  #Number of clusters
  n.clust <- dim(mat.dist.ordered)[2]
  
  #Colors and breaks
  n.col <- 100
  f.col <- colorRampPalette(brewer.pal(9,'YlOrRd'))
  cp.dist <- f.col(n.col)
  break.points <- seq(from = 0, to = max(mat.dist.ordered), length.out = n.col+1)

  #Coloring tree if required
  if(ncolors > 0){
    clus.col <- cutree(cl.tree, ncolors)
    f.col <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Set1'))
    cp <- f.col(ncolors)
    col.id <- clus.col[rownames(mat.dist.ordered)]
  }
  
  #Adding outside box as a "split"
  x.splits <- which(as.logical(diff(col.id) != 0)) + 0.5
  x.splits <- c(0.5,x.splits,n.clust + 0.5)

  #Plotting
  image(mat.dist.ordered, 
        x = 1:n.clust,
        y = 1:n.clust,
        xlim = c(0.5-n.clust/10, n.clust+0.5+n.clust/10),
        ylim = c(0.5-n.clust/10, n.clust+0.5+n.clust/10),
        breaks = break.points,
        col = cp.dist,
        xlab = '',
        ylab = '',
        asp = 1,
        bty = 'n',
        xaxt = 'n',
        yaxt = 'n')
  
  #Adding labels
  text(1:n.clust + 0.2, 0, pos = 2, colnames(mat.dist.ordered), col = cp[col.id], cex = 2, srt = 90)
  text(0, 1:n.clust - 0.1, pos = 2, colnames(mat.dist.ordered), col = cp[col.id], cex = 2)

  #Separating line
  segments(x.splits, rep(-10, length(x.splits)), x.splits, rep(n.clust + 0.5, length(x.splits)), lwd = 0.5, col = 'black')
  segments(rep(-10, length(x.splits)), x.splits, rep(n.clust + 0.5, length(x.splits)), x.splits, lwd = 0.5, col = 'black')  
  
  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
  
}

#Plot the IC load of different clusters for the sake of comparison
plot.ic.cluster.diff <- function(mat.ic.cl.avg, cl.list, df.col = NULL, file.path = 'ic-clusters-difference.pdf'){

  #mat.ic.cl.dist: matrix of IC/cluster distance based on averaged IC expression
  #cl.list(string): list of cluster to plot (clusters.named)
  #df.col: data frame with color information for plotting clusters
  #file.path(string): path or file name, with extension. If null, plots are not saved in a file.
    
  #If saving the plot in a file
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path, width = 8, height = 7) 
  }else{
    plot.new()
  }
  
  #Creating color data frame if required
  if(is.null(df.col))
    df.col <- get.default.df.colors(cl.list)
  
  
  #Selecting sub matrix of interest
  sub.mat <- mat.ic.cl.avg[cl.list,]
  
  #Definning xlim and ylim
  maxx <- max(sub.mat)
  minx <- min(sub.mat)
  xdiff <- maxx-minx
  xlim <- c(minx-0.05*xdiff,maxx+0.05*xdiff)
  n.ic <- dim(sub.mat)[2]
  ylim <- rev(c(0, n.ic+1))
  
  #Sequence of vertical lines
  x.lines <- seq(floor(xlim[2]*2)/2, ceiling(xlim[1]*2)/2, by = -0.5)
  x.lines <- setdiff(x.lines,0)
  
  #Initializing title
  title = 'Comparison of IC loads'
  
  #Margin for the legend
  par(mar= c(3.8, 4.1, 4.1, 10.3),
      mgp = c(1, 1, 0),
      xpd=TRUE)
  
  #Plotting
  plot(1,
       type = 'n',
       main = title,
       xlim = xlim,
       ylim = ylim,
       xlab = 'IC load',
       ylab = '',
       yaxt = 'n',
       xaxt = 'n',
       bty = 'n')
  
  #Grid segments
  segments(xlim[1],1:(dim(sub.mat)[2]),xlim[2],1:(dim(sub.mat)[2]), col = 'lightgray', lty = 2, lwd = 0.5)  
  segments(x.lines,rep(ylim[1],length(x.lines)),x.lines,rep(ylim[2],length(x.lines)), lty = 1, 
           col = 'lightgray', lwd = 0.5, lend = 'butt')
  segments(0,ylim[1],0,ylim[2], lty = 1, col = 'darkgray', lwd = 2.5, lend = 'butt')
  
  #Plotting individuals clusters
  for(i in 1:(dim(sub.mat)[1])){
    points(sub.mat[i,],
           1:(dim(sub.mat)[2]),
           col = df.col[rownames(sub.mat)[i],'color'],
           pch = 16)
  }
  
  #Box
  axis(1, pos = ylim[1], at = xlim, tick = TRUE, labels = NA, lwd.ticks = 0)
  axis(2, pos = xlim[1], at = ylim, tick = TRUE, labels = NA, lwd.ticks = 0)
  axis(3, pos = ylim[2], at = xlim, tick = TRUE, labels = NA, lwd.ticks = 0)
  axis(4, pos = xlim[2], at = ylim, tick = TRUE, labels = NA, lwd.ticks = 0)
  
  #Ticks
  axis(1, pos = ylim[1], at = sort(c(x.lines,0)), cex.axis = 0.5, xlab = 'IC load')
  axis(2, pos = xlim[1], at = 1:(dim(sub.mat)[2]), las = 2, cex.axis = 0.5, labels = colnames(sub.mat))
  
  #Color for the legend
  col.legend <- df.col[cl.list,'color']

  #Plotting the legend
  legend("right", 
         inset = c(-0.41,0), 
         legend = cl.list, 
         pch = rep(16,length(cl.list)), 
         col = col.legend,
         title="Clusters",
         bty = 'n')
  
  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
}

