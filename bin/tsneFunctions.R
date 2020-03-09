## ----------- Functions for manipulating tsne related data ----------- 

#Add the two columns tsne1 and tsne2 to the spots table from a tsne file
append.tsne.to.spots.table <- function(spots.table, tsne.path){
  
  #spots.table
  #tsne.path(char): path to a tsne file
  
  t <- read.table(tsne.path, sep='\t', stringsAsFactors = F, header = TRUE, row.names = 1)
  sel.rows <- dplyr::intersect(rownames(t),rownames(spots.table))
  
  spots.table[sel.rows,'tsne1'] <- t[sel.rows,1]
  spots.table[sel.rows,'tsne2'] <- t[sel.rows,2]
  
  return(spots.table)
}

#Returns a data frame providing colors for every cluster in spots table based on median coordinates in the 3D tsne
get.color.from.3d.tsne <- function(path.to.tsne, spots.table, use.hsv = FALSE, vector.order = c(1, 2, 3)) {
  #path.to.tsne(str): path to a 3d tsne file
  #spots.table
  #vector.order(numeric): gives which coordinate should be used for r,g,b or h,s,v
  
  spots.table$clusters.named <-
    as.character(spots.table$clusters.named)
  
  #Reading the tsne
  t <-
    read.table(path.to.tsne,
               sep =  '\t',
               row.names = 1,
               header = T)
  colnames(t) <- c('tsne1', 'tsne2', 'tsne3')
  
  #Converting into rgb proporitions
  t1 <- t[, vector.order[1]]
  t2 <- t[, vector.order[2]]
  t3 <- t[, vector.order[3]]
  
  #Normalizing vectors between 0 and 1
  vect.1 <- (t1 - min(t1)) / abs(max(t1) - min(t1))
  vect.2 <- (t2 - min(t2)) / abs(max(t2) - min(t2))
  vect.3 <- (t3 - min(t3)) / abs(max(t3) - min(t3))
  
  #Initializing data frame
  df.colors <- unique(spots.table[, c('clusters.named', 'cluster')])
  rownames(df.colors) <- df.colors$clusters.named
  df.colors$clusters.named <- NULL
  colnames(df.colors) <- 'cluster.id'
  rownames(df.colors) <- unique(spots.table$clusters.named)
  
  #Filling it with median coordinates of the cluster
  for (cl in rownames(df.colors)) {
    sel.spots <-
      is.element(rownames(t), rownames(spots.table[spots.table$clusters.named == cl, ]))
    if (use.hsv) {
      df.colors[cl, 'clusters.colors'] <- hsv(median(vect.1[sel.spots]),
                                              median(vect.2[sel.spots]),
                                              median(vect.3[sel.spots]))
    } else{
      df.colors[cl, 'clusters.colors'] <- rgb(median(vect.1[sel.spots]),
                                              median(vect.2[sel.spots]),
                                              median(vect.3[sel.spots]))
    }
  }
  
  return(df.colors)
}

## ----------- Functions for plotting tsne  ----------- 

#Plot a t-SNE. Should not be called directly, unconvenient
.plot.tsne <- function(sp, title ='', list.legend = NULL, file.path = 'tsne.pdf', no.box = FALSE, cex = 0.1){
  
  #sp: spots.table with tsne1 and tsne2 appended + color containing the proper color
  #title(char): plot title
  #list.legend: list with information for plotting legend (cl.list and col.legend for items, colorscale and bin.sep for colorscale). If null, no legend plotted
  #file.path:  #(char): path or file name, with extension. If null, plots are not saved in a file. 
  #no.box(logical): if TRUE; no box neither axis are plotted
  #cex(numeric): cex for plotting
  
  #If saving the plot in a file
  if (!is.null(file.path)){
    file.path <- generate.appropriate.file.name(file.path)
    pdf(file.path, width = 8, height = 8, useDingbats = F) 
  }

  par(mar= c(0.2, 0.2, 0.2, 0.2),
      xpd=TRUE)
  
  #Margin for the legend if required
  if (!is.null(list.legend)){
    
    prev.val <- list(mar = par('mar'),
                     xpd = par('xpd'))

    #Margin for colorscale
    if (!is.null(list.legend$colorscale)){
      par(mar= c(3.8, 4.1, 4.1, 5.3),
          xpd=TRUE)
    
    #Margin for labels
    }else
      par(mar= c(3.8, 4.1, 4.1, 10.3),
          xpd=TRUE)
  }
  
  if(no.box){
    
    #Plotting    
    plot(sp$tsne1, 
         sp$tsne2, 
         col = sp$color, 
         asp = 1, 
         pch = 20, 
         cex = cex, 
         main = title,
         xlab = '',
         ylab = '',
         xaxt = 'n',
         yaxt = 'n',
         bty = 'n')    
  }else{
  
    #Plotting    
    plot(sp$tsne1, 
         sp$tsne2, 
         col = sp$color, 
         asp = 1, 
         pch = 20, 
         cex = cex, 
         main = title,
         xlab = 't-SNE 1',
         ylab = 't-SNE 2')
  }
  
  #Plotting legend if required
  if (!is.null(list.legend)){
    
    #Items case
    if (!is.null(list.legend$cl.list)){
      legend("right", 
             inset = c(-0.41,0), 
             legend = list.legend$cl.list, 
             pch = rep(16,length(list.legend$cl.list)), 
             col = list.legend$col.legend,
             title="Clusters",
             bty = 'n')
    
    #Colorscale case  
    }else{
      
      #x coordinates
      x.left <-  par('usr')[1]
      x.right <- par('usr')[2]
      x.dist <- x.right-x.left
      
      rect.x1 <- x.right + 0.025*x.dist
      rect.x2 <- x.right + 0.075*x.dist
      
      #y coordinates
      y.bot <- par('usr')[3]
      y.top <- par('usr')[4]
      y.dist <- y.top - y.bot
      
      rect.y1 <- y.bot + 0.15 * y.dist
      rect.y2 <- y.top - 0.15 * y.dist
      
      #n bins
      n <- length(list.legend$colorscale)

      #creating vector for all gradient rectangle
      l.x1 <- rep(rect.x1, n)
      l.x2 <- rep(rect.x2, n)
      l.y <- seq(from = rect.y1, to = rect.y2, length.out = (n+1))
      l.y1 <- l.y[1:n]
      l.y2 <- l.y[2:(n+1)]

      #Plotting outside rectangle
      rect(rect.x1, rect.y1, rect.x2, rect.y2, col = NA)
      
      #Plotting gradient rectangles
      rect(l.x1, l.y1, l.x2, l.y2, col = list.legend$colorscale, border = NA)
      
      #Default value for digits
      if (is.null(list.legend$digits))
        list.legend$digits <- 1
      
      y.text <- seq(from = rect.y1, to = rect.y2, length.out = 5)
      str.text <- round(seq(from = list.legend$bin.sep[1], to = list.legend$bin.sep[n+1], length.out = 5), digits = list.legend$digits)

      #Legend labels
      text(rep(rect.x2,5),
           y.text,
           str.text,
           pos = 4)
      
      text(rect.x1 - 0.01*x.dist, mean(c(rect.y1,rect.y2)), 'Expression level', srt = 90, pos = 3)
    }
  }
  
  #Putting back default values
  if (!is.null(list.legend)){
    par(mar = prev.val$mar,
        xpd = prev.val$xpd)
  }
  
  #Closing dev only if writting in pdf
  if (!is.null(file.path))
    dev.off()
}

#Plot all the spots on a t-SNE with different coloring methods
plot.tsne.all <- function(path.to.tsne, #(char): path to the t-SNE file
                          spots.table, #spots.table
                          color.mode = c('ARA','ARA.main','clusters','AP','sections'), #Criteria used for coloring
                          file.path = 'tsne.pdf', #(char): path or file name, with extension. If null, plots are not saved in a file.
                          cex = 0.1, #(float): points cex
                          no.box = FALSE)
{
 
  #Creating a long color vector
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rep(col_vector, 10)
  
  #Appending the tsne
  sp <- append.tsne.to.spots.table(spots.table, path.to.tsne)

  #Selecting only the first color mode
  color.mode <- color.mode[1]
  
  #Extracting the proper color depending on the color mode
  if (color.mode == 'ARA')
    sp$color <- color.from.acronym(sp$acronym)
  else if (color.mode == 'ARA.main')
    sp$color <- color.from.acronym(sp$acronym.parent)
  else if (color.mode == 'clusters')
    sp$color <- mapvalues(as.character(sp$clusters.named), 
                          from = unique(as.character(sp$clusters.named)),
                          to = col_vector[1:length(unique(as.character(sp$clusters.named)))])
  else if (color.mode == 'AP')
    sp$color <- mapvalues(sp$AP, 
                          from = unique(sp$AP), 
                          to = col_vector[1:length(unique(sp$AP))])
  else if (color.mode == 'section')
    sp$color <- mapvalues(sp$slice_index, 
                          from = unique(sp$slice_index), 
                          to = col_vector[1:length(unique(sp$slice_index))])

  #Setting up title
  title <- paste('Color mode:', color.mode, sep=' ')
  
  #Plotting
  .plot.tsne(sp, title, file.path = file.path, cex = cex, no.box = no.box)

}

#Plot only specific spots by cluster on a t-SNE
plot.tsne.clusters <- function(path.to.tsne, #path to the t-SNE file
                               spots.table, #spots.table
                               cl.list, #list of clusters to plot (as clusters.named field)
                               df.col = NULL, #data frame with color information for plotting clusters
                               file.path = 'tsne-clusters.pdf') #(char): path or file name, with extension. If null, plots are not saved in a file.
{

  #Creating color data frame if required
  if(is.null(df.col))
    df.col <- get.default.df.colors(cl.list)
  
  #Appending the tsne
  sp <- append.tsne.to.spots.table(spots.table, path.to.tsne)
  
  #Initializing all color to gray
  sp$color <- 'gray'
  for (cl in cl.list){
    sp[sp$clusters.named == cl, 'color'] <- df.col[cl, 'color']
  }

  #Seting up legend and title
  col.legend <- df.col[cl.list,'color']
  list.legend <- list(cl.list = cl.list, col.legend = col.legend)
  
  title <- paste(cl.list, collapse='/')
  
  #Plotting
  .plot.tsne(sp, title, list.legend = list.legend, file.path = file.path)
}  
    
#Plot the expression level of an IC on the t-SNE
plot.tsne.ic <- function(path.to.tsne, #path to the t-SNE file
                         spots.table, #spots.table
                         mat.ic, #matrix of ic from get.ic.mat
                         ic = 'IC1', #(char): ic name
                         file.path = 'tsne-ic.pdf', #(char): path or file name, with extension. If null, plots are not saved in a file.
                         no.box = FALSE) #(logical): if TRUE, not box is plotted arround the t-SNE
{
 
  #Appending the tsne
  sp <- append.tsne.to.spots.table(spots.table, path.to.tsne)
  
  #Only keeping shared spots between mat.ic and spots.table
  shared.spots <- intersect(rownames(sp), rownames(mat.ic))
  sp <- sp[shared.spots,]
  mat.ic <- mat.ic[shared.spots,]
  
  #Selecting corresponding component
  vect.ic <- mat.ic[,ic]
  
  #Adding the expression
  sp[rownames(mat.ic), 'ic.load'] <- vect.ic
  
  #Binning the expression level
  bin.sep <- seq(from = -max(abs(vect.ic)), to = max(abs(vect.ic)), length.out = 202)
  bin.sep[1] <- bin.sep[1] - 0.1
  bin.sep[length(bin.sep)] <- bin.sep[length(bin.sep)] + 0.1
  sp$ic.bin <- .bincode(sp$ic.load, breaks = bin.sep)
  
  #Creating color palette
  # function.cp <- colorRampPalette(brewer.pal(11, 'RdBu'))
  # cp <- rev(function.cp(201))
  function.cp <- colorRampPalette(c('#0004ff', '#bfbfbf', '#ff0000'))
  cp <- function.cp(201)
  
  list.legend <- list(colorscale = cp,
                      bin.sep = bin.sep)
  
  sp$color <- cp[sp$ic.bin]

  .plot.tsne(sp, ic, file.path = file.path, no.box = no.box, list.legend = list.legend)
}

#Plot the expression level of a gene on the t-SNE
plot.tsne.gene <- function(path.to.tsne, #path to the t-SNE file
                           spots.table, #spots.table
                           st.data, #Gene expression matrix
                           gene, #(char): ic name
                           file.path = 'tsne-gene.pdf', #(char): path or file name, with extension. If null, plots are not saved in a file.
                           no.box = FALSE, #(logical): if TRUE, not box is plotted arround the t-SNE)
                           cex =  0.1) #(numerical): cex for plotting
{
  
  #Appending the tsne
  sp <- append.tsne.to.spots.table(spots.table, path.to.tsne)
  
  #Only keeping shared spots between mat.ic and spots.table
  shared.spots <- intersect(rownames(sp), rownames(st.data))
  sp <- sp[shared.spots,]
  st.data <- st.data[shared.spots,]
  
  #Selecting corresponding component
  vect.gene <- st.data[,gene]
  
  #Adding the expression
  sp[rownames(st.data), 'gene.expr'] <- vect.gene
  
  #Binning the expression level
  bin.sep <- seq(from = min(min(vect.gene),0), to = max(vect.gene), length.out = 202)
  bin.sep[1] <- bin.sep[1] - 0.001
  sp$gene.bin <- .bincode(sp$gene.expr, breaks = bin.sep)
  
  #Ordering to plot top expression on top
  sp <- sp[order(sp$gene.bin),]
  
  #Creating color palette
  # function.cp <- colorRampPalette(brewer.pal(9, 'YlOrRd'))
  function.cp <- colorRampPalette(c('#bfbfbf', '#ff0000'))
  cp <- function.cp(201)
  sp$color <- cp[sp$gene.bin]
  
  list.legend <- list(colorscale = cp,
                      bin.sep = bin.sep,
                      digits = 0)
  
  #Plotting
  .plot.tsne(sp, gene, file.path = file.path, no.box = no.box, list.legend = list.legend, cex = cex)
  
}
