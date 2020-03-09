## ----------- Pixels - stereo conversion ----------- 

#Returns the ratio to convert stereo units into pixels
get.ratio.stereo.to.pixels <- function(section){
  
  vect.mean <- c((max(section$somaX) - min(section$somaX)) / (max(section$ML) - min(section$ML)),
                 (max(section$somaY) - min(section$somaY)) / (max(section$DV) - min(section$DV)))  
  return(mean(vect.mean))
}

#Converts a list of ML coordinates into pixels
ML.to.pixels <- function(ML.list, section){
  return((ML.list - min(section$ML)) * get.ratio.stereo.to.pixels(section) + min(section$somaX))
}

#Converts a list of DV coordinates into pixels
DV.to.pixels <- function(DV.list, section){
  return(-(DV.list - max(section$DV)) * get.ratio.stereo.to.pixels(section) + min(section$somaY))
}

#Converts a list of pixels coordinates into stereo ML
pixels.to.ML <- function(pixels.list, section){
  return((pixels.list - min(section$somaX)) * (1/get.ratio.stereo.to.pixels(section)) + min(section$ML))
}

#Converts a list of pixels coordinates into stereo DV
pixels.to.DV <- function(pixels.list, section){
  return(-(pixels.list - max(section$somaY)) * (1/get.ratio.stereo.to.pixels(section)) + min(section$DV))
}

## ----------- Returning usefeul data for plotting ----------- 

#Returns a data.frame to use with segments function for plotting ARA contour without meaningless lines
get.segments.outlines <- function(outlines, field.name.x = 'xr', field.name.y = 'yr', scale.factor = 1){
  
  #outlines: list of outlines
  #field.name.x(char): name of the x field
  #field.name.y(char): name of the y field
  #scale.factor(numeric): potential scale factor to divide cooridnates by
  
  #Counting the polygon points to initialize data frame
  cnt.n.points <- 0
  for (i in 1:length(outlines))
    cnt.n.points <- cnt.n.points + length(outlines[[i]][[field.name.x]])
  
  
  #Initializing the data frame that will contain all segments
  n.rows.seg <- cnt.n.points - length(outlines)
  df.segments <- data.frame(x0 = as.numeric(rep(NA, n.rows.seg)),
                            y0 = as.numeric(rep(NA, n.rows.seg)),
                            x1 = as.numeric(rep(NA, n.rows.seg)),
                            y1 = as.numeric(rep(NA, n.rows.seg)),
                            stringsAsFactors = FALSE)
  
  cnt.row.seg <- 1
  
  #Going through the paths to plot
  for (i in 1:length(outlines)){
    
    #Detecting splits
    splits <- which(names(outlines[[i]][[field.name.y]]) == 'move')
    
    # Appending length of list as final split
    splits <- c(splits,(length(outlines[[i]][[field.name.y]])+1))
    
    #Going through the splits
    for (j in 1:(length(splits)-1)){
      
      sp.start <- splits[j]
      sp.end <- splits[(j+1)] - 1
      
      range.0 <- sp.start:(sp.end-1)
      range.1 <- (sp.start+1):sp.end
      
      df.segments[cnt.row.seg:(cnt.row.seg + length(outlines[[i]][[field.name.y]][range.1]) -1),
                  c('x0','y0','x1','y1')] <- cbind(
                    outlines[[i]][[field.name.x]][range.0],
                    outlines[[i]][[field.name.y]][range.0],
                    outlines[[i]][[field.name.x]][range.1],
                    outlines[[i]][[field.name.y]][range.1])
      
      #Moving cursor
      cnt.row.seg <- cnt.row.seg + length(outlines[[i]][[field.name.y]][range.1])
      
    }
  }
  
  df.segments <- df.segments[!is.na(df.segments$x0),]
  df.segments <- df.segments / scale.factor
  return(df.segments)
}

## ----------- Returning usefeul data for plotting ----------- 

#OLD 2d mixed
# #Plot both the ARA as right hemisphere and the generated atlas as left hemisphere in 2D
# plot.2d.mixed <- function(grid.3d, id.grid, df.colors, grid.parameters, pdf.name = NULL, right.hemisphere = TRUE, show.saggital = TRUE){
# 
#   #grid.3d: grid 3d object
#   #id.grid(num): coronnal cut in the grid to plot
#   #df.colors: a df.colors object from the tsne type colors of objects
#   #grid.parameters: grid parameters object
#   #pdf.name(char): if not null, the plot is saved in the pdf name (should not include '.pdf')
#   #right.hemisphere(logical): if true, right hemisphere is the ARA. Otherwise, it's the left one
#   #show.saggital(logical): display saggital view
#   
#   if (show.saggital)
#     yl <- c(-8,2)
#   else
#     yl <- c(-8,0)
#   
#   
#   if (right.hemisphere)
#     x.coeff <- -1
#   else
#     x.coeff <- 1
#   
#   #Reversing order  
#   selected.AP <- dim(grid.3d)[3] - id.grid + 1
#   
#   #Finding corresponding AP and closest atlas
#   AP <- grid.parameters$AP[id.grid]
#   id.atlas <- which.min(abs(aligned.atlas$AP - AP))
#   
#   #Vectoring for remapping
#   x.vect.remap <- ML.to.pixels(x.coeff*grid.parameters$ML, aligned.atlas$soma[[id.atlas]])
#   y.vect.remap <- DV.to.pixels(grid.parameters$DV, aligned.atlas$soma[[id.atlas]])
#   
#   #Extracting only 2D infos
#   grid.2d <- grid.3d[,,id.grid]
#   
#   #Detecting edges
#   grid.edges <- matrix(FALSE, nrow = dim(grid.2d)[1], ncol = dim(grid.2d)[2])
#   for (i in 1:(dim(grid.edges)[1]-1)){
#     grid.edges[i,] <- (grid.2d[i+1,] - grid.2d[i,]) != 0
#   }
#   for (i in 1:(dim(grid.edges)[2])-1){
#     grid.edges[,i] <- ((grid.2d[,i+1] - grid.2d[,i]) != 0) | grid.edges[,i]
#   } 
#   pix <- which(grid.edges, arr.ind = T)
#   x.edges <- x.vect.remap[pix[,1]]
#   y.edges <- y.vect.remap[pix[,2]]
#   
#   #Computing distance between pixels
#   d <- dist(pix)
#   d <- as.matrix(d)
#   
#   #Finding close pixels (distance of one, diagonals included)
#   close.p <- which(d <= sqrt(2), arr.ind = T)
#   close.p <- close.p[(close.p[,1] != close.p[,2]),]
#   close.p <- t(apply(close.p, 1, sort))
#   close.p <- unique(close.p)
# 
#   #If a pdf should be open, opening it
#   if (!is.null(pdf.name))
#     pdf(paste(pdf.name,'.pdf',sep=''), width = 10, height = 7)
#   
#   #Plotting the outline on the right
#   plot.brain.outline(aligned.atlas$outlines[[id.atlas]],
#                      aligned.atlas$color[[id.atlas]],
#                      aligned.atlas$soma[[id.atlas]],
#                      'x', 'y', AP = AP, xlim.stereo = c(-6,6), ylim.stereo = yl, right.hemisphere = right.hemisphere)
# 
#   #Plotting all clusters
#   for(cl in setdiff(unique(as.numeric(grid.2d)),0)){
#     pts <- which(grid.2d == cl, arr.ind = T)
#     x <- x.vect.remap[pts[,1]]
#     y <- y.vect.remap[pts[,2]]
#     points(x,y,col = df.colors[df.colors$cluster.id == cl,'clusters.colors'], pch = 15, cex = 0.2)
#   }
#   
#   #Plotting edges
#   # points(x.edges, y.edges, col = 'black', cex = 0.1, pch = 15)
#   segments(x.edges[close.p[,1]], y.edges[close.p[,1]], x.edges[close.p[,2]], y.edges[close.p[,2]], lwd = 0.6)
#   
#   if (show.saggital){
#   
#     #Plotting saggital overivew
#     plot.saggital.view(ML.to.pixels(-5.8, aligned.atlas$soma[[id.atlas]]),
#                        DV.to.pixels(1.8, aligned.atlas$soma[[id.atlas]]),
#                        DV.to.pixels(-0.2, aligned.atlas$soma[[id.atlas]]),
#                        AP)
#     
#   }
#   #If a pdf was opened, closing it
#   if (!is.null(pdf.name))
#     dev.off()
# }


#Plot both the ARA as right hemisphere and the generated atlas as left hemisphere in 2D
plot.2d.mixed <- function(grid.3d, id.grid, df.colors, grid.parameters, pdf.name = NULL, right.hemisphere = TRUE, show.saggital = TRUE, show.title = TRUE, show.box = TRUE, increase.factor = 2, fill.atlas = TRUE){

  #grid.3d: grid 3d object
  #id.grid(num): coronnal cut in the grid to plot
  #df.colors: a df.colors object from the tsne type colors of objects
  #grid.parameters: grid parameters object
  #pdf.name(char): if not null, the plot is saved in the pdf name (should not include '.pdf')
  #right.hemisphere(logical): if true, right hemisphere is the ARA. Otherwise, it's the left one
  #show.saggital(logical): display saggital view

  xlim.stereo <- c(-6,6)
  ylim.stereo <- c(-8,0)
  
  if (show.saggital)
    ylim.stereo <- c(-8,2)

  #Inverting coordinates if right hemisphere
  if (right.hemisphere)
    x.coeff <- -1
  else
    x.coeff <- 1

  #Finding corresponding AP and closest atlas
  AP <- grid.parameters$AP[id.grid]
  id.atlas <- which.min(abs(atlas.stereo$AP - AP))

  outlines <- atlas.stereo$outlines[[id.atlas]]
  
  #Getting the segments to plot
  df.segments <- get.segments.outlines(outlines, 'ML', 'DV', 1)  
  
  #Extracting only 2D infos
  grid.2d <- grid.3d[,,id.grid]

  title <- ''
  
  if(show.title)
    title <- sprintf('Bregma: %.2f mm', AP)
  
   #If a pdf should be open, opening it
  if (!is.null(pdf.name))
    pdf(generate.appropriate.file.name(pdf.name), width = 10, height = 7, useDingbats = F)

  plot(1, type = 'n', asp = 1, ylim = ylim.stereo, xlim = xlim.stereo, 
       main = title, bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='') 

  if(fill.atlas){
    #Plotting the inside filling
    for (x in 1:length(outlines)){
      
      #Current outline
      o <- outlines[[x]]
      
      col <- atlas.stereo$plate.info[[id.atlas]][x,'col']
     
      #Plotting polygons for colors
      polygon(o$ML * x.coeff, o$DV, col = col, border = NA)
      
    }
  }
  #Plotting the contours
  segments(x.coeff*df.segments$x0, df.segments$y0, x.coeff*df.segments$x1, df.segments$y1, lwd = lwd.atlas.segments)

  if(show.box){

    #Displaying the axis
    axis(3, pos = ylim.stereo[2], at = xlim.stereo[1]:xlim.stereo[2], labels = xlim.stereo[1]:xlim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0))
    axis(2, pos = xlim.stereo[1], at = ylim.stereo[1]:ylim.stereo[2], labels = ylim.stereo[1]:ylim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0), las= 1)
    axis(1, pos = ylim.stereo[1], at = xlim.stereo[1]:xlim.stereo[2], labels = xlim.stereo[1]:xlim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0))
    axis(4, pos = xlim.stereo[2], at = ylim.stereo[1]:ylim.stereo[2], labels = ylim.stereo[1]:ylim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0), las= 1)
  }
  
  plot.raster.atlas(grid.3d[,,id.grid], grid.parameters, df.colors, increase.factor = increase.factor, file.path = NULL, add = TRUE, right.hemisphere = right.hemisphere)
  
 
  if (show.saggital){

    #Plotting saggital overivew
    plot.saggital.view(-5.8, 1.8, -0.2, AP)

  }
  #If a pdf was opened, closing it
  if (!is.null(pdf.name))
    dev.off()
}

#Plot a saggital view of the brain, with a line at AP coordinates lines.AP
#The plot should have inverted Y axis
plot.saggital.view <- function(xl, yt, yb, lines.AP = NULL){
  
  #xl(num): leftiest x coordinates desired for the saggital plot
  #yt(num): top coordinates of the desired saggital plot (smaller y value)
  #yb(num): low coordinates of the desired saggital plot (highest y value)  
  #lines.AP(num): AP coordinates of desired lines
  
  #Gives X coordinates of any AP in the saggital coordinate system
  get.x.from.ap <- function(ap.list, ap.bregma, x.bregma){
    
    x.vect <- numeric(length(ap.list))
    cnt <- 0
    for (a in ap.list){
      cnt <- cnt + 1
      prop.a <- (ap.bregma[1] - a) / abs(ap.bregma[1] - ap.bregma[2])
      x.vect[cnt] <- prop.a * abs(x.bregma[1] - x.bregma[2]) + x.bregma[1]
    }
    
    return(x.vect)
  }
  
  #Convert into new x coordinates (based on xl, yt, yb system)
  n.x <- function(x.coord){
    return((x.coord - saggital.view$xscale[1]) / abs(diff(saggital.view$xscale)) * (yb - yt) * abs(diff(saggital.view$xscale)/diff(saggital.view$yscale)) + xl)
  }
  
  #Convert into new x coordinates (based on xl, yt, yb system)  
  n.y <- function(y.coord){
    return((saggital.view$yscale[2] - y.coord) / abs(diff(saggital.view$yscale)) * (yb - yt) + yt)
  }
  
  #Loading the saggital view data
  load(paste(path.bin, 'saggitalView.RData', sep='/'))
  
  df.poly <- saggital.view$df.poly
  df.segments <- saggital.view$df.segments
  
  #Plotting polygons and outlines
  c <- 'white'
  # polygon(n.x(df.poly[df.poly$col == c, 'x']),
  #         n.x(df.poly[df.poly$col == c, 'y']),
  #         col = c, border = NA)
  segments(n.x(df.segments[df.segments$col == c, 'x0']),
           n.y(df.segments[df.segments$col == c, 'y0']), 
           n.x(df.segments[df.segments$col == c, 'x1']),
           n.y(df.segments[df.segments$col == c, 'y1']),
           lend = 'butt')
  
  for (c in setdiff(unique(df.poly$col),'white')){ 
    polygon(n.x(df.poly[df.poly$col == c, 'x']),
            n.y(df.poly[df.poly$col == c, 'y']),
            col = c, border = NA)
    segments(n.x(df.segments[df.segments$col == c, 'x0']),
             n.y(df.segments[df.segments$col == c, 'y0']), 
             n.x(df.segments[df.segments$col == c, 'x1']),
             n.y(df.segments[df.segments$col == c, 'y1']),
             lend = 'butt')
  }
  
  #Plotting AP lines
  if (!is.null(lines.AP)){
    x.coord <- n.x(get.x.from.ap(lines.AP, saggital.view$ap.bregma, saggital.view$x.bregma))
    segments(x.coord, n.y(saggital.view$yscale[1] - 0.1*abs(diff(saggital.view$yscale))),
             x.coord, n.y(saggital.view$yscale[2] + 0.1*abs(diff(saggital.view$yscale))),
             col = 'red', lty = 1, lend = 'butt', lwd = 3)
  }
}

## ----------- Extracting and adding information related to ARA ----------- 

#Returns the depth of an acronym based on the allen ontology
get.acronym.depth <- function(acronyms.vect){
  
  #acronyms.vect(str): an acronym. If vector, applied to all of the vectors.
  
  library(rjson)
  
  #Local recursive function
  loc.function <- function(json, cur.depth = 0, acronym){
    
    #If acronym was found, return current depth
    if (json$acronym == acronym)
      return(cur.depth)
    
    #If leaf, return -1
    if (length(json$children) == 0)
      return(-1)
    
    #Otherwise, initialize vector and call it on every children
    v <- numeric()
    for (i in 1:length(json$children)){
      v <- c(v, loc.function(json$children[[i]], (cur.depth+1), acronym))
    }
    
    #Return the vector
    return(v)
  }
  
  #Loading json data
  json_data <- fromJSON(file= paste(path.bin,'ontology.json',sep='/'))
  json <- json_data$msg[[1]]
  
  #Calling local function and returning maximal depth found
  return(sapply(acronyms.vect,function(x){return(max(loc.function(json, 0,x)))}))
}

#Adds a field to the spots table with all levels from the ontology
add.all.acronyms <- function(spots.table){
  
  #spots.table
  
  #Using a recursive functions, as long as we have not reached the main level.
  rec.func <- function(x){
    if (x == 'root')
      return(x)
    else
      return(c(x,rec.func(get.acronym.parent(x))))
  }
  
  spots.table$all.acronyms <- spots.table$acronym
  spots.table$all.acronyms <- lapply(spots.table$all.acronyms,rec.func)
  
  return(spots.table)
}

#Appends the acronym and the full name of the parent of a spot at the level specificied in the acronym list
add.parent.acronym <- function(spots.table,
                               list.acronym.parents = c('TH','HY','MB','HB',
                                                        'CB','STR','PAL','Isocortex','OLF',
                                                        'HIP','RHP', 'CTXsp','root','fiber tracts','VS')){

  #spots.table
  #list.acronym.parents(str): vector of region acronyms which are the desired parents. Always include root to avoid infinite loops.
  
  library(wholebrain)
  list.full.name.parents <- as.character(name.from.acronym(list.acronym.parents))
  list.full.name.parents[list.full.name.parents == 'root'] <- 'Undefined areas'
  spots.table$acronym.parent <- spots.table$acronym
  
  #Finding appropriate parent for all acronyms
  for (i in 1:length(spots.table$acronym.parent)){
    
    if (!is.na(spots.table$acronym.parent[i])){
      #Getting parents until reaching appropriate level
      while (!is.element(spots.table$acronym.parent[i],list.acronym.parents)){
        spots.table$acronym.parent[i] <- get.acronym.parent(spots.table$acronym.parent[i])
      }
      
      #Converting to full name
      spots.table$full.name.parent[i] <- list.full.name.parents[list.acronym.parents == spots.table$acronym.parent[i]]
    }
  }
  
  #Areas as factor
  spots.table$full.name.parent <- as.factor(spots.table$full.name.parent)
  
  return(spots.table)
}

#Returns the allen color based on the cluster name
get.allen.color.cluster <- function(cl.list, spots.table){
  
  #cl.list(string): list of clusters
  #spots.table: used for getting color from full name
  
  #Initializing data frame
  df <- data.frame(cluster = cl.list, stringsAsFactors = F)
  rownames(df) <- df$cluster
  
  #Extracting main region
  sp <- strsplit(df$cluster, '[-]')
  df$region <- unlist(lapply(sp, function(x){return(x[[1]])}))
  
  #Getting acronym correspondance
  df.remap <- unique(spots.table[,c('acronym.parent','full.name.parent')])
  df.remap$full.name.parent <- as.character(df.remap$full.name.parent)
  df.remap$color <- color.from.acronym(df.remap$acronym.parent)
  df.remap <- rbind(df.remap, c('Mixed','Mixed','#331000'))
  
  #Remapping
  df$color <- mapvalues(df$region, from = df.remap$full.name.parent, to = df.remap$color, warn_missing = FALSE)
  
  return(df)
}


toupper.first.letter <- function(x){
  return(paste0(toupper(substring(x, 1, 1)), 
                substring(x, 2, nchar(x))))
}
