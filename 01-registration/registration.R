# ---------------------- Defining functions ---------------------- 

#Windows only
quartz<-function(width,height){windows(width, height)}

#Returns the alignment data with coordinates of spots and their radius (only for spots on cell)
st.get.cy3 <- function(slice.parameters, Cy3.gray, HE.path,display.segmentation = TRUE, qc.plot = FALSE){
  
  #--------------------- Default parameters and structure definition ------------------
  
  alignment.path <- 'Alignment/selected_adjusted_spots_coordinates.txt'
  
  #Defining structures for the pictures
  HE <- list(path = HE.path,
             grayscale = character(),
             resolution = list(x = 0,y = 0),
             original = numeric()
  )
  
  Cy3 <- list(grayscale = Cy3.gray,
              resolution = list(x = 0,y = 0)
  )
  
  
  
  
  #Default filter working well in order to detect spots from Cy3 picture
  if (length(slice.parameters$feature.filter) == 0){
    slice.parameters$feature.filter <- structure(list(alim = c(120, 1000), threshold.range = c(14L, 255),
                                                      eccentricity = 1000L, Max = 70, Min = 0, brain.threshold = 1L, 
                                                      resize = 0.08*0.25, blur = 4L, downsample = 0.125/2),
                                                 .Names = c("alim","threshold.range", "eccentricity",
                                                            "Max", "Min", "brain.threshold","resize",
                                                            "blur", "downsample"))
  }
  
  #Creating directories
  dir.create('Registration', showWarnings = FALSE)
  dir.create('Registration/QC', showWarnings = FALSE)
  
  #--------------------- Detecting spots from Cy3 picture ------------------
  
  #Reading alignment data
  alignment <- read.table(alignment.path,col.names = c('id.x','id.y','raw.id.x','raw.id.y','x','y'))
  
  #Extracting spots trough wholeBrain segmentation
  spots <- segment(Cy3$grayscale, filter = slice.parameters$feature.filter, get.contour = TRUE, display = display.segmentation)
  
  #--------------------- Computing the equivalent radius of spots ------------------
  
  #Getting trough segmented spots to compute an equivalent radius
  for (i in sort(unique(spots$soma$contour.ID))){
    
    #Extracting the contour of specific spot
    poly.points <- t(rbind(spots$soma$contour.x[spots$soma$contour.ID == i],
                           spots$soma$contour.y[spots$soma$contour.ID == i]))
    
    #Extracting the center of specific spot and making it a matrix
    poly.center <- t(replicate(dim(poly.points)[1], c(spots$soma$x[i+1],spots$soma$y[i+1])))
    
    #Computing the euclidian distance between the points from contour and the center
    poly.dist <- sqrt(rowSums((poly.points - poly.center)^2))
    
    #Defining equivalent radius as the mean distance
    spots$soma$radius[i+1] <- mean(poly.dist)
    
  }
  
  #Creating a matrix to plot the detected contour
  #NA values are inserted between different spots to avoid connecting them together
  matrix.lines = matrix(data = NA, nrow = length(spots$soma$contour.ID) + length(unique(spots$soma$contour.ID)) - 1,ncol = 2)
  
  #Cursor for knowing the current row in the matrix
  curs <- 1
  
  #Going through all contours
  for (i in unique(spots$soma$contour.ID)){
    
    #Finding the x,y coordinates of this contour
    contour.x <- spots$soma$contour.x[spots$soma$contour.ID == i]
    contour.y <- spots$soma$contour.y[spots$soma$contour.ID == i]
    
    #Storing them at the current position in the matrix
    matrix.lines[curs:(curs+length(contour.x)-1),1] <- contour.x
    matrix.lines[curs:(curs+length(contour.x)-1),2] <- contour.y
    
    #Moving the cursor and make sure to let one row of NA values
    curs <- curs + length(contour.x) + 1
  }
  
  if (qc.plot){
    #Plotting real and estimated circle for visual control
    pdf('Registration/QC/circle-approximation.pdf')
    
    #In blue the theoretical circle
    symbols(spots$soma$x,spots$soma$y,circles = spots$soma$radius, inches=FALSE, asp = 1,fg='blue',xlim = c(0,32000),ylim = c(30000,1000),
            lwd = 0.5,xlab = 'x coordinates in Cy3 picture (pixels)',ylab = 'y coordinates in Cy3 picture (pixels)')
    
    #In red the real circles
    lines(matrix.lines,col = 'red',lwd = 0.15)
    
    dev.off()
  }
  
  #--------------------- Remapping the spots center from Cy3 to HE_cropped ------------------
  
  #Estimating dimensions of Cy3 picture
  Cy3$resolution$x <- max(spots$soma$x) - min(spots$soma$x) 
  Cy3$resolution$y <- max(spots$soma$y) - min(spots$soma$y)
  
  #Loading the HE picture to extract resolution
  library(jpeg)
  HE$original <- readJPEG(HE$path)
  HE$resolution$x <- dim(HE$original)[2]
  HE$resolution$y <- dim(HE$original)[1]
  
  #Saving this value
  slice.parameters$HE.resolution$x <- HE$resolution$x
  slice.parameters$HE.resolution$y <- HE$resolution$y
  
  #Estimating coefficient
  coeff.all <- c(HE$resolution$x / Cy3$resolution$x,HE$resolution$y / Cy3$resolution$y)
  coeff <- mean(coeff.all)
  
  #Formatting and printing a string for controlling coefficient values are reasonable.
  str <- sprintf('Estimated scale coefficient: %.3f, %.3f (values should be very similar). Final value: %.3f',
                 coeff.all[1],coeff.all[2],coeff)
  print(str)
  
  #Estimating top left point
  top.left <- numeric(2)
  top.left[1] <- min(spots$soma$x)
  top.left[2] <- min(spots$soma$y)
  
  #Remapping to HE picture
  spots$soma$x.remap  <- (spots$soma$x - top.left[1])*coeff
  spots$soma$y.remap <- (spots$soma$y - top.left[2])*coeff
  
  #Plotting the remapping for visual control
  #Red O: alignment position
  #Black X: remapping position
  if (qc.plot){
    pdf('Registration/QC/cy3-to-he-coordinates.pdf')
    plot(1, type = 'n',xlim = c(0,10000),ylim = c(10000,0),asp = 1,
         xlab = 'x coordinates in HE picture (pixels)',ylab = 'y coordinates in HE picture (pixels)')
    points(alignment$x,alignment$y,pch = 13,col = 'red',lwd = 0.5)
    points(spots$soma$x.remap,spots$soma$y.remap,pch = 4,lwd=0.5)
    dev.off()
  }
  #--------------------- Associating alignment with the closest detected spot ------------------
  
  distance.matrix <- matrix(nrow = dim(alignment)[1],ncol = length(spots$soma$x))
  
  for (i in 1:dim(alignment)[1]){
    
    #Creating column matrix where each row is the x-y coordinates of the current spot
    alignment.pos <- t(replicate(length(spots$soma$x),c(alignment$x[i],alignment$y[i])))
    
    #Computing euclidian distance of all nuclei to the current spot
    dist.align.to.spot <- sqrt(rowSums((cbind(spots$soma$x.remap,spots$soma$y.remap) - alignment.pos)^2))
    
    distance.matrix[i,] <- dist.align.to.spot
    
  }
  
  max.acceptable.dist <- 100
  
  
  #Getting the closest segmented spots for each aligned spots
  min.val <- apply(distance.matrix,1,min)
  min.id <- apply(distance.matrix,1,which.min)
  
  #Splitting aligned spots according to segmented spots being found close enough or not
  align.too.far <- which(min.val > max.acceptable.dist)
  align.valid <- setdiff(1:length(min.val),align.too.far)
  
  #Getting the ids of spots associated with a valid alignment
  spots.valid <- min.id[align.valid]
  
  #Pre-allocating
  #False if the spot could not be segmented
  alignment$segmented <- TRUE   
  
  #Estimated remapped coordinates when spot was segmented
  alignment$x.remap <- NA       
  alignment$y.remap <- NA
  
  #Radius (computed from segmented if possible, mean value of detected radius otherwise)
  alignment$radius <- NA
  
  #Saving the values for valid aligned spots
  alignment$x.remap[align.valid] <- spots$soma$x.remap[spots.valid]
  alignment$y.remap[align.valid] <- spots$soma$y.remap[spots.valid]
  alignment$radius[align.valid] <- spots$soma$radius[spots.valid] * coeff
  
  #Mean radius for not segmented spots and changing boolean
  alignment$radius[align.too.far] <- mean(alignment$radius, na.rm = TRUE)
  alignment$segmented[align.too.far] <- FALSE
  
  #Priting warning message
  if (length(align.too.far > 0)){
    str = ''
    for (i in 1:length(align.too.far)){
      str = paste(str,sprintf('\nxID=%d,yID=%d',alignment$id.x[align.too.far[i]],alignment$id.y[align.too.far[i]]),sep='')    
    }
    warning(sprintf('%d points were not properly segmented:%s.\nRadius was estimated by averaging segmented radius:%.1f.',length(align.too.far),str,mean(alignment$radius, na.rm = TRUE)),immediate. = TRUE)
  }
  
  #Error if segmented spots collide
  if (length(spots.valid) != length(unique(spots.valid))){
    stop('The same segmented spot is associated to multiple aligned spot. Something should be done')
  }
  
  if (qc.plot){
    
    #Plotting re-aligned points with color to check pairs are correct
    pdf('Registration/QC/spots-alignment-pairs.pdf')
    n <- dim(alignment)[1]
    rgbMat <- cbind(sample(1:255,n,replace=TRUE),sample(1:255,n,replace=TRUE),sample(1:255,replace=TRUE,n))  
    rgbMat <- rgbMat / 255
    plot(alignment$x[alignment$segmented == TRUE],alignment$y[alignment$segmented == TRUE],col = rgb(rgbMat[alignment$segmented==TRUE,]),asp=1,pch=13,xlim = c(0,10000),ylim = c(10000,0),
         xlab = 'x coordinates in HE picture (pixels)',ylab = 'y coordinates in HE picture (pixels)')
    points(alignment$x.remap,alignment$y.remap, col = rgb(rgbMat),pch = 3,lwd=0.5)
    points(alignment$x[alignment$segmented == FALSE],alignment$y[alignment$segmented == FALSE],col = 'red',bg = 'red',pch = 25)
    dev.off()
    
    #Plotting the alignement with circles on the picture
    pdf('Registration/QC/spots-on-picture.pdf')
    plot(1,type='n',xlim = c(0,HE$resolution$x), ylim = c(HE$resolution$y,0),asp = 1,
         xlab = 'x coordinates in HE picture (pixels)',ylab = 'y coordinates in HE picture (pixels)')
    rasterImage(as.raster(HE$original), xleft = 1, xright = HE$resolution$x, ybottom =  HE$resolution$y, ytop = 1,angle=0)
    symbols(alignment$x,alignment$y,circles = alignment$radius, inches=FALSE, asp = 1,fg='red',ylim = c(HE$resolution$y,0),add = TRUE, lwd = 0.2)
    dev.off()
    
    #Plotting circles with spots ids
    pdf('Registration/QC/spots-id.pdf')
    plot(alignment$id.x,alignment$id.y,asp=1,ylim=c(36,0),xlab='x spot id',ylab='y spot id',col='lightgray')
    text(alignment$id.x,alignment$id.y,labels=paste(alignment$id.x,alignment$id.y,sep='/'),cex=0.25)
    dev.off()
  }
  
  #Using a list because alignment needs to be returned, and slice.parameters was updated
  return(list(alignment = alignment,slice.parameters = slice.parameters))
  
}

#Returns the segmented nuclei from the cropped picture
st.get.nuclei <- function(slice.parameters, HE.gray, display.segmentation = TRUE){
  
  #Default filter for segmenting cells
  if (length(slice.parameters$cell.body.filter) == 0){
    slice.parameters$cell.body.filter <-structure(list(alim = c(3, 50), threshold.range = c(60,275L),
                                                       eccentricity = 500L, Max = 200, Min = 10, brain.threshold = 50L, 
                                                       resize = 0.05, blur = 4L, downsample = 0.25),
                                                  .Names = c("alim","threshold.range", "eccentricity",
                                                             "Max", "Min", "brain.threshold","resize", "blur", "downsample"))
  }
  
  #Calling the wholebrain function
  nuclei <- segment(HE.gray, filter = slice.parameters$cell.body.filter, display = display.segmentation)  
  
  #Keeping number of segmented cells
  slice.parameters$number.of.cells.segmented <- length(nuclei$soma$x)
  
  #Using a list because nuclei needs to be returned, and slice.parameters was updated
  return(list(nuclei = nuclei, slice.parameters = slice.parameters))
}

#Counts the number of nuclei per spots and add the associated spot to nuclei when relevant
st.add.nuclei.per.spot <- function(alignment, nuclei, qc.plot = FALSE){
  
  #Creating empty columns to fill
  alignment$nuclei <- numeric(dim(alignment)[1])
  nuclei$soma$spot.id.x <- numeric(length(nuclei$soma$x))
  nuclei$soma$spot.id.y <- numeric(length(nuclei$soma$y))
  
  #Going trough all spots
  for (i in 1:(dim(alignment)[1])){
    
    #Creating column matrix where each row is the x-y coordinates of the current spot
    spot.pos <- t(replicate(length(nuclei$soma$x),c(alignment$x[i],alignment$y[i])))
    
    #Radius of current sport
    r <- alignment$radius[i]
    
    #Computing euclidian distance of all nuclei to the current spot
    distance.to.spot <- sqrt(rowSums((cbind(nuclei$soma$x,nuclei$soma$y) - spot.pos)^2))
    
    #Getting id of nuclei in spots
    nuclei.in.spot <- which(distance.to.spot <= r)
    
    #Counting nuclei in spots
    alignment$nuclei[i] <- length(nuclei.in.spot)
    
  }
  
  #Creating QC picture with all nuclei and circles
  if (qc.plot){
    pdf('Registration/QC/spots-and-nuclei.pdf')
    plot(1, type = 'n',xlim = c(0,10000),ylim = c(10000,0),asp = 1,
         xlab = 'x coordinates in HE picture (pixels)',ylab = 'y coordinates in HE picture (pixels)')
    points(nuclei$soma$x,nuclei$soma$y,pch = 19,cex=0.05,col='gray')
    symbols(alignment$x,alignment$y,circles = alignment$radius, inches=FALSE, asp = 1,fg='red',ylim = c(HE$resolution$y,0),add = TRUE, lwd = 0.3)
    dev.off()
  }
  
  #Returning the updated alignment and nuclei files
  # return(list('alignment' = alignment, 'nuclei' = nuclei))  
  return(alignment)
  
}

#To run in order to preprocess the data
st.registration.preprocessing <- function(slice.parameters, HE.path){
  
  #Registering (required to recompute pictures)
  HE.gray.oriented <- rgb2gray(HE.path,rotate = slice.parameters$rotation)
  
  load('Registration/regi.RData')
  regi2 <- registration(HE.gray.oriented, coordinate = slice.parameters$bregma, filter=slice.parameters$cell.body.filter,
                        right.hemisphere=slice.parameters$right.hemisphere, correspondance = regi, display = FALSE)
  regi <- regi2
  save(regi, file='Registration/regi.RData')
  
  if (length(slice.parameters$subset) >= 2){
    
    for (i in 2:length(slice.parameters$subset)){
      regi.path <- sprintf('Registration/regi%d.RData',i)
      load(regi.path)
      regi2 <- registration(HE.gray.oriented, coordinate = slice.parameters$bregma, filter=slice.parameters$cell.body.filter,
                            right.hemisphere=slice.parameters$right.hemisphere, correspondance = regi, display = FALSE)
      regi <- regi2
      save(regi, file=regi.path)
    }
  }  
}

#Registering the spots into the atlasS
st.register.spots <- function(alignment, slice.parameters, qc.plot = FALSE){
  
  #Returns the spots in an appropriate format for using wholeBrain
  #Also, operates rotation to re-orient in the same way as the Cy3 picture
  get.spots.for.inspection <- function(alignment){
    
    if (slice.parameters$rotation == 90){
      spots = list('soma' = list('x' = slice.parameters$HE.resolution$y - alignment$y,'y'=alignment$x, 'intensity' = rep(1,length(alignment$y)),
                                 'area' = rep(1,length(alignment$y)), 'contour.x' = numeric(0),
                                 'contour.y' = numeric(0),'contour.ID' = NULL))
    }else if(slice.parameters$rotation == -90){
      spots = list('soma' = list('x' = alignment$y,'y'= slice.parameters$HE.resolution$x - alignment$x, 'intensity' = rep(1,length(alignment$y)),
                                 'area' = rep(1,length(alignment$y)), 'contour.x' = numeric(0),
                                 'contour.y' = numeric(0),'contour.ID' = NULL))
      
    }else if(slice.parameters$rotation == 0){
      spots = list('soma' = list('x' = alignment$x,'y'= alignment$y, 'intensity' = rep(1,length(alignment$y)),
                                 'area' = rep(1,length(alignment$y)), 'contour.x' = numeric(0),
                                 'contour.y' = numeric(0),'contour.ID' = NULL))   
    } 
    return(spots)
  }
  
  #Simple case
  if (length(slice.parameters$subset) == 0){
    
    #Loading the previsouly defined registration
    load('Registration/regi.Rdata')
    
    #Defining spots from alignment with an appropriate format for using wholebrain
    #Also rotating when registration do not correspond to the alignment
    spots <- get.spots.for.inspection(alignment)
    
    #Registering the spots into the atlas
    dataset <- inspect.registration(regi,spots,forward.warps = TRUE)
    
    if (qc.plot){
      dev.copy2pdf(file = "Registration/dataset.pdf")
    }
    
    dataset$radius <- alignment$radius
    dataset$nuclei <- alignment$nuclei
    dataset$id.x <- alignment$id.x
    dataset$id.y <- alignment$id.y
    dataset$float.id.x <- alignment$raw.id.x
    dataset$float.id.y <- alignment$raw.id.y
    dataset$x.wihtout.rotation <- alignment$x
    dataset$y.without.rotation <- alignment$y
    dataset$segmented <- alignment$segmented
    
  }else{
    
    dataset <- NULL
    cursor <- 1
    
    for (i in 1:length(slice.parameters$subset)){
      
      if (i == 1){
        #Loading the previsouly defined registration
        load('Registration/regi.Rdata')
      }else{
        load(sprintf('Registration/regi%d.RData',i))
      }
      
      cur.alignment <- alignment[slice.parameters$subset[[i]],]
      
      #Defining spots from alignment with an appropriate format for using wholebrain
      #Also rotating when registration do not correspond to the alignment
      spots <- get.spots.for.inspection(cur.alignment)
      
      #Registering the spots into the atlas
      dataset.tmp <- inspect.registration(regi,spots,forward.warps = TRUE)
      
      if (qc.plot){
        if (i == 1){
          #Loading the previsouly defined registration
          dev.copy2pdf(file ='Registration/dataset.pdf')
        }else{
          dev.copy2pdf(file = sprintf('Registration/dataset%d.pdf',i))
        }
      }
      
      dataset.tmp$radius <- cur.alignment$radius
      dataset.tmp$nuclei <- cur.alignment$nuclei
      dataset.tmp$id.x <- cur.alignment$id.x
      dataset.tmp$id.y <- cur.alignment$id.y
      dataset.tmp$float.id.x <- cur.alignment$raw.id.x
      dataset.tmp$float.id.y <- cur.alignment$raw.id.y
      dataset.tmp$x.wihtout.rotation <- cur.alignment$x
      dataset.tmp$y.without.rotation <- cur.alignment$y
      dataset.tmp$segmented <- cur.alignment$segmented
      dataset.tmp$subset <- i
      
      dataset <- rbind(dataset,dataset.tmp)
    }
  }
  
  #Adding ofset is required
  if (length(slice.parameters$ML.offset) != 0){
    dataset$ML <- dataset$ML + slice.parameters$ML.offset
  } 
  
  dataset$animal <- slice.parameters$slice
  
  return(dataset)
}

#Save the variables, and check if new file is identical to previous one and if previous file is lacking)
st.save.variables <- function(alignment.new, slice.parameters.new, dataset.new, df.equality.variables = data.frame(slice = '',stringsAsFactors = FALSE), idx = 1, save.intermediate.r.variables){
  
  al.path <- 'Registration/alignment.RData'
  sp.path <- 'Registration/sliceparameters.RData'
  da.path <- 'Registration/dataset.RData'
  
  #Already an alignment file
  if (file.exists(al.path)){
    
    #Loading it and uploading boolean
    load(al.path)
    df.equality.variables$alignment.exist[idx] <- TRUE
    
    #Checking if identical
    identical.bool <- identical(alignment, alignment.new)
    df.equality.variables$alignment.identical[idx] <- identical.bool
    
    #If not identical, saving old copy
    if (!identical.bool & save.intermediate.r.variables){
      save(alignment, file = 'Registration/alignment-old.RData')  
    }
    
    #Case without file
  }else{
    df.equality.variables$alignment.exist[idx] <- FALSE
  }
  
  alignment <- alignment.new
  
  #Already a dataset file
  if (file.exists(da.path)){
    
    #Loading it and uploading boolean
    load(da.path)
    df.equality.variables$dataset.exist[idx] <- TRUE
    
    #Checking if identical
    identical.bool <- identical(dataset[,2:dim(dataset)[2]], dataset.new[,2:dim(dataset)[2]])
    df.equality.variables$dataset.identical[idx] <- identical.bool
    
    #If not identical, saving old copy
    if (!identical.bool & save.intermediate.r.variables){
      save(dataset, file = 'Registration/dataset-old.RData')  
    }
    
    #Case without file
  }else{
    df.equality.variables$dataset.exist[idx] <- FALSE
  }  
  
  #Saving alignment (eventually similar as before)
  dataset <- dataset.new
  
  #Already an alignment file
  if (file.exists(sp.path)){
    
    #Loading it and uploading boolean
    load(sp.path)
    df.equality.variables$slice.parameters.exist[idx] <- TRUE
    
    #Checking if identical
    identical.bool <- identical(slice.parameters, slice.parameters.new)
    df.equality.variables$slice.parameters.identical[idx] <- identical.bool
    
    #If not identical, saving old copy
    if (!identical.bool & save.intermediate.r.variables){
      save(slice.parameters, file = 'Registration/sliceparameters-old.RData')  
    }
    
    #Case without file
  }else{
    df.equality.variables$slice.parameters.exist[idx] <- FALSE
  }
  
  slice.parameters <- slice.parameters.new
  
  
  df.equality.variables$slice[idx] <- slice.parameters$slice
  
  
  if (save.intermediate.r.variables){
    save(alignment, file = al.path) 
    save(slice.parameters, file = sp.path)
    save(dataset, file = da.path) 
  }
  
  
  
  return(df.equality.variables)
  
}

#Creating spots matrix
st.create.spots.matrix <- function(path.to.slice, dataset, path.to.index){
  
  table.index <- read.table(paste('../',path.to.index,sep=''), sep=';', stringsAsFactors = FALSE, header = TRUE)
  
  spots.table <- data.frame()
  st.data <- data.frame()
  number.spot.discarded <- 0
  st.data.row.names <- list()
  
  #Selecting relevant columns
  spots.table <- dataset[c('float.id.x','float.id.y','ML','DV','AP','acronym','nuclei','name','radius','x','y')]
  
  #Excluding spots that are outside the brain (other hemisphere, random tissue etc.)
  spots.table <- spots.table[!is.na(spots.table$acronym),]
  
  #Experiment number was put as animal name
  spots.table$slice.index <- table.index$id[table.index$slice.name == dataset$animal[1]]
  
  #If this was the left hemisphere, taking the opposite of ML coordinates to have all slices in the same hemisphere  
  if (is.element(FALSE,dataset$right.hemisphere)){
    spots.table$ML <- -spots.table$ML
  }
  
  path1 <- 'Data/stdata_inside.tsv'
  if (file.exists(path1)){
    st.data.file <- read.table(path1)
  }else{
    l.file <- list.files('Data', '-*(_stdata.tsv)$')
    if (length(l.file) != 1){
      stop('Problem in finding the proper expression file')
    }
    path2 <- sprintf('Data/%s',l.file)
    st.data.file <- read.table(path2)
  }
  
  #Parsing the row names from the st.data
  rnames <- row.names(st.data.file)
  rows.split <- strsplit(rnames,'x')
  rows.mat <- matrix(as.numeric(unlist(rows.split)),ncol = 2,byrow = TRUE)
  spots.st.coordinates <- data.frame(x = rows.mat[,1],y=rows.mat[,2])
  
  #Rows to keep in spots.append
  keep.row = logical(dim(spots.table)[1])
  
  #Ordered idx of st.data to match the spot vector
  st.data.idx.order <- NULL
  
  curs <- 0
  
  #Going through rows
  for (j in 1:dim(spots.table)[1]){
    
    #x/y spot id of current sport
    x = spots.table$float.id.x[j]
    y = spots.table$float.id.y[j]
    
    #Looking for the corresponding spot in genetic data (no similar rounding)
    corresponding.spot <- which(abs(spots.st.coordinates$x - x) < 0.1 & abs(spots.st.coordinates$y - y) < 0.1)
    
    #Two spots detected: raise an error
    if (length(corresponding.spot) > 1){
      stop('Only one spot should be detected')
      
      #Zero spot detected: might happen, just discarding the spot
    }else if(length(corresponding.spot) == 0){
      warning('One spot from alignment was not found in st data.')
      number.spot.discarded <- number.spot.discarded + 1
      
      #If everything is fine
    }else{
      keep.row[j] = TRUE
      curs <- curs+1
      st.data.idx.order[[curs]] <- corresponding.spot
    }
  }
  
  #Selecting only appropriate rows
  spots.table <- spots.table[keep.row,]
  st.data <- st.data.file[st.data.idx.order,]
  
  spots.id <- sprintf('%s_%s',table.index$id[table.index$slice.name == dataset$animal[1]],row.names(st.data))
  
  spots.table$spots.id <- spots.id 
  row.names(st.data) <- spots.id
  
  print(sprintf('%d spots discarded because the corresponding st data could not be found.',number.spot.discarded))
  
  rownames(spots.table) <- spots.table$spots.id
  spots.table <- spots.table[c('slice.index','ML','DV','AP','acronym','name','nuclei','radius','x','y')]
  colnames(spots.table)[1] <- 'slice_index'
  
  dir.create('Matrices', showWarnings = FALSE)
  
  write.table(spots.table, file='Matrices/spotstable.tsv', sep='\t', quote=FALSE, col.names = NA)
  write.table(st.data, file='Matrices/exprmat.tsv', sep='\t', quote=FALSE, col.names = NA)
  
}

#Execute all the instructions to register one slice
st.entirely.process.one.slice <- function(slice.dir, display.segmentation = FALSE, qc.plot = TRUE, df.equality.variables = data.frame(slice = '',stringsAsFactors = FALSE), idx = 1, save.intermediate.r.variables = FALSE, path.to.index){
  
  graphics.off()
  
  #Going into slice directory
  setwd(slice.dir)
  
  #Loading the file containing all information about the process
  load('Registration/sliceparameters.RData')
  
  #Defining relevant paths
  HE.path <- paste('HE_Cropped/HE_',slice.parameters$slice,'.jpg',sep='')
  files.Cy3 <- list.files('./Cy3') 
  if (length(files.Cy3) != 1){
    stop('There should be one and only one file in the Cy3 folder.')
  }
  Cy3.path <- paste('Cy3/',files.Cy3,sep='')
  
  #Getting the radius of the spots and aligning them
  Cy3.gray <- rgb2gray(Cy3.path,invert=FALSE, rotate = 180)
  list.output <- st.get.cy3(slice.parameters, Cy3.gray,HE.path, display.segmentation, qc.plot)
  alignment <- list.output$alignment
  slice.parameters <- list.output$slice.parameters
  rm(list.output)
  
  #Segmenting the nuclei
  HE.gray.original <- rgb2gray(HE.path)
  list.output <- st.get.nuclei(slice.parameters, HE.gray.original, display.segmentation)
  nuclei <- list.output$nuclei
  slice.parameters <- list.output$slice.parameters
  rm(list.output)
  
  #Adding labels of spots to nuclei and counting nuclei per spots
  alignment <- st.add.nuclei.per.spot(alignment, nuclei, qc.plot)
  
  #Processing pictures for registration to work in the next step
  st.registration.preprocessing(slice.parameters,HE.path)
  
  #Computing dataset by registering the spots
  dataset <- st.register.spots(alignment, slice.parameters, qc.plot)
  
  #Saving computed files
  df.equality.variables <- st.save.variables(alignment, slice.parameters, dataset, df.equality.variables, idx, save.intermediate.r.variables)
  
  #Creating the spots matrix
  st.create.spots.matrix(slice.dir, dataset, path.to.index)
  
  setwd('..')
  
  #Return a variable to check if the re-computed files are identical to the saved ones
  return(df.equality.variables)
}

#Might be windows only and RStudio only, returns the default path to all slices
st.get.default.path.slices <- function(){
  
  return(list.files(full.names = T, pattern = '^(ID).*'))
}

# ---------------------- Defining parameters ---------------------- 

#Loading required libraries
library(wholebrain)

#Chose the directory with the data for registration
setwd('data')

#Display intermediate windows
display.segmentation <- FALSE

#Saving 
save.intermediate.r.variables <- TRUE

#TRU for re-compuing the QC .pdf plots
qc.plot <- TRUE

#chr list with absolute paths to slices directories
path.to.slices <- st.get.default.path.slices()

#Index to convert slices from old indexing to new one
path.to.index <- 'sections/index-sections.csv'

#Object to check if re-computed data changed
df.equality.variables <- data.frame(slice = '',
                                    alignment.identical = logical(length(path.to.slices)),
                                    dataset.identical = logical(length(path.to.slices)),
                                    slice.parameters.identical = logical(length(path.to.slices)),
                                    alignment.exist = logical(length(path.to.slices)),
                                    dataset.exist = logical(length(path.to.slices)),
                                    slice.parameters.exist = logical(length(path.to.slices)),stringsAsFactors = FALSE)


#Looking through slices
for (idx in 1:length(path.to.slices)){
  df.equality.variables <- st.entirely.process.one.slice(path.to.slices[idx], display.segmentation, qc.plot, df.equality.variables, idx, save.intermediate.r.variables, path.to.index)
}


#------------- Merging matrices together ------------- 

save.as.RData = TRUE
save.as.tsv = TRUE

n.slices <- length(path.to.slices)

#Initializing lists that will contain the spots tabe and expression matrices of all slices
list.spots.tables <- NULL
list.expr.mat <- NULL

#Loading all matrices
for (i in 1:n.slices){
  
  list.spots.tables[[i]] <- read.table(paste(path.to.slices[i],'/Matrices/spotstable.tsv',sep='/'),
                                       sep='\t', header = TRUE, row.names = 1, quote="", stringsAsFactors = FALSE)
  
  list.expr.mat[[i]] <- read.table(paste(path.to.slices[i],'/Matrices/exprmat.tsv',sep='/'),
                                   sep='\t', header = TRUE, row.names = 1, quote="", stringsAsFactors = FALSE)
  
}

#Counting number of spots and genes
n.spots <- 0
n.genes <- 0

for (i in 1:n.slices){
  n.spots <- n.spots + dim(list.expr.mat[[i]])[1]
  n.genes <- n.genes + dim(list.expr.mat[[i]])[2]
}

#Creating vector with all genes (including duplicates) with preallocation
list.genes <- character(n.genes)

curs <- 1
for (i in 1:n.slices){
  list.genes[curs:(curs+dim(list.expr.mat[[i]])[2] -1)] <- colnames(list.expr.mat[[i]])
  curs <- curs + dim(list.expr.mat[[i]])[2]
}

#Keeping only unique genes, will be used as column names
list.genes <- unique(list.genes)
n.genes <- length(list.genes)

#Pre-allocating the spot table
spots.table <- data.frame(slice_index = character(n.spots),
                          ML = double(n.spots),
                          DV = double(n.spots),
                          AP = double(n.spots),
                          acronym = character(n.spots),
                          name = character(n.spots),
                          nuclei = integer(n.spots),
                          radius = double(n.spots),
                          x = integer(n.spots),
                          y = integer(n.spots),
                          stringsAsFactors = FALSE)

#Pre-allocating the expression matrix (st.data)
st.data <- data.frame(matrix(0.0, ncol = n.genes, nrow = n.spots))

#Naming columns by genes
colnames(st.data) <- list.genes

#Filling the st.data and spots.table objects
curs <- 1
for (i in 1:n.slices){
  cur.range <- curs:(curs+dim(list.spots.tables[[i]])[1] -1)
  spots.table[cur.range,] <- list.spots.tables[[i]]
  st.data[cur.range,colnames(list.expr.mat[[i]])] <- list.expr.mat[[i]]
  rownames(spots.table)[cur.range] <- rownames(list.spots.tables[[i]])
  curs <- curs + dim(list.spots.tables[[i]])[1]
}

#Using spots as row names for st.data
rownames(st.data) <- rownames(spots.table)


if (save.as.RData){
  save(st.data, file = 'exprmat.RData')
  save(spots.table, file='metatable.RData')
}

if (save.as.tsv){
  write.table(spots.table, file='metatable.tsv', sep='\t', quote=FALSE, col.names = NA)
  write.table(st.data, file='exprmat.tsv', sep='\t', quote=FALSE, col.names = NA)
}



