#------------------ Librairies ------------------  

library(Rvcg)
library(rjson)
library(nat)

#------------------ Ontology handling ------------------ 

get.all.acronyms.children <- function(acronym){
  
  rec.func <- function(acronym, json, add = FALSE){
    
    if(is.null(json))
      return()
    
    if ((acronym == json$acronym | add))
      return(c(json$acronym, unlist(lapply(json$children, function(x){return(rec.func(acronym, x, TRUE))}))))
    
    if(!add)
      return(unlist(lapply(json$children, function(x){return(rec.func(acronym, x, FALSE))})))
    
  }
  
  #Loading json data
  json_data <- fromJSON(file= paste(path.matrices,'ontology.json',sep='/'))
  json <- json_data$msg[[1]]
  
  rec.func(acronym, json, FALSE)
  
}

get.acronym.from.id <- function(id){
  
  rec.func <- function(id, json){
    
    if(is.null(json))
      return()
    
    if(id == json$id)
      return(json$acronym)
    
    return(unlist(lapply(json$children, function(x){return(rec.func(id, x))})))
  }
  
  #Loading json data
  json_data <- fromJSON(file= paste(path.matrices,'ontology.json',sep='/'))
  json <- json_data$msg[[1]]
  
  return(rec.func(id, json))
  
}

get.name.from.id <- function(id){
  
  rec.func <- function(id, json){
    
    if(is.null(json))
      return()
    
    if(id == json$id)
      return(json$name)
    
    return(unlist(lapply(json$children, function(x){return(rec.func(id, x))})))
  }
  
  #Loading json data
  json_data <- fromJSON(file= paste(path.matrices,'ontology.json',sep='/'))
  json <- json_data$msg[[1]]
  
  return(rec.func(id, json))
  
}

get.id.from.acronym <- function(acronym){
  
  rec.func <- function(acronym, json){
    
    if(is.null(json))
      return()
    
    if(acronym == json$acronym)
      return(json$id)
    
    return(unlist(lapply(json$children, function(x){return(rec.func(acronym, x))})))
  }
  
  #Loading json data
  json_data <- fromJSON(file= paste(path.matrices,'ontology.json',sep='/'))
  json <- json_data$msg[[1]]
  
  return(rec.func(acronym, json))
}

#------------------ 3D handling ------------------ 

mesh3d.allen.annot.from.id <- function(id, max.nrows = 349696, bregma = c(5200, 650, 5700)/1000, no.normals = F){ #5400, 650, 5700 In mm, AP/DV/ML
  
  fpath <- sprintf('%s/%d.obj', allen.annot.path, id)
  
  vert.start <- 0
  vert.end <- 0
  norm.start <- 0
  norm.end <- 0
  face.start <- 0
  face.end <- 0
  i <- 0
  
  con <- file(fpath, 'r')
  while(TRUE){
    
    i <- i + 1
    line = readLines(con, n = 1)
    
    if (length(line) == 0)
      break
    
    if(startsWith(line, 'v ') & vert.start == 0)
      vert.start <- i
    
    if(vert.start != 0 & vert.end == 0 & !startsWith(line, 'v '))
      vert.end <- i
    
    if(startsWith(line, 'vn') & norm.start == 0)
      norm.start <- i
    
    if(norm.start !=0 & norm.end == 0 & !startsWith(line, 'vn'))
      norm.end <- i
    
    if(startsWith(line, 'f') & face.start == 0)
      face.start <- i
    
    if(face.start !=0 & face.end == 0 & !startsWith(line, 'f'))
      face.end <- i
    
  }
  
  if(face.end == 0)
    face.end <- i - 1
  
  close(con)
  
  con <- file(fpath, 'r')
  line = readLines(con, n = (vert.start - 1))
  line = readLines(con, n = vert.end - vert.start)
  vertices <- strsplit(line, ' ')
  vertices.ul <- unlist(vertices)
  df.vertices <- data.frame(x = as.numeric(vertices.ul[seq(from = 2, to = 4*length(vertices), by = 4)]),
                            y = as.numeric(vertices.ul[seq(from = 3, to = 4*length(vertices), by = 4)]),
                            z = as.numeric(vertices.ul[seq(from = 4, to = 4*length(vertices), by = 4)]))
  close(con)
  
  con <- file(fpath, 'r')
  line = readLines(con, n = (norm.start - 1))
  line = readLines(con, n = norm.end - norm.start)
  norm <- strsplit(line, ' ')
  norm.ul <- unlist(vertices)
  df.norm <-  data.frame(   x = as.numeric(norm.ul[seq(from = 2, to = 4*length(norm), by = 4)]),
                            y = as.numeric(norm.ul[seq(from = 3, to = 4*length(norm), by = 4)]),
                            z = as.numeric(norm.ul[seq(from = 4, to = 4*length(norm), by = 4)]))
  close(con)
  
  con <- file(fpath, 'r')
  line = readLines(con, n = (face.start - 1))
  line = readLines(con, n = face.end - face.start)
  face <- strsplit(line, ' ')
  face <- lapply(face, function(x){return(unlist(strsplit(x, '//')))})
  face.ul <- unlist(face)
  
  df.face <-  data.frame(   x = as.numeric(face.ul[seq(from = 2, to = 7*length(face), by = 7)]),
                            y = as.numeric(face.ul[seq(from = 4, to = 7*length(face), by = 7)]),
                            z = as.numeric(face.ul[seq(from = 6, to = 7*length(face), by = 7)]))
  close(con)
  
  if(!no.normals){
    mesh <- tmesh3d(vertices = rbind(t(df.vertices),1),
                    indices = t(df.face),
                    normals = t(df.norm)[c(3,2,1),])
  }else{
    mesh <- tmesh3d(vertices = rbind(t(df.vertices),1),
                    indices = t(df.face))
  }
  
  mesh$vb[1:3,] <- mesh$vb[1:3,] / 1000 
  
  mesh$vb[1,] <- -mesh$vb[1,] + bregma[1]
  mesh$vb[2,] <- -mesh$vb[2,] + 0.05 #Maybe not the +0.05
  mesh$vb[3,] <- -mesh$vb[3,] + bregma[3]
  
  tmp <- mesh$vb[1,]  
  mesh$vb[1,] <- mesh$vb[3,]
  mesh$vb[3,] <- tmp
  
  return(mesh)
  
}

mesh3d.cluster.from.id <- function(cl.list, df.colors){
  
  if(!exists('grid.3d') | !exists('grid.parameters'))
    load(grid.svm.path)
  
  adding.points.data <- NULL
  adding.points.color <- NULL
  cl.added <- NULL
  
  
  for(i in 1:length(cl.list)){
    
    cl <- cl.list[i]
    t <- array(data = 0, dim = c(dim(grid.3d)[1], dim(grid.3d)[2], dim(grid.3d)[3]))
    t[grid.3d == cl] <- 1
    
    #Not enough points
    if(sum(as.numeric(t)) == 0)
      next
    
    # browser()
    
    mesh <- NULL
    
    tryCatch({
      mesh <- vcgIsosurface(t,from = 0.5, to = 1.5, spacing = c(mean(abs(diff(grid.parameters$ML))),
                                                                -mean(abs(diff(grid.parameters$DV))),
                                                                -mean(abs(diff(grid.parameters$AP)))))
      },
      error = function(x){
        message("Mesh is null")
    })
    
    if(is.null(mesh)){
      next
    }

    mesh$vb[1,] <- mesh$vb[1,] + 1
    mesh$vb[2,] <- mesh$vb[2,] - 8
    mesh$vb[3,] <- mesh$vb[3,] + grid.parameters$AP[1]
    adding.points.data[[i]] <- mesh
    
    adding.points.color[[i]] <- df.colors[df.colors$cluster.id == cl, 'clusters.colors']
    cl.added <- c(cl.added, cl)
  }
  
  return(list(data = adding.points.data,
              color = adding.points.color,
              clusters = cl.added))
}

mesh3d.show.outline <- function(col = 'lightgray', alpha = 0.2){
  wire3d(mesh3d.allen.annot.from.id(get.id.from.acronym('root')), col = col, alpha = alpha)
}

mesh3d.new.window <- function(show.outline = TRUE){
  
  rgl.open()
  rgl.bg(color = 'white')
  par3d(windowRect = c(0, 0, 1920, 1080))
  
  if(show.outline)
    mesh3d.show.outline()
}




# #Brain outline
# get.brain.outline <- function(){
# 
#   indices.path <- 'C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/190211-postTreatment3D/indices.csv'
#   vertices.path <- 'C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/190211-postTreatment3D/vertices.csv'
#   
#   indices <- as.matrix(t(read.table(indices.path, sep = ',')))
#   indices <- indices[c(2,3,1),]
#   
#   vertices <- read.table(vertices.path, sep = ',')/100
#   vertices <- vertices[,c(2,3,1)]
#   vertices <- rbind(t(as.matrix(vertices)),1)
#   vertices[3,] <- - vertices[3,]
#   vertices[2,] <- vertices[2,] - 0.5
#   vertices[1,] <- vertices[1,] 
#   
#   tmesh.outline <- tmesh3d(vertices, indices)
#   
#   return(tmesh.outline)
# }
