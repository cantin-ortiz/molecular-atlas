library(wholebrain)
library(plotly)
library(RColorBrewer)
## ----------- Glassbrain ----------- 

#Modified version of glassbrain to avoid having random AP coordinates
.glassbrain.modified <- function (dataset, high.res = FALSE, dim = c(0, 0, 720, 1080), device = TRUE,
                                 col = "region", cex = 0.5, hemisphere = "right", spheres = FALSE,
                                 laterality = TRUE, plane = "coronal")
{
  attach(loadNamespace('wholebrain'), warn.conflicts = FALSE)
  
  if (sum(dataset$color == "#000000") > 0) {
    dataset <- dataset[-which(dataset$color == "#000000"),
                       ]
  }
  if (device) {
    open3d(windowRect = dim)
  }
  par3d(persp)
  if (high.res) {
    drawScene.rgl(list(VOLUME))
  }
  else {
    drawScene.rgl(list(VOLUMESMALL))
  }
  
  if(is.null(dataset))
    return()
  
  if (col == "region") {
    color <- as.character(dataset$color)
  }
  else {
    color <- col
  }
  if (hemisphere == "right") {
    hemisphere <- 2
  }
  else {
    hemisphere <- 1
  }
  
  smp.AP <- 0
  smp.ML <- 0
  
  if (laterality) {
    if (length(unique(dataset$AP)) > 1) {
      laterality <- table(dataset$AP, dataset$right.hemisphere)
      for (i in 1:nrow(laterality)) {
        if (hemisphere == "right") {
          if (laterality[i, 1] > laterality[i, 2]) {
            dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])] <- -dataset$ML[which(dataset$AP ==
                                                                                                         as.numeric(row.names(laterality))[i])]
          }
        }
        else {
          if (laterality[i, 2] > laterality[i, 1]) {
            dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])] <- -dataset$ML[which(dataset$AP ==
                                                                                                         as.numeric(row.names(laterality))[i])]
          }
        }
      }
    }
  }
  if (spheres) {
    spheres3d(paxTOallen(dataset$AP) - 530/2 + smp.AP, -dataset$DV *
                1000/25 - 320/2, dataset$ML * 1000/25 + smp.ML, col = color,
              radius = cex, alpha = dataset$alpha)
  }
  else {
    points3d(paxTOallen(dataset$AP) - 530/2 + smp.AP, (-dataset$DV *
                                                         1000/25 * 0.95) - 320/2, dataset$ML * 1000/25 + smp.ML,
             col = color, size = cex, dataset$alpha)
  }
}

#Function to plot a spot table with 3D with glassbrain outline, only right hemisphere
.plot.3d.glassbrain.col <- function(spots.table, cex = 1.5, HD = FALSE, dim = c(0, 0, 720, 1080)){
  
  #spots.table
  #cex(num): sphere cex
  #HD(logical): resolution of the outline
  #dim: dimensions of the windon
  .glassbrain.modified(spots.table, laterality = FALSE, spheres = TRUE, cex = cex, high.res = HD, dim = dim)
  
}

plot.3d.glassbrain <- function(spots.table, clusters = NULL, view = '3d', cex = 1.5, HD = FALSE, col.vect = NULL, DV.offset = 0, ML.offset = 0){
  library(wholebrain)
  
  if (is.null(clusters))
    clusters <- levels(spots.table$clusters.named)
  
  if (is.null(col.vect)){
    col.tab.paired <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff00','#b15928')
    col.vect <- rep(col.tab.paired,40)
    col.vect <- col.vect[1:length(clusters)]    
  }
  spots.plot <- spots.table[is.element(spots.table$clusters.named,clusters),]
  spots.plot$color <- mapvalues(spots.table$clusters.named,from = levels(spots.table$clusters.named), to = col.vect)
  spots.plot$right.hemisphere <- TRUE
  spots.plot$DV <- spots.plot$DV - DV.offset
  spots.plot$ML <- spots.plot$ML - ML.offset
  glassbrain(spots.plot,laterality = FALSE, spheres = TRUE, cex = cex, high.res = HD, dim = dim)
  plot.3d.glassbrain.setview(view)
  
}

plot.3d.expr.glassbrain <- function(gene, spots.plot, st.data, view = '3d', cex = 1.5, HD = FALSE, 
                                    min.expr = 0.8, log.scale = FALSE, DV.offset = 0, ML.offset = 0,
                                    dim = c(0, 0, 720, 1080)){

  
  spots.plot$expr <- st.data[rownames(spots.plot),gene]
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('white','red'))
  #This adds a column of color values
  # based on the y values
  if (min.expr >= 1)
    spots.plot <- spots.plot[spots.plot$expr >= min.expr,]
  else{
    spots.plot <- spots.plot[(spots.plot$expr > as.numeric(quantile(spots.plot$expr,min.expr))) == TRUE,]
  }
  
  if (log.scale){
    spots.plot$expr[spots.plot$expr == 0] <- 0.98
    spots.plot$expr <- spots.plot$expr + 0.01
    spots.plot$expr <- ceiling(5*log(spots.plot$expr))
  }
  
  spots.plot$color <- rbPal(max(spots.plot$expr)+1)[as.numeric(spots.plot$expr)+1]  
  spots.plot$right.hemisphere <- TRUE
  spots.plot$DV <- spots.plot$DV - DV.offset
  spots.plot$ML <- spots.plot$ML - ML.offset

  # glassbrain(spots.plot,laterality = FALSE, spheres = TRUE, cex = cex, high.res = HD, dim = dim)
  .plot.3d.glassbrain.col(spots.plot, dim = dim, HD = HD, cex = cex)
  plot.3d.glassbrain.setview(view)
  
}

plot.3d.glassbrain.setview <- function(view = 'dorsal', zoom = 0.6){
  
  if (view == 'dorsal')
    rgl.viewpoint(userMatrix = rotationMatrix(-pi/2, 1, 0, 0), zoom = zoom)
  else if(view =='ventral')
    rgl.viewpoint(userMatrix = rotationMatrix(pi/2, 1, 0, 0), zoom = zoom)
  else if(view == 'posterior')
    rgl.viewpoint(userMatrix = rotationMatrix(pi, 1, 0, 1), zoom = zoom)
  else if(view == 'anterior')
    rgl.viewpoint(userMatrix = rbind(c(0,0,-1,0),c(0,-1,0,0),c(-1,0,0,0),c(0,0,0,1)), zoom = zoom)
  else if(view == 'medial')
    rgl.viewpoint(userMatrix = rotationMatrix(pi, 1, 0, 0), zoom = zoom)
  else if(view == 'lateral')
    rgl.viewpoint(userMatrix = rbind(c(-1,0,0,0),c(0,-1,0,0),c(0,0,1,0),c(0,0,0,1)), zoom = zoom)
  else if(view == '3d')
    rgl.viewpoint(userMatrix = rbind(c(0.563,0.15,-0.821,0),c(0.307,-0.948,0,0),c(-0.77,-0.27,-0.58,0),c(0,0,0,1)), zoom = zoom)
  else if(view == '3d2')
    rgl.viewpoint(userMatrix = rbind(c(0.434,0.05,0.899,0),c(-0.05,-0.992,0.04,0),c(0.90,-0.08,-0.45,0),c(0,0,0,1)), zoom = zoom)
}

## ----------- 3D Plotly ----------- 

plot.3d.all <- function(spots.table, col.vect = NULL, title = ''){
  library(plotly)
  library(RColorBrewer)
  
  if (is.null(col.vect)){
    col.base <- colorRampPalette(brewer.pal(9, "Set1"))(18)
    col.base <- rep(col.base,40)
    col.vect <- col.base[1:length(levels(spots.table$cluster))]
    # col.vect <- colorRampPalette(brewer.pal(9, "Set1"))(length(levels(spots.table$clusters.named)))
  }
  
  updatemenus = list(list(active = -1, type = 'buttons', x = 0.9, y = 1, buttons = list()))
  updatemenus[[1]]$buttons[[1]] <- list(label = 'Hide all', method = 'update',args = list(list(visible = rep('legendonly',length(levels(spots.table$clusters.named))))))
  updatemenus[[1]]$buttons[[2]] <- list(label = 'Show all', method = 'update', args = list(list(visible = rep(TRUE,length(levels(spots.table$clusters.named))))))
  
  
  p <- plot_ly(type = 'scatter3d', mode='markers',data = spots.table, color=~clusters.named, x=~ML, y=~AP, z=~DV, colors = col.vect,
               marker=list(size=5), visible = 'true', showlegend = TRUE, hoverinfo='text', 
               text = ~sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.3f<br />Region: %s/%s<br />Cluster: %s',
                               rownames(spots.table),DV,ML,AP,full.name.parent,name,clusters.named)) %>%

    
        layout(scene = list(yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3)),
                            zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0), scaleanchor = "y"),
                            xaxis = list(title = 'x = Medial-lateral (mm)', range = c(0,5), scaleanchor = "y"),
                            aspectratio = list(x=(5/8.9), y=1, z=(7.9/8.9))),             
               title = title,
               legend = list(traceorder = 'normal'),
               updatemenus = updatemenus)
  
  return(p)
}

#Generate a plotly 3d plot from a clustering file
plot.3d.all.clust.file <- function(cluster.file, spots.table = NULL){
  
  #Reading the cluster file
  t <- read.table(cluster.file)
  
  #First column: spot id
  t$V1 <- as.character(t$V1)
  
  #Loading the spots table if not provided
  if (is.null(spots.table)){
    spots.table <- load.spots.table()
    spots.table <- add.parent.acronym(spots.table)
  }
    
  #Only selecting clustered spots
  spots.table.cur <- spots.table[t$V1,]

  #Appending the cluster
  spots.table.cur[t$V1,'cluster'] <- factor(t$V2)
  
  #Reording clusters for proper renaming
  df <- plyr::count(spots.table.cur$cluster)
  df <- df[order(df$freq,decreasing = TRUE),]
  df$new.clust <- 0:(dim(df)[1]-1)

  spots.table.cur$cluster <- mapvalues(spots.table.cur$cluster, from = df$x, to = df$new.clust)
  spots.table.cur$cluster <- factor(spots.table.cur$cluster, levels = df$new.clust)
  spots.table.cur <- name.clusters(spots.table.cur)

  #Defining vector color
  col.base <- colorRampPalette(brewer.pal(9, "Set1"))(18)
  col.base <- rep(col.base,40)
  col.vect <- col.base[1:length(levels(spots.table.cur$clusters.named))]
  
  #Plotting
  p <- plot.3d.all(spots.table.cur,col.vect, title = paste(length(levels(spots.table.cur$cluster)),'clusters',sep=' '))

  return(p)
}

plot.3d.all.shaded <- function(spots.table){
  
  library(plotly)
  load('C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/180329-boundaries/vertices.RData')
  
  l.vertices[[3]] <- NULL
  l.vertices[[13]] <- NULL
  
  boolean.list <- rep(TRUE,length(l.vertices) + length(levels(spots.table$clusters.named)))
  
  updatemenus = list(list(active = -1, type = 'buttons', buttons = list()))
  updatemenus[[1]]$buttons[[1]] <- list(label = 'All', method = 'update',args = list(list(visible = rep(TRUE,length(l.vertices)+3))))
  updatemenus[[1]]$buttons[[2]] <- list(label = 'None', method = 'update', args = list(list(visible = boolean.list)))
  
  
  for (i in 1:length(l.vertices)){
    bool.update <- boolean.list
    bool.update[i] <- TRUE
    
    updatemenus[[1]]$buttons[[i+2]] <- list(label = as.character(name.from.acronym(l.vertices[[i]]$acornym)),
                                            method = "update",
                                            args = list(list(visible = bool.update)))
  }
  
  
  p <- plot_ly(type = 'mesh3d')
  
  for (i in 1:length(l.vertices)){
    p <- add_trace(p,type = 'mesh3d', visible = TRUE, data = l.vertices[[i]],x=~vertices$ML,y=~vertices$AP,z = ~vertices$DV, name = ~name.from.acronym(acornym),
                   i = ~indices$V1-1, j=~indices$V2-1, k=~indices$V3-1, inherit = FALSE, opacity = 0.1,facecolor = ~rep(toRGB(color),1,length(indices$V1)), hoverinfo = "name", showlegend = FALSE)
  }
  
  
  #Plotting
  p <- add_markers(p,type = 'scatter3d', mode='markers',data = spots.table, color=~clusters.named, x=~ML, y=~AP, z=~DV, colors = col.vect,
                   marker=list(size=5), visible = 'true', showlegend = TRUE, hoverinfo='text', 
                   text = ~sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.3f<br />Region: %s/%s<br />Cluster: %s',
                                   rownames(spots.table),DV,ML,AP,full.name.parent,name,clusters.named))
  
  
  p <- layout(p, updatemenus=updatemenus,
              scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                           zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                           yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))))
  
  return(p)
  
}

plot.3d.specific.clusters <- function(spots.table, clusters.list){
  
  library(RColorBrewer)
  spots.cur <- spots.table[is.element(spots.table$clusters.named,clusters.list), ]
  
  cols <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(spots.cur$clusters.named)))
  
  spots.cur$clusters.named <- factor(spots.cur$clusters.named)
  
  spots.cur$col <- mapvalues(as.character(spots.cur$clusters.named),from = levels(spots.cur$clusters.named),to = cols)
  
  p <- plot_ly(type = 'scatter3d', mode='markers',data = spots.cur, color=~clusters.named, x=~ML, y=~AP, z=~DV, colors = cols,
               marker=list(size=5), visible = 'true', showlegend = TRUE, hoverinfo='text', 
               text = ~sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.3f<br />Region: %s/%s<br />Cluster: %s',
                               rownames(spots.cur),DV,ML,AP,full.name.parent,name,clusters.named)) %>%
    
    layout(scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                        zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                        yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))),
           legend = list(traceorder = 'normal'))
  
  list.output<- list(p = p, spots.cur = spots.cur)
  return(list.output)
  
}

plot.3d.specific.area <- function(spots.table, acronym.list, min.spots = 0){
  
  clust <- NULL
  for (acr in acronym.list)
    clust <- c(clust,as.character(spots.table[(spots.table$acronym == acr) | (spots.table$acronym.parent == acr),"clusters.named"]))
  
  spots.table$clusters.named[spots.table$acronym == 'MOp1']
  
  clust <- unique(clust)
  plot.3d.specific.clusters(spots.table,clust)
  
}

plot.3d.expr <- function(gene, spots.table, st.data, normalize.by.nuclei = FALSE, change.opacity = TRUE){
  
  spots.table$expr.level <- st.data[,gene]
  
  spots.plot <- spots.table
  
  # #Discarding outside of brain spots
  # print(sprintf('%d spots discarded because they are outside the brain contour.',sum(is.na(spots.plot$acronym))))
  # spots.plot <- spots.table[!is.na(spots.plot$acronym),]
  # 
  # #Discarding spots for which gene was not sequenced
  # print(sprintf('%d spots discarded because the selected gene was not sequenced.',sum(is.na(spots.plot$expr.level))))
  # spots.plot <- spots.plot[!is.na(spots.plot$expr.level),]
  # 
  
  str.title <- ''
  if (normalize.by.nuclei == TRUE){
    spots.plot$nuclei[spots.plot$nuclei == 0] <- 1
    spots.plot$expr.level <- spots.plot$expr.level / spots.plot$nuclei
    str.title <- '(Normalized by number of nuclei)'
  }
  
  #Top and lower quantile
  if (is.null(change.opacity)){
    
    spots.plot.1 <- spots.plot[spots.plot$expr.level >= quantile(spots.plot$expr.level,0.9),]
    spots.plot.2 <- spots.plot[spots.plot$expr.level <= quantile(spots.plot$expr.level, 0.1),]
    spots.plot.3 <- spots.plot[setdiff(setdiff(rownames(spots.plot),rownames(spots.plot.1)),rownames(spots.plot.2)),]

    
    #Plotting
    p <- plot_ly(spots.plot.1,x=~ML,y=~AP,z=~DV,type = 'scatter3d',mode = 'markers',name='High expression',
                 marker=list(size=2, color = ~expr.level, colorscale = c('#3300ff', '#683531'), showscale = TRUE,
                             opacity = 1, cmin = min(spots.plot$expr.level), cmax = max(spots.plot$expr.level)),
                 hoverinfo = 'text',
                 text = ~sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                 DV,ML,AP,full.name.parent,name,expr.level)) %>%
      
      add_markers(x=spots.plot.2$ML,y=spots.plot.2$AP,z=spots.plot.2$DV,name=sprintf('Low expression'),
                  marker=list(size=2, color = spots.plot.2$expr.level,colorscale =  c('#3300ff', '#683531'), showscale = TRUE,
                              opacity = 1,cmin = min(spots.plot$expr.level),cmax = max(spots.plot$expr.level),name='Low expression') ,
                  hoverinfo = 'text',
                  text = sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                 spots.plot.2$DV,spots.plot.2$ML,spots.plot.2$AP,spots.plot.2$full.name.parent,spots.plot.2$name,spots.plot.2$expr.level))%>%
      
      add_markers(x=spots.plot.3$ML,y=spots.plot.3$AP,z=spots.plot.3$DV,name=sprintf('Average expression'),
                  marker=list(size=2, color = spots.plot.3$expr.level, colorscale =  c('#3300ff', '#683531'), showscale = TRUE,
                              opacity = 0.5,cmin = min(spots.plot$expr.level),cmax = max(spots.plot$expr.level)) ,
                  hoverinfo = 'text',
                  text = sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                 spots.plot.3$DV,spots.plot.3$ML,spots.plot.3$AP,spots.plot.3$full.name.parent,spots.plot.3$name,spots.plot.3$expr.level))%>%
      
      
      layout(title=sprintf("Gene: %s %s",gene,str.title),
             scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                          zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                          yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))))    
    
    
  }else if (change.opacity){
  
    change.opacity <- 0.1*max(spots.plot$expr.level)
    
    spots.plot.1 <- spots.plot[spots.plot$expr.level >= change.opacity,]
    spots.plot.2 <- spots.plot[spots.plot$expr.level < change.opacity,]
    
    #Plotting
    p <- plot_ly(spots.plot.1,x=~ML,y=~AP,z=~DV,type = 'scatter3d',mode = 'markers',name='High expression',
                 marker=list(size=2, color = ~expr.level,colorscale = c('#FFE1A1', '#683531'), showscale = TRUE,
                             opacity = 1,cmin = min(spots.plot$expr.level),cmax = max(spots.plot$expr.level)),
                 hoverinfo = 'text',
                 text = ~sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                 DV,ML,AP,full.name.parent,name,expr.level)) %>%
      
      add_markers(x=spots.plot.2$ML,y=spots.plot.2$AP,z=spots.plot.2$DV,name=sprintf('Low expression (<%0.f)',change.opacity),
                  marker=list(size=2, color = spots.plot.2$expr.level,colorscale =  c('#FFE1A1', '#683531'), showscale = TRUE,
                              opacity = 0.05,cmin = min(spots.plot$expr.level),cmax = max(spots.plot$expr.level),name='Low expression') ,
                  hoverinfo = 'text',
                  text = sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                 spots.plot.2$DV,spots.plot.2$ML,spots.plot.2$AP,spots.plot.2$full.name.parent,spots.plot.2$name,spots.plot.2$expr.level))%>%
      
      layout(title=sprintf("Gene: %s %s",gene,str.title),
             scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                          zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                          yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))))
  }else{
    p <- plot_ly(spots.plot,x=~ML,y=~AP,z=~DV,type = 'scatter3d',mode = 'markers',name='High expression',
                 marker=list(size=2, color = ~expr.level,colorscale = c('#FFE1A1', '#683531'), showscale = TRUE,
                             opacity = 0.5,cmin = min(spots.plot$expr.level),cmax = max(spots.plot$expr.level)),
                 hoverinfo = 'text',
                 text = ~sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                 DV,ML,AP,full.name.parent,name,expr.level)) %>%
      
       layout(title=sprintf("Gene: %s %s",gene,str.title),
             scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                          zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                          yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))))    
  }
  
  p
  return(p)
  
}

plot.3d.ica <- function(spots.table, ica.matrix, ica.vect, quantile.ratio = 0.1, absolute.quantile = TRUE){
  
  spots.plot <- spots.table[intersect(rownames(spots.table),rownames(ica.matrix)),]
  ica.matrix <- ica.matrix[intersect(rownames(spots.table),rownames(ica.matrix)),]
  spots.plot$expr.level <- ica.matrix[,ica.vect]
  
  # #Discarding outside of brain spots
  # print(sprintf('%d spots discarded because they are outside the brain contour.',sum(is.na(spots.plot$acronym))))
  # spots.plot <- spots.table[!is.na(spots.plot$acronym),]
  # 
  # #Discarding spots for which gene was not sequenced
  # print(sprintf('%d spots discarded because the selected gene was not sequenced.',sum(is.na(spots.plot$expr.level))))
  # spots.plot <- spots.plot[!is.na(spots.plot$expr.level),]
  # 
  
  str.title <- ''
  
  if (absolute.quantile == FALSE){

    lower.quantile <- quantile(spots.plot$expr.level,quantile.ratio)
    higher.quantile <- quantile(spots.plot$expr.level,1-quantile.ratio)
  
    
    ma <- max(max(spots.plot$expr.level),-min(spots.plot$expr.level))
    
    spots.plot.1 <- spots.plot[spots.plot$expr.level >= higher.quantile,]
    spots.plot.2 <- spots.plot[spots.plot$expr.level <= lower.quantile,]
    spots.plot.3 <- spots.plot[spots.plot$expr.level < higher.quantile & spots.plot$expr.level > lower.quantile,]
    
    #Plotting
    p <- plot_ly(spots.plot.1,x=~ML,y=~AP,z=~DV,type = 'scatter3d',mode = 'markers',name='Top quantile',
                 marker=list(size=2, color = ~expr.level,colorscale = c('#FFE1A1', '#683531'), showscale = TRUE,
                             opacity = 1,cmin = -ma,cmax = ma),
                 hoverinfo = 'text',
                 text = ~sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.2f',
                                 rownames(spots.plot.1),DV,ML,AP,full.name.parent,name,expr.level)) %>%
      
      
      add_markers(x=spots.plot.2$ML,y=spots.plot.2$AP,z=spots.plot.2$DV,name='Lower quantile',
                  marker=list(size=2, showscale = TRUE, color = spots.plot.2$expr.level,colorscale = c('#FFE1A1', '#683531'), 
                              opacity = 1,cmin = -ma,cmax = ma),
                  hoverinfo = 'text',
                  text = sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.2f',
                                 rownames(spots.plot.2),spots.plot.2$DV,spots.plot.2$ML,spots.plot.2$AP,spots.plot.2$full.name.parent,spots.plot.2$name,spots.plot.2$expr.level))%>%
      
      
      add_markers(x=spots.plot.3$ML,y=spots.plot.3$AP,z=spots.plot.3$DV,name='Average expression',
                  marker=list(size=2, showscale = TRUE, color = spots.plot.3$expr.level,colorscale = c('#FFE1A1', '#683531'), 
                              opacity = 0.1,cmin = -ma,cmax = ma),
                  hoverinfo = 'text',
                  text = sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.2f',
                                 rownames(spots.plot.3),spots.plot.3$DV,spots.plot.3$ML,spots.plot.3$AP,spots.plot.3$full.name.parent,spots.plot.3$name,spots.plot.3$expr.level))%>%
  
             layout(title=sprintf("Ica vector: %s %s",ica.vect,str.title),
                    scene = list(yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3)),
                                 zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0), scaleanchor = "y"),
                                 xaxis = list(title = 'x = Medial-lateral (mm)', range = c(0,5), scaleanchor = "y"),
                                 aspectratio = list(x=(5/8.9), y=1, z=(7.9/8.9))))       
  }
  else{
    
    higher.quantile <- quantile(abs(spots.plot$expr.level),1-quantile.ratio)
    
    
    ma <- max(max(spots.plot$expr.level),-min(spots.plot$expr.level))
    
    spots.plot.1 <- spots.plot[spots.plot$expr.level >= higher.quantile | spots.plot$expr.level <= -higher.quantile,]
    spots.plot.3 <- setdiff(spots.plot,spots.plot.1)
    
    #Plotting
    p <- plot_ly(spots.plot.1,x=~ML,y=~AP,z=~DV,type = 'scatter3d',mode = 'markers',name='Top quantile',
                 marker=list(size=2, color = ~expr.level,colorscale = c('#FFE1A1', '#683531'), showscale = TRUE,
                             opacity = 1,cmin = -ma,cmax = ma),
                 hoverinfo = 'text',
                 text = ~sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.2f',
                                 rownames(spots.plot.1),DV,ML,AP,full.name.parent,name,expr.level)) %>%
      
      
      add_markers(x=spots.plot.3$ML,y=spots.plot.3$AP,z=spots.plot.3$DV,name='Average expression',
                  marker=list(size=2, showscale = TRUE, color = spots.plot.3$expr.level,colorscale = c('#FFE1A1', '#683531'), 
                              opacity = 0.1,cmin = -ma,cmax = ma),
                  hoverinfo = 'text',
                  text = sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.2f',
                                 rownames(spots.plot.3),spots.plot.3$DV,spots.plot.3$ML,spots.plot.3$AP,spots.plot.3$full.name.parent,spots.plot.3$name,spots.plot.3$expr.level))%>%
      
      layout(title=sprintf("Ica vector: %s %s",ica.vect,str.title),
             scene = list(yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3)),
                          zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0), scaleanchor = "y"),
                          xaxis = list(title = 'x = Medial-lateral (mm)', range = c(0,5), scaleanchor = "y"),
                          aspectratio = list(x=(5/8.9), y=1, z=(7.9/8.9))))           
    
    
  }
 
  return(p)
  
}

plot.3d.expr.contour <- function(gene, spots.table, st.data, normalize.by.nuclei = FALSE){
  load('C:/Users/MatLab/Desktop/transcripBrainAtlas/exp/180329-boundaries/vertices.RData')
  
  l.vertices[[3]] <- NULL
  l.vertices[[13]] <- NULL
  
  boolean.list <- rep(FALSE,length(l.vertices)+3)
  boolean.list[14:15] <- TRUE
  
  updatemenus = list(list(active = -1, type = 'buttons', buttons = list()))
  updatemenus[[1]]$buttons[[1]] <- list(label = 'All', method = 'update',args = list(list(visible = rep(TRUE,length(l.vertices)+3))))
  updatemenus[[1]]$buttons[[2]] <- list(label = 'None', method = 'update', args = list(list(visible = boolean.list)))
  
  
  for (i in 1:length(l.vertices)){
    bool.update <- boolean.list
    bool.update[i+1] <- TRUE
    
    updatemenus[[1]]$buttons[[i+2]] <- list(label = as.character(name.from.acronym(l.vertices[[i]]$acornym)),
                                            method = "update",
                                            args = list(list(visible = bool.update)))
  }
  
  
  p <- plot_ly(type = 'mesh3d')
  
  for (i in 1:length(l.vertices)){
    p <- add_trace(p,type = 'mesh3d', visible = TRUE, data = l.vertices[[i]],x=~vertices$ML,y=~vertices$AP,z = ~vertices$DV, name = ~name.from.acronym(acornym),
                   i = ~indices$V1-1, j=~indices$V2-1, k=~indices$V3-1, inherit = FALSE, opacity = 0.1,facecolor = ~rep(toRGB(color),1,length(indices$V1)), hoverinfo = "name", showlegend = FALSE)
  }
  
  
  spots.table$expr.level <- st.data[,gene]
  
  spots.plot <- spots.table
  
  #Discarding outside of brain spots
  print(sprintf('%d spots discarded because they are outside the brain contour.',sum(is.na(spots.plot$acronym))))
  spots.plot <- spots.table[!is.na(spots.plot$acronym),]
  
  #Discarding spots for which gene was not sequenced
  print(sprintf('%d spots discarded because the selected gene was not sequenced.',sum(is.na(spots.plot$expr.level))))
  spots.plot <- spots.plot[!is.na(spots.plot$expr.level),]
  
  
  str.title <- ''
  if (normalize.by.nuclei == TRUE){
    spots.plot$nuclei[spots.plot$nuclei == 0] <- 1
    spots.plot$expr.level <- spots.plot$expr.level / spots.plot$nuclei
    str.title <- '(Normalized by number of nuclei)'
  }
  
  change.opacity <- 0.1*max(spots.plot$expr.level)
  
  spots.plot.1 <- spots.plot[spots.plot$expr.level >= change.opacity,]
  spots.plot.2 <- spots.plot[spots.plot$expr.level < change.opacity,]
  
  #Plotting
  p <- add_markers(p, data=spots.plot.1,x=~ML,y=~AP,z=~DV,type = 'scatter3d',mode = 'markers',name='High expression',
                   marker=list(size=2, color = ~expr.level,colorscale = c('#FFE1A1', '#683531'), showscale = TRUE,
                               opacity = 1,cmin = 0,cmax = max(spots.plot$expr.level)),
                   hoverinfo = 'text',
                   text = ~sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                   DV,ML,AP,full.name.parent,name,expr.level)) 
  
  p <- add_markers(p, x=spots.plot.2$ML,y=spots.plot.2$AP,z=spots.plot.2$DV,name=sprintf('Low expression (<%0.f)',change.opacity),
                   marker=list(size=2, color = ~expr.level,colorscale =  c('#FFE1A1', '#683531'), showscale = TRUE,
                               opacity = 0.05,cmin = 0,cmax = max(spots.plot$expr.level),name='Low expression') ,
                   hoverinfo = 'text',
                   text = sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %.0f',
                                  spots.plot.2$DV,spots.plot.2$ML,spots.plot.2$AP,spots.plot.2$full.name.parent,spots.plot.2$name,spots.plot.2$expr.level))
  
  
  p <- layout(p,title=sprintf("Gene: %s %s",gene,str.title), updatemenus=updatemenus,
              scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                           zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                           yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))))
  
  return(p)
  
}

plot.3d.diff.expr <- function(gene.list, spots.table, st.data, clust, percentile = 0.85){
  
  p <- plot_ly()
  show.bool <- TRUE
  updatemenus = list(list(active = -1, type = 'buttons', buttons = list()))
  boolean.list <- rep(FALSE,2*length(gene.list))
  for (i in 1:length(gene.list)){
    
    spots.table$expr.level <- as.numeric(t(st.data[gene.list[i],]))
    
    spots.plot <- spots.table
    
    change.opacity <- quantile(spots.plot$expr.level,percentile)
    
    spots.plot.1 <- spots.plot[spots.plot$expr.level > change.opacity,]
    spots.plot.2 <- spots.plot[spots.plot$expr.level <= change.opacity,]
    
    #Plotting
    p <- add_markers(p,data = spots.plot.1, x=~ML, y=~AP, z=~DV, type = 'scatter3d', mode = 'markers', name=sprintf('%s, High expression',gene.list[i]), visible = show.bool,
                     marker=list(size=5, color = ~expr.level,colorscale = c('#FFE1A1', '#683531'), showscale = FALSE,
                                 opacity = 1, cmin = 0, cmax = max(spots.plot$expr.level)),
                     hoverinfo = 'text',
                     text = ~sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %d',
                                     DV,ML,AP,full.name.parent,name,expr.level))
    
    p <- add_markers(p, x=spots.plot.2$ML, y=spots.plot.2$AP, z=spots.plot.2$DV, name=sprintf('%s, Low expression (%0.2f percentile)',gene.list[i],percentile), visible = show.bool,
                     marker=list(size=5, color = spots.plot.2$expr.level, colorscale = c('#FFE1A1', '#683531'), showscale = FALSE,
                                 opacity = 0.5,cmin = 0, cmax = max(spots.plot$expr.level),name='Low expression') ,
                     hoverinfo = 'text',
                     text = sprintf('DV: %0.2f<br />ML: %0.2f<br />AP: %0.2f<br />Region: %s/%s<br />Expression level: %d',
                                    spots.plot.2$DV,spots.plot.2$ML,spots.plot.2$AP,spots.plot.2$full.name.parent,spots.plot.2$name,spots.plot.2$expr.level))
    
    show.bool <- FALSE
    bool.update <- boolean.list
    bool.update[(2*i-1):(2*i)] <- TRUE
    updatemenus[[1]]$buttons[[i]] <- list(label = gene.list[[i]],
                                          method = "update",
                                          args = list(list(visible = bool.update))) 
    
  }
  
  p <- layout(p, title=sprintf("Cluster: %s",clust),
              scene = list(xaxis = list(title = 'x = Medial-lateral (mm)',range = c(0,5)),
                           zaxis = list(title = 'z = Dorsal-ventral (mm)',range = c(-7.9,0)),
                           yaxis = list(title = 'y = Anterior-posterior (mm)',range = c(-5.9,3))),
              updatemenus=updatemenus)
  
  return(p)
  
}

open.3d.plot <- function(p){
  
  path <- paste(tempfile(),'html',sep='.')
  p <- htmlwidgets::saveWidget(p, path, selfcontained = FALSE)
  browseURL(path)
}

## ----------- 2D Plotly ----------- 

#Plotting expression levels in the atlas
plot.2d.atlas.expr <- function(gene, spots.table, st.data, AP = NULL, slice_index = NULL, min.expr = 0, cex = 1, new.window = TRUE, text.plot = TRUE, log.scale = FALSE, normalized.scale = FALSE, win.dim = c(14,14)){
  
  #gene(str): gene to plot
  #spots.table(df)
  #st.data(df)
  #AP(float): AP coordinates. If does not match any sections, takes the closest. If two values, plot all sections in-between.
  #slice_index(str): If specifcied, replaces AP coordinates
  #min.expr(float): Expression threhsold to plot spot. If between 0 and 1, defines a percentile. Otherwise, defines a number of reads (not log normalized)
  #cex(float): dimension of text
  #new.window(bool): force new window
  #text.plot(bool): show the legend and title 
  #log.scale(bool): use a logarithmic scale instead
  #normalized.scale(bool): all expression are a value between 0 and 100 (maximum expression across the whole brain)
  #win.dim(numeric vector): dimension of the graphical window
  
  library(data.table)
  
  #Case with two coordinates
  if(length(AP) == 2){
    AP.list <- unique(spots.table$AP)
    AP.list <- sort(AP.list[AP.list >= min(AP) & AP.list <= max(AP)], decreasing = TRUE)
    for (i in AP.list)
      plot.2d.atlas.expr(gene,spots.table,st.data,i,NULL,min.expr,cex,TRUE,text.plot,log.scale, normalized.scale)
  }else{
    
    if (!is.null(AP))
      AP <- select.closest.AP(AP, spots.table)
    
    #Case where plotting virtual slice (both dupplicates at a specific coordinates)
    if (is.null(slice_index))
      slice_index <- unique(spots.table[spots.table$AP == AP,'slice_index'])
    
    #Loading the slice table for correspondance between old and new indexing
    slices.table <- load.slice.table()
    
    old.slice.index <- slices.table[slice_index,'slice_old_id']

    #Create a function to generate a continuous color palette
    rbPal <- colorRampPalette(c('white','red'))
    
    dataset <- get.expression.slice(gene, spots.table, st.data, min.expr, log.scale, normalized.scale, slice_index, rbPal)
    AP <- dataset[[1]]$AP[1]
    
    #Crashing case
    if (length(old.slice.index) > 1){
      for(j in 2:length(old.slice.index)){
        if(dataset[[j]]$AP[1] != AP)
          stop('Unxpected multiple AP coordinates')
      }
    }
      
    m.val <- dataset[[length(dataset)]]
    
    #Selecting proper infos from the atlas
    AP.infos <- atlas.spots$outlines$AP[atlas.spots$outlines$AP$AP == AP,]
    outlines <- atlas.spots$outlines$outlines[[AP.infos$i]]
    
    width <- AP.infos$xmax -  AP.infos$xmin
    height <- AP.infos$ymax - AP.infos$ymin
    margin.x <- 0.4 * width
    margin.y <- 0.05 * height
    xMin <- AP.infos$xmin - margin.x
    xMax <- AP.infos$xmax + margin.x
    yMin <- AP.infos$ymin - margin.y
    yMax <- AP.infos$ymax + margin.y

    #Initializing the plot
    if (new.window)
      win.graph(width = win.dim[1],height = win.dim[2])
    
    if (text.plot)
      title.str <- paste("Bregma:", AP,"mm\nSlices: ",paste(slice_index, collapse = ' & '),'\nGene:',gene)
    else
      title.str <- ''
    
    plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin),
         xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "",
         col = 0, main = title.str, font.main = 1)
    
    if (text.plot){
      mtext("Dorso-ventral (mm)", side = 2, line = 1.5)
      mtext("Medio-lateral (mm)", side = 1, line = -1.5)
    } 
    
    #Plotting outlines
    for (j in 1:length(outlines)){
      o <- outlines[[j]]
      
      #Colors for fiber and ventricle
      if (as.character(o$col) == '#aaaaaa')
        col = 'black'
      else if (as.character(o$col) == '#cccccc')
        col = 'darkgray'
      else
        col = NA
      
      polygon(o$x,o$y,border = "black",col = col)
    }
    
    #Plotting the spots
    for(i in 1:length(old.slice.index)){
      
      #Selecting the spots registered into the atlas
      cur.soma <- atlas.spots$spots[rownames(dataset[[i]]),]
      points(cur.soma$somaX, cur.soma$somaY, pch=21, bg = dataset[[i]]$col,cex = cex)
    }

    #Legend plotting
    leg.ys <- (yMin+yMax)/2-0.25*height
    leg.ye <- (yMin+yMax)/2+0.25*height
    y.list <- seq(from = leg.ys, to = leg.ye, length.out = 101)
    col.vect <- rev(rbPal(100))
    
    for(i in 1:(length(y.list)-1))
      rect(xMax-0.25*width,y.list[i],xMax-0.2*width,y.list[i+1],col=col.vect[i], border = NA)

    rect(xMax-0.25*width,leg.ys,xMax-0.2*width,leg.ye)
    text(xMax-0.2*width,leg.ye,'0',cex = cex,pos = 4)
    text(xMax-0.2*width,leg.ys,m.val,cex = cex,pos = 4)
    if(log.scale)
      text(xMax-0.225*width,leg.ys,'Log scale',cex = cex,pos=3)
    if(normalized.scale)
      str.s <- 'Normalized scale'
    else
      str.s <- 'Number of reads'
    text(xMax-0.15*width,(yMin + yMax)/2,str.s,cex = cex,pos=3,srt = -90)
  }
}

#Plotting expression levels on a slice
plot.2d.slice.expr <- function(gene,spots.table, st.data, slice_index,  min.expr = 0, cex = 1, text.plot = TRUE, new.window = TRUE, log.scale = FALSE, normalized.scale = FALSE, win.dim = c(14,14), res.ds = 10){

  #gene(str): gene to plot
  #spots.table(df)
  #st.data(df)
  #slice_index(str): Section to plot
  #min.expr(float): Expression threhsold to plot spot. If between 0 and 1, defines a percentile. Otherwise, defines a number of reads (not log normalized)
  #cex(float): dimension of text
  #new.window(bool): force new window
  #text.plot(bool): show the legend and title 
  #log.scale(bool): use a logarithmic scale instead
  #normalized.scale(bool): all expression are a value between 0 and 100 (maximum expression across the whole brain)
  #win.dim(numeric vector): dimension of the graphical window
  #res.ds(integer): downsample factor for plotting the HE picture faster
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('white','red'))
  
  dataset <- get.expression.slice(gene, spots.table, st.data, min.expr, log.scale, normalized.scale, slice_index, rbPal)
  spots.cur <- dataset[[1]]
  spots.cur$radius <- spots.table[rownames(spots.cur),'radius']
  m.val <- dataset[[2]]
  
  slices.table <- load.slice.table()
  old.slice.index <- slices.table[slice_index,'slice_old_id']
  rotation <- slices.table[slice_index,'rotation']
  
  HE.path <- sprintf('%s/%s/HE_cropped/HE_%s.jpg',path.slices,old.slice.index,old.slice.index)
  library(jpeg)
  HE.original <- readJPEG(HE.path)
  HE.resolution.x <- dim(HE.original)[2]
  HE.resolution.y <- dim(HE.original)[1]

  HE.downscale <- HE.original[seq(1,HE.resolution.y,by=res.ds),seq(1,HE.resolution.x,by=res.ds),]
  
  if(new.window)
    win.graph(width = win.dim[1],height = win.dim[2])
  
  if(rotation == 90){
    xlim <- c(0,HE.resolution.x)
    ylim <- c(2*HE.resolution.y,HE.resolution.y)
    angle <- -rotation
    spots.cur$y <- spots.cur$y + HE.resolution.y 
  }
  if(rotation == -90){
    xlim <- c(-HE.resolution.x,0)
    ylim <- c(HE.resolution.y,0)
    angle <- -rotation
    spots.cur$x <- spots.cur$x - HE.resolution.y 
    spots.cur$y <- spots.cur$y + (abs(HE.resolution.x - HE.resolution.y))
  }
  if(rotation == 0){
    xlim <- c(0,HE.resolution.x)
    ylim <- c(HE.resolution.y,0)
    angle <- 0
  }
  
  if (text.plot)
    title.str <- paste("Bregma:", spots.cur$AP[1],"mm\nSlices: ",paste(slice_index, collapse = ' & '),'\nGene:',gene)
  else
    title.str <- ''

  width <- abs(xlim[1] - xlim[2])
  
  plot(c(xlim[1],xlim[2] + 0.2*width), ylim, xlim = c(xlim[1],xlim[2] + 0.4*width), ylim = ylim,
       asp = 1, axes = F, xlab = "", ylab = "",
       col = 0, main = title.str, font.main = 1)
  
  
  rasterImage(as.raster(HE.downscale), xleft = 1, xright = HE.resolution.x, ybottom =  HE.resolution.y, ytop = 1, angle=angle, asp = 1)
  symbols(spots.cur$x,spots.cur$y,circles = spots.cur$radius, bg = spots.cur$col, inches=FALSE, asp = 1,fg='black',add = TRUE, lwd = 0.2, asp = 1)
  
  leg.ys <- mean(ylim)-(0.2*abs(ylim[1]-ylim[2]))
  leg.ye <- mean(ylim)+(0.2*abs(ylim[1]-ylim[2]))
  y.list <- seq(from = leg.ys, to = leg.ye, length.out = 101)
  col.vect <- rev(rbPal(100))
  
  xMax <- max(xlim)

  for(i in 1:(length(y.list)-1))
    rect(xMax+0.1*width,y.list[i],xMax+0.15*width,y.list[i+1],col=col.vect[i], border = NA)
  
  rect(xMax+0.1*width,leg.ys,xMax+0.15*width,leg.ye)
  text(xMax+0.15*width,leg.ye,'0',cex = cex,pos = 4)
  text(xMax+0.15*width,leg.ys,m.val,cex = cex,pos = 4)
  if(log.scale)
    text(xMax+0.125*width,leg.ys,'Log scale',cex = cex,pos=3)
  if(normalized.scale)
    str.s <- 'Normalized scale'
  else
    str.s <- 'Number of reads'
  text(xMax+0.17*width,mean(ylim),str.s,cex = cex,pos=3,srt = -90)
  
  
}

plot.2d.specific.coordinates <- function(spots.table, coordinate, order.by.cluster.size = FALSE, use.distinct.shapes = TRUE){
  
  spots.plot <- spots.table[spots.table$AP == coordinate,]
  
  if(use.distinct.shapes ==  TRUE){
    pch.list <- 15:20
    symb.list <- numeric(2*length(ordered.col.vect))
    curs <- 1
    cur.pch <- 1
    while(curs <= length(ordered.col.vect)){
      symb.list[curs:(curs+length(col.tab.paired)-1)] <- pch.list[cur.pch]
      cur.pch <- cur.pch + 1
      if (cur.pch > length(pch.list))
        cur.pch = 1
      curs <- curs + length(col.tab.paired)
    }
    symb.list <- symb.list[1:length(ordered.col.vect)]
  }else
    symb.list <- rep(19,length(ordered.col.vect))
  
  if (order.by.cluster.size){
    library(plyr)
    df.count <- as.data.frame(table(spots.plot$clusters.named))
    ord <- order(df.count$Freq,decreasing = TRUE)
    spots.plot$clusters.named <- factor(spots.plot$clusters.named, levels = levels(spots.plot$clusters.named)[ord])
    ordered.col.vect <- col.vect[ord]
    symb.list <- symb.list[ord]
  }else
    ordered.col.vect <- col.vect
  
  plot_ly(type = 'scatter', mode = 'markers', data = spots.plot, color=~clusters.named,
          x = ~ML, y = ~DV, colors = ordered.col.vect, marker = list(size=8), symbol = ~clusters.named, symbols = symb.list, hoverinfo='text', xaxis = "x", yaxis = "y",
          text = ~sprintf('Spot: %s<br />DV: %0.2f<br />ML: %0.2f<br />AP: %0.3f<br />Region: %s/%s<br />Cluster: %s',
                          rownames(spots.plot),DV,ML,AP,full.name.parent,name,clusters.named)) %>%
    layout(yaxis = list(scaleanchor = "x", range = c(-7.5,0.5)),
           xaxis = list(range = c(0,7)))
  
}

plot.2d.on.slice <- function(spots.cur, slice){
  
  path.slices <- 'C:/Users/MatLab/Desktop/transcripBrainAtlas/data-registration-final/data-raw'
  index.slices <- read.table('C:/Users/MatLab/Desktop/transcripBrainAtlas/data-registration-final/scripts/index-slices.csv',sep=';', stringsAsFactors = FALSE, header = TRUE)
  cur.slice <- index.slices$slice.name[index.slices$id == slice] 
  
  spots.cur <- spots.cur[spots.cur$slice.index == slice,]
  
  HE.path <- sprintf('%s/%s/HE_cropped/HE_%s.jpg',path.slices,cur.slice,cur.slice)
  library(jpeg)
  HE.original <- readJPEG(HE.path)
  HE.resolution.x <- dim(HE.original)[2]
  HE.resolution.y <- dim(HE.original)[1] 
  
  spots.cur$without.slice.index <- transpose(as.data.frame(strsplit(spots.cur$spots.id,c('_')), stringsAsFactors = FALSE)[2,])
  spots.cur$id.x <- as.numeric(transpose(as.data.frame(strsplit((spots.cur$without.slice.index$V1),'x'), stringsAsFactors = FALSE))[,1])
  spots.cur$id.y <- as.numeric(transpose(as.data.frame(strsplit((spots.cur$without.slice.index$V1),'x'), stringsAsFactors = FALSE))[,2])
  spots.cur$pixels.x <- (spots.cur$id.x - 1) * (HE.resolution.x /32)
  spots.cur$pixels.y <- (spots.cur$id.y - 1) * (HE.resolution.y /34)
  
  count <- NULL
  curs <- 0
  for (j in unique(spots.cur$clusters.named)){
    curs <- curs + 1
    count[[curs]] <- sum(spots.cur$clusters.named == j)
  }
  
  #Plotting the alignement with circles on the picture
  win.graph()
  plot(1,type='n',xlim = c(0,HE.resolution.x), ylim = c(HE.resolution.y,0), main = sprintf('Slice: %s - AP: %0.3f', slice, spots.cur$AP[1]),
       xlab = 'x coordinates in HE picture (pixels)',ylab = 'y coordinates in HE picture (pixels)', asp = 1)
  rasterImage(as.raster(HE.original), xleft = 1, xright = HE.resolution.x, ybottom =  HE.resolution.y, ytop = 1,angle=0, asp = 1)
  symbols(spots.cur$pixels.x,spots.cur$pixels.y,circles = spots.cur$radius, bg = spots.cur$col, inches=FALSE, asp = 1,fg='black',add = TRUE, lwd = 0.2, asp = 1)
  
  legend.x <- rep(-3500,length(unique(spots.cur$col)))
  legend.y <- seq(2000,8000,length.out = length(unique(spots.cur$col)))
  
  legend.matrix <- NULL
  legend.matrix$name.clust <- (unique(spots.cur$clusters.named))
  legend.matrix$str <- paste(unique(spots.cur$clusters.named),'(N = ',count,')',sep='')
  legend.matrix <- data.frame(legend.matrix, stringsAsFactors = FALSE)
  legend.matrix <- legend.matrix[order(legend.matrix$name.clust),]
  for (i in 1:length(legend.matrix$name.clust)){
    legend.matrix$col[i] <- spots.cur$col[spots.cur$clusters.named == legend.matrix$name.clust[i]][1]
  }
  
  symbols(legend.x,legend.y, bg = legend.matrix$col, circles = rep(70,length(unique(spots.cur$col))), inches = FALSE, add = TRUE, asp = 1)
  text(legend.x,legend.y,legend.matrix$str, pos = 4)
  
}

#Returns a plotly 2d atlas figure with the spots and their cluster
plot.2d.atlas.clust <- function(spots.table, AP, slice_index = NULL, col.vect = NULL, order.by.cluster.size = TRUE, minimum.spots = 0){
  
  #spots.table
  #AP(float): AP coordinates of the virtual section to plot (potential two sections at the same coordinate)
  #slice_index(char): list of section to plot (all in the same AP coordinates, useful to plot one section only)
  #col.vect(df char*char): data frame that associates a cluster name and a color for plotting. Default color is used if null
  #order.by.cluster.size(bool): if TRUE, clusters in the legned are ordered by their number of elements on the section
  #minimum.spots(int): clusters with less element on the section than this value will not be displayed
  
  library(plyr)
  
  #Case where plotting virtual slice (both dupplicates at a specific coordinates)
  if (is.null(slice_index)){
    AP <- select.closest.AP(AP,spots.table)
    slice_index <- unique(spots.table[spots.table$AP == AP,'slice_index'])
  }else{
    AP <- spots.table[spots.table$slice_index == slice_index, 'AP'][1]
  }
  
  #Selecting only appropriate spots
  selection.bool <- sapply(spots.table$slice_index,function(x){return(is.element(x,slice_index))})
  spots.selected <- spots.table[selection.bool,]

  #Counting spots in each cluster
  cnt <- plyr::count(spots.selected$clusters.named)
  
  if (minimum.spots > 0){
    cnt <- cnt[cnt$freq >= minimum.spots,]
    spots.selected <- spots.selected[sapply(spots.selected$clusters.named, function(x){return(is.element(x,cnt$x))}),]
  }
  
  #If no color vector provided, generating it automatically
  if (is.null(col.vect)){
    
    #Case with reordered colors
    if (order.by.cluster.size){
      cnt <- cnt[order(cnt$freq,decreasing = TRUE),]
      col.vect <- data.frame(clusters.named = cnt$x,
                             color = colorRampPalette(brewer.pal(9, "Set1"))(dim(cnt)[1]),
                             stringsAsFactors = FALSE)
      
    #Basic case
    }else
      col.vect <- data.frame(clusters.named = as.character(unique(spots.selected$clusters.named)),
                             color = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(spots.selected$clusters.named))),
                            stringsAsFactors = FALSE)
  
  #Otherwise, ordering it if required  
  }else{
    col.vect <- col.vect[sapply(col.vect$clusters.named, function(x){return(is.element(x,unique(as.character(spots.selected$clusters.named))))}),]
    if (order.by.cluster.size){
      col.vect <- col.vect[order(col.vect$clusters.named),]
      col.vect <- col.vect[order(cnt$freq,decreasing = TRUE),]
    }   
  }
  
  #Selecting proper infos from the atlas
  AP.infos <- atlas.spots$outlines$AP[atlas.spots$outlines$AP$AP == AP,]
  outlines <- atlas.spots$outlines$outlines[[AP.infos$i]]

  #Initializing the plot
  p <- plot_ly()
  first.bool <- TRUE
  
  #Plotting outlines
  for (j in 1:length(outlines)){
    o <- outlines[[j]]
    
    #Colors for fiber and ventricle
    if (as.character(o$col) == '#aaaaaa')
      col = 'rgb(0,0,0)'
    else if (as.character(o$col) == '#cccccc')
      col = 'I(darkgray)'
    else
      col = 'transparent'
    
    #Adding the traces
    p <- add_trace(p, type = "scatter", mode = 'lines', x = o$x, y = o$y, showlegend = first.bool, fill = 'toself',fillcolor = col,
                   hoverinfo = 'skip',line = list(color = 'rgb(0,0,0)'), hoveron = 'points', legendgroup = 'Brain outline', name = 'Brain outline')
    
    #Only legend for the first trace (do not want > 30 "Brain outline" legends)
    if (first.bool)
     first.bool <- FALSE
  }
  
  cnt.trace <- 0
  
  #Going through the clusters to plot
  for (k in 1:dim(col.vect)[1]){

    #Selecting the spots registered into the atlas
    cur.soma <- atlas.spots$spots[rownames(spots.selected[spots.selected$clusters.named == col.vect[k,'clusters.named'],]),]
    
    #If there are spots to plot
    if(dim(cur.soma)[1] > 0){
      cur.soma$full.name <- sapply(spots.selected[rownames(cur.soma),'acronym'],name.from.acronym)
  
      #Plotting them
      p <- add_trace(p, type = "scatter", mode = 'markers', data = cur.soma, x = ~somaX, y = ~somaY, hoverinfo='text',
                     showlegend = TRUE, marker = list(size=10,color=col.vect[k, 'color'], opacity = 1), name = col.vect[k,'clusters.named'],
                     text = sprintf('Spot: %s<br />Cluster: %s<br />%s', rownames(cur.soma),col.vect[k,'clusters.named'],cur.soma$full.name))
      
      cnt.trace <- cnt.trace + 1
    }
  }

  
  updatemenus = list(list(active = -1, type = 'buttons', x = 0.9, y = 1, buttons = list()))
  updatemenus[[1]]$buttons[[1]] <- list(label = 'Hide all', method = 'update',args = list(list(visible = c(rep(TRUE,length(outlines)),rep('legendonly',cnt.trace)))))
  updatemenus[[1]]$buttons[[2]] <- list(label = 'Show all', method = 'update', args = list(list(visible =  c(rep(TRUE,length(outlines)),rep(TRUE,cnt.trace)))))
  
  #Adding appropriate layout
  p <- layout(p, yaxis = list(scaleanchor = "x", range = c(AP.infos$ymax,AP.infos$ymin), title = 'Dorsal-Ventral', showticklabels = FALSE,
                              zeroline = FALSE, showline = FALSE, showgrid = FALSE),
              xaxis = list(range = c(AP.infos$xmin, AP.infos$xmax), title = 'Medial-Lateral', showticklabels = FALSE,
                           zeroline = FALSE, showline = FALSE, showgrid = FALSE),
              showlegend = TRUE,
              title = sprintf('Slices: %s<br />AP: %.3f',paste(slice_index,collapse = ' & '),AP),
              legend = list(traceorder = 'normal+grouped'),
              updatemenus = updatemenus)
  
  return(p)
}

plot.2d.ica <- function(spots.table, ica.matrix, ica.vect, by.section = FALSE){
  
  
  spots.table <- spots.table[intersect(rownames(spots.table),rownames(ica.matrix)),]
  ica.matrix <- ica.matrix[intersect(rownames(spots.table),rownames(ica.matrix)),]
  
  slice.table <- load.slice.table()
  if (by.section)
    x.val <- sort(unique(rownames(slice.table)),decreasing = FALSE)
  else
    x.val <- sort(unique(slice.table$AP),decreasing = TRUE)
  
  ica <- ica.matrix[,ica.vect]
  
  y.val <- numeric(length = length(x.val))

  idx <- 0
  for(i in x.val){
    idx <- idx + 1
    if (by.section)
      l <- spots.table[spots.table$slice_index == i,]
    else
      l <- spots.table[spots.table$AP == i,]
    
    y.val[idx] <- mean(ica[rownames(l)])
  }
  
  plot(x.val, y.val)
}

## ----------- AP axis plots----------- 

#Plot the mean expression along the AP axis of genes
plot.ap.expr <- function(spots.table, st.data, gene = NULL, limits = NULL, title = NULL, draw.line.min.expr = FALSE){
  
  #Case where spots in spots table and matrix do not match
  if (dim(st.data)[1] != dim(spots.table)[1]){
    r <- intersect(rownames(spots.table),rownames(st.data))
    spots.table <- spots.table[r,]
    st.data <- st.data[r,]
  }
  
  #If no gene is specified, all of them are selected
  if (is.null(gene))
    gene <- colnames(st.data)

  #Default title in case no gene is specified
  if (is.null(title))
    title = gene
  
  #List of all APs
  ap.list <- sort(unique(spots.table$AP), decreasing = TRUE)
  
  #Only TRUE during first loop
  first.loop <- TRUE
  
  #Going through all genes (possibly only one)
  for (i in 1:length(gene)){
    
    #Selecting gene
    g <- gene[[i]]
    
    #Will contain Y value (means)
    y <- NULL
    e <- NULL
    
    #Looping  through AP and computing mean value
    for (ap in ap.list){
      y <- c(y,mean(st.data[spots.table$AP == ap,g]))
      e <- c(e,sd(st.data[spots.table$AP == ap,g]))
    }
    
    #Default ylim
    if (first.loop & is.null(limits)){
      mi.v <- min(y-e)
      ma.v <- max(y+e)
      diff <- ma.v - mi.v
      limits <- c(mi.v - 0.1*diff, ma.v + 0.1*diff)
    }
    
    #Initial plot
    if (first.loop)
      plot(ap.list, y, type = 'o', ylim = limits, xlim = c(ap.list[1]-0.2,ap.list[length(ap.list)]+0.2), xlab= 'AP', ylab = 'Mean reads per spot', main = title)
    
    #Not first loop, only adding lines
    else
      lines(ap.list,y)
    
    
    arrows(ap.list, y-e, ap.list, y+e, length=0.05, angle=90, code=3)
    
    #Changing value of boolean
    first.loop <- FALSE
    
    #Drawing line if required
    if (draw.line.min.expr)
      segments(-1000,min(y),1000,min(y), lty=2, col = 'lightgray')
    
  }
}

## ----------- Clusters -----------  

find.clusters.composition <- function(spots.table){
  
  library(ggplot2)
  library(RColorBrewer)
  library(plyr)
  
  spots.table <- spots.table[!is.na(spots.table$clusters.named),]
  list.parents <- unique(spots.table$acronym.parent)
  clusters <- unique(spots.table$clusters.named)
  
  for(cl in clusters){
    spots.table.loop <- spots.table[spots.table$clusters.named == cl,]
    df <- plyr::count(spots.table.loop$acronym.parent)
    df <- df[order(df$freq, decreasing = TRUE),]
    df$x <- factor(df$x, levels = df$x[order(df$freq,decreasing = TRUE)])
    df$col.plot <- color.from.acronym(as.character(df$x))
    df$col.plot <- factor(df$col.plot,levels = df$col.plot)
    df$full.name <- name.from.acronym(as.character(df$x))
    df$txt <- sprintf('%.1f%%',(df$freq / sum(df$freq)*100))
    win.graph(14,14)
    g <- ggplot(df,aes(x,freq, fill = col.plot))
    g <- g + geom_bar(stat='identity',color = 'black')  
    g <- g + scale_fill_manual(values = as.character(df$col.plot), labels = df$full.name)
    g <- g + geom_text(aes(x,freq,label = txt),size=3.5,color = 'black',vjust = 'bottom',hjust = 'middle', nudge_y = df$freq[1]*0.01)
    g <- g + labs(y = 'Number of spots',x = '', title = paste('Cluster: ',cl,sep=''))
    g <- g + scale_y_continuous(expand = c(0,0), limits = c(0,1.1*df$freq[1])) 
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
    print(g)
  }
  
  col.vect <- colorRampPalette(brewer.pal(9, "Set1"))(length(levels(spots.table$clusters.named)))
  df.col <- data.frame(cluster = levels(spots.table$clusters.named), color = col.vect, stringsAsFactors = FALSE)
  for (reg in unique(spots.table$acronym.parent)){
    spots.table.loop <- spots.table[spots.table$acronym.parent == reg,]
    df <- plyr::count(spots.table.loop$clusters.named)
    df <- df[order(df$freq, decreasing = TRUE),]
    df$x <- factor(df$x, levels = df$x[order(df$freq,decreasing = TRUE)])
    df$col.plot <- mapvalues(df$x,from = df.col$cluster, to = df.col$color, warn_missing = FALSE)
    df$col.plot <- factor(df$col.plot,levels = df$col.plot)
    # df$full.name <- (as.character(df$x))
    df$txt <- sprintf('%.1f%%',(df$freq / sum(df$freq)*100))
    win.graph(14,14)
    g <- ggplot(df,aes(x,freq, fill = col.plot))
    g <- g + geom_bar(stat='identity',color = 'black')  
    g <- g + scale_fill_manual(values = as.character(df$col.plot), labels = df$x)
    g <- g + geom_text(aes(x,freq,label = txt),size=3.5,color = 'black',vjust = 'bottom',hjust = 'middle', nudge_y = df$freq[1]*0.01)
    g <- g + labs(y = 'Number of spots',x = '', title = paste('Brain region: ',reg,sep=''))
    g <- g + scale_y_continuous(expand = c(0,0), limits = c(0,1.1*df$freq[1])) 
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
    print(g)
  }
}

#Appends in the field 'color' the kellys color by cluster
append.kellys.color.spots.table <- function(spots.table, mode = 'kellys'){
  
  if(mode == 'kellys')
	kellys.color <- c("gray95", "gray13", "gold2", "plum4", "darkorange1", "lightskyblue2", "firebrick", "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2", "mediumpurple4", "orange", "maroon", "yellow3", "brown4", "yellow4", "sienna4", "chocolate", "gray19")
 
 #Vivid
 else
	kellys.color <- c('#E58606','#5D69B1','#52BCA3','#99C945','#CC61B0','#24796C','#DAA51B','#2F8AC4','#764E9F','#ED645A','#CC3A8E','#A5AA99')
  
 
 k <- rep(kellys.color,100)
  spots.table$color <- as.character(mapvalues(spots.table$clusters.named, from = levels(spots.table$clusters.named), to = k[1:length(levels(spots.table$clusters.named))]))
  return(spots.table)
}

append.common.vivid.colors <- function(spots.table){
  
  load(paste(path.bin, 'vivid-colors.RData', sep='/'))
  
  spots.table$color <- as.character(mapvalues(spots.table$clusters.named, 
                                              from = rownames(df.colors.vivid),
                                              to = df.colors.vivid$clusters.colors))
  return(spots.table)
}

#----------- Brain outline ----------- 

#Plot in a nice way the coronal cut of brain outline with color inside
plot.brain.outline <- function(outlines, color, soma, name.field.x = 'x', name.field.y = 'y', scale.factor = 1,
                               open.new.plot = TRUE, xlim.stereo = c(-1,6), ylim.stereo = c(-8,0), AP = NULL, 
                               right.hemisphere = TRUE){
  
  #outlines: list of outlines
  #color(char): list of colors that are paired with the polygons from outlines list
  #soma(df): points with pixels and ML/DV coordinates for proper remaping
  #name.field.x(char): name of the field in outlines that contains the X values
  #name.field.y(char): name of the field in outlines that contains the y values
  #scale.factor(num): scale factor (if = 1, no scaling)
  #open.new.plot(logical): force opening of a new plot. If FALSE, limits are useless
  #xlim.stereo(num): 2 values that define the x limits in stereotaxic coordinates
  #ylim.stereo(num): 2 values that define the y limits in stereotaxic coordinates
  #AP(string or num): title
  #right.hemisphere(logical): if true, right hemisphere is the ARA. Otherwise, it's the left one
  
  if(right.hemisphere){
    x.coeff <- 1
    ML.0 <- 0
  }else{
    x.coeff <- -1
    ML.0 <- ML.to.pixels(0, soma)    
  }

  #Converting into pixels limit
  xlim.pixels <- ML.to.pixels(xlim.stereo, soma)
  ylim.pixels <- DV.to.pixels(ylim.stereo, soma)  
  
  #Getting the segments to plot
  df.segments <- get.segments.outlines(outlines, name.field.x, name.field.y, scale.factor)
  
  #Defining title
  if (is.null(AP))
    title <- ''
  else
    title <- sprintf('Bregma: %.2f mm', AP)
  
  #Opening new plot
  if (open.new.plot){
    plot(1, type = 'n', asp = 1, ylim = ylim.pixels, xlim = xlim.pixels, 
         main = title, bty = "n", xaxt = 'n', yaxt = 'n', xlab = '', ylab='') 
  }
  
  #Plotting the inside filling
  for (x in 1:length(outlines)){
    
    #Current outline
    o <- outlines[[x]]
    
    # print(o)
    
    #X and Y scaled
    x.scaled <- (x.coeff*(o[[name.field.x]]-ML.0) / scale.factor) + ML.0
    y.scaled <- o[[name.field.y]] / scale.factor
    
    #Plotting polygons for colors
    polygon(x.scaled, y.scaled, col = color[x], border = NA)
    
  }

  #Plotting the contours
  segments(x.coeff*(df.segments$x0-ML.0) + ML.0, df.segments$y0, x.coeff*(df.segments$x1-ML.0) + ML.0, df.segments$y1, lwd = 0.8)
  
  #Displaying the axis
  axis(3, pos = DV.to.pixels(ylim.stereo[2], soma), at = ML.to.pixels(xlim.stereo[1]:xlim.stereo[2], soma),labels = xlim.stereo[1]:xlim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0))
  axis(2, pos = ML.to.pixels(xlim.stereo[1], soma), at = DV.to.pixels(ylim.stereo[1]:ylim.stereo[2], soma),labels = ylim.stereo[1]:ylim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0), las= 1)
  axis(1, pos = DV.to.pixels(ylim.stereo[1], soma), at = ML.to.pixels(xlim.stereo[1]:xlim.stereo[2], soma),labels = xlim.stereo[1]:xlim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0))
  axis(4, pos = ML.to.pixels(xlim.stereo[2], soma), at = DV.to.pixels(ylim.stereo[1]:ylim.stereo[2], soma),labels = ylim.stereo[1]:ylim.stereo[2], tck  = 0.02, mgp=c(3, .5, 0), las= 1)
  
}

#Plot the atlas with a "double plate" layout: on the left, the spots with the actual color. On the right, the ARA
.plot.2d.atlas.double.plate.col <- function(spots.table, cex = 1, lwd = 1){

  #spots.table: spots.table. Must include only one AP coordinates. Color is used to color the spots
  #cex: points cex
  #lwd: points lwd
  
  #Selecting unique AP
  AP <- unique(spots.table$AP)[1]
  id.atlas <- which.min(abs(aligned.atlas$AP - AP))
  
  #Plotting brain outline
  plot.brain.outline(aligned.atlas$outlines[[id.atlas]],
                     aligned.atlas$color[[id.atlas]],
                     aligned.atlas$soma[[id.atlas]],
                     AP = AP, xlim.stereo = c(-6,6))


  #Getting the segments to plot
  df.segments <- get.segments.outlines(aligned.atlas$outlines[[id.atlas]], 'x', 'y', 1)
  
  #Remapping arround 0
  df.segments$x0 <- ML.to.pixels(0, aligned.atlas$soma[[id.atlas]]) - (df.segments$x0 - ML.to.pixels(0, aligned.atlas$soma[[id.atlas]]))
  df.segments$x1 <- ML.to.pixels(0, aligned.atlas$soma[[id.atlas]]) - (df.segments$x1 - ML.to.pixels(0, aligned.atlas$soma[[id.atlas]]))
  
  #Plotting the contours
  segments(df.segments$x0, df.segments$y0, df.segments$x1, df.segments$y1, lwd = 0.8)
  
  #Converting the ML/DV into coordinates
  somaX <- ML.to.pixels(0, aligned.atlas$soma[[id.atlas]]) - 
          (ML.to.pixels(spots.table$ML, aligned.atlas$soma[[id.atlas]]) - ML.to.pixels(0, aligned.atlas$soma[[id.atlas]]))
  somaY <- DV.to.pixels(spots.table$DV, aligned.atlas$soma[[id.atlas]])
  
  #Plotting
  points(somaX, somaY, col = 'black', bg = spots.table$color, pch = 21, cex = cex, lwd = lwd)
  
}  

#
# for(cl in setdiff(unique(as.numeric(grid.2d)),0)){
#   pts <- which(grid.2d == cl, arr.ind = T)
#   x <- x.vect.remap[pts[,1]]
#   y <- y.vect.remap[pts[,2]]
#   points(x,y,col = df.colors[df.colors$cluster.id == cl,'clusters.colors'], pch = 15, cex = 0.2)
# }
# points(x.edges, y.edges, col = 'black', cex = 0.1, pch = 15)
# dev.off()
# 
# }
## ---- OLD -------

# plot.2d.atlas.clust <- function(spots.table, AP, slice_index = NULL, col.vect = NULL){
#   
#   library(data.table)
#   library(plotly)
#   setwd(path.slices)
#   
#   #Case where plotting virtual slice (both dupplicates at a specific coordinates)
#   if (is.null(slice_index)){
#     slice_index <- unique(spots.table[spots.table$AP == AP,'slice_index'])
#   }else{
#     AP <- spots.table[spots.table$slice_index == slice_index, 'AP'][1]
#   }
#   
#   if (is.null(col.vect))
#     col.vect <- colorRampPalette(brewer.pal(9, "Set1"))(length(levels(spots.table$clusters.named)))
#   
#   #Loading the slice table for correspondance between old and new indexing
#   slices.table <- load.slice.table()
#   
#   old.slice.index <- slices.table[slice_index,'slice_old_id']
#   
#   #Initializing lists that will contain extreme values
#   list.max.x <- NULL
#   list.max.y <- NULL
#   list.min.x <- NULL
#   list.min.y <- NULL
#   scale.factor <- NULL
#   
#   #Keeping registration, dataset and slice.parameters as lists.
#   registration <- NULL
#   dataset <- NULL
#   slice.parameters.list <- NULL
#   
#   #Plotting one or two slices
#   for (i in 1:length(old.slice.index)){
#     
#     #Loading the appropriates files
#     load(paste(old.slice.index[[i]],'Registration/sliceparameters.RData',sep = '/'))
#     load(paste(old.slice.index[[i]],'Registration/regi.RData',sep = '/'))
#     
#     #Saving in list
#     registration[[i]] <- get.forward.warpRCPP(regi)
#     slice.parameters.list[[i]] <- slice.parameters
#     dataset[[i]] <- spots.table[spots.table$slice_index == slice_index[i],c('ML','DV','AP','acronym','x','y','clusters.named')]
#     
#     #Coordinate of the slice
#     coordinate <- dataset[[i]]$AP[1]
#     
#     #Scale factor from downsizing picture
#     scale.factor[[i]] <- mean(dim(registration[[i]]$transformationgrid$mx)/c(registration[[i]]$transformationgrid$height,
#                                                                              registration[[i]]$transformationgrid$width))
#     
#     #Getting the regions from the atlas
#     numPaths <- registration[[i]]$atlas$numRegions
#     outlines <- registration[[i]]$atlas$outlines
#     
#     #Initializing values to find extreme points
#     list.max.x[[i]] <- 0
#     list.max.y[[i]] <- 0
#     list.min.x[[i]] <- Inf
#     list.min.y[[i]] <- Inf
#     
#     #Finding the extrema
#     for (x in(1:numPaths)) {
#       
#       if (slice.parameters.list[[i]]$right.hemisphere){
#         if (max(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE) > list.max.x[[i]])
#           list.max.x[[i]] <- max(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE)
#         if (max(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE) > list.max.y[[i]])
#           list.max.y[[i]] <- max(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE)  
#         if (min(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE) < list.min.y[[i]])
#           list.min.y[[i]] <- min(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE)
#         if (min(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE) < list.min.x[[i]])
#           list.min.x[[i]] <- min(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE)    
#       }else{
#         if (max(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE) > list.max.x[[i]])
#           list.max.x[[i]] <- max(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE)
#         if (max(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE) > list.max.y[[i]])
#           list.max.y[[i]] <- max(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE)  
#         if (min(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE) < list.min.y[[i]])
#           list.min.y[[i]] <- min(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE)
#         if (min(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE) < list.min.x[[i]])
#           list.min.x[[i]] <- min(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE)         
#       }
#     }
#     
#     if (!slice.parameters.list[[i]]$right.hemisphere){
#       tmp <- list.min.x[[i]]
#       list.min.x[[i]] <- -list.max.x[[i]]
#       list.max.x[[i]] <- -tmp
#     }
#   }
#   
#   #Dimensions and limits of the plot, a margin of 5% fron the extrema is added
#   width <- list.max.x[[1]] - list.min.x[[1]]
#   height <- list.max.y[[1]] - list.min.y[[1]]
#   margin.x <- 0.05 * width
#   margin.y <- 0.05 * height
#   xMin <- list.min.x[[1]] - margin.x
#   xMax <- list.max.x[[1]] + margin.x
#   yMin <- list.min.y[[1]] - margin.y
#   yMax <- list.max.y[[1]] + margin.y
#   
#   #Computing offset to align both slices
#   offset.x <- list.max.x[1] - list.max.x[2]
#   offset.y <- list.max.y[1] - list.max.y[2]
#   
#   soma.merged <- NULL
#   
#   #Plotting for both slices
#   for (i in 1:length(old.slice.index)){
#     
#     #Using the remapped space
#     index <- round(scale.factor[[i]] * cbind(dataset[[i]]$y, dataset[[i]]$x))
#     
#     #Getting the spots coordinates
#     if (slice.parameters.list[[i]]$right.hemisphere){
#       if (i == 1){
#         somaX <- registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]]
#         somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]]
#       }else{
#         somaX <- registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]] + offset.x
#         somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]] + offset.y
#       }
#     }else{
#       if (i == 1){
#         somaX <- -registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]]
#         somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]]
#       }else{
#         somaX <- -registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]] + offset.x
#         somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]] + offset.y
#       }
#     }
#     
#     df <- data.frame(somaX = somaX, somaY = somaY, clusters.named = dataset[[i]][,'clusters.named'], full.name = name.from.acronym(dataset[[i]][,'acronym']))
#     rownames(df) <- rownames(dataset[[i]])
#     
#     soma.merged <- rbind(soma.merged, df)
#     
#     #Plotting the regions
#     numPaths <- registration[[i]]$atlas$numRegions
#     outlines <- registration[[i]]$atlas$outlines
#     
#     
#     p <- plot_ly()
#     first.bool <- TRUE
#     
#     for (x in 1:numPaths){
#       if (as.character(registration[[i]]$atlas$col[x]) == '#aaaaaa')
#         col = 'rgb(0,0,0)'
#       else if (as.character(registration[[i]]$atlas$col[x]) == '#cccccc')
#         col = 'I(darkgray)'
#       else
#         col = 'transparent'
#       
#       
#       if (slice.parameters.list[[i]]$right.hemisphere){
#         #Adding offset to align second slice with the first one
#         if (i == 1){
#           x.poly <- outlines[[x]]$xr/scale.factor[[i]]
#           y.poly <-outlines[[x]]$yr/scale.factor[[i]]
#         }else{
#           x.poly <- outlines[[x]]$xr/scale.factor[[i]]+ offset.x
#           y.poly <- outlines[[x]]$yr/scale.factor[[i]]+ offset.y
#         }
#       }else{
#         if (i == 1){
#           x.poly <- -outlines[[x]]$xl/scale.factor[[i]]
#           y.poly <- outlines[[x]]$yl/scale.factor[[i]]
#         }else{
#           x.poly <- -outlines[[x]]$xl/scale.factor[[i]] + offset.x
#           y.poly <- outlines[[x]]$yl/scale.factor[[i]] + offset.y
#         }
#       }
#       p <- add_trace(p, type = "scatter", mode = 'lines', x = x.poly, y = y.poly, showlegend = first.bool, fill = 'toself',fillcolor = col,
#                      hoverinfo = 'skip',line = list(color = 'rgb(0,0,0)'), hoveron = 'points', legendgroup = 'Brain outline', name = 'Brain outline')
#       if (first.bool)
#         first.bool <- FALSE
#     }
#   }
#   
#   for (k in 1:length(col.vect)){
#     cur.soma <- soma.merged[soma.merged$clusters.named == levels(soma.merged$clusters.named)[k],]
#     if (dim(cur.soma)[1] != 0){
#       p <- add_trace(p, type = "scatter", mode = 'markers', data = cur.soma, x = ~somaX, y = ~somaY, hoverinfo='text',
#                      showlegend = TRUE, marker = list(size=10,color=col.vect[k], opacity = 1), name = ~unique(clusters.named),
#                      text = sprintf('Spot: %s<br />Cluster: %s<br />%s', rownames(cur.soma),cur.soma$clusters.named,cur.soma$full.name))
#     }
#   }
#   
#   p <- layout(p, yaxis = list(scaleanchor = "x", range = c(yMax,yMin), title = 'Dorsal-Ventral', showticklabels = FALSE,
#                               zeroline = FALSE, showline = FALSE, shotickslabels = FALSE, showgrid = FALSE),
#               xaxis = list(range = c(xMin,xMax),title = 'Medial-Lateral', showticklabels = FALSE,
#                            zeroline = FALSE, showline = FALSE, shotickslabels = FALSE, showgrid = FALSE),
#               showlegend = TRUE,
#               title = sprintf('Slices: %s<br />AP: %.3f',paste(slice_index,collapse = ' & '),AP),
#               legend = list(traceorder = 'normal+grouped'))
#   
#   return(p)
# }



# plot.2d.atlas.expr <- function(gene, spots.table, st.data, AP = NULL, slice_index = NULL, min.expr = 0, cex = 1, new.window = TRUE, text.plot = TRUE, log.scale = FALSE, normalized.scale = FALSE, win.dim = c(14,14)){
#   
#   #gene(str): gene to plot
#   #spots.table(df)
#   #st.data(df)
#   #AP(float): AP coordinates. If does not match any sections, takes the closest. If two values, plot all sections in-between.
#   #slice_index(str): If specifcied, replaces AP coordinates
#   #min.expr(float): Expression threhsold to plot spot. If between 0 and 1, defines a percentile. Otherwise, defines a number of reads (not log normalized)
#   #cex(float): dimension of text
#   #new.window(bool): force new window
#   #text.plot(bool): show the legend and title 
#   #log.scale(bool): use a logarithmic scale instead
#   #normalized.scale(bool): all expression are a value between 0 and 100 (maximum expression across the whole brain)
#   #win.dim(numeric vector): dimension of the graphical window
#   
#   library(data.table)
#   
#   #Case with two coordinates
#   if(length(AP) == 2){
#     AP.list <- unique(spots.table$AP)
#     AP.list <- sort(AP.list[AP.list >= min(AP) & AP.list <= max(AP)], decreasing = TRUE)
#     for (i in AP.list)
#       plot.2d.atlas.expr(gene,spots.table,st.data,i,NULL,min.expr,cex,TRUE,text.plot,log.scale, normalized.scale)
#   }else{
#     
#     if (!is.null(AP))
#       AP <- select.closest.AP(AP, spots.table)
#     
#     #Case where plotting virtual slice (both dupplicates at a specific coordinates)
#     if (is.null(slice_index))
#       slice_index <- unique(spots.table[spots.table$AP == AP,'slice_index'])
#     
#     #Loading the slice table for correspondance between old and new indexing
#     slices.table <- load.slice.table()
#     
#     old.slice.index <- slices.table[slice_index,'slice_old_id']
#     
#     #Initializing lists that will contain extreme values
#     list.max.x <- NULL
#     list.max.y <- NULL
#     list.min.x <- NULL
#     list.min.y <- NULL
#     scale.factor <- NULL
#     
#     #Keeping registration, dataset and slice.parameters as lists.
#     registration <- NULL
#     dataset <- NULL
#     slice.parameters.list <- NULL
#     
#     #Create a function to generate a continuous color palette
#     rbPal <- colorRampPalette(c('white','red'))
#     
#     dataset <- get.expression.slice(gene, spots.table, st.data, min.expr, log.scale, normalized.scale, slice_index, rbPal)
#     
#     m.val <- dataset[[length(dataset)]]
#     
#     #Plotting one or two slices
#     for (i in 1:length(old.slice.index)){
#       
#       #Loading the appropriates files
#       load(paste(path.slices,old.slice.index[[i]],'Registration/sliceparameters.RData',sep = '/'))
#       load(paste(path.slices,old.slice.index[[i]],'Registration/regi.RData',sep = '/'))
#       
#       #Saving in list
#       registration[[i]] <- get.forward.warpRCPP(regi)
#       slice.parameters.list[[i]] <- slice.parameters
#       
#       #Coordinate of the slice
#       coordinate <- dataset[[i]]$AP[1]
#       
#       #Scale factor from downsizing picture
#       scale.factor[[i]] <- mean(dim(registration[[i]]$transformationgrid$mx)/c(registration[[i]]$transformationgrid$height,
#                                                                                registration[[i]]$transformationgrid$width))
#       
#       #Getting the regions from the atlas
#       numPaths <- registration[[i]]$atlas$numRegions
#       outlines <- registration[[i]]$atlas$outlines
#       
#       #Initializing values to find extreme points
#       list.max.x[[i]] <- 0
#       list.max.y[[i]] <- 0
#       list.min.x[[i]] <- Inf
#       list.min.y[[i]] <- Inf
#       
#       #Finding the extrema
#       for (x in(1:numPaths)) {
#         
#         if (slice.parameters.list[[i]]$right.hemisphere){
#           if (max(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE) > list.max.x[[i]])
#             list.max.x[[i]] <- max(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE)
#           if (max(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE) > list.max.y[[i]])
#             list.max.y[[i]] <- max(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE)  
#           if (min(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE) < list.min.y[[i]])
#             list.min.y[[i]] <- min(outlines[[x]]$yr/scale.factor[[i]],na.rm = TRUE)
#           if (min(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE) < list.min.x[[i]])
#             list.min.x[[i]] <- min(outlines[[x]]$xr/scale.factor[[i]],na.rm = TRUE)    
#         }else{
#           if (max(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE) > list.max.x[[i]])
#             list.max.x[[i]] <- max(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE)
#           if (max(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE) > list.max.y[[i]])
#             list.max.y[[i]] <- max(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE)  
#           if (min(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE) < list.min.y[[i]])
#             list.min.y[[i]] <- min(outlines[[x]]$yl/scale.factor[[i]],na.rm = TRUE)
#           if (min(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE) < list.min.x[[i]])
#             list.min.x[[i]] <- min(outlines[[x]]$xl/scale.factor[[i]],na.rm = TRUE)         
#         }
#       }
#       
#       if (!slice.parameters.list[[i]]$right.hemisphere){
#         tmp <- list.min.x[[i]]
#         list.min.x[[i]] <- -list.max.x[[i]]
#         list.max.x[[i]] <- -tmp
#       }
#       
#     }
#     
#     #Dimensions and limits of the plot, a margin of 5% fron the extrema is added
#     width <- list.max.x[[1]] - list.min.x[[1]]
#     height <- list.max.y[[1]] - list.min.y[[1]]
#     margin.x <- 0.4 * width
#     margin.y <- 0.05 * height
#     xMin <- list.min.x[[1]] - margin.x
#     xMax <- list.max.x[[1]] + margin.x
#     yMin <- list.min.y[[1]] - margin.y
#     yMax <- list.max.y[[1]] + margin.y
#     
#     #Computing offset to align both slices
#     offset.x <- list.max.x[1] - list.max.x[2]
#     offset.y <- list.max.y[1] - list.max.y[2]
#     
#     # #Initializing plotting
#     if (new.window)
#       win.graph(width = win.dim[1],height = win.dim[2])
#     
#     if (text.plot)
#       title.str <- paste("Bregma:", registration[[1]]$coordinate,"mm\nSlices: ",paste(slice_index, collapse = ' & '),'\nGene:',gene)
#     else
#       title.str <- ''
#     
#     plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin),
#          xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "",
#          col = 0, main = title.str, font.main = 1)
#     
#     if (text.plot){
#       mtext("Dorso-ventral (mm)", side = 2, line = 1.5)
#       mtext("Medio-lateral (mm)", side = 1, line = -1.5)
#     } 
#     soma.merged <- NULL
#     
#     #Plotting for both slices
#     for (i in 1:length(old.slice.index)){
#       
#       #Using the remapped space
#       index <- round(scale.factor[[i]] * cbind(dataset[[i]]$y, dataset[[i]]$x))
#       
#       #Getting the spots coordinates
#       if (slice.parameters.list[[i]]$right.hemisphere){
#         if (i == 1){
#           somaX <- registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]]
#           somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]]
#         }else{
#           somaX <- registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]] + offset.x
#           somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]] + offset.y
#         }
#       }else{
#         if (i == 1){
#           somaX <- -registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]]
#           somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]]
#         }else{
#           somaX <- -registration[[i]]$transformationgrid$mxF[index]/scale.factor[[i]] + offset.x
#           somaY <- registration[[i]]$transformationgrid$myF[index]/scale.factor[[i]] + offset.y
#         }      
#       }
#       
#       soma.merged <- rbind(soma.merged, data.frame(somaX = somaX, somaY = somaY, expr = dataset[[i]]$expr, col = dataset[[i]]$col,stringsAsFactors = FALSE))
#       
#       if (i == 1){
#         
#         #Plotting the regions
#         numPaths <- registration[[i]]$atlas$numRegions
#         outlines <- registration[[i]]$atlas$outlines
#         
#         lapply(1:numPaths, function(x){
#           if (as.character(registration[[i]]$atlas$col[x]) == '#aaaaaa')
#             col = 'black'
#           else if (as.character(registration[[i]]$atlas$col[x]) == '#cccccc')
#             col = 'darkgray'
#           else
#             col = NA
#           
#           
#           if (slice.parameters.list[[i]]$right.hemisphere){
#             #Adding offset to align second slice with the first one
#             if (i == 1)
#               polygon(outlines[[x]]$xr/scale.factor[[i]], outlines[[x]]$yr/scale.factor[[i]],border = "black",col = col)
#             else
#               polygon((outlines[[x]]$xr)/scale.factor[[i]]+ offset.x, (outlines[[x]]$yr)/scale.factor[[i]]+offset.y,border = "black",col = col)
#           }else{
#             if (i == 1)
#               polygon(-outlines[[x]]$xl/scale.factor[[i]], outlines[[x]]$yl/scale.factor[[i]],border = "black",col = col)
#             else
#               polygon((-outlines[[x]]$xl)/scale.factor[[i]]+ offset.x, (outlines[[x]]$yl)/scale.factor[[i]]+offset.y,border = "black",col = col)
#           }
#         })
#       }
#     }
#     
#     points(soma.merged$somaX,soma.merged$somaY, pch=21, bg = soma.merged$col,cex = cex)
#     
#     leg.ys <- (yMin+yMax)/2-0.25*height
#     leg.ye <- (yMin+yMax)/2+0.25*height
#     y.list <- seq(from = leg.ys, to = leg.ye, length.out = 101)
#     col.vect <- rev(rbPal(100))
#     for(i in 1:(length(y.list)-1))
#       rect(xMax-0.25*width,y.list[i],xMax-0.2*width,y.list[i+1],col=col.vect[i], border = NA)
#     
#     rect(xMax-0.25*width,leg.ys,xMax-0.2*width,leg.ye)
#     text(xMax-0.2*width,leg.ye,'0',cex = cex,pos = 4)
#     text(xMax-0.2*width,leg.ys,m.val,cex = cex,pos = 4)
#     if(log.scale)
#       text(xMax-0.225*width,leg.ys,'Log scale',cex = cex,pos=3)
#     if(normalized.scale)
#       str.s <- 'Normalized scale'
#     else
#       str.s <- 'Number of reads'
#     text(xMax-0.15*width,(yMin + yMax)/2,str.s,cex = cex,pos=3,srt = -90)
#   }
# }
