# -------------- INCLUDES AND LOADINGS  --------------

source('bin/includes.R')
setwd('figures/fig1/A')
load("~/R/win-library/3.4/wholebrain/data/SAGITTALatlas.RData")
load("~/R/win-library/3.4/wholebrain/data/EPSatlas.RData")

spots.table <- load.spots.table()

# -------------- FUNCTIONS --------------

#Returns x from ap
get.x.from.ap <- function(ap.list){
  
  x.vect <- numeric(length(ap.list))
  cnt <- 0
  for (a in ap.list){
    
    cnt <- cnt + 1
    prop.a <- (ap.bregma.1 - a) / abs(ap.bregma.1 - ap.bregma.2)
    
    x.vect[cnt] <- prop.a * abs(x.bregma.1 - x.bregma.2) + x.bregma.1
    
  }
  
  return(x.vect)
}

# -------------- COMPUTING OBJECTS FOR PLOTTING --------------

#Selecting desired plate
o <- SAGITTALatlas$plates[[7]]

#Correspondance points for AP
x.bregma.1 <- -5250 #x pixels
ap.bregma.1 <- 6  #mm bregma

x.bregma.2 <- 123200 #x pixels
ap.bregma.2 <- -9 #mm bregma

#Range of paths that are fiber / ventricle and everything else
range.fiber <- 157:173
range.ventricle <- 174:175
range.rest <- setdiff(1:length(o@paths),union(range.fiber, range.ventricle))

#Coloring the ranges
col <- rep('white', length(o@paths))
col[157:173] <- 'darkgrey'
col[174:175] <- 'black'

#Counting the polygon points to initialize data frame
cnt.n.points <- 0
for (i in 1:length(o@paths))
  cnt.n.points <- cnt.n.points + length(o@paths[[i]]@x)

n.rows <- cnt.n.points + length(o@paths) - 1

df.poly <- data.frame(x = as.numeric(rep(NA, n.rows)),
                      y = as.numeric(rep(NA, n.rows)),
                      col = as.character(rep(NA, n.rows)),
                      stringsAsFactors = FALSE)

n.rows.seg <- cnt.n.points - length(o@paths)
df.segments <- data.frame(x0 = as.numeric(rep(NA, n.rows.seg)),
                          y0 = as.numeric(rep(NA, n.rows.seg)),
                          x1 = as.numeric(rep(NA, n.rows.seg)),
                          y1 = as.numeric(rep(NA, n.rows.seg)),
                          col = as.character(rep(NA, n.rows.seg)),
                          stringsAsFactors = FALSE)


cnt.row <- 1
cnt.row.seg <- 1

#Going through the paths to plot
for (i in 1:length(o@paths)){
  
  #Detecting splits
  splits <- which(names(o@paths[[i]]@y) == 'move')
  
  #If not fiber or ventricle, just plot outline
  if (is.na(col[i]))
    border <- 'black'
  
  #Otherwise, fill the inside but do not plot outline
  else
    border <- NA
  
  #Adding points to the polygon plotting list
  df.poly[cnt.row:(cnt.row + length(o@paths[[i]]@x) - 1),c('x','y')] <- cbind(as.numeric(o@paths[[i]]@x),
                                                                              as.numeric(o@paths[[i]]@y))
  df.poly[cnt.row:(cnt.row + length(o@paths[[i]]@x) - 1),'col'] <- col[i]
  
  #Moving cursor
  cnt.row <- cnt.row + length(o@paths[[i]]@x) + 1
  
  # Appending length of list as final split
  splits <- c(splits,(length(o@paths[[i]]@y)+1))
  
  #Going through the splits
  for (j in 1:(length(splits)-1)){
    
    sp.start <- splits[j]
    sp.end <- splits[(j+1)] - 1
    
    range.0 <- sp.start:(sp.end-1)
    range.1 <- (sp.start+1):sp.end
    
    df.segments[cnt.row.seg:(cnt.row.seg + length(o@paths[[i]]@y[range.1]) -1),
                c('x0','y0','x1','y1')] <- cbind(
                  o@paths[[i]]@x[range.0],
                  o@paths[[i]]@y[range.0],
                  o@paths[[i]]@x[range.1],
                  o@paths[[i]]@y[range.1])
    df.segments[cnt.row.seg:(cnt.row.seg + length(o@paths[[i]]@y[range.1]) -1),'col'] <- col[i]
    
    cnt.row.seg <- cnt.row.seg + length(o@paths[[i]]@y[range.1])
    
  }
}

df.segments <- df.segments[!is.na(df.segments$x0),]


#Counting spots per section
sp.clusters <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 0)
cnt <- count(sp.clusters$slice_index)
slice.table <- load.slice.table()
cnt$AP <- mapvalues(cnt$x, from = rownames(slice.table), to = as.numeric(slice.table$AP))
cnt$AP <- as.numeric(as.character(cnt$AP))
cnt$start.from <- 0 
for (i in 2:(dim(cnt)[1])){
  
  #Finding spots with same AP before
  pos <-  (which(cnt[1:(i-1),'AP'] == cnt[i,'AP']))

  #If any, starting from it
  if (length(pos) == 1)
    cnt[i,'start.from'] <- cnt[pos,'freq']
}

max.range <- max(cnt$freq + cnt$start.from)

y.range <- abs(diff(as.numeric(o@summary@yscale)))

height.spots.prop.brain <- 0.5
margin.prop.brain <- 0.1

top.y <-as.numeric(o@summary@yscale['ymin']) - margin.prop.brain*y.range
bottom.y <- top.y - height.spots.prop.brain*y.range
y.range.spots <- top.y - bottom.y

cnt$pixels.start <- (cnt$start.from/max.range) * y.range.spots + bottom.y
cnt$pixels.end <- (cnt$start.from + cnt$freq)/max.range * y.range.spots + bottom.y

cnt$color <- 'firebrick4'
cnt$color[cnt$start.from != 0] <- 'firebrick1'
cnt$x.coord <- get.x.from.ap(cnt$AP)

cnt.pixels.per.spots <- (cnt[which.max(cnt$freq),('pixels.end')] - cnt[which.max(cnt$freq),('pixels.start')]) /
                         cnt[which.max(cnt$freq),('freq')]

#Showing which animal corresponds to which section
animal.rect <- unique(cnt[,c('AP','x.coord')])
n <- dim(animal.rect)[1]

slice.table <- load.slice.table()
slice.table.unique.AP <- unique(slice.table[,c('AP','animal')])
animal.rect$animal <- mapvalues(animal.rect$AP, from = slice.table.unique.AP$AP, to = slice.table.unique.AP$animal)
animal.rect$x.start <- animal.rect$x.coord
animal.rect$x.end <- animal.rect$x.coord
animal.rect[1:(n-1), 'x.end'] <- diff(animal.rect$x.coord)/2 + animal.rect[1:(n-1), 'x.coord']
animal.rect[2:n, 'x.start'] <- -diff(animal.rect$x.coord)/2 + animal.rect[2:n, 'x.coord']
animal.rect$color <- mapvalues(animal.rect$animal, from = c('A1','A2','A3'), to = c('red','green','blue'))
animal.rect$x.start[1] <- animal.rect$x.coord[1] - (animal.rect$x.end[1]-animal.rect$x.coord[1])
animal.rect$x.end[n] <- animal.rect$x.coord[n] - (animal.rect$x.start[n]-animal.rect$x.coord[n])
# -------------- GENERATING PLOT --------------

box.col <- 'gray42'
section.color <- 'gray20'
grid.col <- 'gray'
lwd.out <- 1

#Opening plot
pdf('AP-distribution.pdf')
par(xpd = TRUE, mgp = c(3,0.5,0))

#Initializing window
plot(1, type = 'n', xlim = o@summary@xscale, ylim = c(bottom.y, o@summary@yscale['ymax']), asp = 1, bty = "n", 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab='')

#Round AP to plot lines for
ap.round.plot <- -9:6

#Plotting the desired AP
x.round.plot <- get.x.from.ap(ap.round.plot)

#Plotting the grid lines
#First the AP
segments(x0 = x.round.plot, y0 = rep(par("usr")[4], length(x.round.plot)),
         x1 = x.round.plot, y1 = rep(min(cnt$pixels.start),length(x.round.plot)),
         col = grid.col, lwd = 0.5, lend = 'butt')

#Then the spots
spots.scale <- seq(from = 500, to = 1000, by = 500)
segments(x0 = rep(get.x.from.ap(6.75), length(spots.scale)), y0 = spots.scale * cnt.pixels.per.spots + min(cnt$pixels.start),
         x1 = rep(get.x.from.ap(-9.75), length(spots.scale)), y1 = spots.scale * cnt.pixels.per.spots + min(cnt$pixels.start),
         col = grid.col, lwd = 0.5, lend = 'butt')

#Rectangle animals
rect.y.middle <-  mean(c(par("usr")[4],o@summary@yscale['ymax']))
rect.height <- (par("usr")[4] - o@summary@yscale['ymax'])/4
rect.y <- rect.y.middle + c(-rect.height, rect.height)/2
rect(animal.rect$x.start,rect.y[1],animal.rect$x.end,rect.y[2],
     border = NA,
     col = animal.rect$color)


c <- 'white'
polygon(df.poly[df.poly$col == c, 'x'],
        df.poly[df.poly$col == c, 'y'],
        col = c, border = NA)
segments(df.segments[df.segments$col == c, 'x0'],
         df.segments[df.segments$col == c, 'y0'], 
         df.segments[df.segments$col == c, 'x1'],
         df.segments[df.segments$col == c, 'y1'],
         lwd = lwd.out,
         lty = 1,
         lend = 'butt')


for (c in setdiff(unique(df.poly$col),'white')){ 
  polygon(df.poly[df.poly$col == c, 'x'],
          df.poly[df.poly$col == c, 'y'],
          col = c, border = NA)
  segments(df.segments[df.segments$col == c, 'x0'],
           df.segments[df.segments$col == c, 'y0'], 
           df.segments[df.segments$col == c, 'x1'],
           df.segments[df.segments$col == c, 'y1'],
           lend = 'butt',
           lty = 1,
           lwd = lwd.out)
}

#Plotting our sections
segments(cnt$x.coord, rep(par("usr")[4],dim(cnt)[1]), cnt$x.coord, cnt$pixels.end, col = section.color, lty = 5, lend = 'butt')

#Ploting number of spots
segments(cnt$x.coord, cnt$pixels.start, cnt$x.coord, cnt$pixels.end, lwd = 2, col = cnt$color, lend = 'butt')

axis(3, at = x.round.plot, col = box.col, col.ticks = box.col, labels = ap.round.plot, tck  = 0.02, mgp=c(3, .25, 0))
axis(3, at = get.x.from.ap(6.5:-9.5), col = box.col, col.ticks = box.col, labels = FALSE, tck  = 0.015)
axis(3, at = get.x.from.ap(6.75:-9.25), col = box.col, col.ticks = box.col, labels = FALSE, tck  = 0.01)
axis(3, at = get.x.from.ap(6.25:-9.75), col = box.col, col.ticks = box.col, labels = FALSE, tck  = 0.01)

axis(1, at = x.round.plot, pos = min(cnt$pixels.start), col = box.col, col.ticks = box.col, labels = ap.round.plot, tck  = 0.02, mgp=c(3, .25, 0))
axis(1, at = get.x.from.ap(6.5:-9.5), pos = min(cnt$pixels.start), col = box.col, col.ticks = box.col, labels = FALSE, tck  = 0.015)
axis(1, at = get.x.from.ap(6.75:-9.25), pos = min(cnt$pixels.start), col = box.col, col.ticks = box.col, labels = FALSE, tck  = 0.01)
axis(1, at = get.x.from.ap(6.25:-9.75), pos = min(cnt$pixels.start), col = box.col, col.ticks = box.col, labels = FALSE, tck  = 0.01)

axis(2, at = c(min(cnt$pixels.start),par("usr")[4]), pos = get.x.from.ap(6.75), col = box.col, tck = 0, labels = FALSE)
axis(2, at = c(min(cnt$pixels.start) + c(0,spots.scale) * cnt.pixels.per.spots), pos = get.x.from.ap(6.75), 
     col.ticks = box.col, labels = c(0,spots.scale), tck  = 0.02, las = 1)

axis(4, at = c(min(cnt$pixels.start),par("usr")[4]), pos = get.x.from.ap(-9.75), col = box.col, tck = 0, labels = FALSE)

text(labels = 'Number of spots',
     x = get.x.from.ap(8.5),
     y = min(cnt$pixels.start) + max(spots.scale) / 2 * cnt.pixels.per.spots,
     srt = 90,
     cex = 1)

text(labels = 'Animal',
     x = get.x.from.ap(9),
     y = rect.y.middle,
     srt = 0,
     cex = 1,
     pos = 4)

dev.off()

