#--------------- Loading functions --------------- 

#Loads the meta file in meta.path
sc.load.meta.file <- function(meta.path){
  
  #meta.path(string): path to meta file
  
  meta <- read.table(meta.path, header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_name
  
  return(meta)
}

#Loads a prediction file. Returns meta with predicted cluster and the probability table
sc.load.prediction.file <- function(prediction.path, meta, df.cl.name){
  
  t <- read.table(prediction.path, row.names = 1)
  colnames(t) <- c('predicted',1:(dim(t)[2]-1))
  
  meta <- meta[intersect(rownames(t),rownames(meta)),]
  t <- t[intersect(rownames(t),rownames(meta)),]
  
  meta[rownames(t),'predicted'] <- t$predicted 
  
  t <- t[2:(dim(t)[2])]
  
  colnames(t) <- mapvalues(colnames(t), 
                           from = as.character(df.cl.name$cluster), 
                           to = df.cl.name$clusters.named,
                           warn_missing = F)
  
  meta$predicted.name <- mapvalues(meta$predicted,
                                   from = df.cl.name$cluster,
                                   to = df.cl.name$clusters.named,
                                   warn_missing = F)
  return(list(meta = meta, prob = t))
  
}

#--------------- Accuracy functions --------------- 

#Appends a column "correct.mapping" to a meta containing a "prediction" field 
#based on the prediction being a in "clX" field
sc.append.correct.logical <- function(meta.pred, meta.ground.truth){
  
  r <- intersect(rownames(meta.pred), rownames(meta.ground.truth))
  meta.pred <- meta.pred[r,]
  meta.ground.truth <- meta.ground.truth[r,]
  
  #Columns containing ground truth information
  cols.ground.truth <- colnames(meta.ground.truth)[grep('cl\\d', colnames(meta.ground.truth))]
  
  meta.cur <- meta.ground.truth
  meta.cur$predicted <- meta.pred$predicted
  meta.cur$predicted.name <- meta.pred$predicted.name
  
  meta.cur$correct.mapping <- apply(meta.cur,1,function(x){is.element(as.numeric(x['predicted']),
                                                             as.numeric(x[cols.ground.truth]))})
  
  return(meta.cur)
}

#Prints the acuracy of a prediction based on a ground truth
sc.get.prediction.accuracy <- function(path.prediction, meta, meta.ground.truth, text = NULL, df.cl.name){
  
  #Loading a prediction file
  l.output <- sc.load.prediction.file(path.prediction, meta, df.cl.name)
  meta.pred <- l.output$meta
  rm(l.output)
  
  #Getting the column with correct prediction appended
  meta.correct.logical <- sc.append.correct.logical(meta.pred, meta.ground.truth)
  
  n.correct <- sum(meta.correct.logical$correct.mapping)
  n.total <- dim(meta.correct.logical)[1]
  
  #Get prediction accuracy
  cat(sprintf('File: %s\nAccuracy: %.1f%% (%d/%d correctly predicted cells)\n%s\n\n', 
              path.prediction, n.correct/n.total*100, n.correct, n.total, text))
  
}

#Print the repartition of the prediction from a meta data frame with prediction
sc.plot.prediction <- function(meta.pred, df.cl.name, title = '', new.plot = TRUE, ys = 0.35, ye = 0.65, label.cutoff = 0.03, main = NULL, cex = 0.6){
  
  cnt <- count(meta.pred$predicted.name)
  cnt <- cnt[order(cnt$x),]
  cnt$ratio <- cnt$freq / sum(cnt$freq)
  cnt$x <- as.character(cnt$x)
  
  if(new.plot){
    plot(1, type = 'n', xlim = c(-0.3,1), ylim = c(0,1), #used to be asp = 1 
         xlab = '', ylab = '', xaxt = 'n', yaxt ='n', bty = 'n', main = main)
  }
  
  for(i in 2:(dim(cnt)[1])){
    cnt[i,'start'] <- sum(cnt$ratio[1:(i-1)])
    cnt[i,'end'] <- sum(cnt$ratio[1:i])
  }
  cnt[1,'start'] <- 0
  cnt[1,'end'] <- cnt$ratio[1]
  
  rect(cnt$start, ys, cnt$end, ye, col = df.colors[cnt$x,'clusters.colors'], lwd = 0.25)
  
  labs <- subset(cnt, cnt$ratio > label.cutoff) 
  labs$x.coord <- rowMeans(cbind(labs$start, labs$end))
  labs$y.coord <- mean(c(ys,ye))
  
  #Renaming labels that are too long
  to.rename.id <- grep('Hippocampal', labs$x)
  rename.vect <- unlist(strsplit(labs$x[to.rename.id],'[-]'))
  if(!is.null(rename.vect))
    rename.vect <- paste('Hipp. form-',rename.vect[seq(from = 2, to = length(rename.vect), by = 2)],sep='')
  labs$x[to.rename.id] <- rename.vect  
  
  text(labs$x.coord, 
       labs$y.coord, 
       sprintf('%s (%.1f)%%',labs$x, labs$ratio*100),
       srt = 90, 
       cex = cex)
  text(-0.05, labs$y.coord[1], title, pos = 2, cex = 0.6)
}

#Print the repartition of the prediction from a meta data frame with prediction
sc.plot.multiple.prediction <- function(prediction.paths, meta.ground.truth, df.cl.name, label.cutoff = 0.03, main = NULL){

  p.list <- NULL
  vect.acc <- numeric(length(prediction.paths))
  
  for(i in 1:length(prediction.paths)){
    l <- sc.load.prediction.file(prediction.paths[i],
                                 meta.ground.truth,
                                 df.cl.name)
    p.list[[i]] <- sc.append.correct.logical(l$meta,meta.ground.truth)
    vect.acc[i] <- sum(p.list[[i]]$correct.mapping) / dim(p.list[[i]])[1]
    
  }
  
  yl <- seq(from = 0, to = 1, length.out = (length(p.list) + 1))
  ys <- yl[1:(length(yl)-1)]
  ye <- yl[2:length(yl)]
  
  ord.plot <- order(vect.acc)
  
  first.loop.bool <- T
  cnt <- 1
  
  title.list <- sapply(strsplit(prediction.paths,'/'), function(x){return(x[length(x)-1])})
  
  for(o in ord.plot){
    
    title <- sprintf('%s\nAccuracy: %.1f%%', title.list[o], vect.acc[o]*100)
    
    sc.plot.prediction(p.list[[o]], df.cl.name, new.plot = first.loop.bool, title = title,
                       ys = ys[cnt], ye = ye[cnt], label.cutoff = label.cutoff, main = main)
    
    if(first.loop.bool)
      first.loop.bool <- F
    
    cnt <- cnt + 1
  }

}

#--------------- Alluvial plots --------------- 

#Returns a data frame suitable for plotting a Sankey diagram. Returns a list with df.sankey and sankey.colors for plotting.
sc.get.df.sankey <- function(cl.sc, meta, meta.ground.truth = NULL, cutoff.display.flow = 3, cutoff.display.ST = 3){
  
  #cl.sc(char): single cells cluster that should be included in the diagram
  #meta: a meta file that contains the cells with the ST predicted clusters
  #meta.ground.truth: a ground truth file for adding if a remapping is correct/incorrect. If null, "correct" column is not added
  #cutoff.display.flow(int): links smaller than this will be discarded
  #cutoff.display.ST(int): ST clusters smaller than this(and related links) will be deleted
  #use.dissected.region(logical): if TRUE, dissected region is use to decide the brain region. If false, cluster name is used instead. If NULL, uses cluster name if possible
  
  #Preallocating
  df.sankey.flow <- data.frame(sc.cluster.name = character(0),
                               sc.cluster.sankey.id = numeric(0),
                               st.cluster.name = character(0),
                               st.cluster.sankey.id = numeric(0),
                               flow = numeric(0),
                               stringsAsFactors = F)
  
  #Filling
  row <- 0
  
  for(i in 1:length(cl.sc)){
    
    cur.cl <- cl.sc[i]
    cnt <- count(subset(meta, cell_cluster == cur.cl, 'predicted.name'))
    
    for(j in 1:dim(cnt)[1]){
      
      row <- row + 1
      df.sankey.flow[row,'sc.cluster.name'] <- cur.cl
      df.sankey.flow[row,'sc.cluster.sankey.id'] <- (i-1)
      df.sankey.flow[row,'st.cluster.name'] <- cnt[j,'predicted.name']
      df.sankey.flow[row,'flow'] <- cnt[j,'freq']
      
    }
  }
  
  #Vector with correspondance nodes names and node IDs
  st.clusters.sankey <- sort(unique(df.sankey.flow$st.cluster.name))
  sankey.nodes <- c(cl.sc, st.clusters.sankey)
  sankey.nodes.id <- 0:(length(sankey.nodes)-1)
  
  #Adding ID
  df.sankey.flow$st.cluster.sankey.id <- as.numeric(mapvalues( df.sankey.flow$st.cluster.name, 
                                                               from = sankey.nodes,
                                                               to = sankey.nodes.id,
                                                               warn_missing = FALSE))
  
  #Couting the flow for every connection
  cnt <- sapply(unique(df.sankey.flow$st.cluster.name), function(x){return(sum(subset(df.sankey.flow,
                                                                                      st.cluster.name == x,
                                                                                      'flow')))})
  cnt <- as.data.frame(cnt)
  cnt$discarded <- cnt$cnt < cutoff.display.ST
  
  #Deleting the ST that are too small
  df.sankey.flow <- df.sankey.flow[!is.element(df.sankey.flow$st.cluster.name,
                                               rownames(cnt)[cnt$discarded]),]
  
  #Deleting the flows that are too small
  df.sankey.flow <- subset(df.sankey.flow, flow >= cutoff.display.flow)

  #Correct remapping or not
  if(!is.null(meta.ground.truth)){

    meta.correct.bool <- sc.append.correct.logical(meta, meta.ground.truth)
    df.remappings.correct <- unique(meta.correct.bool[,c('cell_cluster','predicted.name','correct.mapping')])

    
    
  #   #Case where each cluster is either correctly either wrongly remapped (glutamatergic)
  #   if(dim(unique(df.remappings.correct[,c('cell_cluster','predicted.name')]))[1] == dim(df.remappings.correct)[1]){
  # 
  #     df.sankey.flow$correct <- unlist(
  #       apply(df.sankey.flow,1,function(x){return(subset(df.remappings.correct, 
  #                                                        cell_cluster == x['sc.cluster.name'] & predicted.name == x['st.cluster.name'],
  #                                                        'correct.mapping'))}))
  #   
  #   #Case where looking at each individual cell dissected_region
  #   }else{
      
    
    
    #Subdividing the links depending on correct / incorrect
    meta.correct.bool <- sc.append.correct.logical(meta, meta.ground.truth)
    df.sankey.flow.2 <- df.sankey.flow
    
    #Preallocating
    df.sankey.flow <- data.frame(sc.cluster.name = character(0),
                                 sc.cluster.sankey.id = numeric(0),
                                 st.cluster.name = character(0),
                                 st.cluster.sankey.id = numeric(0),
                                 flow = numeric(0),
                                 correct = logical(0),
                                 stringsAsFactors = F)
    
    cnt <- 0
    
    #Looping through links
    for(i in 1:(dim(df.sankey.flow.2)[1])){
      
      cur.cells <- subset(meta.correct.bool, 
                          cell_cluster == df.sankey.flow.2[i,'sc.cluster.name']  & predicted.name == df.sankey.flow.2[i,'st.cluster.name'])
      n.corr <- sum(cur.cells$correct.mapping)
      n.wrong <- sum(!cur.cells$correct.mapping)
      
      if(n.corr > 0){
        cnt <- cnt + 1
        df.sankey.flow[cnt,1:4] <- df.sankey.flow.2[i,1:4]
        df.sankey.flow[cnt,'flow'] <- n.corr
        df.sankey.flow[cnt,'correct'] <- TRUE
      }
      
      if(n.wrong > 0){
        cnt <- cnt + 1
        df.sankey.flow[cnt,1:4] <- df.sankey.flow.2[i,1:4]
        df.sankey.flow[cnt,'flow'] <- n.wrong
        df.sankey.flow[cnt,'correct'] <- FALSE
      }
    }
  }
  
  
# }

  #Coloring
  cp <- brewer.pal(12,'Set3')
  cp <- rep(cp,10)
  sankey.color <- c(cp[1:length(cl.sc)],
                    mapvalues(sankey.nodes[(length(cl.sc)+1):length(sankey.nodes)],
                              from = rownames(df.colors),
                              to = df.colors$clusters.colors,
                              warn_missing = F))
  
  #Naming the colors
  names(sankey.color) <- sankey.nodes
  
  return(list(df.sankey = df.sankey.flow, sankey.color = sankey.color))
  
}

#Plot and alluvial diagrams. 
sc.plot.alluvial <- function(df.sankey, sankey.color, col.correct = 'green', col.incorrect = 'red', cex = 0.3, main = NULL, cw = 0.15, show.proba = T, axis_labels = c('SC cluster', 'ST cluster')){

  internal.sc.plot.alluvial <-  function (..., freq, col = "gray", border = 0, layer, hide = FALSE, 
                                  alpha = 0.5, gap.width = 0.05, xw = 0.2, cw = 0.15, blocks = TRUE, 
                                  ordering = NULL, axis_labels = NULL, cex = par("cex"), cex.axis = par("cex.axis"), 
                                  rect.col = NULL, main = NULL, prob = NULL) 
  {
    p <- data.frame(..., freq = freq, col, alpha, border, hide, 
                    stringsAsFactors = FALSE)
    np <- ncol(p) - 5
    if (!is.null(ordering)) {
      stopifnot(is.list(ordering))
      if (length(ordering) != np) 
        stop("'ordering' argument should have ", np, " components, has ", 
             length(ordering))
    }
    n <- nrow(p)
    if (missing(layer)) {
      layer <- 1:n
    }
    p$layer <- layer
    d <- p[, 1:np, drop = FALSE]
    p <- p[, -c(1:np), drop = FALSE]
    p$freq <- with(p, freq/sum(freq))
    col <- col2rgb(p$col, alpha = TRUE)
    if (!identical(alpha, FALSE)) {
      col["alpha", ] <- p$alpha * 256
    }
    p$col <- apply(col, 2, function(x) do.call(rgb, c(as.list(x), 
                                                      maxColorValue = 256)))
    isch <- sapply(d, is.character)
    d[isch] <- lapply(d[isch], as.factor)
    if (length(blocks) == 1) {
      blocks <- if (!is.na(as.logical(blocks))) {
        rep(blocks, np)
      }
      else if (blocks == "bookends") {
        c(TRUE, rep(FALSE, np - 2), TRUE)
      }
    }
    if (is.null(axis_labels)) {
      axis_labels <- names(d)
    }
    else {
      if (length(axis_labels) != ncol(d)) 
        stop("`axis_labels` should have length ", names(d), 
             ", has ", length(axis_labels))
    }
    
    getp <- function(i, d, f, w = gap.width) {
      a <- c(i, (1:ncol(d))[-i])
      if (is.null(ordering[[i]])) {
        o <- do.call(order, d[a])
      }
      else {
        d2 <- d
        d2[1] <- ordering[[i]]
        o <- do.call(order, d2[a])
      }
      x <- c(0, cumsum(f[o])) * (1 - w)
      x <- cbind(x[-length(x)], x[-1])
      gap <- cumsum(c(0L, diff(as.numeric(d[o, i])) != 0))
      mx <- max(gap)
      if (mx == 0) 
        mx <- 1
      gap <- gap/mx * w
      (x + gap)[order(o), ]
    }
    dd <- lapply(seq_along(d), getp, d = d, f = p$freq)
    rval <- list(endpoints = dd)
    op <- par(mar = c(2, 1, 1, 1))
    plot(NULL, type = "n", xlim = c(1 - cw, np + cw), ylim = c(0,1), xaxt = "n", yaxt = "n", 
         xaxs = "i", yaxs = "i", xlab = "", 
         ylab = "", frame = FALSE, main = main)
    ind <- which(!p$hide)[rev(order(p[!p$hide, ]$layer))]
    for (i in ind) {
      for (j in 1:(np - 1)) {
        xspline(c(j, j, j + xw, j + 1 - xw, j + 1, j + 1, j + 1 - xw, j + xw, j) + rep(c(cw, -cw, cw), c(3, 4, 2)),
                c(dd[[j]][i, c(1, 2, 2)], rev(dd[[j + 1]][i, c(1, 1, 2, 2)]), dd[[j]][i, c(1, 1)]),
                shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
                open = FALSE,
                col = p$col[i],
                border = p$border[i])
        # text(j + 2.5*xw, mean(dd[[1]][i,]), paste(i, p[i,'layer'],d[i,1],d[i,2],prob[i],sep='-'), cex = 0.2)
        
        if(show.proba)
          text(j + 0.8*cw, mean(dd[[1]][i,]), sprintf('Prob: %.1f%%',100*prob[i]), cex = cex, pos = 4) 
      }
    }
    
    for (j in seq_along(dd)) {
      
      ax <- lapply(split(dd[[j]], d[, j]), range)
      
      if (blocks[j]) {
        for (k in seq_along(ax)) {
          if(!is.null(rect.col))
            col <- rect.col[names(ax)[k]]
          else
            col <- NA
          rect(j - cw, ax[[k]][1], j + cw, ax[[k]][2], col = col)
        }
      }
      else {
        for (i in ind) {
          
          print(i)
          x <- j + c(-1, 1) * cw
          y <- t(dd[[j]][c(i, i), ])
          w <- xw * (x[2] - x[1])
          xspline(x = c(x[1], x[1], x[1] + w, x[2] - w,
                        x[2], x[2], x[2] - w, x[1] + w,
                        x[1]),
                  y = c(y[c(1,2, 2), 1],
                        y[c(2, 2, 1, 1), 2],
                        y[c(1, 1),1]),
                  shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
                  open = FALSE,
                  col = p$col[i],
                  border = p$border[i])
        }
      }
      for (k in seq_along(ax)) {
        text(j, mean(ax[[k]]), labels = names(ax)[k], cex = cex)
      }
    }
    axis(1, at = rep(c(-cw, cw), ncol(d)) + rep(seq_along(d), 
                                                each = 2), line = 0.5, col = "white", col.ticks = "black", 
         labels = FALSE)
    axis(1, at = seq_along(d), tick = FALSE, labels = axis_labels, 
         cex.axis = cex.axis)
    par(op)
    invisible(rval)
  }
  
  if(!is.null(df.sankey$correct))
    col <- ifelse(df.sankey$correct, col.correct, col.incorrect)
  else
    col <- 'gray'
  
  if(!is.null(df.sankey$prob))
    prob <- df.sankey$prob
  else
    prob <- NULL
  
  internal.sc.plot.alluvial(df.sankey[,c(1,3)], freq = df.sankey$flow, 
                            col = col, cex = cex, main = main,
                            axis_labels = axis_labels, rect.col = sankey.color, 
                            prob = prob, cw = cw)
  
}
#--------------- Genes functions --------------- 

sc.plot.top.genes <- function(data, cluster.id = 1, main = NULL, n.top = 30, show.ylabs = TRUE){
  
  #n.top: number of top genes to plot
  
  
  if(is.null(main))
    main <- sprintf('Cluster: %d', cluster.id)
  
  #Keeping only top genes (and sorting them)
  top.genes <- data[cluster.id,order(abs(data[cluster.id,]), decreasing = TRUE)]
  top.genes <- sort(top.genes[1:n.top])

    #Defining xlim
  diff.lim <- max(top.genes) - min(top.genes)
  xlim = c(min(top.genes) - 0.4*diff.lim, max(top.genes) + 0.05*diff.lim)
  left.axis <- min(top.genes) - 0.05*diff.lim
  
  #Potential step roundings.
  desired.step.roundings <- c(100,250,500,1000,2000,2500,5000,10000)/10000
  
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
  plot(as.numeric(top.genes), 1:n.top, xaxt = 'n', yaxt = 'n', main = main,
       ylim = c(-5,n.top+1), xlim = xlim, xlab = '', ylab = '', bty = 'n', pch = 19)
  par(xpd = TRUE)
  axis(1, pos = 0, at = c(left.axis,xlim[2]), labels = NA, lwd.ticks = 0)
  axis(2, pos = left.axis, at = c(0,n.top+1), labels = NA, lwd.ticks = 0)
  axis(3, pos = n.top+1, at = c(left.axis,xlim[2]), labels = NA, lwd.ticks = 0)
  axis(4, pos = xlim[2], at = c(0,n.top+1), labels = NA, lwd.ticks = 0)
  
  axis(1, pos = 0, at = ticks)
  
  if(show.ylabs){
    axis(2, pos = left.axis, at = 1:n.top, labels = colnames(top.genes), las = 2, cex.axis = 0.65)
  }else{
    text(left.axis, n.top/2, 'Genes', srt = 0, pos = 2)
  }
  
  #0 line
  segments(0,0,0,n.top+1,lty=2)
  
  #xlab
  text(mean(c(left.axis,xlim[2])),-1000,'Gene loading')
}


