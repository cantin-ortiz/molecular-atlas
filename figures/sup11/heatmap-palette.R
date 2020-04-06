#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/sup11')


#------------- FUNCTIONS -----------------

plot.phylo.modified <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
                                 show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
                                 edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
                                 adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
                                 label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
                                 direction = "rightwards", lab4ut = NULL, tip.color = "black", 
                                 plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
                                 align.tip.label = FALSE, new = FALSE, xx.offset = 0, yy.offset = 0, ...) 
{
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[, 
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label) 
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.ultrametric(x)) 
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
      is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) 
      yy <- .nodeHeight(z$edge, Nedge, yy)
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                                 z$edge.length)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin)) 
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) 
        y.lim <- c(1, Ntip)
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  
  if(new){
    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
                 ylab = "", axes = FALSE, asp = asp, ...)
  }
  if (plot) {
    if (is.null(adj)) 
      adj <- if (phyloORclado && direction == "leftwards") 
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt 
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      phylogram.plot(x$edge, Ntip, Nnode, xx+xx.offset, yy+yy.offset, horizontal, 
                     edge.color, edge.width, edge.lty)
    }
    else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                          edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1) 
        edge.color
      else "black"
      rootw <- if (length(edge.width) == 1) 
        edge.width
      else 1
      rootlty <- if (length(edge.lty) == 1) 
        edge.lty
      else 1
      if (type == "fan") {
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw, 
                 lty = rootlty)
      }
      else {
        switch(direction, rightwards = segments(0, yy[ROOT], 
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, 
                                                lty = rootlty), leftwards = segments(xx[ROOT], 
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], 
                                                                                     col = rootcol, lwd = rootw, lty = rootlty), 
               upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, 
                                  col = rootcol, lwd = rootw, lty = rootlty), 
               downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                                    yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw, 
                                    lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label)) 
        underscore <- TRUE
      if (!underscore) 
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]), 
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip], 
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip], 
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]), 
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, 
                   lty = align.tip.label.lty)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label, 
             adj = adj, font = font, srt = srt, cex = cex, 
             col = tip.color)
      }
      else {
        angle <- if (type == "unrooted") 
          XY$axe
        else atan2(yy[1:Ntip], xx[1:Ntip])
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted") 
            "horizontal"
          else "axial"
        }
        else match.arg(lab4ut, c("horizontal", "axial"))
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 * 
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] * 
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj * 
                 cex, x$tip.label, adj = c(adj, 0), font = font, 
               srt = srt, cex = cex, col = tip.color)
        }
        else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips, 
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                                 x$tip.label[i], font = font[i], cex = cex[i], 
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label) 
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
           x$node.label, adj = adj, font = font, srt = srt, 
           cex = cex)
  }
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  invisible(L)
}


#------------- PARAMETERS -----------------

path.cl.266  <- paste(path.matrices, '181-clusters-266-genes.tsv', sep = '/')
path.ic.mat.266 <- paste(path.matrices, 'ic-matrix-266.tsv', sep = '/')
method = 'ward.D2'
breaks <- seq(from = -0.01, to = 1.01, by = 0.01)
col <- colorRampPalette(brewer.pal(9,'YlOrRd'))(length(breaks)-1)
length.colors <- 4

col.hipp <- '#ffa6f8'
col.block <- '#f3ffad'
x.offset <- 38.5

#------------- LOADINGS -----------------

raw.cluster.266 <- read.table(path.cl.266, row.names = 1)
spots.table <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), cl.file, min.cluster.size = 10)
spots.table.266 <- append.cluster.to.spots.table(add.parent.acronym(load.spots.table()), path.cl.266, min.cluster.size = 10)

mat.ic.all <- as.matrix(read.table(paste(path.matrices, 'ic-matrix-scores.tsv', sep = '/'), sep = '\t', row.names = 1, header = T))
mat.ic.all.2 <- mat.ic.all[,ic.kept]
mat.ic.cl.avg.all <- get.ic.cluster.average(spots.table, mat.ic.all.2)
mat.ic.cl.dist <- get.ic.cluster.dist(mat.ic.cl.avg.all)
df.cl.id <- load.df.cl.id.name(cl.id.names.path)

ic.matrix.266 <- as.matrix(read.table(path.ic.mat.266, sep = '\t', row.names = 1, header = T ))
mat.ic.cl.avg.266 <- get.ic.cluster.average(spots.table.266, ic.matrix.266)
mat.ic.cl.dist.266 <- get.ic.cluster.dist(mat.ic.cl.avg.266)
df.cl.id.266 <- unique(spots.table.266[,c('cluster','clusters.named')])

#------------- Initial computation -----------------

common.rows <- intersect(rownames(spots.table), rownames(spots.table.266))

sp <- spots.table[common.rows,]
sp$cluster.all <- sp$cluster
sp$cluster <- NULL
sp$cluster.266 <- spots.table.266[common.rows, 'cluster']
sp$clusters.named.266 <- spots.table.266[common.rows, 'clusters.named']

df.colors.all <- get.color.from.3d.tsne(tsne.3d.path, spots.table)
df.colors.266 <- get.color.from.3d.tsne(tsne.3d.path, spots.table.266)

#Getting clustering tree
cl.tree.all <- get.ic.hiearchical.tree(mat.ic.cl.dist, spots.table, method)
ordered.labels.all <- cl.tree.all$labels[cl.tree.all$order]
ordered.labels.all <- paste('original-', ordered.labels.all, sep = '')

cl.tree.266 <- get.ic.hiearchical.tree(mat.ic.cl.dist.266, spots.table.266, method)
ordered.labels.266 <- cl.tree.266$labels[cl.tree.266$order]
ordered.labels.266 <- paste('palette-', ordered.labels.266, sep = '')

#------------- Computation -----------------

df.count <- as.data.frame(matrix(0, nrow = length(unique(sp$cluster.all)), ncol = length(unique(sp$cluster.266))))
colnames(df.count) <- sprintf('palette-%03d',sort(unique(as.numeric(as.character(sp$cluster.266)))))
rownames(df.count) <- sprintf('original-%03d',sort(unique(as.numeric(as.character(sp$cluster.all)))))

for(cl in unique(as.numeric(as.character(sp$cluster.all)))){
  
  tmp <- subset(sp, cluster.all == cl)
  tmp$cluster.266 <- as.numeric(as.character(tmp$cluster.266))
  cnt <- plyr::count(tmp$cluster.266)
  df.count[sprintf('original-%03d',cl),
           sprintf('palette-%03d',cnt$x)] <- cnt$freq
}

tmp <- unlist(strsplit(rownames(df.count), '-'))
tmp <- as.numeric(tmp[seq(from = 2, to = length(tmp), by = 2)])
df.cl.id <- unique(spots.table[,c('cluster','clusters.named')])
rownames(df.count) <- paste('original-', 
                            mapvalues(tmp, from = as.numeric(as.character(df.cl.id$cluster)), to = as.character(df.cl.id$clusters.named)),
                            sep = '')

tmp <- unlist(strsplit(colnames(df.count), '-'))
tmp <- as.numeric(tmp[seq(from = 2, to = length(tmp), by = 2)])
df.cl.id <- unique(spots.table.266[,c('cluster','clusters.named')])
colnames(df.count) <- paste('palette-', 
                            mapvalues(tmp, from = as.numeric(as.character(df.cl.id$cluster)), to = as.character(df.cl.id$clusters.named)),
                            sep = '')

norm.original <- apply(df.count, 1, function(x){return(x/sum(x))})
norm.palette <- apply(df.count, 2, function(x){return(x/sum(x))})


#------------- Preparing plot ----------------

cur.norm <- t(norm.palette)
cur.norm <- cur.norm[,ordered.labels.all]

#Group by peak value
cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
max.val <- unique(apply(cur.norm, 1, which.max))

for(val in max.val){
  sel <- apply(cur.norm, 1, which.max) == val
  reorder <- order(cur.norm[sel, val], decreasing = T)
  cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
  rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
}

labels.palette <- unlist(strsplit(rownames(cur.norm), 'palette-'))[seq(from = 2, to = 2*nrow(cur.norm), by = 2)]
labels.original <- unlist(strsplit(colnames(cur.norm), 'original-'))[seq(from = 2, to = 2*ncol(cur.norm), by = 2)]

cut.palette.1 <- which(rownames(cur.norm) == 'palette-Isocortex-31')
cut.palette.2 <- which(rownames(cur.norm) == 'palette-Olfactory areas-07')

cut.original.1 <- which(colnames(cur.norm) == 'original-Olfactory areas-10')
cut.original.2 <- which(colnames(cur.norm) == 'original-Isocortex-29')

clusters.original <- c(labels.original[cut.original.1:cut.original.2], labels.original[startsWith(labels.original, 'Hippocampal')]) 
clusters.palette <- c(labels.palette[cut.palette.1:cut.palette.2], labels.palette[startsWith(labels.palette, 'Hippocampal')]) 

all.col <- unique(df.colors.vivid$clusters.colors)
hash.col <- 1:length(all.col)-1000.5
col.used <- c(NA, all.col, col)
breaks.used <- c(-3000,-1500,hash.col,breaks[2:length(breaks)])

#Adding original colors
row.append <- rep(-2000, ncol(cur.norm))
names(row.append) <- colnames(cur.norm)
row.append[paste('original-', clusters.original, sep = '')] <- mapvalues(df.colors.vivid[clusters.original, 'clusters.colors'],
                                                                         from = all.col,
                                                                         to = hash.col,
                                                                         warn_missing = F)
row.append <- as.numeric(row.append)   
plot.norm <- rbind(t(replicate(length.colors, row.append)),cur.norm)

col.append <- rep(-2000, nrow(plot.norm))
names(col.append) <- rownames(plot.norm)
for(cl in clusters.palette){
  cnt <- plyr::count(spots.table[rownames(spots.table.266[spots.table.266$clusters.named == cl,]),'clusters.named'])
  col.append[paste('palette-',cl,sep='')] <- mapvalues(df.colors.vivid[as.character(cnt$x)[which.max(cnt$freq)],'clusters.colors'],
                                                       from = all.col,
                                                       to = hash.col,
                                                       warn_missing = F)
}

col.append <- as.numeric(col.append)
plot.norm <- cbind(replicate(length.colors, col.append),plot.norm)

labels.palette <- c(rep('',nrow(plot.norm) - length(labels.palette)),labels.palette)
labels.original <-c(rep('',ncol(plot.norm) - length(labels.original)),labels.original)

x.rect.palette <- which(is.element(labels.palette, clusters.palette)) + x.offset - 0.5
cols.rect.palette <- ifelse(startsWith(labels.palette[x.rect.palette+0.5-x.offset], 'Hippocampal'),
                            col.hipp,
                            col.block)

y.rect.original <- which(is.element(labels.original, clusters.original)) - 0.5
cols.rect.original <- ifelse(startsWith(labels.original[y.rect.original+0.5], 'Hippocampal'),
                            col.hipp,
                            col.block)
y.rect.original <- y.rect.original - length.colors
#---------- Plot ---------- 

pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram.pdf'),
    useDingbats = F,
    paper = 'a4',
    width = 8.27,
    height = 11.69)
par(mar= c(0.2, 0.2, 0.2, 0.2),
    xpd=TRUE)
plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
     xlim = c(-0.5,nrow(plot.norm)+0.5+x.offset), ylim = c(-0.5-length.colors, ncol(plot.norm) + 0.5), asp = 1)

rect(x.rect.palette,0.5-length.colors, x.rect.palette+1, -25-length.colors, border = NA, col = cols.rect.palette)
rect(0.5+x.offset, y.rect.original, -25+x.offset, y.rect.original + 1, border = NA, col = cols.rect.original)

image(x = 1:nrow(plot.norm)+x.offset,
      y = (1:ncol(plot.norm)) - length.colors,
      z = plot.norm, 
      breaks = breaks.used, 
      col = col.used,
      xaxt= 'n',
      xlab = '',
      yaxt ='n',
      ylab = '',
      add = T)

axis(1,at = 1:nrow(plot.norm)+x.offset, labels = labels.palette, las = 2, lwd = 0, cex.axis = 0.25, pos = 6 - length.colors)
axis(2,at = (1:ncol(plot.norm)) - length.colors, labels = labels.original, las = 2, lwd = 0, cex.axis = 0.25, pos = x.offset+6)

plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.5)

rect(cut.palette.1+x.offset-0.5+length.colors,
     cut.original.1-0.5,
     cut.palette.2+x.offset+0.5+length.colors, 
     cut.original.2+0.5)

# mtext('Clusters from reduced gene palette', side = 3, line = -5)
# mtext('Clusters from all genes', side = 4, line = -1)
text(x.offset+nrow(plot.norm)+2, (ncol(plot.norm)+1)/2, 'Clusters from all genes', srt = -90, pos= 3)
text(x.offset+(nrow(plot.norm)+1)/2, ncol(plot.norm) - length.colors, 'Clusters from reduced gene palette', pos= 3)
dev.off()


# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,sort(colnames(cur.norm), decreasing = T)]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# # for(val in max.val){
# #   sel <- apply(cur.norm, 1, which.max) == val
# #   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][order(cur.norm[sel, val], decreasing = T)],]
# # }
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-name.pdf'),
#     useDingbats = F)
# image(x = 1:nrow(cur.norm),
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '')
# 
# axis(1,at = 1:nrow(cur.norm), labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2)
# dev.off()

# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# #Coloring tree if required
# clus.col <- cutree(cl.tree.all, 13)
# names(clus.col) <- paste('original', names(clus.col), sep = '-')
# clus.col <- clus.col[ordered.labels.all]
# x.splits <- which(as.logical(diff(clus.col) != 0)) + 0.5
# 

# 
# x.offset <- 38.5
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram-lines.pdf'),
#     useDingbats = F)
# 
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(cur.norm)+0.5+x.offset), ylim = c(-0.5, ncol(cur.norm) + 0.5), asp = 1)
# 
# image(x = 1:nrow(cur.norm)+x.offset,
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(cur.norm)+x.offset, labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = 4)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = x.offset+6)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.5)
# 
# segments(x.offset+0.5, x.splits, x.offset+nrow(cur.norm)+0.5, x.splits, lwd = 0.5, col = 'black', lend='butt')
# segments(x.offset+0.5, x.splits, -50, x.splits, lwd = 0.5, lty = 2, col = 'black', lend='butt')
# 
# dev.off()
# 

# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# #Group by peak value
# cur.norm <- cur.norm[order(apply(cur.norm, 1, which.max), decreasing = T),]
# max.val <- unique(apply(cur.norm, 1, which.max))
# for(val in max.val){
#   sel <- apply(cur.norm, 1, which.max) == val
#   reorder <- order(cur.norm[sel, val], decreasing = T)
#   cur.norm[which(sel),] <- cur.norm[rownames(cur.norm)[sel][reorder],]
#   rownames(cur.norm)[which(sel)] <- rownames(cur.norm)[which(sel)][reorder]
# }
# 
# #Coloring tree if required
# clus.col <- cutree(cl.tree.all, 13)
# names(clus.col) <- paste('original', names(clus.col), sep = '-')
# clus.col <- clus.col[ordered.labels.all]
# x.splits <- which(as.logical(diff(clus.col) != 0)) + 0.5
# m.all.val <- apply(cur.norm, 1, which.max)
# y.splits <- sapply(x.splits, function(x){return(which(x > m.all.val)[1])}) - 0.5
# 

# 
# x.offset <- 38.5
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram-lines-vertical.pdf'),
#     useDingbats = F)
# 
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(cur.norm)+0.5+x.offset), ylim = c(-0.5, ncol(cur.norm) + 0.5), asp = 1)
# 
# image(x = 1:nrow(cur.norm)+x.offset,
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(cur.norm)+x.offset, labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = 8)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = x.offset+6)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.5)
# 
# segments(x.offset+0.5, x.splits, x.offset+nrow(cur.norm)+0.5, x.splits, lwd = 0.25, col = 'black', lend='butt')
# segments(x.offset+0.5, x.splits, -50, x.splits, lwd = 0.25, lty = 2, col = 'black', lend='butt')
# 
# segments(y.splits+x.offset,0.5, y.splits+x.offset, ncol(cur.norm)+0.5, lwd = 0.25, lty = 1, col = 'black', lend='butt')
# segments(y.splits+x.offset,0.5, y.splits+x.offset, -60, lwd = 0.25, lty = 2, col = 'black', lend='butt')
# 
# dev.off()
# 

# 
# cur.norm <- t(norm.palette)
# cur.norm <- cur.norm[,ordered.labels.all]
# 
# cur.norm <- cur.norm[ordered.labels.266,]
# 

# 
# x.offset <- 38.5
# y.offset <- -39
# 
# pdf(generate.appropriate.file.name('heatmap-palette-ordered-by-dendogram-both.pdf'),
#     useDingbats = F)
# 
# plot(1, type = 'n',xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
#      xlim = c(-0.5,nrow(cur.norm)+0.5+x.offset), ylim = c(-34, ncol(cur.norm) + 0.5), asp = 1)
# 
# image(x = 1:nrow(cur.norm)+x.offset,
#       y = 1:ncol(cur.norm),
#       z = cur.norm, 
#       breaks = breaks, 
#       col = col,
#       xaxt= 'n',
#       xlab = '',
#       yaxt ='n',
#       ylab = '',
#       add = T)
# 
# axis(1,at = 1:nrow(cur.norm)+x.offset, labels = rownames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = 8)
# axis(2,at = 1:ncol(cur.norm), labels = colnames(cur.norm), las = 2, lwd = 0, cex.axis = 0.2, pos = x.offset+7)
# 
# plot.phylo.modified(as.phylo(cl.tree.all), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.25,
#                     xx.offset = 3)
# plot.phylo.modified(as.phylo(cl.tree.266), type = 'phylogram', cex = 1, show.tip.label = FALSE, new  = FALSE, edge.width = 0.25, direction = 'upwards',
#                     xx.offset = x.offset, yy.offset = y.offset)
# dev.off()
# 
# ##------------ LEGEND ------------ 
# 
# breaks.legend <- (breaks + 0.01) / 1.02
# 
# pdf('legend.pdf', useDingbats = F)
# plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
# rect(breaks.legend[1:(length(breaks.legend)-1)],0.4,breaks.legend[2:length(breaks.legend)],0.6, col = col, border = NA)
# rect(0,0.4,1,0.6)
# axis(1, pos = 0.3)
# dev.off()
# 
# 
