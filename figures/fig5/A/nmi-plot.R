#------------------- Includes  ------------------- 

source('bin/includes.R')
setwd('figures/fig5/A')

#------------------- Parameters -------------------

library(NMI)
main.dir <- paste(path.matrices, 'nmi-plot', sep = '/')
dirs.to.explore <- c('svm','ic')

#------------------- Initial loadings -------------------

spots.table <- add.parent.acronym(load.spots.table())

#-------------- Computing NMI to compare the clustering quality ------------

l.vector.mi <- NULL
cnt <- 0
cur.dir <- getwd()

for(d in dirs.to.explore){

  sub.dirs <- list.files(path = paste(main.dir,d, sep = '/'), no.. = TRUE, full.names = TRUE)
  
  for(sd in sub.dirs){
    
    setwd(sd)
    
    clustering.files <- list.files(pattern = '^(clusters).*(.tsv)$')
    
    sp.2 <- append.cluster.to.spots.table(spots.table, cl.file, min.cluster.size = 10)
  
    vector.nmi <- numeric(length(clustering.files))
    
    for(i in 1:length(clustering.files)){
      
      sp.1 <- append.cluster.to.spots.table(spots.table, clustering.files[i], min.cluster.size = 10)
      
      sp.1.intersect <- sp.1[intersect(rownames(sp.1), rownames(sp.2)),]
      sp.2.intersect <- sp.2[intersect(rownames(sp.1), rownames(sp.2)),]
      
      m1 <- data.frame(id = rownames(sp.1.intersect),
                       label = as.numeric(as.character(sp.1.intersect$cluster)),
                       stringsAsFactors = F)
      
      m2 <- data.frame(id = rownames(sp.2.intersect),
                       label = as.numeric(as.character(sp.2.intersect$cluster)),
                       stringsAsFactors = F)
      
      n.val <- NMI(m1,m2)
      vector.nmi[i] <- n.val
      
      tmp.name <- strsplit(clustering.files[i], 'genes-')[[1]][2]
      tmp.name <- strsplit(tmp.name, '[.]')[[1]][1]
      names(vector.nmi)[i] <- tmp.name
    }
    
    cnt <- cnt + 1
    l.vector.mi[[cnt]] <- vector.nmi
    names(l.vector.mi)[[cnt]] <- sd

    setwd(cur.dir)
  }
  
  
  
  vector.ngenes <- numeric(length(vector.nmi))
  
  lf <- list.files(main.dir, pattern = sprintf('^(top-genes-%s-)1*', d), full.names = T)
  
  for(i in 1:length(lf)){
    
    f <- lf[i]
    t <- read.table(f)
    vector.ngenes[i] <- dim(t)[1]
    names(vector.ngenes)[i] <- strsplit(strsplit(lf[i],'genes-')[[1]][2], '[.]')[[1]][1]
    
  }
  
  for(i in 1:length(l.vector.mi)){
    
    names.list <- unlist(lapply(strsplit(names(vector.ngenes), 'genes-'), function(x){return(x[2])}))
    names.list <- unlist(lapply(strsplit(names.list, '[.]'), function(x){return(x[1])}))
    
    # if(typeof(l.vector.mi[[i]]) == 'list')
    if(is.null(dim(l.vector.mi[[i]]))){
      l.vector.mi[[i]] <- data.frame(file = names(l.vector.mi[[i]]),
                                     nmu = as.numeric(l.vector.mi[[i]]),
                                     n.genes = vector.ngenes[names(l.vector.mi[[i]])])
    }
  }
}

for(i in 1:length(l.vector.mi)){
  genes.list.type <- strsplit(names(l.vector.mi)[i], '/')[[1]][1]
}


labs <- c(100,500,1000,5000)

fname <- 'NMI.pdf'
pdf(generate.appropriate.file.name(fname), useDingbats = F)

plot(log10(l.vector.mi[[1]][,3]), l.vector.mi[[1]][,2], pch = 16, 
     xlab = 'Number of genes', ylab = 'Normalized mutual information', xlim = c(2,4), ylim = c(0.4,1), 
     cex.axis = 1.5, cex.lab = 1.5, xaxt = 'n')

points(log10(l.vector.mi[[2]][,3]), l.vector.mi[[2]][,2], pch = 16, col = 'red')

x.sel <- 266
id.seg <- which(l.vector.mi[[2]][,3] == x.sel)
y.seg <- l.vector.mi[[2]][id.seg,2]

legend('topleft', c('SVM weights', 'ICS loads'), pch = 16, col = c('black','red'), cex = 1.5,
       bty = 'n')

segments(log10(266), 0, log10(266), y.seg, lty=2)
segments(0, y.seg, log10(266), y.seg, lty=2)

axis(1, at  = log10(labs), lab = labs, cex.axis = 1.5)
axis(1, at  = log10(seq(from = 100, to = 1000, by = 100)), lab = NA)
axis(1, at  = log10(seq(from = 1000, to = 10000, by = 1000)), lab = NA)
dev.off()

