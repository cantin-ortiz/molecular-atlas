#------------- INCLUDES -----------------

source('bin/includes.R')
setwd('figures/fig2/J')

#------------- PARAMETERS -----------------

N <- c(0.935,-0.192,0.299)
d <- 0.397303260

AP.list <- c(-2.88)
color.list <- c('red', 'blue', 'green')
list.allen <- c('root', 'TH','HY','MB','HB',
                'CB','STR','PAL','Isocortex','OLF',
                'HIP','RHP', 'CTXsp','fiber tracts','VS')

obs <- list(x = 0, y = -0.4, z = 12)
zoom <- 0.6
M <- rbind(c(-1,0,0,0),c(0,0.296854227781296,0.954922616481781,0),c(0,0.954922795295715,-0.296854168176651,0),c(0,0,0,1))
M2 <- rbind(c(0.561396062374115,0,-0.827547252178192,0),c(-0.662099361419678,0.599900841712952,-0.449158608913422,0),c(0.496446281671524,0.800074398517609,0.336781978607178,0),c(0,0,0,1))
zoom2 <- 0.45

M3 <- rbind(c(0.960396707057953,0,0.278636187314987,0),c(0.0686225295066833,0.969198703765869,-0.236526533961296,0),c(-0.270053833723068,0.246280029416084,0.93081521987915,0),c(0,0,0,1))

#------------- Plotting a smoothed atlas -----------------

allen.id <- sapply(list.allen, get.id.from.acronym)
allen.color <- color.from.acronym(list.allen)
allen.color[1] <- 'black'

outlines.list <- list()
mesh.outline <- mesh3d.allen.annot.from.id(get.id.from.acronym('root'))
l.plane <- get.plane.equation(-2, 'horizontal')
outlines.list[[1]] <- get.meshes.list.plane.intersection(list(mesh.outline), 0, l.plane$N, l.plane$d)
l.plane <- get.plane.equation(1.750, 'saggital')
outlines.list[[2]] <- get.meshes.list.plane.intersection(list(mesh.outline), 0, l.plane$N, l.plane$d)


#----------- V1 ----------- 

mesh3d.new.window(T)
par3d(userMatrix = M3)
par3d(zoom = 0.45)
observer3d(obs$x, obs$y, obs$z)

for(i in 1:length(outlines.list)){
  for(j in 1:length(outlines.list[[i]])){

    df.cur.outline <- data.frame(x = numeric(),
                                 y = numeric(),
                                 z = numeric())

    o <- outlines.list[[i]][[j]]
    if(is.null(o))
      next
    if(length(o) == 0)
      next

    for(k in 1:length(o)){

      if(is.null(o[[k]]))
        next

      df.cur.outline <- rbind(df.cur.outline,
                              o[[k]])
    }

    segments3d(df.cur.outline, col = 'black', lwd = 4)

  }
}

rgl.snapshot('cut-position.png')
