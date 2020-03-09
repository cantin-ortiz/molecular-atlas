source('bin/includes.R')
setwd('figures/fig1/E')

# --------- Loading tables ---------

spots.table <- add.all.acronyms(load.spots.table())

spots.table.main.group <- add.parent.acronym(spots.table, list.acronym.parents = c('CH','BS','CB','grey','fiber tracts','VS','root'))

# --------- Grouping ---------

cort.list <- c('Isocortex','OLF','HPF','CTXsp')
spots.table.cort <- spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,cort.list)) != 0)}),]
spots.table.cort <- add.parent.acronym(spots.table.cort, list.acronym.parents = c(cort.list,'grey','fiber tracts','VS','root'))

cnu.list <- c('STR','PAL')
spots.table.cnu <- spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,cnu.list)) != 0)}),]
spots.table.cnu <- add.parent.acronym(spots.table.cnu, list.acronym.parents = c(cnu.list,'grey','fiber tracts','VS','root'))

ib.list <- c('TH','HY')
spots.table.ib <- spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,ib.list)) != 0)}),]
spots.table.ib <- add.parent.acronym(spots.table.ib, list.acronym.parents = c(ib.list,'grey','fiber tracts','VS','root'))

iso.list <- c("RSP","PTLp","VISC","ACA","AUD","PERI","TEa","AI","ECT","GU","MO","SS","ILA","FRP","ORB","VIS","PL")
spots.table.iso <- spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,iso.list)) != 0)}),]
spots.table.iso <- add.parent.acronym(spots.table.iso, list.acronym.parents = c(iso.list,'grey','fiber tracts','VS','root'))

hipp.list <- c('CA','DG','FC','IG','ENT','PAR','POST','PRE','SUB','ProS','HATA','Apr')
spots.table.hipp <- spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,hipp.list)) != 0)}),]
spots.table.hipp <- add.parent.acronym(spots.table.hipp, list.acronym.parents = c(hipp.list,'grey','fiber tracts','VS','root'))

n.cort <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'CTX')) != 0)}),])[1]
n.cnu <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'CNU')) != 0)}),])[1]
n.ib <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'IB')) != 0)}),])[1]
n.mb <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'MB')) != 0)}),])[1]
n.hb <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'HB')) != 0)}),])[1]
n.isocortex <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'Isocortex')) != 0)}),])[1]
n.hipp <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'HPF')) != 0)}),])[1]
n.hipp.acc <- dim(spots.table[sapply(spots.table$all.acronyms,function(x){return(length(intersect(x,'HIP')) != 0)}),])[1]

# --------- Plotting ---------

pdf('spot-distribution.pdf', width = 20, height = 14)
plot(1,type='n',xlim = c(0,40000), ylim = c(-100,100),
     xlab = '',ylab = '',xaxt = 'n',yaxt= 'n', bty = 'n')

cex.scale <- 1.1

segments(0,-20,0,30,col = 'gray', lty = 2)

segments(n.isocortex*2,0,n.isocortex,30,col = 'gray', lty = 2)
segments(n.isocortex*2,0,n.isocortex*2,-20,col = 'gray', lty = 2)
text(n.isocortex,-20,'Isocortex', col = 'darkgray', cex = 1.5, pos = 3)

segments((n.isocortex + n.hipp)*2,0,n.isocortex + n.hipp,30,col = 'gray', lty = 2)
segments((n.isocortex + n.hipp)*2,0,(n.isocortex + n.hipp)*2,-20,col = 'gray', lty = 2)
text(n.isocortex*2 + n.hipp.acc,-20,'Hippocampal', col = 'darkgray', cex = 1.5, pos = 3)

segments((n.isocortex + n.hipp.acc)*2,0,(n.isocortex + n.hipp.acc)*2,-20,col = 'gray', lty = 2)
text(n.isocortex*2 + n.hipp.acc*2 + n.hipp - n.hipp.acc,-20,'Retro-\nhippocampal', col = 'darkgray', cex = 1.5,pos = 3)

df <- plyr::count(spots.table.main.group$acronym.parent)
df <- df[order(df$freq, decreasing= TRUE),]
df$x <- as.character(df$x)
# spots.table <- add.all.acronyms(spots.table)


cnt <- 0
for(idx in 1:dim(df)[1]){
  df[idx,'start'] <- cnt 
  df[idx,'end'] <- cnt + df[idx,'freq']
  df[idx,'col'] <- color.from.acronym(df[idx,'x'])
  cnt <- cnt + df[idx,'freq']
}

df$full.name <- as.character(name.from.acronym(df$x))
df[df$full.name == 'fiber tracts','full.name']  <- 'Fiber tracts'
df[df$x == 'grey','full.name']  <- 'Unspecifc grey matter'
df[df$x == 'VS','full.name']  <- 'Ventricular systems'
df$count.str <- format(df$freq, nsmall=0, big.mark=" ", trim = TRUE)

offset <- -200 + 100
for (idx in 1:dim(df)[1]){
  rect(df[idx,'start'],80,df[idx,'end'],100,border = 'black',lty = 1, col = df[idx,'col'])
  if (idx <= 3)
    text(mean(c(df[idx,'start'],df[idx,'end'])),90,paste(df[idx,'full.name'],'\n(',df[idx,'count.str'],')',sep=''),cex = cex.scale, col='black')
  else{
    text((df[idx,'end'] + offset),79, paste(df[idx,'full.name'],' (',df[idx,'count.str'],')',sep=''),srt = 90, pos=2, cex = cex.scale, col='black')
    offset <- offset + 200
  }
}

## ------------------

x.offset.bot <- +100 
y.bot <- 28

df.cort <- plyr::count(spots.table.cort$acronym.parent)
df.cort <- df.cort[order(df.cort$freq, decreasing= TRUE),]
df.cort$x <- as.character(df.cort$x)

cnt <- 0
for(idx in 1:dim(df.cort)[1]){
  df.cort[idx,'start'] <- cnt 
  df.cort[idx,'end'] <- cnt + df.cort[idx,'freq']
  df.cort[idx,'col'] <- color.from.acronym(df.cort[idx,'x'])
  cnt <- cnt + df.cort[idx,'freq']
}

df.cort$full.name <- as.character(name.from.acronym(df.cort$x))
df.cort$count.str <- format(df.cort$freq, nsmall=0, big.mark=" ", trim = TRUE)

for (idx in 1:dim(df.cort)[1]){
  rect(df.cort[idx,'start'],30,df.cort[idx,'end'],50,border = 'black',lty = 1, col = df.cort[idx,'col'])
  text(mean(c(df.cort[idx,'start'],df.cort[idx,'end']))+x.offset.bot,y.bot, paste(df.cort[idx,'full.name'],' (',df.cort[idx,'count.str'],')',sep=''),srt = 45, pos=2, cex = cex.scale, col = 'black')
}

rect(df.cort[idx,'end'],30,n.cort,50,border = 'black',lty = 1, col = color.from.acronym('CTX'))
text(mean(c(df.cort[idx,'end']+200,n.cort)) + x.offset.bot,y.bot, paste('Unspecific cortex',' (',-df.cort[idx,'end']+n.cort,')',sep=''),srt = 45, pos=2, cex = cex.scale, col = 'black')
## ------------------

df.cnu <- plyr::count(spots.table.cnu$acronym.parent)
df.cnu <- df.cnu[order(df.cnu$freq, decreasing= TRUE),]
df.cnu$x <- as.character(df.cnu$x)

cnt <- n.cort
for(idx in 1:dim(df.cnu)[1]){
  df.cnu[idx,'start'] <- cnt
  df.cnu[idx,'end'] <- cnt + df.cnu[idx,'freq']
  df.cnu[idx,'col'] <- color.from.acronym(df.cnu[idx,'x'])
  cnt <- cnt + df.cnu[idx,'freq']
}

df.cnu$full.name <- as.character(name.from.acronym(df.cnu$x))
df.cnu$count.str <- format(df.cnu$freq, nsmall=0, big.mark=" ", trim = TRUE)

for (idx in 1:dim(df.cnu)[1]){
  rect(df.cnu[idx,'start'],30,df.cnu[idx,'end'],50,border = 'black',lty = 1, col = df.cnu[idx,'col'])
  text(mean(c(df.cnu[idx,'start'],df.cnu[idx,'end'])) + x.offset.bot,y.bot, paste(df.cnu[idx,'full.name'],' (',df.cnu[idx,'count.str'],')',sep=''),srt = 45, pos=2, cex = cex.scale, col='black')
}
rect(cnt,30,df[df$x == 'CH','freq'],50, border = 'black',lty = 1, col = color.from.acronym('CH'))

## --------------------------


df.ib <- plyr::count(spots.table.ib$acronym.parent)
df.ib <- df.ib[order(df.ib$freq, decreasing= TRUE),]
df.ib$x <- as.character(df.ib$x)

cnt <- n.cort + n.cnu
for(idx in 1:dim(df.ib)[1]){
  df.ib[idx,'start'] <- cnt
  df.ib[idx,'end'] <- cnt + df.ib[idx,'freq']
  df.ib[idx,'col'] <- color.from.acronym(df.ib[idx,'x'])
  cnt <- cnt + df.ib[idx,'freq']
}

df.ib$full.name <- as.character(name.from.acronym(df.ib$x))
df.ib$count.str <- format(df.ib$freq, nsmall=0, big.mark=" ", trim = TRUE)

for (idx in 1:dim(df.ib)[1]){
  rect(df.ib[idx,'start'],30,df.ib[idx,'end'],50,border = 'black',lty = 1, col = df.ib[idx,'col'])
  text(mean(c(df.ib[idx,'start'],df.ib[idx,'end'])) + x.offset.bot,y.bot, paste(df.ib[idx,'full.name'],' (',df.ib[idx,'count.str'],')',sep=''),srt = 45, pos=2, cex = cex.scale, col='black')
}


## ------------------
x.offset.bot <- 200
y.bot <- -40

df.iso <- plyr::count(spots.table.iso$acronym.parent)
df.iso <- df.iso[order(df.iso$freq, decreasing= TRUE),]
df.iso$x <- as.character(df.iso$x)

cnt <- 0
for(idx in 1:dim(df.iso)[1]){
  df.iso[idx,'start'] <- cnt
  df.iso[idx,'end'] <- cnt + df.iso[idx,'freq']
  df.iso[idx,'col'] <- color.from.acronym(df.iso[idx,'x'])
  cnt <- cnt + df.iso[idx,'freq']
}

df.iso$full.name <- as.character(name.from.acronym(df.iso$x))
df.iso$count.str <- format(df.iso$freq, nsmall=0, big.mark=" ", trim = TRUE)

df.iso$start <- df.iso$start * 2
df.iso$end <- df.iso$end * 2

offset <- -100
for (idx in 1:dim(df.iso)[1]){
  rect(df.iso[idx,'start'],y.bot,df.iso[idx,'end'],y.bot+20,border = 'black',lty = 1, col = df.iso[idx,'col'])
  if (idx <= 11)
    text(mean(c(df.iso[idx,'start'],df.iso[idx,'end'])) + x.offset.bot,y.bot-2, paste(df.iso[idx,'full.name'],' (',df.iso[idx,'count.str'],')',sep=''),srt = 45, pos=2, cex = cex.scale, col='black')
  else{
    text(mean(c(df.iso[idx,'start'],df.iso[idx,'end'])) + x.offset.bot + offset,y.bot-2, paste(df.iso[idx,'full.name'],' (',df.iso[idx,'count.str'],')',sep=''),srt = 90, pos=2, cex = 1, col='black')
    offset <- offset + 50
  }
}

## ------------------


df.hipp <- plyr::count(spots.table.hipp$acronym.parent)
rownames(df.hipp) <- df.hipp$x
df.hipp <- df.hipp[c('CA','DG','FC','IG','ENT','SUB','POST','PRE','PAR'),]


df.hipp$x <- as.character(df.hipp$x)



cnt <- n.isocortex
for(idx in 1:dim(df.hipp)[1]){
  df.hipp[idx,'start'] <- cnt
  df.hipp[idx,'end'] <- cnt + df.hipp[idx,'freq']
  df.hipp[idx,'col'] <- color.from.acronym(df.hipp[idx,'x'])
  cnt <- cnt + df.hipp[idx,'freq']
}

df.hipp$full.name <- as.character(name.from.acronym(df.hipp$x))
df.hipp$count.str <- format(df.hipp$freq, nsmall=0, big.mark=" ", trim = TRUE)

df.hipp$start <- df.hipp$start * 2
df.hipp$end <- df.hipp$end * 2

df.hipp$offset <- 0
df.hipp['FC','offset'] <- -150
df.hipp['IG','offset'] <- 150
df.hipp['POST','offset'] <- -175
df.hipp['PAR','offset'] <- 175


for (idx in 1:dim(df.hipp)[1]){
  rect(df.hipp[idx,'start'],y.bot,df.hipp[idx,'end'],y.bot+20,border = 'black',lty = 1, col = df.hipp[idx,'col'])
  # if (idx <= 11)
    # text(mean(c(df.hipp[idx,'start'],df.hipp[idx,'end'])) + x.offset.bot,y.bot-2, paste(df.hipp[idx,'full.name'],' (',df.hipp[idx,'count.str'],')',sep=''),srt = 45, pos=2, cex = cex.scale, col='black')
  # else
    text(mean(c(df.hipp[idx,'start'],df.hipp[idx,'end'])) + x.offset.bot + df.hipp[idx,'offset'],y.bot-2, paste(df.hipp[idx,'full.name'],' (',df.hipp[idx,'count.str'],')',sep=''),srt = 90, pos=2, cex = 1, col='black')
}


## -----------------

s.ib <- n.cort+n.cnu
e.ib <- n.cort+n.cnu+n.ib
s.mb <- e.ib
e.mb <- s.mb + n.mb
s.hb <- e.mb
e.hb <- s.hb + n.hb

rect(1,55,n.cort,75,border = 'black',lty = 1, col = color.from.acronym('CTX'))
rect(n.cort,55,n.cort + n.cnu,75,border = 'black',lty = 1, col = color.from.acronym('CNU'))
rect(s.ib,55,e.ib,75,border = 'black',lty = 1, col = color.from.acronym('IB'))
rect(s.mb,55,e.mb,75,border = 'black',lty = 1, col = color.from.acronym('MB'))
rect(s.hb,55,e.hb,75,border = 'black',lty = 1, col = color.from.acronym('HB'))

text(mean(c(1,n.cort)),65,paste('Cerebral cortex','\n(',format(n.cort, nsmall=0, big.mark=" ", trim = TRUE),')',sep=''),cex = cex.scale, col='black')
text(mean(c(n.cort,s.ib)),65,paste('Cerebral nuclei','\n(',format(n.cnu, nsmall=0, big.mark=" ", trim = TRUE),')',sep=''),cex = cex.scale, col='black')
text(mean(c(s.ib,e.ib)),65,paste('Interbrain','\n(',format(n.ib, nsmall=0, big.mark=" ", trim = TRUE),')',sep=''),cex = cex.scale, col='black')
text(mean(c(s.mb,e.mb)),65,paste('Midbrain','\n(',format(n.mb, nsmall=0, big.mark=" ", trim = TRUE),')',sep=''),cex = cex.scale, col='black')
text(mean(c(s.hb,e.hb)),65,paste('Hindbrain','\n(',format(n.hb, nsmall=0, big.mark=" ", trim = TRUE),')',sep=''),cex = cex.scale, col='black')

dev.off()
