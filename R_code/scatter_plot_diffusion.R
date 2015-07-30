library(diffusionMap)
library(gclus)
library(R.utils)
library(randomForest)
library(ggplot2)
library(gtable)
library(gridExtra)

source('R_code/get_training.R')
source('R_code/get_test.R')
source('R_code/util.R')
source('R_code/figure_of_merits.R')

par.defaults <- par(no.readonly=TRUE)
homeDir <- '~/sandbox/miniature-adventure/'
plotDir <- paste(homeDir, 'results/Plots/', sep='')

colors <- c('darkorange', 'seagreen3', 'blue2')
dmap.x.coord <- 10
dmap.y.coord <- 15
dmap.z.coord <- 3


path <- 'results/SIMGEN_PUBLIC_FIT/RBF_test-length/'
path.to.dump <- 'train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP'
fileNameRoot <- 'RBF_test-length'
if (!exists('trainingSet')){
    trainingSet <- get_training(path, fileNameRoot)
}
if (!exists('testSet')){
    testSet <- get_test(path, fileNameRoot, path.to.dump)
}
if (!exists('dmap')){
    dmap <- loadObject(paste(path, 'distance_matrix/',  'diffusion_map.RData', sep=''))
}
if (!exists('sn.rf')){
    sn.rf <- loadObject(paste(homeDir, 'results/SIMGEN_PUBLIC_FIT/', 'RBF_test-lengthrandomForestModel.RData', sep=''))
}
if (exists('dump.file')){
    rm(dump.file)
}
dump.file <- get_dump('train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP')
if (exists('trainingSet.with.predictions')){
    rm(trainingSet.with.predictions)
}
trainingSet.with.predictions <- cbind(trainingSet, sn.rf$predicted, dump.file$GENZ[match(trainingSet$snid, dump.file$CID)])
names(trainingSet.with.predictions) <- c(names(trainingSet.with.predictions)[1:4], 'type.pred', 'GENZ')

if (exists('testSet.with.predictions')){
    rm(testSet.with.predictions)
}
testSet.with.predictions <- cbind(testSet, sn.rf$test$predicted, dump.file$GENZ[match(testSet$snid, dump.file$CID)])
names(testSet.with.predictions) <- c(names(testSet.with.predictions)[1:4], 'type.pred', 'GENZ')

varImpPlot(sn.rf, type=1, pch=21, main='', cex=1, bg='black')

W.Ia.false <- 3

N.Ia.true <- c()
N.Ia.tot <- c()
N.Ia.false <- c()

fom.Ia <- c()
eff.Ia <- c()
pseudo.purity.Ia <- c()
purity.Ia <- c()
criteria.true <- c()
criteria.false <- c()
criteria.false <- c()
# z.seq <- seq(0, max(dump.file$GENZ), by=0.1)
z.seq <- seq(0, 1.25, by=0.25)
for (i in seq(length(z.seq))){
  criteria.true <- which((trainingSet.with.predictions$type == trainingSet.with.predictions$type.pred) &
                           (z.seq[i] <= trainingSet.with.predictions$GENZ) &
                           (trainingSet.with.predictions$GENZ < z.seq[i+1]))
  criteria.false <- which(trainingSet.with.predictions$type != trainingSet.with.predictions$type.pred &
                            (z.seq[i] <= trainingSet.with.predictions$GENZ) &
                            (trainingSet.with.predictions$GENZ < z.seq[i+1]))
  criteria.tot <- which((z.seq[i] <= trainingSet.with.predictions$GENZ) &
                          (trainingSet.with.predictions$GENZ < z.seq[i+1]))
  
  N.Ia.tot <- c(N.Ia.tot, length(which(trainingSet$type[criteria.tot] == 'snIa')))
  N.Ia.true <- c(N.Ia.true, length(which(trainingSet.with.predictions$type.pred[criteria.true] == 'snIa')))
  N.Ia.false <- c(N.Ia.false, length(which(trainingSet.with.predictions$type.pred[criteria.false] == 'snIa')))
 
  fom.Ia <- c(fom.Ia, CFoM.Ia(tail(N.Ia.tot, n=1), tail(N.Ia.true, n=1), tail(N.Ia.false, n=1)))
  eff.Ia <- c(eff.Ia, tail(N.Ia.true, n=1)/tail(N.Ia.tot, n=1))
  pseudo.purity.Ia <- c(pseudo.purity.Ia, tail(N.Ia.true, n=1)/(tail(N.Ia.true, n=1)+W.Ia.false*tail(N.Ia.false, n=1)))
  purity.Ia <- c(purity.Ia, tail(N.Ia.true, n=1)/(tail(N.Ia.true, n=1)+tail(N.Ia.false, n=1)))
}
# plot(z.seq[1:length(z.seq)-1], fom.Ia[1:length(fom.Ia)-1], type='l', lwd=3, 
#      xlab='redshift', ylab='FoM-Ia', ylim=c(0,1))
plot(z.seq[1:length(z.seq)-1], eff.Ia[1:length(eff.Ia)-1], lwd=3, col='seagreen', type='l',
     xlab='redshift', ylab='Ia efficiency', ylim=c(0,1))
# plot(z.seq[1:length(z.seq)-1], pseudo.purity.Ia[1:length(pseudo.purity.Ia)-1], type='l',
#      lwd=3, xlab='redshift', ylab='Ia purity', ylim=c(0,1), col='red')
# plot(z.seq[1:length(z.seq)-1], purity.Ia[1:length(purity.Ia)-1], lwd=3, col='blue')


N.Ia.true.test <- c()
N.Ia.tot.test <- c()
N.Ia.false.test <- c()

fom.Ia.test <- c()
eff.Ia.test <- c()
pseudo.purity.Ia.test <- c()
purity.Ia.test <- c()
for (i in seq(length(z.seq))){
  criteria.true <- which((testSet.with.predictions$type == trainingSet.with.predictions$type.pred) &
                           (z.seq[i] <= trainingSet.with.predictions$GENZ) &
                           (trainingSet.with.predictions$GENZ < z.seq[i+1]))
  criteria.false <- which(trainingSet.with.predictions$type != trainingSet.with.predictions$type.pred &
                            (z.seq[i] <= trainingSet.with.predictions$GENZ) &
                            (trainingSet.with.predictions$GENZ < z.seq[i+1]))
  criteria.tot <- which((z.seq[i] <= trainingSet.with.predictions$GENZ) &
                          (trainingSet.with.predictions$GENZ < z.seq[i+1]))
  
  N.Ia.tot.test <- c(N.Ia.tot.test, length(which(trainingSet$type[criteria.tot] == 'snIa')))
  N.Ia.true.test <- c(N.Ia.true.test, length(which(trainingSet.with.predictions$type.pred[criteria.true] == 'snIa')))
  N.Ia.false.test <- c(N.Ia.false.test, length(which(trainingSet.with.predictions$type.pred[criteria.false] == 'snIa')))
  
  fom.Ia.test <- c(fom.Ia.test, CFoM.Ia(tail(N.Ia.tot.test, n=1), tail(N.Ia.true.test, n=1), tail(N.Ia.false.test, n=1)))
  eff.Ia.test <- c(eff.Ia.test, tail(N.Ia.true.test, n=1)/tail(N.Ia.tot.test, n=1))
  pseudo.purity.Ia.test <- c(pseudo.purity.Ia.test, tail(N.Ia.true.test, n=1)/(tail(N.Ia.true.test, n=1)+tail(W.Ia.false, n=1)*tail(N.Ia.false.test, n=1)))
  purity.Ia.test <- c(purity.Ia, tail(N.Ia.true.test, n=1)/(tail(N.Ia.true.test, n=1)+tail(N.Ia.false.test, n=1)))
}
#lines(z.seq[1:length(z.seq)-1], fom.Ia.test[1:length(fom.Ia.test)-1], lty=2, lwd=3)
lines(z.seq[1:length(z.seq)-1], eff.Ia.test[1:length(eff.Ia.test)-1], lty=2, lwd=3, col='seagreen')
# lines(z.seq[1:length(z.seq)-1], pseudo.purity.Ia.test[1:length(pseudo.purity.Ia.test)-1], 
#       lwd=3, lty=2, col='red')
# lines(z.seq[1:length(z.seq)-1], purity.Ia[1:length(purity.Ia)-1], lwd=3, col='blue')

legend('bottomleft',c('training', 'test'), lty=c(1,2), lwd=c(3,3), col=c('seagreen', 'seagreen'))
# 
# CFoM.Ia(length(which(trainingSet$type[criteria.tot] == 'snIa')), 
#         length(which(trainingSet.with.predictions$type.pred[criteria.true] == 'snIa')), 
#         length(which(trainingSet.with.predictions$type.pred[criteria.false] == 'snIa')))
# 

## Dividing diffusion map coordinates by SN type
trainingSet.Ia <- trainingSet[which(trainingSet$type == 'snIa'),]
trainingSet.Ibc <- trainingSet[which(trainingSet$type == 'snIbc'),]
trainingSet.II <- trainingSet[which(trainingSet$type == 'snII'),]

testSet.Ia <- testSet[which(testSet$type == 'snIa'),]
testSet.Ibc <- testSet[which(testSet$type == 'snIbc'),]
testSet.II <- testSet[which(testSet$type == 'snII'),] 

dmapX.trainingSet.Ia <- dmap$X[trainingSet.Ia$idx+1,]
dmapX.trainingSet.Ibc <- dmap$X[trainingSet.Ibc$idx+1,]
dmapX.trainingSet.II <- dmap$X[trainingSet.II$idx+1,]

dmapX.testSet.Ia <- dmap$X[testSet.Ia$idx+1,]
dmapX.testSet.Ibc <- dmap$X[testSet.Ibc$idx+1,]
dmapX.testSet.II <- dmap$X[testSet.II$idx+1,]

## Prepare data frames for ggplot 
df.training.Ia <- data.frame(cbind(dmapX.trainingSet.Ia[,dmap.x.coord], dmapX.trainingSet.Ia[,dmap.y.coord]))
df.training.II <- data.frame(cbind(dmapX.trainingSet.II[,dmap.x.coord], dmapX.trainingSet.II[,dmap.y.coord]))
df.training.Ibc <- data.frame(cbind(dmapX.trainingSet.Ibc[,dmap.x.coord], dmapX.trainingSet.Ibc[,dmap.y.coord]))

type <- rep('SN Ia',times=dimension(dmapX.trainingSet.Ia)[1])
df.training.Ia <- cbind(df.training.Ia, type)

type <- rep('SN II',times=dimension(dmapX.trainingSet.II)[1])
df.training.II <- cbind(df.training.II, type)

type <- rep('SN Ibc',times=dimension(dmapX.trainingSet.Ibc)[1])
df.training.Ibc <- cbind(df.training.Ibc, type)

xlabel <- paste('Diffusion Coordinate', dmap.x.coord)
ylabel <- paste('Diffusion Coordinate', dmap.y.coord)
zlabel <- paste('Diffusion Coordinate', dmap.z.coord)

names(df.training.Ia) <-  c('x', 'y', 'type')
names(df.training.II) <-  c('x', 'y', 'type')
names(df.training.Ibc) <- c('x', 'y', 'type')

df.training <- rbind(df.training.Ia, df.training.II, df.training.Ibc)

## Test set data frames
df.test.Ia <- data.frame(cbind(dmapX.testSet.Ia[,dmap.x.coord], dmapX.testSet.Ia[,dmap.y.coord]))
df.test.II <- data.frame(cbind(dmapX.testSet.II[,dmap.x.coord], dmapX.testSet.II[,dmap.y.coord]))
df.test.Ibc <- data.frame(cbind(dmapX.testSet.Ibc[,dmap.x.coord], dmapX.testSet.Ibc[,dmap.y.coord]))

type <- rep('SN Ia',times=dimension(dmapX.testSet.Ia)[1])
df.test.Ia <- cbind(df.test.Ia, type)

type <- rep('SN II',times=dimension(dmapX.testSet.II)[1])
df.test.II <- cbind(df.test.II, type)

type <- rep('SN Ibc',times=dimension(dmapX.testSet.Ibc)[1])
df.test.Ibc <- cbind(df.test.Ibc, type)

names(df.test.Ia) <-  c('x', 'y', 'type')
names(df.test.II) <-  c('x', 'y', 'type')
names(df.test.Ibc) <- c('x', 'y', 'type')

df.test <- rbind(df.test.Ia, df.test.II, df.test.Ibc)

## ggplot
plot.xlim <- c(-2, 1) #c(-0.5, 1.)
plot.ylim <- c(-1,1) #c(-1, 1)
marginalTop.theme <- theme(legend.position = "none",          
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           plot.margin = unit(c(4,-5.5,4,7), "mm"))

marginalRight.theme <- theme(legend.position = "none",          
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             plot.margin = unit(c(3,5,27,3), "mm"))

scatter.theme <- theme(legend.position = "bottom",
                       plot.margin = unit(c(3,-5.5,4,3), "mm"))

empty.theme <- theme(panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.title.x = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks = element_blank(), 
                     axis.text.x = element_blank(), 
                     axis.text.y = element_blank())



scalex.align <- scale_x_continuous(breaks = seq(from=min(plot.xlim), to=max(plot.xlim), by=0.5),
                                   limits = plot.xlim)#,
#                                    expand = c(.05,.05))

scaley.align <- scale_x_continuous(breaks = seq(from=min(plot.ylim), to=max(plot.ylim), by=0.5),
                                   limits = plot.ylim)#,
#                                    expand = c(.05,.05))

empty <- ggplot() + geom_point(aes(1,1), colour="white") + empty.theme


sp.training.densTop <- ggplot(df.training, fill=factor(type)) + 
    #    geom_histogram(aes(x=x, colour=factor(type)),position='identity', alpha=0.3) + 
    geom_line(stat='density', bw='SJ', adjust=0.5, aes(x=x, y=..scaled.., colour=factor(type), linetype=factor(type))) + 
    marginalTop.theme  + scalex.align + scale_colour_manual(values=c('darkorange', 'seagreen3', 'blue2'))

sp.training <- ggplot(df.training, aes(x=x, y=y))
sp.training.point <- sp.training + geom_point(alpha=0.5, aes(colour=factor(type), shape=factor(type))) + scatter.theme + #scalex.align + scaley.align +
    xlab(xlabel) + ylab(ylabel) + labs(colour='type', shape='type') + xlim(plot.xlim) + ylim(plot.ylim) + scale_colour_manual(values=c('darkorange', 'seagreen3', 'blue2'))

sp.training.densRight <- ggplot(df.training) + 
    geom_line(stat='density', bw='SJ', adjust=0.5, aes(x=y, y=..scaled.., colour=factor(type), linetype=factor(type))) +
    coord_flip() + marginalRight.theme + scaley.align + scale_colour_manual(values=c('darkorange', 'seagreen3', 'blue2'))

save.plot <- FALSE
if (save.plot){
    ppi <- 300
    png('results/Plots/diffusion_space_proj_coo-3-10.png', width=10*ppi, height=10*ppi, res=ppi)
}
grid.arrange(sp.training.densTop, empty, sp.training.point, sp.training.densRight, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

if(save.plot){
    dev.off()
}


plot.ylim <- c(-1, 1)
scaley.align <- scale_x_continuous(breaks = seq(from=min(plot.ylim), to=max(plot.ylim), by=0.5),
                                   limits = plot.ylim)
sp.test.densTop <- ggplot(df.test) + 
    geom_line(stat='density', bw='SJ', adjust=0.5, aes(x=x, y=..scaled.., colour=factor(type), linetype=factor(type))) + 
    marginalTop.theme  + scalex.align + scale_colour_manual(values=c('darkorange', 'seagreen3', 'blue2'))

sp.test <- ggplot(df.test, aes(x=x, y=y))
sp.test.point <- sp.test + geom_point(alpha=0.5, aes(colour=factor(type), shape=factor(type))) + scatter.theme +
    xlab(xlabel) + ylab(ylabel) + labs(colour='type', shape='type') + xlim(plot.xlim) + ylim(plot.ylim) + 
    scale_colour_manual(values=c('darkorange', 'seagreen3', 'blue2'))

sp.test.densRight <- ggplot(df.test) + 
    geom_line(stat='density', bw='SJ', adjust=0.8, aes(x=y, y=..scaled.., colour=factor(type), linetype=factor(type))) +
    coord_flip() + marginalRight.theme + scaley.align + scale_colour_manual(values=c('darkorange', 'seagreen3', 'blue2'))

save.plot <- FALSE
if (save.plot){
    ppi <- 300
    png('results/Plots/diffusion_space_proj_coo-3-10_TESTSET.png', width=10*ppi, height=10*ppi, res=ppi)
}
grid.arrange(sp.test.densTop, empty, sp.test.point, sp.test.densRight, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

if(save.plot){
    dev.off()
}

## par(par.defaults)
print(sp.training.point)
print(sp.test.point)
# message(paste(plotDir, fileTraining, sep=''))
# message(paste(plotDir, fileTest, sep=''))





# Variables for rgl plot
if (TRUE){
    library(rgl)

    plot3d(dmapX.trainingSet.Ia[,dmap.x.coord], dmapX.trainingSet.Ia[,dmap.y.coord], 
           dmapX.trainingSet.Ia[,dmap.z.coord], col=colors[1], type='p',
           xlab='', ylab='', zlab='',
           xlim=c(-1, 1), ylim=c(-1, 0.5), zlim=c(-1,3))
    
    plot3d(dmapX.trainingSet.II[,dmap.x.coord], dmapX.trainingSet.II[,dmap.y.coord], 
           dmapX.trainingSet.II[,dmap.z.coord], col=colors[2], type='p', 
           xlim=c(-1,1), ylim=c(-1, 0.5), zlim=c(-1,3),
           add=TRUE)
    
    plot3d(dmapX.trainingSet.Ibc[,dmap.x.coord], dmapX.trainingSet.Ibc[,dmap.y.coord], 
           dmapX.trainingSet.Ibc[,dmap.z.coord], col=colors[3], type='p', 
           xlim=c(-1,1), ylim=c(-1, 0.5), zlim=c(-1,3),
           add=TRUE)

    rgl.bbox(color="grey50",          # grey60 surface and black text
             emission="grey50",       # emission color is grey50
             xlen=0, ylen=0, zlen=0)  # Don't add tick marks
    axes3d(edges=c("x--", "y+-", "z--"), ntick=6, cex=1)
    rgl::mtext3d(xlabel, edge="x--", line=2)
    rgl::mtext3d(ylabel, edge="y+-", line=3)
    rgl::mtext3d(zlabel, edge="z--", line=3)
    
    
    # 3D plot TEST
    zlim.3d <- c(-2,3.5)
    
    plot3d(dmapX.testSet.Ia[,dmap.x.coord], dmapX.testSet.Ia[,dmap.y.coord], 
           dmapX.testSet.Ia[,dmap.z.coord], col=colors[1], type='p', size=2,
           xlab='', ylab='', zlab='',
           xlim=c(-1,1), ylim=c(-4, 1), zlim=zlim.3d)
#            xlim=c(-0.6, 1.1), ylim=c(-4, 1), zlim=c(-1,1),
    
    plot3d(dmapX.testSet.II[,dmap.x.coord], dmapX.testSet.II[,dmap.y.coord], 
           dmapX.testSet.II[,dmap.z.coord], col=colors[2], type='p', size=2,
           xlim=c(-1,1), ylim=c(-4, 1), zlim=zlim.3d,
           add=TRUE)
    
    plot3d(dmapX.testSet.Ibc[,dmap.x.coord], dmapX.testSet.Ibc[,dmap.y.coord], 
           dmapX.testSet.Ibc[,dmap.z.coord], col=colors[3], type='p', size=2,
           xlim=c(-1,1), ylim=c(-4, 1), zlim=zlim.3d,
           add=TRUE)
#            xlim=c(-0.6, 1.1), ylim=c(-4, 1), zlim=c(-1,1))
    
    rgl.bbox(color="grey50",          # grey60 surface and black text
             emission="grey50",       # emission color is grey50
             xlen=0, ylen=0, zlen=0)  # Don't add tick marks
    axes3d(edges=c("x--", "y+-", "z--"), ntick=6, cex=1)
    rgl::mtext3d(xlabel, edge="x--", line=2)
    rgl::mtext3d(ylabel, edge="y+-", line=3)
    rgl::mtext3d(zlabel, edge="z--", line=3)
    
}