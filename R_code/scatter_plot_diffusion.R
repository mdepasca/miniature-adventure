library(diffusionMap)
library(gclus)
library(R.utils)
library(randomForest)
library(ggplot2)

source('R_code/get_training.R')
source('R_code/get_test.R')
source('R_code/util.R')
source('R_code/figure_of_merits.R')

par.defaults <- par(no.readonly=TRUE)
homeDir <- '~/sandbox/miniature-adventure/'
plotDir <- paste(homeDir, 'results/Plots/', sep='')

colors <- c('darkorange', 'seagreen3', 'blue2')
dmap.x.coord <- 10
dmap.y.coord <- 3
dmap.z.coord <- 15

path <- 'results/SIMGEN_PUBLIC_FIT/RBF_test-length/'
fileNameRoot <- 'RBF_test-length'
if (!exists('trainingSet')){
    trainingSet <- get_training(path, fileNameRoot)
}
if (!exists('testSet')){
    testSet <- get_test(path, fileNameRoot)
}
if (!exists('dmap')){
    dmap <- loadObject(paste(path, 'distance_matrix/',  'diffusion_map_t_2.RData', sep=''))
}
if (!exists('sn.rf')){
    sn.rf <- loadObject(paste(homeDir, 'results/SIMGEN_PUBLIC_FIT/', 'RBF_test-lengthrandomForestModel.RData', sep=''))
}
dump.file <- get_dump('train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP')

trainingSet.with.predictions <- cbind(trainingSet, sn.rf$predicted, dump.file$GENZ[match(trainingSet$snid, dump.file$CID)])
names(trainingSet.with.predictions) <- c(names(trainingSet.with.predictions)[1:4], 'type.pred', 'GENZ')

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

dmapX.trainingSet.Ia <- dmap$X[trainingSet.Ia$idx,]
dmapX.trainingSet.Ibc <- dmap$X[trainingSet.Ibc$idx,]
dmapX.trainingSet.II <- dmap$X[trainingSet.II$idx,]

dmapX.testSet.Ia <- dmap$X[testSet.Ia$idx,]
dmapX.testSet.Ibc <- dmap$X[testSet.Ibc$idx,]
dmapX.testSet.II <- dmap$X[testSet.II$idx,]

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

sp.training <- ggplot(df.training, aes(x=x, y=y))
sp.training.point <- sp.training + geom_point(alpha=0.5, aes(colour=factor(type), shape=factor(type))) +
    xlab(xlabel) + ylab(ylabel) + labs(colour='type', shape='type')

sp.test <- ggplot(df.test, aes(x=x, y=y))
sp.test.point <- sp.test + geom_point(alpha=0.5, aes(colour=factor(type), shape=factor(type))) +
    xlab(xlabel) + ylab(ylabel) + labs(colour='type', shape='type') + xlim(-0.6, 1.1) + ylim(-4, 1)


## par(par.defaults)
print(sp.training.point)
print(sp.test.point)
# message(paste(plotDir, fileTraining, sep=''))
# message(paste(plotDir, fileTest, sep=''))





# Variables for rgl plot
plot3d(dmapX.trainingSet.Ia[,dmap.x.coord], dmapX.trainingSet.Ia[,dmap.y.coord], 
       dmapX.trainingSet.Ia[,dmap.z.coord], col=colors[1], type='p',
       xlab=xlabel, ylab=ylabel, zlab=zlabel)

plot3d(dmapX.trainingSet.II[,dmap.x.coord], dmapX.trainingSet.II[,dmap.y.coord], 
       dmapX.trainingSet.II[,dmap.z.coord], col=colors[2], type='p',
       #xlab=xlabel, ylab=ylabel, zlab=zlabel,
       add=TRUE)

plot3d(dmapX.trainingSet.Ibc[,dmap.x.coord], dmapX.trainingSet.Ibc[,dmap.y.coord], 
       dmapX.trainingSet.Ibc[,dmap.z.coord], col=colors[3], type='p',
       #xlab=xlabel, ylab=ylabel, zlab=zlabel,
       add=TRUE)


# 3D plot TEST

plot3d(dmapX.testSet.Ia[,dmap.x.coord], dmapX.testSet.Ia[,dmap.y.coord], 
       dmapX.testSet.Ia[,dmap.z.coord], col=colors[1], type='p', size=1.4,
       xlab=xlabel, ylab=ylabel, zlab=zlabel,
       xlim=c(-0.6, 1.1), ylim=c(-4, 1), zlim=c(-1,1),
       forceClipregion=TRUE)

plot3d(dmapX.testSet.II[,dmap.x.coord], dmapX.testSet.II[,dmap.y.coord], 
       dmapX.testSet.II[,dmap.z.coord], col=colors[2], type='p', size=1.4,
       #xlab=xlabel, ylab=ylabel, zlab=zlabel,
       add=TRUE,
       xlim=c(-0.6, 1.1), ylim=c(-4, 1), zlim=c(-1,1))

plot3d(dmapX.testSet.Ibc[,dmap.x.coord], dmapX.testSet.Ibc[,dmap.y.coord], 
       dmapX.testSet.Ibc[,dmap.z.coord], col=colors[3], type='p', size=1.4,
       #xlab=xlabel, ylab=ylabel, zlab=zlabel,
       add=TRUE,
       xlim=c(-0.6, 1.1), ylim=c(-4, 1), zlim=c(-1,1))
