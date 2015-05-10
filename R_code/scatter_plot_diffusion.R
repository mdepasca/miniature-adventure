library(diffusionMap)
library(gclus)
library(R.utils)
library(randomForest)
library(ggplot2)

source('R_code/get_training.R')
source('R_code/get_test.R')

par.defaults <- par(no.readonly=TRUE)
homeDir <- '~/sandbox/miniature-adventure/'
plotDir <- paste(homeDir, 'results/Plots/', sep='')

colors <- c('darkorange', 'seagreen3', 'blue2')
dmap.x.coord <- 10#4
dmap.y.coord <- 3#3

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

names(df.training.Ia) <- c('x', 'y', 'type')
names(df.training.II) <- c('x', 'y', 'type')
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

names(df.test.Ia) <- c('x', 'y', 'type')
names(df.test.II) <- c('x', 'y', 'type')
names(df.test.Ibc) <- c('x', 'y', 'type')

df.test <- rbind(df.test.Ia, df.test.II, df.test.Ibc)

## ggplot

sp.training <- ggplot(df.training, aes(x=x, y=y))
sp.training.point <- sp.training + geom_point(alpha=0.1, aes(colour=factor(type), shape=factor(type)))

sp.test <- ggplot(df.test, aes(x=x, y=y))
sp.test.point <- sp.test + geom_point(alpha=0.1, aes(colour=factor(type), shape=factor(type)))

X11()
dev.new()
sp.training.point

dev.new()
sp.test.point


par(omi=c(0,0,0,0))
## dev.new()
fileTraining <- paste('trainingSet', dmap.x.coord, dmap.y.coord, '.png', sep='_')
## png(paste(plotDir, fileTraining, sep=''), width=7, height=7, units='in', res=300)
dev.new()
plot(dmapX.trainingSet.Ia[,dmap.x.coord], dmapX.trainingSet.Ia[,dmap.y.coord], pch=20, col=colors[1], cex=.5, xlab=paste('Diffusion Coordinate', dmap.x.coord), ylab=paste('Diffusion Coordinate', dmap.y.coord), main='SNPhotCC Training Set')
points(dmapX.trainingSet.II[,dmap.x.coord], dmapX.trainingSet.II[,dmap.y.coord], pch=15, col=colors[2], cex=.5, xlab=paste('Diffusion Coordinate', dmap.x.coord), ylab=paste('Diffusion Coordinate', dmap.y.coord))
points(dmapX.trainingSet.Ibc[,dmap.x.coord], dmapX.trainingSet.Ibc[,dmap.y.coord], pch=17, col=colors[3], cex=.5, xlab=paste('Diffusion Coordinate', dmap.x.coord), ylab=paste('Diffusion Coordinate', dmap.y.coord))

legend('topleft', c('Ia', 'II', 'Ibc'), col=colors, pch=c(20,15,17))
## dev.off()

dev.new()
fileTest <- paste('testSet', dmap.x.coord, dmap.y.coord, '.png', sep='_')
## png(paste(plotDir, fileTest, sep=''), width=7, height=7, units='in', res=300)

plot(dmapX.testSet.Ia[,dmap.x.coord], dmapX.testSet.Ia[,dmap.y.coord], pch=20, col=colors[1], cex=.3, xlab=paste('Diffusion Coordinate', dmap.x.coord), ylab=paste('Diffusion Coordinate', dmap.y.coord), main='SNPhotCC data set')
points(dmapX.testSet.II[,dmap.x.coord], dmapX.testSet.II[,dmap.y.coord], pch=15, col=colors[2], cex=.3, xlab=paste('Diffusion Coordinate', dmap.x.coord), ylab=paste('Diffusion Coordinate', dmap.y.coord))
points(dmapX.testSet.Ibc[,dmap.x.coord], dmapX.testSet.Ibc[,dmap.y.coord], pch=17, col=colors[3], cex=.3, xlab=paste('Diffusion Coordinate', dmap.x.coord), ylab=paste('Diffusion Coordinate', dmap.y.coord))
##dev.off()
par(par.defaults)

message(paste(plotDir, fileTraining, sep=''))
message(paste(plotDir, fileTest, sep=''))
