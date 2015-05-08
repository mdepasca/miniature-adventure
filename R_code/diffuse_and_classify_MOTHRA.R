## File managing diffusion map calculation and successive classification with
## Random Forest (Breiman 2001).

## NOTE: to open distance matrix file of 20,000x20,000 needs at least 48 GB
## of RAM.
library(diffusionMap)
library(randomForest)
library(bootstrap)
library(snowfall)
library(bigmemory)
library(verification)
## Initialize nodes
## message("Initializing nodes...")
## sfInit(parallel=TRUE, cpus=8, type='SOCK', socketHosts=rep('localhost', 8))

diffuse.and.classify <- function(dmap.eps.val, distance.matrix, training.set, test.set){

    print(dmap.eps.val)
    dmap <- diffuse(distance.matrix, eps.val=dmap.eps.val, neigen=120)

    ## extracting training and test set diffusion coordinate
    dmapCoord.train <- data.frame(dmap$X[training.set$idx,])
    dmapCoord.test <- data.frame(dmap$X[test.set$idx,])
    ## adding `type' column
    dmapCoord.train <- cbind(dmapCoord.train, training.set$type)
    dmapCoord.test <- cbind(dmapCoord.test, test.set$type)

    ncols <- dim(dmapCoord.train)[2]
    names(dmapCoord.train)[ncols] <- 'type'
    names(dmapCoord.test)[ncols] <- 'type'

    types <- c('snIa', 'snIbc', 'snII')
    ## training the random forest classifier
    ## sn.rf <- randomForest(type ~ ., data=dmapCoord.train, xtest=dmapCoord.test[,1:ncols], ytest=dmapCoord.test[,ncols+1], importance=TRUE)
    sn.rf <- randomForest(x=dmapCoord.train[,1:(ncols-1)],
                    y=as.factor(dmapCoord.train[,ncols]),
                    xtest=dmapCoord.test[,1:(ncols-1)],
                    ytest=as.factor(dmapCoord.test[,ncols]), importance=TRUE)

    penalty <- 3
    N.total.Ia <- sum(sn.rf$confusion['snIa', 1:3])
    N.true.Ia <- sn.rf$confusion['snIa', 'snIa']
    N.false.Ia <- sum(sn.rf$confusion[2:3, 'snIa'])
    fom <- (1/N.total.Ia)*(N.true.Ia**2/(N.true.Ia+penalty*N.false.Ia))
    result <- list(dmap=dmap, sn.rf=sn.rf, fom=fom)
    return(result)
}

## Load the data
if (!exists('dump.table')){
    message("Reading simgen dump file...")
    dump.table <- read.table('train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP', header=TRUE, skip=1)
}

if (!exists('training.set')){
    message("Building training set...")
    training.set.snIa <- read.table('products/SIMGEN_PUBLIC_DES_FIT.Ia.TRAIN', col.names=c('idx', 'snid', 'path', 'type'))
    training.set.snIbc <- read.table('products/SIMGEN_PUBLIC_DES_FIT.Ibc.TRAIN', col.names=c('idx', 'snid', 'path', 'type'))
    training.set.snII <- read.table('products/SIMGEN_PUBLIC_DES_FIT.II.TRAIN', col.names=c('idx', 'snid', 'path', 'type'))

    training.set <- rbind(training.set.snIa, training.set.snIbc, training.set.snII)
    ## training.set <- cbind(training.set, rbind(training.set.snIa[3], training.set.snIbc[3], training.set.snII[3]))
    ## adding 1 to idx values: R indexes from 1, the read indeces comes from Python, which indexes from 0
    training.set$idx <- training.set$idx + 1
}

if (!exists('test.set')){
    message("Building test set...")
    test.set <- read.table('products/SIMGEN_PUBLIC_DES_FIT.TEST', col.names=c('idx', 'snid', 'path', 'type'))
    match.idx <- match(test.set$snid,dump.table$CID)
    test.set$type <- dump.table$GENTYPE[match.idx]
    test.set$type[test.set$type == 1] <- 'snIa'
    test.set$type[test.set$type == 21] <- 'snII'
    test.set$type[test.set$type == 22] <- 'snII'
    test.set$type[test.set$type == 23] <- 'snII'
    test.set$type[test.set$type == 3] <- 'snIbc'
    test.set$type[test.set$type == 32] <- 'snIbc'
    test.set$type[test.set$type == 33] <- 'snIbc'
    test.set$type<-as.factor(test.set$type)
    test.set$idx <- test.set$idx + 1
    rm(match.idx)
}
## Reading distance matrix
if (!exists('distMat')){
    message("Reading distance matrix...")
    filePath <- 'products/distance_matrix/dist_matrix_Sum.txt'
    distMat <- matrix(scan(filePath, comment.char='#'), ncol=21304, nrow=21304, byrow=TRUE)
} else{
    if (!is.big.matrix(distMat) & !exists('big.distMat')){
        big.distMat <- as.big.matrix(distMat)
    }
}

## reducing training set to exclude LCs not computed
if (!exists('training.set.red')){
    training.set.red <- training.set[training.set$idx < dim(big.distMat)[1], ]
}
if (!exists('big.training.set.red')){
    big.training.set.red <- as.big.matrix(training.set.red)
}


eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

## message("Loading packages  on nodes...")
## sfLibrary(diffusionMap)
## sfLibrary(randomForest)
## sfLibrary(bigmemory)

## message("Exporting data to nodes...")
## sfExport("eps.grid")
## sfExport("big.distMat")
## sfExport("training.set.red")
## sfExport("test.set")

## sfClusterSetupRNG()
## message("Calculating!")
## result <- sfLapply(eps.grid, diffuse.and.classify, big.distMat, training.set.red, test.set)#, (distMat, trainSet))
## print(eps.grid[3])


dmap <- diffuse(distMat, eps.val=eps.grid[2], neigen=120)

## extracting training and test set diffusion coordinate
dmapCoord.train <- data.frame(dmap$X[training.set$idx,])
dmapCoord.test <- data.frame(dmap$X[test.set$idx,])
## adding `type' column
dmapCoord.train <- cbind(training.set$type, dmapCoord.train)
dmapCoord.test <- cbind(test.set$type, dmapCoord.test)

names(dmapCoord.train)[1] <- 'type'
names(dmapCoord.test)[1] <- 'type'

types <- c('snIa', 'snIbc', 'snII')

## Partition the data

## Setting number of folds
## k <- 10

## ## Size of each fold
## n <- floor(nrow(training.set.red)/k)
## ## Vector which will store the errors
## err.vect <- rep(NA, k)

## ## CV of random forest

## ##looping over folds
## for (i in 1:k){
##     ## Partitioning of i-fold
##     s1 <- ((i - 1) * n + 1) # start index of subset
##     s2 <- (i * n)           # end index of subset
##     subset <- s1:s2         # range of subset

##     cv.train <- dmapCoord.train[-subset,]
##     cv.test <- dmapCoord.train[subset,]

##     ## training the random forest classifier
##     fit <- randomForest(x=cv.train[,-1], y=as.factor(cv.train[,1]), importance=TRUE)
##     prediction <- predict(fit, newdata=cv.test[,-1], type='prob')[,2]

##     ## err.vect[i] <- roc.area(cv.test[,1], prediction)$A
##     ## print(paste("AUC for fold", i, ":", err.vect[i]))
## }
## print(paste("Average AUC:", mean(err.vect)))
##     penalty <- 3
##     N.total.Ia <- sum(sn.rf$confusion['snIa', 1:3])
##     N.true.Ia <- sn.rf$confusion['snIa', 'snIa']
##     N.false.Ia <- sum(sn.rf$confusion[2:3, 'snIa'])
##     fom <- (1/N.total.Ia)*(N.true.Ia**2/(N.true.Ia+penalty*N.false.Ia))


result <- diffuse.and.classify(eps.grid[3], distMat, training.set.red, test.set)

## stop snowfall
## sfStop()
