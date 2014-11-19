## File managing diffusion map calculation and successive classification with
## Random Forest (Breiman 2001).

## NOTE: to open distance matrix file of 20,000x20,000 needs at least 48 GB
## of RAM.
library(diffusionMap)
library(randomForest)
library(bootstrap)
library(snowfall)

## Initialize nodes
## message("Initializing nodes...")
## sfInit(parallel=TRUE, cpus=8, type='SOCK', socketHosts=rep('localhost', 8))

diffuse.and.classify <- function(dmap.esp.val, distance.matrix, training.set){

    print(dmap.esp.val)
    dmap <- diffuse(distance.matrix, eps.val=dmap.eps.val, neigen=120)
    
    ## extracting training set diffusion coordinate
    dmapCoord.train <- data.frame(dmap$X[training.set$idx,])

    ## adding column of types
    dmapCoord.train <- cbind(dmapCoord.train, training.set$type)
    ncols <- dim(dmapCoord.train)[2]
    names(dmapCoord.train)[ncols] <- 'type'

    ## training the random forest classifier
    sn.rf <- randomForest(type ~ ., data=dmapCoord.train, importance=TRUE)

    penalty <- 3
    N.total.Ia <- sum(sn.rf$confusion['snIa', 1:3])
    N.true.Ia <- sn.rf$confusion['snIa', 'snIa']
    N.false.Ia <- sum(sn.rf$confusion[2:3, 'snIa'])
    fom <- (1/N.total.Ia)*(N.true.Ia**2/(N.true.Ia+penalty*N.false.Ia))
    result <- list(dmap=dmap, sn.rf=sn.rf, fom=fom)
    return(result)
}

## Load the data

## Reading indeces of training set SNe
if (!exists('trainSet')){
    message("Building training set...")
    training.set.snIa <- read.table('products/SIMGEN_PUBLIC_DES_FIT.Ia.TRAIN', col.names=c('idx', 'path', 'type'))
    training.set.snIbc <- read.table('products/SIMGEN_PUBLIC_DES_FIT.Ibc.TRAIN', col.names=c('idx', 'path', 'type'))
    training.set.snII <- read.table('products/SIMGEN_PUBLIC_DES_FIT.II.TRAIN', col.names=c('idx', 'path', 'type'))
 
    training.set <- rbind(training.set.snIa[1], training.set.snIbc[1], training.set.snII[1])
    training.set <- cbind(training.set, rbind(training.set.snIa[3], training.set.snIbc[3], training.set.snII[3]))
    ## adding 1 to idx values: R indexes from 1, the read indeces comes from Python, which indexes from 0
    training.set$idx <- training.set$idx + 1
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
if (!exists('trainSetRed')){
    trainSetRed <- trainSet[trainSet$idx < dim(big.distMat)[1], ]
}
## if (!exists('big.trainSetRed')){
##     big.trainSetRed <- as.big.matrix(trainSetRed)
## }


eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

## message("Loading packages  on nodes...")
## sfLibrary(diffusionMap)
## sfLibrary(randomForest)
## sfLibrary(bigmemory)

## message("Exporting data to nodes...")
## sfExport("eps.grid")
## sfExport("big.distMat")
## sfExport("trainSetRed")


## sfClusterSetupRNG()
## message("Calculating!")
## result <- sfLapply(eps.grid, diffuse.and.classify, big.distMat, trainSetRed)#, (distMat, trainSet))
print(eps.grid[3])

result <- diffuse.and.classify(eps.grid[3], big.distMat, trainSetRed)

# stop snowfall
## sfStop()
