source('calc_diffusion_map.R')
source('get_training.R')
source('get_test.R')
source('classify_RF.R')

path <- 'results/SIMGEN_PUBLIC_FIT/'
## specific.path <- 'RATQUAD-with_prior'
specific.path <- 'RBF-with_prior'
## specific.path <- 'RBF_test-length'

path.distance <- paste(path, specific.path, '/distance_matrix/', sep='')

neigen <- 50
eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

message(paste('Data from directory:', paste(path, specific.path, sep='')))

## eps.val should be near the value for which the Euclidean distance is considered too big
## for the 2 light curves to be similar. In this case is 5
eps.val <- 10
if (exists('dmap')){
    ## to avoid any possiblity of duplicating the size of dmap
    rm(dmap)
}
dmap <- calc_diffusion_map(path.distance, eps.val=eps.val, neigen=neigen)

if (exists('trainingSet')){
    ## if the variable exists, assignment as below would only
    ## append a new version of it
    rm(trainingSet)
}
message(paste('Building training set from', paste(path, specific.path, '/', sep='')))
trainingSet <- get_training(path=paste(path, specific.path, '/', sep=''), fileNameRoot=specific.path)

if (exists('testSet')){ rm(testSet) }
message(paste('Building test set from', paste(path, specific.path, '/', sep='')))
testSet <- get_test(path=paste(path, specific.path, '/', sep=''), fileNameRoot=specific.path)

if (exists('sn.rf')){ rm(sn.rf) }
sn.rf <- classify_RF(dmap, trainingSet, testSet)
