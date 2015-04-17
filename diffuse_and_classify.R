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
dmap <- calc_diffusion_map(path.distance, eps.val=5, neigen=neigen)

message(paste('Building training set from', paste(path, specific.path, '/', sep='')))
trainingSet <- get_training(path=paste(path, specific.path, '/', sep=''), fileNameRoot=specific.path)

message(paste('Building test set from', paste(path, specific.path, '/', sep='')))
testSet <- get_test(path=paste(path, specific.path, '/', sep=''), fileNameRoot=specific.path)

sn.rf <- classify_RF(dmap, trainingSet, testSet)
