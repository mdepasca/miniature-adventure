source('calc_diffusion_map.R')
source('get_training.R')
source('get_test.R')
source('classify_RF.R')

path <- 'results/SIMGEN_PUBLIC_FIT/'
## specific.path <- 'RATQUAD-with_prior'
## specific.path <- 'RBF-with_prior'
specific.path <- 'RBF_test-length'

path.distance <- paste(path, specific.path, '/distance_matrix/', sep='')

neigen <- 120
eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

dmap <- calc_diffusion_map(path.distance, eps.val=eps.grid[2], neigen=neigen)

message(paste('Building training set from', paste(path, specific.path, '/', sep='')))
trainingSet <- get_training(path=paste(path, specific.path, '/', sep=''), fileNameRoot=specific.path)

message(paste('Building test set from', paste(path, specific.path, '/', sep='')))
testSet <- get_test(path=paste(path, specific.path, '/', sep=''), fileNameRoot=specific.path)

sn.rf <- classify_RF(dmap, trainingSet, testSet)
