source('diffusion.R')
source('get_training.R')
source('get_test.R')
source('classify_RF.R')

path <- 'results/SIMGEN_PUBLIC_FIT/'
specific.path <- 'RATQUAD-with_prior'

path <- paste(path, specific.path, '/distance_matrix/', sep='')

neigen <- 120
eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

dmap <- calc_diffusion_map(path, eps.val=eps.grid[2], neigen=niegen)

trainingSet <- get_training(path=paste(path,specific.path,'/', sep=''), 
	fileNameRoot=specific.path)

testSet <- get_test(path=paste(path,specific.path,'/', sep=''), 
	fileNameRoot=specific.path)

sn.rf <- classify_RF(dmap, trainingSet, testSet)