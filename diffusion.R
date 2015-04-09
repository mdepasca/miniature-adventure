# File managing diffusion map calculation
library(diffusionMap)
library(R.utils)

filePath <- 'results/RATQUAD-with_prior/distance-matrix/dist_matrix_Sum.txt'
fileLines <- countLines(filePath)
distMat <- matrix(scan(filePath, comment.char='#'),
    ncol=fileLines, nrow=fileLines, byrow=TRUE)

dmap <- diffuse(distMat, eps.val=eps.grid[2], neigen=120)
