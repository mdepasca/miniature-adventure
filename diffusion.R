# File managing diffusion map calculation
library(diffusionMap)
library(R.utils)

filePath <- 'results/SIMGEN_PUBLIC_FIT/RATQUAD-with_prior/distance_matrix/dist_matrix_Sum.txt'
if (!exists('fileLines')){
    fileLines <- countLines(filePath)
}

eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

if (!exists('distMat')){
    distMat <- matrix(scan(filePath, comment.char='#'),
        ncol=fileLines, nrow=fileLines, byrow=TRUE)
}

dmap <- diffuse(distMat, eps.val=eps.grid[2], neigen=120)
