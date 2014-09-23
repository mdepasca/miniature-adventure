# File managing diffusion map calculation
library(diffusionMap)
# distMatFrame <- read.table('Pymatrix.txt')
filePath <- 'products/distance_matrix/Pymatrix.txt'
distMat <- matrix(scan(filePath, comment.char='#'), 
    ncol=18080, nrow=18080, byrow=TRUE)
# distMat <- data.matrix(distMatFrame)
# rm(distMatFrame)

# matrix(scan("matrix.dat", n = 200*2000), 200, 2000, byrow = TRUE) 
