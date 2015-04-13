# File managing diffusion map calculation
library(diffusionMap)
library(R.utils)

# filePath <- 'results/SIMGEN_PUBLIC_FIT/RATQUAD-with_prior/distance_matrix/dist_matrix_Sum.txt'
# filePath <- 'results/SIMGEN_PUBLIC_FIT/RBF-with_prior/distance_matrix/dist_matrix_Sum.txt'
filePath <- 'results/SIMGEN_PUBLIC_FIT/RBF_test-length/distance_matrix/dist_matrix_Sum.txt'
if (!exists('fileLines')){
    fileLines <- countLines(filePath)
    message(paste('Number of lines in input file = ', fileLines))
}

eps.grid <- as.vector(seq(from=2, to=5, by=0.2))

if (!exists('distMat')){
    message(paste('Reading distance matrix from ', filePath, ' ...'))
    distMat <- matrix(scan(filePath, comment.char='#'),
        ncol=fileLines, nrow=fileLines, byrow=TRUE)
}

message('Creating diffusion map...')
dmap <- diffuse(distMat, eps.val=eps.grid[2], neigen=120)

dmap.corr <- abs(cor(dmap$X[,1:5]))
myorder <- order.single(dmap.corr)
mycolors <- dmat.color(dmap.corr)
cpairs(dmap$X[,1:5], myorder, panel.colors=mycolors, gap=.5)
# pairs(~dmap$X[,1]+dmap$X[,2]+dmap$X[,3]+dmap$X[,4]+dmap$X[,5]+dmap$X[,6])
