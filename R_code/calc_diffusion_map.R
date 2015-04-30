# File managing diffusion map calculation
library(diffusionMap)
library(R.utils)

# filePath <- 'results/SIMGEN_PUBLIC_FIT/RATQUAD-with_prior/distance_matrix/dist_matrix_Sum.txt'
# filePath <- 'results/SIMGEN_PUBLIC_FIT/RBF-with_prior/distance_matrix/dist_matrix_Sum.txt'
# filePath <- 'results/SIMGEN_PUBLIC_FIT/RBF_test-length/distance_matrix/dist_matrix_Sum.txt'
# neigen <- 120
calc_diffusion_map <- function(filePath, eps.val, neigen, old.dmap=FALSE){
    message(paste('Reading diffusion map from directory:', filePath))
    print(eps.val)
    fileLines <- countLines(paste(filePath,'dist_matrix_Sum.txt', sep=''))
    message(paste('Number of lines in input file = ', fileLines))

    if (!file.exists(paste(filePath,'distance_matrix.RData',sep=''))){
        message(paste('Reading distance matrix from ', filePath,
                      'dist_matrix_Sum.txt ...', sep=''))
        distMat <- matrix(scan(paste(filePath, 'dist_matrix_Sum.txt', sep=''),
                               comment.char='#'),
                          ncol=fileLines, nrow=fileLines, byrow=TRUE)
        saveObject(distMat, paste(filePath, 'distance_matrix.RData', sep=''),
                   compress=TRUE, safe=TRUE)
    }else{
         message(paste('Loading distance matrix from file: ', filePath,
                        'distance_matrix.RData', sep=''))
         distMat <- loadObject(paste(filePath, 'distance_matrix.RData', sep=''))
    }

    if (old.dmap){
        message('Loading diffusion map from file...')
        dmap <- loadObject(paste(filePath, 'diffusion_map.RData', sep=''))
    }else{
         message('Creating diffusion map...')
         dmap <- diffuse(distMat, eps.val=eps.val, neigen=neigen)
         rm(distMat, fileLines)
         saveObject(dmap, paste(filePath, 'diffusion_map.RData', sep=''),
                    compress=TRUE, safe=TRUE)
     }
    return(dmap)
}

# dmap.corr <- abs(cor(dmap$X[,1:5]))
# myorder <- order.single(dmap.corr)
# mycolors <- dmat.color(dmap.corr)
# cpairs(dmap$X[,1:5], myorder, panel.colors=mycolors, gap=.5)
# pairs(~dmap$X[,1]+dmap$X[,2]+dmap$X[,3]+dmap$X[,4]+dmap$X[,5]+dmap$X[,6])
