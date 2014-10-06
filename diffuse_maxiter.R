# needs to import
library(diffusionMap)
library(Matrix)
library(igraph)

eps.val <- 0.06
delta <- 10^-3
max.iter <- 10000
neigen <- 120
cat('eps.val = ', eps.val, '\n')
cat('delta = ', delta, '\n')
cat('max.iter = ', max.iter, '\n')
cat('neigen = ', neigen, '\n')
diffuse.maxiter <- function(D, eps.val=epsilonCompute(D), neigen=NULL, t=0, maxdim=50, delta=10^-5, maxiter=NULL){
    start <- proc.time()[3]

    if ((is.matrix(D)) == FALSE){
        D <- as.matrix(D)
    }
    n <- dim(D)[1]
    K <- exp(-D^2/(eps.val))
    v <- sqrt(apply(K, 1, sum))
    A <- K/(v%*%t(v))                   #symmetric graph Laplacian
    
    ind <- which(A>delta, arr.ind=TRUE)
    Asp <- sparseMatrix(i=ind[,1], j=ind[,2], x=A[ind], dims=c(n, n))
    
    f <- function(x, a=NULL){
        as.matrix(A%*%x)
    }
    
    cat('Performing eigendecomposition\n')
    if (is.null(neigen)){
        neff <- min(maxdim+1, n)
    }else{
        neff <- min(neigen+1, n)
    }

    # eigeindecomposition using ARPACK
    if (is.null(maxiter)){
        arpack.options <- list(which='LA', nev=neff, n=n, ncv=max(30, 2 * neff))
    }else{
        arpack.options <- list(which='LA', nev=neff, n=n, ncv=max(30, 2 * neff), maxiter=maxiter)
    }

    decomp <- arpack(f, extra=Asp, sym=TRUE, options=arpack.options)
    psi <- decomp$vectors / (decomp$vectors[, 1]%*%matrix(1, 1, neff)) # right eigenvectors
    phi <- decomp$vectors * (decomp$vectors[, 1]%*%matrix(1, 1, neff)) # left eigenvectors
    eigenvals <- decomp$values          # eigenvalues

    cat('Computing diffusion coordinates\n')
    if (t<=0){                          # use multi-scale geometry
        lambda <- eigenvals[-1] / (1-eigenvals[-1])
        lambda <- rep(1, n) %*% t(lambda)
    }
    else{                               #use fixed scale t
        lambda <- eigenvals[-1]^t
        lambda <- rep(1, n) %*% t(lambda)
    }
    if (is.null(neigen)){# use no. of dimensions corresponding to 95% dropoff
        lam <- lambda[1, ]/lambda[1, 1]
        neigen <- min(which(lam<.05)) # default no. of eigenvalues
        neigen <- min(neigen, maxdim)
        eigenvals <- eigenvals[1:(neigen+1)]
        cat('Used default value:', neigen, 'dimensions\n')
    }
    X <- psi[, 2:(neigen+1)] * lambda[, 1:neigen] # diffusion coordinates X
    
    cat('Elapsed time:', signif(proc.time()[3]-start, digits=4), 'seconds\n')

    y <- list(X=X,
              phi0=phi[, 1],
              eigenvals=eigenvals[-1],
              eigenmult=lambda[1, 1:neigen],
              psi=psi,
              phi=phi,
              neigen=neigen,
              epsilon=eps.val)
    class(y) <- "dmap"
    return(y)
}


#load(file='distMat4520.RData')
filePath <- 'products/distance_matrix/Pymatrix_pc017261.ads.eso.org_1411167030.143.txt'
distMat <- matrix(scan(filePath, comment.char='#'), ncol=4520, nrow=4520, byrow=TRUE)

#distMat2 <- distMat[1:1000, 1:1000]

dmap <- diffuse.maxiter(distMat, eps.val=eps.val, neigen=neigen, delta=delta, maxiter=max.iter)
