library(randomForest)

classify_RF <- function(dmap, trainingSet, testSet){
	## extracting training and test set diffusion coordinate
    dmapCoord.train <- data.frame(dmap$X[trainingSet$idx,])
    dmapCoord.test <- data.frame(dmap$X[testSet$idx,])
    ## adding `type' column
    dmapCoord.train <- cbind(dmapCoord.train, trainingSet$type)
    dmapCoord.test <- cbind(dmapCoord.test, testSet$type)
    
    ncols <- dim(dmapCoord.train)[2]
    names(dmapCoord.train)[ncols] <- 'type'
    names(dmapCoord.test)[ncols] <- 'type'

    types <- c('snIa', 'snIbc', 'snII')

    ## training the random forest classifier and testing it
    sn.rf <- randomForest(x=dmapCoord.train[,1:(ncols-1)], 
    	y=as.factor(dmapCoord.train[,ncols]), 
    	xtest=dmapCoord.test[,1:(ncols-1)], 
    	ytest=as.factor(dmapCoord.test[,ncols]), importance=TRUE)
 	
 	return(sn.rf)
}