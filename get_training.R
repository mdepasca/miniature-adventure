get_training <- function(path, fileNameRoot){
	fNames <- dir(path, pattern=paste(fileNameRoot, '.I*.TRAIN', sep=''))
    for (f in fNames){
    	if (!exist('trainingSet')){
    		trainingSet <- read.table(paste(path,f, sep=''), 
    			col.names=c('idx', 'snid', 'path', 'type'))
    	}else{
    		trainingSet <- rbind(trainingSet, read.table(paste(path,f, sep=''), 
    			col.names=c('idx', 'snid', 'path', 'type')))
    	}
    }
    trainingSet$idx <- trainingSet$idx + 1

    return(trainingSet)
}
