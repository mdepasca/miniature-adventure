get_training <- function(path, fileNameRoot){

    fNames <- list.files(path, pattern=paste(fileNameRoot, '\\.I[[:alnum:]]{1,2}\\.TRAIN', sep=''))

    for (f in fNames){
        message(paste('Reading file', f))
    	if (!exists('trainingSet')){
            trainingSet <- read.table(paste(path,f, sep=''), 
                                      col.names=c('idx', 'snid', 'path', 'type'))
    	}else{
             trainingSet <- rbind(trainingSet, read.table(paste(path,f, sep=''), 
                                                          col.names=c('idx', 'snid', 'path', 'type')))
         }
    }
    message('Increasing by 1 idx column')
    trainingSet$idx <- trainingSet$idx + 1

    return(trainingSet)
}
