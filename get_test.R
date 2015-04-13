source('util.R')

get_test <- function(path, fileNameRoot){
	testSet <- read.table(paste(path, fileNameRoot, '.TEST', sep=''), 
		col.names=c('idx', 'snid', 'path', 'type'))
	dump <- get_dump('train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP')

	match.idx <- match(testSet$snid,dump$CID)
    testSet$type <- dump$GENTYPE[match.idx]

    testSet$type[testSet$type == 1] <- 'snIa'
    testSet$type[testSet$type == 21] <- 'snII'
    testSet$type[testSet$type == 22] <- 'snII'
    testSet$type[testSet$type == 23] <- 'snII'
    testSet$type[testSet$type == 3] <- 'snIbc'
    testSet$type[testSet$type == 32] <- 'snIbc'
    testSet$type[testSet$type == 33] <- 'snIbc'
    testSet$type<-as.factor(testSet$type)
    testSet$idx <- test.set$idx + 1
    rm(match.idx, dump)

    return(testSet)
}