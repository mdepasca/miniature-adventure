source('R_code/util.R')

get_test <- function(path, fileNameRoot){
    ## column `type` is set to 'test' for each SN
    testSet <- read.table(paste(path, fileNameRoot, '.TEST', sep=''),
                          stringsAsFactors=FALSE,
                          col.names=c('idx', 'snid', 'path'))
    
    dump <- get_dump('train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP')

    match.idx <- match(testSet$snid, dump$CID)

    testSet <- cbind(testSet, dump$GENTYPE[match.idx])
    colnames(testSet)[length(colnames(testSet))] <- 'type'

    testSet$type <- sapply(testSet$type, as.numeric)
        
    testSet$type[which(testSet$type == 1 )] <- 'snIa'
    testSet$type[which(testSet$type == 2 )] <- 'snII'
    testSet$type[which(testSet$type == 21)] <- 'snII'
    testSet$type[which(testSet$type == 22)] <- 'snII'
    testSet$type[which(testSet$type == 23)] <- 'snII'
    testSet$type[which(testSet$type == 3 )] <- 'snIbc'
    testSet$type[which(testSet$type == 32)] <- 'snIbc'
    testSet$type[which(testSet$type == 33)] <- 'snIbc'

    testSet$type <- as.factor(testSet$type)
    testSet$idx <- testSet$idx + 1
    
    rm(match.idx, dump)

    return(testSet)
}
