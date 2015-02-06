library(Hmisc)

redshift.distrib <- function(dump.file,save=FALSE){
    dump <- read.table(dump.file, header=TRUE, skip=1)

    if (is.name(dump$Z)){
        d <- density(dump$Z)
    }else {
        d <- density(dump$GENZ)
    }
    if (save){
        message('saving plot on file redshift_distribution.png...')
        png(filename='redshift_distribution.png', width=4, height= 4, units='in', res=300)
        plot(d, main='Redshift distribution')
        dev.off()
        message('done!')
    }else{
        plot(d, main='Redshift distribution')
    }
    return(dump)
}

explore.lightcurves <- function(dir.path, band='r'){
    file.list <- list.files(path=dir.path, pattern='\\.DAT$', full.names=TRUE)
    for (f in file.list){
        print(f)
        lc <- read.table(f, header=TRUE, skip=63, fill=TRUE)
        del.rows <- which(is.na(lc$FLUXCALERR))
        lc <- lc[-del.rows,]
        lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
        print(class(lc$FLUXCALERR))
        print(lc$FLUXCALERR)
        ## plot(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)])
        errbar(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
               lc$FLUXCAL[which(lc$FLT==band)]+lc$FLUXCALERR[which(lc$FLT==band)],
               lc$FLUXCAL[which(lc$FLT==band)]-lc$FLUXCALERR[which(lc$FLT==band)])
        cin <- readline('Press Enter to next plot: ')
    }
    return(file.list)
}

explore.minmax <- function(dir.path, band='r'){
    min.vec <- c()
    max.vec <- c()
    file.list <- list.files(path=dir.path, pattern='\\.DAT$', full.names=TRUE)
    for (f in file.list){
        lc <- read.table(f, header=TRUE, skip=63, fill=TRUE)
        del.rows <- which(is.na(lc$FLUXCALERR))
        lc <- lc[-del.rows,]
        lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
        min.vec <- c(min.vec, min(lc$MJD))
        max.vec <- c(max.vec, max(lc$MJD))
    }
    return(c(min(min.vec), max(max.vec)))
}
