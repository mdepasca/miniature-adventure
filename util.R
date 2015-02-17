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

explore.lightcurves <- function(dir.path, band='r', skip=63, err.bars=FALSE, multi.plot=FALSE, save=FALSE){
    file.list <- list.files(path=dir.path, pattern='\\.DAT$', full.names=TRUE)
    lc.n <- 50                          #number of light curves to plot
    ## print(class(file.list))
    if (multi.plot){
        pchs <- c(0:25)
        ## pch.col <- heat.colors(length(file.list))
        pch.col <- gray(0:(lc.n-1) / lc.n + 0.01)
        i <- 0
        j <- 1
    }
    if (save){
        png(filename='light_curves.png', width=5, height= 5, units='in', res=300)
    }
    title <- paste(lc.n, "SN non-1a simulated by SNANA", sep=' ')
    for (f in file.list[1:lc.n]){
        ## print(file.list[1])
        
        lc <- read.table(f, header=TRUE, skip=skip, fill=TRUE)
        del.rows <- which(is.na(lc$FLUXCALERR))
        lc <- lc[-del.rows,]
        lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
        ## print(class(lc$FLUXCALERR))
        ## print(lc$FLUXCALERR)
        if (err.bars){
        ## 
            errbar(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
                   lc$FLUXCAL[which(lc$FLT==band)]+lc$FLUXCALERR[which(lc$FLT==band)],
                   lc$FLUXCAL[which(lc$FLT==band)]-lc$FLUXCALERR[which(lc$FLT==band)])
        }else{
            if (!multi.plot){
                plot(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)])
            }
        }
        if (multi.plot){
            if (f == file.list[1]){
                plot(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
                     col=pch.col[j],
                     xlim=c(56170,56360), ylim=c(0,50000), type='l',
                     xlab='epoch [MJD]', ylab='flux [ADU]',
                     main=title)
            }else{
                lines(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)], pch=1, col=pch.col[j])
            }
            ## print(pch.col[j])
            if (i == length(pchs)) i <- 0 else i <- i+1
            if (j == length(pch.col)) j <- 1 else j <- j+1
        }else{
            cin <- readline('Press Enter to next plot: ')
        }
    }
    if (save){
        dev.off()
    }
    ## return(file.list)
}

explore.minmax <- function(dir.path, band='r', save=FALSE){
    mjd.min.vec <- c()
    mjd.max.vec <- c()
    mjd.flux.max.vec <- c()
    idx.flux.max.vec <- c()
    ext.mjd <- c()
    file.list <- list.files(path=dir.path, pattern='\\.DAT$', full.names=TRUE)
    for (f in file.list){
        # message(f)
        r <- readLines(f)
        skip <- grep('NVAR:', r)
        lc <- read.table(f, header=TRUE, skip=skip, fill=TRUE)
        if (class(lc$FLUXCAL) != 'numeric'){
            del.rows <- which(is.na(lc$FLUXCALERR))
            lc <- lc[-del.rows,]
            lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
        }
        mjd.min.vec <- c(mjd.min.vec, min(lc$MJD))
        mjd.max.vec <- c(mjd.max.vec, max(lc$MJD))
        ext.mjd <- c(ext.mjd, (max(lc$MJD)-min(lc$MJD)))
        mjd.flux.max.vec <- c(mjd.flux.max.vec, lc$MJD[which(lc$FLUXCAL == max(lc$FLUXCAL), arr.ind=TRUE)])
        idx.flux.max.vec <- c(idx.flux.max.vec, which(lc$FLUXCAL == max(lc$FLUXCAL), arr.ind=TRUE))
    }
    ext.mjd <- ext.mjd[!is.na(ext.mjd)]
    d.mjd <- density(mjd.flux.max.vec)
    d.idx <- density(idx.flux.max.vec)
    d.ext <- density(ext.mjd)
    if (save){
        png(filename='results/plots/mjd_distributions.png', width=7, height= 7, units='in', res=300)
    }
    par(mfrow=c(3,1))
    plot(d.mjd, main='Distribution of epoch of maximum light', xlab='epoch [MJD]')
    rug(c(min(mjd.min.vec), max(mjd.max.vec), median(mjd.flux.max.vec)))
    plot(d.idx, main='Distribution of index of maximum light')
    plot(d.ext, main='Distribution of light curves duration', xlab='epoch [MJD]')

    par(mfrow=c(1,1))
    if (save){dev.off()}
    message('MIN MJD   MAX MJD    MIN MJD MAX    MAX MJD MAX ')
    return(c(min(mjd.min.vec), max(mjd.max.vec), min(mjd.flux.max.vec), max(mjd.flux.max.vec), ext.mjd))
}
