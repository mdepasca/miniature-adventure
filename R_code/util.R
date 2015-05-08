library(Hmisc)
## library(sm) # for comparative density plots.
library(plotrix) # boxed labels

vec.nObs <- function(path){
    file.list <- list.files(path=path, pattern='\\.DAT$', full.names=TRUE)
    nObs <- c()
    for (f in file.list){
        print(f)
        r <- readLines(f)
        skip <- grep('OBS:', r)
        tmp <- read.table(f, header=TRUE, skip=skip[2]-1, fill=TRUE) 
        nObs <- c(nObs, dimension(tmp)[1])
        rm(tmp)
    }
    return(nObs)
}

get_dump <- function(dump.file){
    dump <- read.table(dump.file, header=TRUE, skip=1)
    dump$GENTYPE[which(dump$GENTYPE == 21)] <- 2
    dump$GENTYPE[which(dump$GENTYPE == 22)] <- 2
    dump$GENTYPE[which(dump$GENTYPE == 23)] <- 2

    dump$GENTYPE[which(dump$GENTYPE == 32)] <- 3
    dump$GENTYPE[which(dump$GENTYPE == 33)] <- 3

    dump$GENTYPE <- as.factor(dump$GENTYPE)

    return(dump)
}

snphotcc.param.distrib <- function(dump.df, param='GENZ', save=FALSE){
    if (class(param) == 'character'){
        param.numeric <- which(names(dump.df) == param)
    }else{
        param.numeric <- param
    }
    d.param <- density(as.numeric(dump.df[param.numeric]))
    plot(d, xlab=paste(param), main='')
}

redshift.distrib <- function(dump.df, save=FALSE){
    # dump <- read.table(dump.file, header=TRUE, skip=1)
    # par(mfrow=c(1,1))

    d <- density(dump.df$GENZ, bw=0.03)
    z.ref <- d$x[which(d$y==max(d$y))]
    z.thr <- 0.17
    if (save){
        message('saving plot on file...')
        pdf(file='redshift_distribution.pdf', width=7, height= 7)
    }
    # plot.new()
    plot(d, xlab='Redshift', main='')#, sub=paste('N =', length(dump.df$CID),
        # ' Bandwidth =', d$bw, sep=' '))
    ## polygon(c(0,d$x[which(d$x<=z.thr)],z.thr), c(0,d$y[which(d$x<=z.thr)],0),
    ##     density=10, angle=45)
    ## boxed.labels(0.1, 0.1, labels=paste(length(
    ##     dump.df$GENZ[which(dump.df$GENZ<z.thr)])), cex=1, border=FALSE)

    ## polygon(c(z.ref-z.thr/2, d$x[which((z.ref-z.thr/2)<d$x &
    ##     d$x<(z.ref+z.thr/2))], z.ref+z.thr/2),
    ## c(0,d$y[which((z.ref-z.thr/2)<d$x & d$x<(z.ref+z.thr/2))],0),
    ## density=10, angle=135, lty=1)
    ## boxed.labels(z.ref, 0.7, labels=paste(length(
    ##     dump.df$GENZ[which((z.ref-z.thr/2)<dump.df$GENZ &
    ##     dump.df$GENZ<(z.ref+z.thr/2))])), cex=1, border=FALSE, ypad=1.5)
        # col='darkgray', border='darkgray')
    # lines(d, xlab='redshift', main='', sub=paste('N =', length(dump.df$CID),
    #     ' Bandwidth =', d$bw, sep=' '))

    rug(d$x[which(d$y==max(d$y))])
    if (save){
        dev.off()
        message('done!')
    }
    return(d)
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


compare.densities <- function(dump.df,save=FALSE, train=FALSE){
    if (save){
        pdf(file='compare_redshift_distribution.pdf', width=7, height= 7)
    }
    gentype.f <- factor(dump.df$GENTYPE, levels=c(1,2,3),
        labels=c("type Ia", "type II", "type Ib/c"))
    d <- density(dump.df$GENZ, bw=0.03)
    if (train){
        train.df <- read.table('results/SIMGEN_PUBLIC_DES.TRAIN',
            col.names=c('idx', 'snid', 'path', 'type'))
        dump.df <- dump.df[match(train.df$snid, dump.df$CID),]
    }

    d1 <- density(dump.df$GENZ[which(dump.df$GENTYPE == 1)], bw=0.03)
    d2 <- density(dump.df$GENZ[which(dump.df$GENTYPE == 2)], bw=0.03)
    d3 <- density(dump.df$GENZ[which(dump.df$GENTYPE == 3)], bw=0.03)

    z.ref <- d$x[which(d$y==max(d$y))]
    z.thr <- 0.17/2

    colfill <- c('darkorange', 'seagreen3', 'blue2')
    par.default <- par(no.readonly=TRUE)
    par(mfrow=c(3,1), lty=3, ann=FALSE, mar=c(0,0,0,0), oma=c(5,6,4,2))
    plot(d3, col=colfill[3], main='', xlab='Redshift', xlim=c(0,max(d1$x)),
        ylim=c(0,3.5), axes=FALSE, type='l')
    # plot(d3, col=colfill[3], main='', xlab='Redshift', xlim=c(0,max(d1$x)),
    #     ylim=c(0,3.5), axes=FALSE)
    axis(1, at=seq(0,max(d1$x),0.2), label=FALSE, cex.axis=1.5)
    axis(2, at=seq(0,3.5,1), cex.axis=1.5)

    box()
    polygon(c(z.ref-z.thr,d3$x[which((z.ref-z.thr)<d3$x &
        d3$x<(z.ref+z.thr))],z.ref+z.thr),
    c(0,d3$y[which((z.ref-z.thr)<d3$x & d3$x<(z.ref+z.thr))],0),
    col=colfill[3], border=colfill[3], density=10, angle=45, lty=1)
    mtext("type Ib/c", side=3, line=-2, adj=0.9)
    mtext(paste("N = ", length(which(dump.df$GENTYPE == 3))), line=-4, adj=0.9)
    boxed.labels(z.ref-.01, 0.5, labels=paste(length(dump.df$GENZ[which((dump.df$GENTYPE == 3)
        & (z.ref-z.thr)<dump.df$GENZ & dump.df$GENZ<(z.ref+z.thr))])), cex=1, border=FALSE)

    par(lty=2)
    plot(d2,col=colfill[2], main='', xlab='Redshift', xlim=c(0,max(d1$x)),
        ylim=c(0,3.5), axes=FALSE)
    axis(1, at=seq(0,max(d1$x),0.2), label=FALSE, cex.axis=1.5)
    axis(2, at=seq(0,3.5,1), cex.axis=1.5)
    box()
    polygon(c(z.ref-z.thr,d2$x[which((z.ref-z.thr)<d2$x &
        d2$x<(z.ref+z.thr))],z.ref+z.thr),
    c(0,d2$y[which((z.ref-z.thr)<d2$x & d2$x<(z.ref+z.thr))],0),
    col=colfill[2], border=colfill[2], density=10, angle=45, lty=1)
    mtext("type II", side=3, line=-2, adj=0.9)
    mtext(paste("N = ", length(which(dump.df$GENTYPE == 2))), line=-4, adj=0.93)
    boxed.labels(z.ref, 0.5, labels=paste(length(dump.df$GENZ[which((dump.df$GENTYPE == 2)
        & (z.ref-z.thr)<dump.df$GENZ & dump.df$GENZ<(z.ref+z.thr))])), cex=1, border=FALSE)

    par(lty=1)
    plot(d1, col=colfill[1], main='', xlab='Redshift', xlim=c(0,max(d1$x)),
        ylim=c(0,3.5), axes=FALSE)
    axis(1, at=seq(0,max(d1$x),0.2), cex.axis=1.5)
    axis(2, at=seq(0,3.5,1), cex.axis=1.5)
    box()
    polygon(c(z.ref-z.thr,d1$x[which((z.ref-z.thr)<d1$x &
        d1$x<(z.ref+z.thr))],z.ref+z.thr),
    c(0,d1$y[which((z.ref-z.thr)<d1$x & d1$x<(z.ref+z.thr))],0),
    col=colfill[1], border=colfill[1], density=10, angle=45, lty=1)

    rug(z.ref)

    boxed.labels(z.ref, 0.5, labels=paste(length(dump.df$GENZ[which((dump.df$GENTYPE == 1)
        & (z.ref-z.thr)<dump.df$GENZ & dump.df$GENZ<(z.ref+z.thr))])), cex=1, border=FALSE)
    mtext("type Ia", side=3, line=-2, adj=0.9)
    mtext(paste("N = ", length(which(dump.df$GENTYPE == 1))), line=-4, adj=0.92)
    mtext("Redshift", side=1, outer=TRUE, line=3.2)
    mtext("Density", side=2, outer=TRUE, line=3.5)
    # sm.density.compare(dump.df$GENZ, dump.df$GENTYPE, h=0.03, xlab="redshift")
    # if (train){
    #     legend(0.8, 3, levels(gentype.f), fill=colfill)
    # }
    # else{
    #     legend(0, 1.8, levels(gentype.f), fill=colfill)
    # }
    if (save){
        dev.off()
    }
    par(par.default)
}

plot.simulated.sample <- function(dir.path, lc.num=20, band='r', save=FALSE){
    file.list <- list.files(path=dir.path, pattern='\\.DAT$', full.names=TRUE)
    # plot.new()
    # colors <- grey(0:(lc.num-1)/(lc.num+1))
    colors <- heat.colors(lc.num)
    ran <- sample(1:length(file.list), lc.num, replace=FALSE)
    for (i in c(1:lc.num)){
    # for (f in file.list[1:20]){
        r <- readLines(file.list[ran[i]])
        # next line gives the position of the row right before the head of the
        #
        # table containing the data. In this way, all is written before can be
        #
        # skipped.
        skip <- grep('NVAR:', r)
        lc <- read.table(file.list[ran[i]], header=TRUE, skip=skip, fill=TRUE)
        if (class(lc$FLUXCAL) != 'numeric'){
            del.rows <- which(is.na(lc$FLUXCALERR))
            lc <- lc[-del.rows,]
            lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
        }
        if (i == 1){
            plot(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
                type='l',
                xlim=c(56171, 56351), ylim=c(0,200000),
                xlab='epoch [mjd]', ylab='flux [adu]', col=colors[i])
        }else{
            lines(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
                col=colors[i])
        }
        errbar(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
            yplus=lc$FLUXCAL[which(lc$FLT==band)]+lc$FLUXCALERR[which(lc$FLT==band)],
            yminus=lc$FLUXCAL[which(lc$FLT==band)]-lc$FLUXCALERR[which(lc$FLT==band)], add=T)
    }
    return (file.list)
}


plot.lightcurve <- function(file.path, save=FALSE){
    r <- readLines(file.path)
    # gsub("[[:space:]]", "", x)
    sn.type <- gsub('[[:space:]]', '',
        unlist(strsplit(r[grep('SNTYPE:', r)], ':'))[2])

    if (sn.type == 1){
        sn.type.lab <- 'type Ia'
    }else if (sn.type %in% c(2,21,22,23)){
        sn.type.lab <- 'type II'
    }else{
        sn.type.lab <- 'type Ib/c'
    }
    message(sn.type.lab)

    skip <- grep('NVAR:', r)
    lc <- read.table(file.path, header=TRUE, skip=skip, fill=TRUE)
    if (class(lc$FLUXCAL) != 'numeric'){
        del.rows <- which(is.na(lc$FLUXCALERR))
        lc <- lc[-del.rows,]
        lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
    }
    if (save){
        pdf(file='four_bands_lightcurve.pdf', width=7, height= 7)
    }
    par.default <- par(no.readonly=TRUE)
    par(mfrow=c(4,1), ann=FALSE, mar=c(0,0,0,0), oma=c(5,6,4,2))
    for (band in c('g','r','i','z')){
        plot(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
            axes=FALSE, pch=20,
            xlim=c(floor(min(lc$MJD)), floor(max(lc$MJD))),
            ylim=c(floor(min(lc$FLUXCAL)), floor(max(lc$FLUXCAL))))
        errbar(lc$MJD[which(lc$FLT==band)], lc$FLUXCAL[which(lc$FLT==band)],
            yplus=lc$FLUXCAL[which(lc$FLT==band)]+lc$FLUXCALERR[which(lc$FLT==band)],
            yminus=lc$FLUXCAL[which(lc$FLT==band)]-lc$FLUXCALERR[which(lc$FLT==band)],
            add=T, cap=0)
        axis(2, at=seq(0, max(lc$FLUXCAL), 30000), cex.axis=1.5)
        axis(1, at=seq(floor(min(lc$MJD)), floor(max(lc$MJD)), 10), label=FALSE,
            cex.axis=1.5, tick=TRUE)
        box()
        mtext(band, side=3, line=-2, adj=0.9)
    }
    axis(1, at=seq(floor(min(lc$MJD)), floor(max(lc$MJD)), 10), cex.axis=1.5)
    mtext("Epoch [mjd]", side=1, outer=TRUE, line=3.2)
    mtext("Flux [adu]", side=2, outer=TRUE, line=3.5)
    if (save){
        dev.off()
    }
    par(par.default)
    return(lc)
}

plot.screen <- function(dir.path, save=FALSE){
    par.default <- par(no.readonly=TRUE)
    file.list <- list.files(path=dir.path, pattern='\\.DAT$', full.names=TRUE)
    ran <- sample(1:length(file.list), length(file.list), replace=FALSE)
    m <- rbind(
        c(0.1-0.035, 0.5-0.035, 0.9-0.035, 1-0.035),
        c(0.1-0.035, 0.5-0.035, 0.8-0.035, 0.9-0.035),
        c(0.1-0.035, 0.5-0.035, 0.7-0.035, 0.8-0.035),
        c(0.1-0.035, 0.5-0.035, 0.6-0.035, 0.7-0.035),
        #-----
        c(0.6-0.035, 1-0.035, 0.9-0.035, 1-0.035),
        c(0.6-0.035, 1-0.035, 0.8-0.035, 0.9-0.035),
        c(0.6-0.035, 1-0.035, 0.7-0.035, 0.8-0.035),
        c(0.6-0.035, 1-0.035, 0.6-0.035, 0.7-0.035),
        #-----
        c(0.1-0.035, 0.5-0.035, 0.4-0.035, 0.5-0.035),
        c(0.1-0.035, 0.5-0.035, 0.3-0.035, 0.4-0.035),
        c(0.1-0.035, 0.5-0.035, 0.2-0.035, 0.3-0.035),
        c(0.1-0.035, 0.5-0.035, 0.1-0.035, 0.2-0.035),
        #-----
        c(0.6-0.035, 1-0.035, 0.4-0.035, 0.5-0.035),
        c(0.6-0.035, 1-0.035, 0.3-0.035, 0.4-0.035),
        c(0.6-0.035, 1-0.035, 0.2-0.035, 0.3-0.035),
        c(0.6-0.035, 1-0.035, 0.1-0.035, 0.2-0.035)
    )
    if (save){
        pdf(file='four_bands_four_lightcurve.pdf', width=8.7, height=11.7)
    }
    split.screen(m)
    band <- c('g','r','i','z')
    for (i in (1:16)){
        screen(i)
        par(mar = c(0, 0, 0, 0), new=TRUE, tcl=-0.3, mgp=c(3,0.3,0))
        if ((i %% 4) == 1){ # open new file
            while (TRUE){
                r <- readLines(file.list[ran[1]])

                skip <- grep('NVAR:', r)
                lc <- read.table(file.list[ran[1]], header=TRUE, skip=skip, fill=TRUE)

                if (class(lc$FLUXCAL) != 'numeric'){
                    del.rows <- which(is.na(lc$FLUXCALERR))
                    lc <- lc[-del.rows,]
                    lc$FLUXCAL <- as.numeric(levels(lc$FLUXCAL))[lc$FLUXCAL]
                }
                ran <- ran[-1]
                if (!is.na(floor(min(lc$MJD))) && !is.na(floor(max(lc$MJD)))) {
                    break
                }
            }

            r.max <- max(lc$FLUXCAL[which(lc$FLT=='r')])
            flux.max <- max(lc$FLUXCAL)
            mjd.r.max <- round(lc$MJD[which(lc$FLUXCAL==r.max)])
            message(mjd.r.max)
            lc$MJD <- lc$MJD - mjd.r.max
            lc$FLUXCAL <- lc$FLUXCAL/r.max*100
            lc$FLUXCALERR <- lc$FLUXCALERR/r.max*100
        }

        if (i%%4 == 0){ax <- TRUE}
        plot(lc$MJD[which(lc$FLT==band[i%%4+1])],
            lc$FLUXCAL[which(lc$FLT==band[i%%4+1])],
            axes=FALSE, pch=20, cex=0.5,
            xlim=c(floor(min(lc$MJD)), floor(max(lc$MJD))),
            ylim=c(floor(min(lc$FLUXCAL)), floor(max(lc$FLUXCAL))))
        errbar(lc$MJD[which(lc$FLT==band[i%%4+1])],
            lc$FLUXCAL[which(lc$FLT==band[i%%4+1])],
            yplus=(lc$FLUXCAL[which(lc$FLT==band[i%%4+1])]+
                lc$FLUXCALERR[which(lc$FLT==band[i%%4+1])]),
            yminus=(lc$FLUXCAL[which(lc$FLT==band[i%%4+1])]-
                lc$FLUXCALERR[which(lc$FLT==band[i%%4+1])]),
            add=T, cap=0, cex=0.5)
        box()
        axis(2, at=seq(0, max(lc$FLUXCAL), 40), cex.axis=0.8)

        axis(1, at=seq(floor(min(lc$MJD)), floor(max(lc$MJD)), 10), label=FALSE,
            tick=TRUE, cex.axis=1)

        if (i %in% c(2,6,10,14)){
            ylab.coord <- c(0,floor(min(lc$FLUXCALERR[which(lc$FLT==band[i%%4+1])])))
            mtext(paste("Flux / ", r.max/100," [adu]", sep=''), side=2, at=ylab.coord, line=1.5, cex=0.8)
        }

        if (i %in% c(1,5,9,13)){mtext(band[1], side=3, line=-1.5, adj=0.95)}
        if (i %in% c(2,6,10,14)){mtext(band[2], side=3, line=-1.5, adj=0.95)}
        if (i %in% c(3,7,11,15)){mtext(band[3], side=3, line=-1.5, adj=0.95)}
        if (i %in% c(4,8,12,16)){mtext(band[4], side=3, line=-1.5, adj=0.95)}

        if ((i %% 4) == 0){
            axis(1, at=seq(floor(min(lc$MJD)), floor(max(lc$MJD)), 10), cex.axis=0.8)
            mtext(paste("Epoch - ", mjd.r.max," [mjd]", sep=''), side=1, line=1.5, cex=0.8)
        }
        if ((i%%4) == 1){
            mtext("SN type Ia - MLCS2kUV", side=3, line=0.2, cex=0.8)
        }
    }
    par(par.default)
    close.screen(all.screens = TRUE)
    if (save){
        dev.off()
    }
}
