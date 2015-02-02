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
