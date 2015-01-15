redshift.distrib <- function(dump.file,save=FALSE){
    dump <- read.table(dump.file, header=TRUE, skip=1)
    message(is.name(dump$Z))
    if (is.name(dump$Z)){
        d <- density(dump$Z)
    }else {
        d <- density(dump$GENZ)
    }
    if (save){
        png(filename='redshift_distribution.png', width=1024, height= 1024, res=300)
        plot(d, main='Redshift distribution')
        dev.off()
    }else{
        plot(d, main='Redshift distribution')
    }
    return(dump)
}
