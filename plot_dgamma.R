plot_dgamma <- function(xlim=c(0,1000,1), gamma.shape, gamma.rate=1,
                        plot.xlim=c(0, 20), xlab='x', ylab='Gamma p.d.f.'){
    x <- seq(from=xlim[1], to=xlim[2], by=xlim[3])
    y <- dgamma(x, gamma.shape, gamma.rate)
    print(max(y))
    plot(x, y, xlim=plot.xlim, ylim=c(0,max(y)), type='l', xlab=xlab, ylab=ylab,
        col='#7f0000', lwd=2)
}
