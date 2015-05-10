library(diffusionMap)
library(gclus)
library(R.utils)
library(randomForest)

# To display diffusionMap coordinates correlation

if (!exists('dmap')){
    message('diffusionMap object called `dmap` does not exists!')
}else{
     par.defaults <- par(no.readonly=TRUE)
     par(pch='.')
     min <- 2
     max <- 21
     homeDir <- '~/sandbox/miniature-adventure/'
     plotDir <- paste(homeDir, 'results/Plots/', sep='')
     fileName <- paste('dmap_scatter_plot_matrix', min, max, '.png', sep='_')
     mydata <- dmap$X[,c(min:max)]
     mydata.corr <- abs(cor(mydata))
     mycolors <- dmat.color(mydata.corr)
     myorder <- order.single(mydata.corr)
     png(file=paste(plotDir, fileName, sep=''), width=7, height=7, units='in', res=300)
     cpairs(mydata,
            myorder,
            panel.colors=mycolors,
            gap=.1,
            pch='.',
            upper.panel=NULL,
            labels=as.character(myorder+min-1),
            main="Variables Ordered and Colored by Correlation"
            )
     dev.off()
     message(paste(plotDir, fileName, sep=''))
     par(par.defaults)
}
