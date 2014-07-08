import numpy as np
import pylab as pb
pb.ion()
import GPy as GPy
import pickle as pkl
import gzip as gzip
from matplotlib import pyplot as plt
import snphotcc_infile_reformat as snphot
import time as time

# The dataset is a simulated one coming from the SuperNova Classification Challange (Kessler+ 2010).
# The ASCII files are read using code from Newling+ 2011. A catalog of LCs is produced.

# import utility
# 
# snCatalog = utility.create_sn_catalog('../DES_BLIND+HOSTZ')
# utility.pickle_sn_catalog(snCatalog, 'snCatalog.pkl')



# Opening pickle file containing the catalog


t1 = time.mktime(time.gmtime())
dirTrainData = 'train_data/'
fileZip = 'snCatalog.gz'
# f = open('snCatalog.pkl', 'rb')
print '  Unzipping traning catalog ...'
f = gzip.open(dirTrainData + fileZip, 'rb')
print '  Unpickling training catalog ...'
snCatalog = pkl.load(f)
f.close()

len(snCatalog.SNID)

gLimMag = 25.2
rLimMag = 25.4
iLimMag = 25.1
zLimMag = 24.9


# Looking for SN with ID 142 (whose LC as already been fitted with Faraway's R code)


#snIdx = np.where(snCatalog.SNID == 142)[0][0]
snIdx = np.random.random_integers(low=0, high=len(snCatalog.SNID))
print snCatalog.SNID[snIdx]



# Band Selection


offset = 3
nLines = 6
nCols = 8

for band in ('g', 'r', 'i', 'z'):
# <codecell>
# Plot some of the light curves on the same page (48 per page)
    fig, ax = plt.subplots(nrows=nLines, ncols=nCols)
    fig.suptitle('band %s' %band)
    plt.subplots_adjust(wspace=0.05, hspace=0.05, top=0.95, bottom=0.05, left=0.05, right=0.95)
    j = 0 # row index in tuple ax
    k = 0 # col index in tuple ax

    if band == 'g': LimMag = gLimMag
    if band == 'r': LimMag = rLimMag 
    if band == 'i': LimMag = iLimMag
    if band == 'z': LimMag = zLimMag       
    
    rangeLim = nLines * nCols

    for i in (range(rangeLim)):
        snIdx = i + offset * (nLines * nCols)
        numObs = len(snCatalog.sne[snIdx].lightCurvesDict[band].mjd)


#numObs = len(snCatalog.sne[snIdx].g.mjd)
        X = np.reshape(snCatalog.sne[snIdx].lightCurvesDict[band].mjd, (numObs, 1))
        X = X - np.min(X) # to avoid problems when producing the model
        Y = np.reshape(snCatalog.sne[snIdx].lightCurvesDict[band].flux, (numObs, 1))
        errY = np.reshape(snCatalog.sne[snIdx].lightCurvesDict[band].fluxErr, (numObs, 1))
        
        Ymag = snphot.flux_to_mag(np.squeeze(Y), LimMag, False)
        errYmag = snphot.error_to_mag(np.squeeze(Y), np.squeeze(errY))
        
        Ymag = np.reshape(Ymag, (numObs, 1))
        errYmag = np.reshape(errYmag, (numObs, 1))
        errYmag


        # 
        # Setting the kernel function from which depends the cross-validation between different inputs
        # 
        kern = GPy.kern.Bias(1) + GPy.kern.RatQuad(1)#GPy.kern.RBF(1)


        # 
        # Creating the GP model
        # 

        # 
        # Model 1
        # 
        m = GPy.models.GPHeteroscedasticRegression(X, Ymag, kern)
        m['.*Gaussian_noise'] = errYmag.flatten() #Set the noise parameters to the error in Y
        [m['.*Gaussian_noise_%s' %i].constrain_fixed() for i in range(numObs)] #Fix the noise parameters, we know its value so we don't need to learn it
        m.checkgrad(verbose=0) # 1
        m.optimize_restarts(num_restarts=10)
    #ax[j][k].set_title('sn_id%(snid)i_band%(photband)s' % {"snid": snCatalog.SNID[snIdx], "photband": band})
        m.plot_f(fignum=1, ax=ax[j][k])
        ax[j][k].errorbar(X.flatten(), Ymag.flatten(), yerr=errYmag.flatten(), fmt=None, ecolor='r', zorder=1)
    # ax[j][k].ylim(20, 28)
        ax[j][k].set_ylim(28, 20)
        ax[j][k].set_xticklabels([])
        ax[j][k].set_yticklabels([])

        fillX = np.linspace(ax[j][k].get_xlim()[0], ax[j][k].get_xlim()[1])
        fillY = fillX * 0 + LimMag
        fillY2 = fillX * 0 + ax[j][k].get_ylim()[0]
        ax[j][k].fill_between(fillX, fillY, fillY2, facecolor='gray', alpha=0.3)
        ax[j][k].text(ax[j][k].get_xlim()[0]+0.3*ax[j][k].get_xlim()[1], ax[j][k].get_ylim()[1]+0.05*ax[j][k].get_ylim()[0], 'idx %d' %snIdx, size='medium')
    # ax[j][k].set_title(snCatalog.sne[snIdx])
    # pb.title('Model 1')
    # ax[j][k].gca().invert_yaxis()

        if (k < nCols - 1):
            k += 1
        else:
            k = 0
            j += 1
    
    # figFile = 'img/samples/sn_sample%(block)03d_%(band)s.pdf' % {"block": offset+1, "band": band}
    # pb.savefig(figFile, format='pdf', dpi=300)

t2 = time.mktime(time.gmtime())
deltaT = time.localtime(t2 - t1)
print 'It took %(mins)d mins and %(secs)d secs.' %{"mins": deltaT[4], "secs": deltaT[5]}
