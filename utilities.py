"""
Implementation of general use functions.
"""
import numpy as np
import numpy.ma as ma
import pickle as pkl
import gzip as gzip
import GPy as GPy
import classes as classes
import time as time
import argparse
import matplotlib.pyplot as plt

if __name__ == '__main__':
    """
    Parsing input of parameters for test.
    """
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description="Test of general functions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-s", "--sn-candidate", dest="candidate", 
        type=np.int64, default=None, 
        help="Candidate id")
    
    parser.add_argument("-b", "--band", dest="band", default='r', 
                        help="Photometric band.")
    
    parser.add_argument(
        "-c", "--catalog", dest="catalog", default=None,
        help="SN catalog.")
    
    parser.add_argument(
        "-p", "--catalog_path", dest="path", 
        default="train_data/snCatalog.gz",
        help="Path to SN catalog.")

    parser.add_argument(
        "-t", "--test-lengthscale", dest="testLength",
        default=False, help="Flag to randomize the lengthscale parameter.")

    parser.add_argument(
        "-m", "--magnitudes", dest="mag",
        default=False, help="Flag to switch to magnitudes.")

    args = parser.parse_args()


def flux_to_mag(flux, limFlux=False):
    """
    Converts fluxes to magnitudes using the following law (from Kessler+2010):
        flux = 10^(-0.4 * m + 11) => m = -2.5 * (log_10(flux) - 11)
    If the flux is below the limit of the instrument, returns the 
    corresponding limiting magnitude.

    INPUT:
        flux: numpy array of fluxes
        limFlux: specifies the limiting flux
        
    OUTPUT:
        mag: magnitude-converted fluxes
    """

    if limFlux:
        limMag = -2.5 * (-11 + np.log10(limFlux))

        # applying the mask to detection below the limiting flux
        maFlux = ma.masked_where(flux < limFlux, flux)
        
        # to avoid warnings due to values passed to np.log10
        # fluxMask = maFlux.mask
        # maMag = -2.5 * (-11.0 + np.log10(ma.filled(maFlux,1)))
        maMag = -2.5 * (-11.0 + np.log10(maFlux))
        
        mag = ma.filled(maMag, limMag)
    else:
        if flux > 0:
            mag = -2.5 * (-11 + np.log10(flux))
        else:
            mag = None

    return mag


def flux_error_to_mag_error(fluxErr, flux):
    # error prop. zErr = abs(dF(x)/dx) xErr
    magErr = np.multiply(np.abs(np.divide(-2.5, np.log(10) * flux)), 
                         np.abs(fluxErr)) 

    return magErr


def mag_to_flux(mag, limMag=False):
    """
    Converts magnitudes to fluxes using the following law (from Kessler+2010):
        flux = 10^(-0.4 * m + 11)
    If the magnitude is above the limit of the instrument, returns the 
    corresponding limiting flux.

    INPUT:
        mag: numpy array of magnitude
        limMag: specifies the limiting magnitude

    OUTPUT:
        flux: flux-converted magnitudes
    """

    if limMag:
        exponent = (-0.4 * limMag) + 11
        limFlux = np.power(10, exponent)

        # masked arrays if magnitude is grater then limit
        maMag = ma.masked_where(mag > limMag, mag)
        maExponent = (-0.4 * maMag) + 11
        maFlux =  np.power(10, exponent)
        maFlux.set_fill_value(limFlux)

        flux = ma.filled(maFlux)
    else:
        exponent = (-0.4 * mag) + 11
        flux = 10**exponent

    return flux


def open_gzip_pkl_catalog(path):
    f = gzip.open(path, 'rb')
    catalog = pkl.load(f)
    f.close()

    return catalog


def pick_random_sn(catalog, band):
    """
    Extract random observation in specified band from catalog. 
    Returns phase, flux, flux errors arrays and index in the catalog.
    Phase has zeropoint on the maximum flux in r band
    """
    idx = np.random.random_integers(low=0, high=len(catalog.SNID))
    numObs = len(catalog.sne[idx].lightCurvesDict[band].mjd)

    phase = catalog.sne[idx].lightCurvesDict[band].mjd
    # phase = phase - phase[catalog.sne[idx].lightCurvesDict['r'].flux.argmax()]
    phase = phase - phase.min()

    flux = catalog.sne[idx].lightCurvesDict[band].flux
    
    errFlux = catalog.sne[idx].lightCurvesDict[band].fluxErr
    
    return phase, flux, errFlux, idx

def get_sn(catalog, band, idx):
    """
    Extract specified supernova observation in specified band from catalog.
    Returns time, flux, flux errors arrays.
    Time has zeropoint on the maximum flux in r band
    """

    numObs = len(catalog.sne[idx].lightCurvesDict[band].mjd)

    phase = catalog.sne[idx].lightCurvesDict[band].mjd
    # phase = phase - phase[catalog.sne[idx].lightCurvesDict['r'].flux.argmax()]
    phase = phase - phase.min()

    flux = catalog.sne[idx].lightCurvesDict[band].flux
    
    errFlux = catalog.sne[idx].lightCurvesDict[band].fluxErr
    
    return phase, flux, errFlux

def get_sn_from_file(pathToSN):
    """Reads photometric data of SN from file formatted by SNPhotCC"""

    sn = classes.Supernova(pathToSN)
    return sn

def reshape_for_GPy(vec):
    return np.reshape(vec, (len(vec), 1))


def gp_fit(X, Y, errY, kernel, n_restarts=0, test_length=False):
    """
    Performs gaussian process regression
    NOTE
    check on shap of input should be added
    """
 
    medPhaseStep = np.median(np.abs(np.subtract(phase[0:-2], phase[1:-1])))
    maxPhaseStep = phase[-1] - phase[0]
    print "  Median phase step {:<5.3f}".format(medPhaseStep)
    print "  Time range {:<5.3f}".format(maxPhaseStep)

    rsX = reshape_for_GPy(X)
    rsY = reshape_for_GPy(Y)

    gpModel = GPy.models.GPHeteroscedasticRegression(rsX, rsY, kernel)
    gpModel['.*Gaussian_noise'] = errY

    [gpModel['.*Gaussian_noise_%s' %i].constrain_fixed(warning=False) 
     for i in range(X.size)]

    if test_length:
        np.random.RandomState
        length = np.random.uniform(low=medPhaseStep, high=maxPhaseStep)
        print "  Randomized lengthscale {:<5.3f}".format(length)
        gpModel['.*lengthscale'].constrain_fixed(length, warning=True)
    else:
        pass
        # gpModel['.*lengthscale'].constrain_bounded(
        #     2*medPhaseStep, maxPhaseStep,
        #     warning=True)

    if n_restarts > 0:
        gpModel.optimize_restarts(num_restarts=n_restarts, 
                                  verbose=False,
                                  robust=True,
                                  messages=False)

    predPhase = reshape_for_GPy(np.arange(X.min(), X.max(), 1))

    meanY, var = gpModel._raw_predict(predPhase, full_cov=True)
    return meanY, var, gpModel

#
#
# Module Testing
#
#
if __name__ == '__main__':
    # these limiting magnitudes are for simulations of DES survey. They should
    # be changed.
    limMag = {'g': 25.2, 
              'r': 25.4, 
              'i': 25.1, 
              'z': 24.9}
    
    dataSN = "../DES_BLIND+HOSTZ/"

    if args.catalog is None:
        # args.catalog = open_gzip_pkl_catalog(args.path)
        if args.candidate is None:
            args.candidate = np.random.random_integers(
                                low=0, high=18321)# num SN in SNPhotCC

            # phase, flux, errFlux, args.candidate = pick_random_sn(
            #                                         args.catalog, args.band)
        pathToSN = dataSN + "DES_SN" + "{:>06}".format(args.candidate) + ".DAT"
        sn = get_sn_from_file(pathToSN)

        phase = sn.lightCurvesDict[args.band].mjd
        flux = sn.lightCurvesDict[args.band].flux
        errFlux = sn.lightCurvesDict[args.band].fluxErr
    else:
        pass
        # phase, flux, errFlux = get_sn(
        #                         args.catalog, args.band, args.candidate)
    
    print "  Candidate ID {:>06}".format(args.candidate)
    print "  Testing lengthscale {:<5}".format(args.testLength)
    print "  Magnitudes {:<5}".format(args.mag)
    
    # tranforming to magnitudes
    if args.mag:
        limFlux = mag_to_flux(limMag[args.band])
        mag = flux_to_mag(flux, limFlux)
        errMag = flux_error_to_mag_error(errFlux, flux)

    # check the bias. Neale: should be added to the mean of the rational 
    # quadratic
    #kern = GPy.kern.Bias(1) + GPy.kern.RatQuad(1)

    kern = GPy.kern.RatQuad(1)


    if args.mag:
        mu, var, GPModel = gp_fit(
            phase, mag, errMag, 
            kern, n_restarts=10, test_length=args.testLength)
    else:    
        mu, var, GPModel = gp_fit(
            phase, flux, errFlux, 
            kern, n_restarts=10, test_length=args.testLength)
    
    print GPModel['.*lengthscale']
    print GPModel['.*power']

    if plt.get_fignums():
        figNum = plt.get_fignums()[-1]+1
    else:
        figNum = 1

    GPModel.plot_f(fignum=figNum)

    if args.mag:
        ylim = plt.ylim()
        plt.ylim(ylim[1], ylim[0])
        plt.errorbar(phase, mag, 
                     yerr=errMag, fmt=None, ecolor='black', zorder=1)
    else:
        plt.errorbar(phase, flux, 
                     yerr=errFlux, fmt=None, ecolor='black', zorder=1)

    plt.text(
        plt.xlim()[0], 
        plt.ylim()[1], "{:>06}".format(args.candidate))
    print "The process took {:5.3f} secs.".format(time.time()-start_time)
    plt.show()