"""
Implementation of general use functions.
"""
import numpy as np
import numpy.ma as ma
import pickle 
import gzip 
import GPy 
import classes 
import time 
import argparse
import matplotlib.pyplot as plt
from cStringIO import StringIO

# code from 
#
# http://stackoverflow.com/questions/16571150/
#        how-to-capture-stdout-output-from-a-python-function-call
import sys
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

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
        "-d", "--data-catalog", dest="catalog", default=None,
        help="SN catalog.")
    
    parser.add_argument(
        "-c", "--catalog-path", dest="path", 
        default="train_data/snCatalog.gz",
        help="Path to SN catalog.")

    parser.add_argument(
        "-t", "--test-lengthscale", dest="testLength",
        action="store_true", 
        help="Flag to randomize the lengthscale parameter.")

    parser.add_argument(
        "-m", "--magnitudes", dest="mag",
        action="store_true", help="Flag to switch to magnitudes.")

    parser.add_argument(
        "-p", "--prior-test", dest="testPrior",
        action="store_true", help="Flag to test prior on lengthscale")

    args = parser.parse_args()
else:
    # the file has been imported as a module
    pass

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


def gp_fit(
    X, Y, errY, kernel, n_restarts=0, 
    test_length=False,
    test_prior=False,
    verbose=False):
    """
    Performs gaussian process regression
    NOTE
    check on shape of input should be added
    """
 
    medXStep = np.median(np.abs(np.subtract(X[0:-2], X[1:-1])))
    maxXStep = X[-1] - X[0]
    if verbose:
        print "  Median X step {:<5.3f}".format(medXStep)
        print "  X range {:<5.3f}".format(maxXStep)

    rsX = reshape_for_GPy(X)
    rsY = reshape_for_GPy(Y)

    gpModel = GPy.models.GPHeteroscedasticRegression(rsX, rsY, kernel)
    gpModel['.*Gaussian_noise'] = errY

    # Block to capture unwanted output from .constrain_fixed() function
    with Capturing() as output:
        [gpModel['.*Gaussian_noise_%s' %i].constrain_fixed(warning=False) 
         for i in range(X.size)]

    if test_length:
        np.random.RandomState
        length = np.random.uniform(low=medXStep, high=maxXStep)
        if verbose:
            print "  Randomized lengthscale {:<5.3f}".format(length)
        with Capturing() as output:
            gpModel['.*lengthscale'].constrain_fixed(length, warning=False)
    elif test_prior:
        prior = GPy.core.parameterization.priors.Gamma(1, 20)
        gpModel['.*lengthscale'].set_prior(prior, warning=False)
    else:
        pass

    if n_restarts > 0:
        
        gpModel.optimize_restarts(num_restarts=n_restarts, 
                                  verbose=False,
                                  robust=True,
                                  messages=False)
    else:
        gpModel.optimize(optimizer='scg')

    predX = reshape_for_GPy(np.arange(X.min(), X.max(), 1))

    meanY, var = gpModel._raw_predict(predX, full_cov=True)
    return meanY, var, gpModel


def save_fit(
    candidate,
    fitX, 
    fit, errFit):
    pass
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
        '''
        The use of SN catalog as from Newling has been deprecated. Its 
        handling is too timecosuming. Getting data from SN file seems much
        faster
        '''

        if args.candidate is None:
            # Picking random candidate
            #
            # high set max number of SN in SNPhotCC 
            candidateIdx = np.random.random_integers(
                low=0, high=18321)
            print candidateIdx
            args.candidate = np.genfromtxt(
                dataSN+"DES_BLIND+HOSTZ.LIST", dtype=None)[candidateIdx]

            # Setting path and getting data
            pathToSN = dataSN + args.candidate
        else:
            pathToSN = dataSN + \
                        "DES_SN" + "{:>06}".format(args.candidate) + ".DAT"
        sn = get_sn_from_file(pathToSN)

        phase = sn.lightCurvesDict[args.band].mjd
        flux = sn.lightCurvesDict[args.band].flux
        errFlux = sn.lightCurvesDict[args.band].fluxErr
    else:
        pass
    
    print "  Candidate ID          {:<06}".format(args.candidate)
    print "  Testing lengthscale ? {:<5}".format(args.testLength)
    print "  Magnitudes ?          {:<5}".format(args.mag)
    
    # Tranforming to magnitudes
    if args.mag:
        limFlux = mag_to_flux(limMag[args.band])
        mag = flux_to_mag(flux, limFlux)
        errMag = flux_error_to_mag_error(errFlux, flux)

    # check the bias. Neale: should be added to the mean of the rational 
    # quadratic
    #kern = GPy.kern.Bias(1) + GPy.kern.RatQuad(1)

    kern = GPy.kern.RatQuad(1)
    # kern = GPy.kern.RBF(1)

    # Fitting the data points
    # 
    # TBD: the fit is OK if passes the model validation procedure (which has 
    # 
    # to be done)
    if sn.badCurve:
        if args.mag:
            mu, var, GPModel = gp_fit(
                phase, mag, errMag, 
                kern, n_restarts=10, 
                test_length=args.testLength, 
                test_prior=args.testPrior)
        else:
            mu, var, GPModel = gp_fit(
                phase, flux, errFlux, 
                kern, n_restarts=0, 
                test_length=args.testLength,
                test_prior=args.testPrior)

        
        print GPModel['.*lengthscale|.*power']
        
        print "  Model log likelihood = {: <6}".format(GPModel.log_likelihood())

        print "  Fit to data:", mu, "\n"
        print "  Normalised fit to data:", mu/mu.max(), "\n"
        # 
        # Plot
        # 
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
        print "  The process took {:5.3f} secs.".format(time.time()-start_time)
        plt.show()
    else:
        print "  This was flagged as BAD CURVE!"