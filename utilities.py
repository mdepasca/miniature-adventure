"""
Implementation of general use functions.
"""
import numpy as np
import numpy.ma as ma
import cPickle 
import gzip 
import GPy 
import classes 
import time
from os import path
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

class bcolors:
    txtblk='\033[0;30m' # Black - Regular
    txtred='\033[0;31m' # Red
    txtgrn='\033[0;32m' # Green
    txtylw='\033[0;33m' # Yellow
    txtblu='\033[0;34m' # Blue
    txtpur='\033[0;35m' # Purple
    txtcyn='\033[0;36m' # Cyan
    txtwht='\033[0;37m' # White

    bldblk='\033[1;30m' # Black - Bold
    bldred='\033[1;31m' # Red
    bldgrn='\033[1;32m' # Green
    bldylw='\033[1;33m' # Yellow
    bldblu='\033[1;34m' # Blue
    bldpur='\033[1;35m' # Purple
    bldcyn='\033[1;36m' # Cyan
    bldwht='\033[1;37m' # White

    unkblk='\033[4;30m' # Black - Underline
    undred='\033[4;31m' # Red
    undgrn='\033[4;32m' # Green
    undylw='\033[4;33m' # Yellow
    undblu='\033[4;34m' # Blue
    undpur='\033[4;35m' # Purple
    undcyn='\033[4;36m' # Cyan
    undwht='\033[4;37m' # White

    bakblk='\033[40m'   # Black - Background
    bakred='\033[41m'   # Red
    badgrn='\033[42m'   # Green
    bakylw='\033[43m'   # Yellow
    bakblu='\033[44m'   # Blue
    bakpur='\033[45m'   # Purple
    bakcyn='\033[46m'   # Cyan
    bakwht='\033[47m'   # White
    txtrst='\033[0m'    # Text Reset
    
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


def open_pkl(filePath):
    # print '>>> Loading pkl file from: ' + filePath
    # elTime = time.time()
    
    fileHandler = open(filePath,'r')
    data = cPickle.load(fileHandler)
    fileHandler.close()

    # elTime -= time.time()
    # print '>>> Done in ' + str(abs(elTime)) + ' secs'
    return data

def dump_pkl(filePath, dataStruct):
    # print '>>> Dumping data structure into: ' + filePath
    # ERROR_STATE = 0
    # whileOn = True
    # i = 0
    # while whileOn:
    #     if path.exists(filePath):
    #             i += 1
    #             pklIdx = filePath.rfind('.pkl')
    #             if i > 1: pklIdx -= 3
    #             filePath = filePath[0:pklIdx] + '({:<1}).pkl'.format(i)
    #     else:
    #         whileOn = False            

    fileHandler = open(filePath,'w')
        
    cPickle.dump(dataStruct, fileHandler)
    fileHandler.close()


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
        action="store_true", help="Flag to test prior on lengthscale.")

    parser.add_argument(
        "-v", "--verbose", dest="verbose",
        action="store_true", help="Enable task progress report.")

    parser.add_argument(
        "-k", "--kernel", dest="kern",
        default="RatQuad", help="Kernel to use in GP regression. Accepts: "+\
        "RatQuad (Rational Quadratic) or RBF (Radial Basis Functions).")
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

def rewrite_file(fileName):
    """Rewrites files produced after fit of old code version.
    It removes `#` at beginning of the first 10 lines, leaving as it is the
    first line.
    If necessary it adds `MJD_MAX_FLUX-CCF:      0.000`.
    Adds trivial column `OBS`.
    """

    inFile = file(fileName, 'r')
    lines = inFile.readlines()
    inFile.close()

    outFile = open(fileName, 'w')

    line = ''
    for i in range(len(lines)):
        if i == 0:
            outFile.write(lines[i])
            continue
        if i > 0 and i < 10:
            if i == 7 and ('REDSHIFT_SPEC:' not in lines[i]):
                line = 'REDSHIFT_SPEC: -9.0000 +- 9.0000\n'
                outFile.write(line)    
            
            line = lines[i]
            line = line[2:]
            outFile.write(line)
            continue

        if i == 10:
            if lines[i] != '\n': 
                outFile.write(lines[i])
            else:
                line = 'MJD_MAX_FLUX-CCF:      0.000\n'
                outFile.write(line)
                outFile.write(lines[i])
                print i, lines[i]
            continue

        # empty space between details and table
        if lines[i] == '\n' and i > 10:
            outFile.write(lines[i])
            print i, lines[i]
            continue

        if lines[i][0] == '#':
            outFile.write(lines[i])
            continue
        
        if 'MJD' in lines[i]:
            if 'FIELD' not in lines[i]:
                line = lines[i][0:15] + '   FIELD' + lines[i][15:]
            if 'OBS' not in lines[i]:
                line = lines[i]
                line =  ' OBS  ' + line
            outFile.write(line)
            continue

        if '----' in lines[i]:
            line = lines[i][0:15] + '  ------' + lines[i][15:]
            line = '----  ' + line
            outFile.write(line)
            continue


        line = lines[i][0:15] + '  NULL  ' + lines[i][15:]
        line = 'OBS:  ' + line
        outFile.write(line)

    outFile.close()






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
    X, Y, errY, kernel, n_restarts=0, parallel=True,
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
    # with Capturing() as output:
    [gpModel['.*Gaussian_noise_%s' %i].constrain_fixed(warning=False) 
         for i in range(X.size)
         ]

    if test_length:
        np.random.RandomState
        length = np.random.uniform(low=medXStep, high=maxXStep)
        if verbose:
            print "  Randomized lengthscale {:<5.3f}".format(length)
        # with Capturing() as output:
        gpModel['.*lengthscale'].constrain_fixed(length, warning=False)
    elif test_prior:
        prior = GPy.core.parameterization.priors.Gamma(1, 20)
        gpModel['.*lengthscale'].set_prior(prior, warning=False)
    else:
        pass

    if n_restarts > 0:
        # optimize_restart is from GPy/core/model.py
        gpModel.optimize_restarts(num_restarts=n_restarts,
                                    parallel=parallel,
                                    verbose=False,
                                    robust=True,
                                    messages=False)
    else:
        gpModel.optimize(optimizer='scg')

    predX = reshape_for_GPy(np.arange(X.min(), X.max(), 1.))

    # _raw_predict is from GPy/core/gp.py
    meanY, var = gpModel._raw_predict(predX, full_cov=False)
    return predX, meanY, var, gpModel

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
    
    print "  Candidate ID          {:>06}".format(args.candidate)
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

    
    if args.kern == "RatQuad":
        kern = GPy.kern.RatQuad(1)
    elif args.kern == "RBF":
        kern = GPy.kern.RBF(1)
    else:
        raise Exception

    # Fitting the data points
    # 
    # TBD: the fit is OK if passes the model validation procedure (which has 
    # 
    # to be done)
    if not sn.lightCurvesDict[args.band].badCurve:
        if args.mag:
            predPhase, mu, var, GPModel = gp_fit(
                                            phase, mag, errMag, 
                                            kern, n_restarts=10, 
                                            test_length=args.testLength, 
                                            test_prior=args.testPrior,
                                            verbose=args.verbose)
        else:
            predPhase, mu, var, GPModel = gp_fit(
                                            phase, flux, errFlux, 
                                            kern, n_restarts=10, 
                                            test_length=args.testLength,
                                            test_prior=args.testPrior,
                                            verbose=args.verbose)

        
        print GPModel['.*lengthscale|.*power']
        
        print "  Model log likelihood = {: <6}".format(GPModel.log_likelihood())

        # print "  Fit to data:", mu, "\n"
        # print "  Normalised fit to data:", mu/mu.max(), "\n"
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