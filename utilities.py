"""
Implementation of general use functions.
"""
import numpy as np
import numpy.ma as ma                     # Masked arrays
import pickle as pkl
import gzip as gzip
import GPy as GPy

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
    snIdx = np.random.random_integers(low=0, high=len(catalog.SNID))
    numObs = len(catalog.sne[snIdx].lightCurvesDict[band].mjd)

    phase = catalog.sne[snIdx].lightCurvesDict[band].mjd
    phase = phase - phase[catalog.sne[snIdx].lightCurvesDict['r'].flux.argmax()]

    flux = catalog.sne[snIdx].lightCurvesDict[band].flux
    
    errFlux = catalog.sne[snIdx].lightCurvesDict[band].fluxErr
    
    return phase, flux, errFlux, snIdx

def get_sn(catalog, band, idx):
    """
    Extract specified supernova observation in specified band from catalog.
    Returns time, flux, flux errors arrays.
    Time has zeropoint on the maximum flux in r band
    """

    numObs = len(catalog.sne[idx].lightCurvesDict[band].mjd)

    phase = catalog.sne[idx].lightCurvesDict[band].mjd
    phase = phase - phase[catalog.sne[idx].lightCurvesDict['r'].flux.argmax()]
    # t = t - np.min(t)

    flux = catalog.sne[idx].lightCurvesDict[band].flux
    
    errFlux = catalog.sne[idx].lightCurvesDict[band].fluxErr
    
    return phase, flux, errFlux

def reshape_for_GPy(vec):
    return np.reshape(vec, (len(vec), 1))


def gp_fit(phase, mag, errMag, kernel, n_restarts=0):
    """
    Performs gaussian process regression
    NOTE
    check on shap of input should be added
    """
    rsPhase = reshape_for_GPy(phase)
    rsMag = reshape_for_GPy(mag)

    gpModel = GPy.models.GPHeteroscedasticRegression(rsPhase, rsMag, kern)
    gpModel['.*Gaussian_noise'] = errMag
    [gpModel['.*Gaussian_noise_%s' %i].constrain_fixed() 
     for i in range(phase.size)]
    if n_restarts > 0:
        gpModel.optimize_restarts(num_restarts=n_restarts, 
                                  verbose=False,
                                  robust=True,
                                  messages=False)