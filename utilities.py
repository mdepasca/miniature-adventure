"""
Implementation of general use functions.
"""
import numpy as np
import numpy.ma as ma                     # Masked arrays
import pickle as pkl
import gzip as gzip

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
    Returns time, flux and flux errors arrays.
    The output is formatted for GPy regression.
    """
    snIdx = np.random.random_integers(low=0, high=len(catalog.SNID))
    numObs = len(snCatalog.sne[snIdx].lightCurvesDict[band].mjd)

    t = np.reshape(snCatalog.sne[snIdx].lightCurvesDict[band].mjd, (numObs, 1))
    t = t - np.min(t)

    flux = np.reshape(snCatalog.sne[snIdx].lightCurvesDict[band].flux, 
                      (numObs, 1))
    errFlux = np.reshape(snCatalog.sne[snIdx].lightCurvesDict[band].fluxErr, 
                      (numObs, 1))
    
    return t, flux, errFlux