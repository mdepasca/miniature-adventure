import numpy as np
import os
import sys
import glob
import cPickle
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn, vstack, hstack
import time
import argparse
import warnings
import subprocess

warnings.filterwarnings(
    'error', 
    message=".*divide by zero encountered in double_scalars.*", 
    category=RuntimeWarning
    )
from math import sqrt
from scipy import interpolate, polyfit

if __name__ == '__main__':
    import utilities as util
    import GPy
    from matplotlib import pyplot as plt

    

    parser = argparse.ArgumentParser(
        description="Test of general functions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    actionGroup = parser.add_argument_group('ACTION')
    inputGroup = parser.add_argument_group('INPUT')

    """

    ACTION OPTIONS
    
    """
    actionGroup.add_argument(
        "--k-correction", dest='kcor', 
        action='store_true', help='Switch on k correction.'
        )

    actionGroup.add_argument(
        '--distance', dest='distance',
        action='store_true', help='Calculate distance between fitted lightcurves \
        in same band.'
        )

    actionGroup.add_argument(
        '--test-prior', dest='testPrior',
        action='store_true', help='Test prior in GP regression.'
        )

    """

    INPUT OPTIONS
    
    """
    inputGroup.add_argument(
        "--data-directory", dest="dirData",
        default="train_data" + os.sep + "SIMGEN_PUBLIC_DES",
        help="Path to directory containing training data.")

    inputGroup.add_argument(
        "-b", "--band", dest="band", default='r', 
        help="Photometric band.")

    inputGroup.add_argument(
        "-c1", "--candidate1", dest="candidate1", 
        type=np.int32, default=None, 
        help="First candidate idx")

    inputGroup.add_argument(
        "-c2", "--candidate2", dest="candidate2", 
        type=np.int32, default=None, 
        help="Second candidate idx")    

    inputGroup.add_argument(
        "--mag", dest="mag",
        action="store_true",
        help="Reads in magnitudes from file."
        )

    args = parser.parse_args()



class photoBand():
    def __init__(self, name, effLambda):
        """
        Initialize the photometric filter with `name' and
        effective lambda, `effLambda', in nanometers.
        """
        self.name = name
        self.effLambda = effLambda

class LightCurve():
    """
    Once fully initiated, instances of this class have the following important 
    properties
    band            (string) 'g','r','i' or 'z'
    lim             (float) brightness threshold for the band (in correct units
                            flux or mag)
    mjd             (array) modified julian dates of observations
    flux            (array) the observed flux
    fluxErr         (array) the error in the observed flux
    shifted_mjd     (array) mjd shifted such that the peak has mjd = 0. To be 
                            modified
    magFlag         (boolean) if the flux is expressed in magnitudes
    """
    badCurve = False
    shifted_mjd = np.zeros(0)
    normFlux = list()
    normErr = list()
    """
    TRY TO USE __slots__
    """
    __slots__ = ['band', 'lim', 'mjd', 'shiftedMjd', 'flux', 'fluxErr', 'badCurve',
        'shifted_mjd', 'normFlux', 'normErr', 'magFlag']
    def __init__(self, band, magFlag=False, lim=0):
        self.band = band
        self.mjd = list()#np.zeros(0, dtype=float)
        self.shiftedMjd = list()#np.zeros(0, dtype=float)
        self.flux = list()#np.zeros(0, dtype=float)
        self.fluxErr = list()#np.zeros(0, dtype=float)
        self.magFlag = magFlag
        self.lim = lim

    @classmethod
    def data(band, mjd, flux, fluxErr):
        self.band = band
        self.mjd = mjd
        self.flux = flux
        self.fluxErr = fluxErr
        self.set_badCurve()

    def set_badCurve(self): 

        # this could be generalised as "if min(self.flux) <= self.lim:"
        if len(self.flux) == 0:
            self.badCurve = True
        elif self.magFlag:
            if min(self.flux) <= self.lim:
                self.badCurve = True
        elif max(self.flux) == 0:
                self.badCurve = True

    def set_shifted_mjd(self, distance):
        """
        Construct shiftedMjd, by subtracting 'distance' from 'self.flux'
        """
        self.shiftedMjd = [self.mjd[i]-distance for i in range(len(self.mjd))]
    
    @property
    def max_flux(self):
        if not self.badCurve:
            result = min(self.flux) if self.magFlag else max(self.flux)
        else:
            result = self.lim

        return result

    @property
    def max_error(self):
        if not self.badCurve:
            result = max(self.fluxErr)
        else:
            result = 0

        return result

    @property
    def max_flux_index(self):
        """
        Return the index of the maximum flux
        """
        # return np.argmax(self.flux)
        return self.flux.index(min(self.flux)) if self.magFlag else self.flux.index(max(self.flux))

    @property
    def size(self):
        # return self.mjd.size
        return len(self.mjd)

class Supernova():
    """
    Has the following properties

    g               (LightCurve)
    r               (LightCurve)
    i               (LightCurve)
    z               (LightCurve)
    lightCurvesDict (dictionary) 4 entries, g,r,i,z returning the corresponding 
                                 LightCurves
    SNID            (int)   supernova ID
    SNTypeInt       (int)   supernova type integer (see relationship between 
                                                    number and type)
    zSpec           (float) If known via spectroscope, otherwise None
    hostGalaxyID    (int)   THe host galaxy ID (all supernovae (in +zPhotHost) 
                                                catalog have this)
    zPhotHost       (float) The redshift of the host galaxy (all supernovae in 
                                                        the catalog have this)
    zPhotHostErr    (float) Error in zPhotHost
    """
    """
    TRY TO USE __slots__
    """
    def __init__(self, inFileName, magFlag=False):
        """
        Parses all the light curve data in inFileName into a Supernova object.
        """

        inFile = file(inFileName, "r")
        lines = inFile.readlines()
        inFile.close()

        self.g = LightCurve("g", magFlag=magFlag)
        self.r = LightCurve("r", magFlag=magFlag)
        self.i = LightCurve("i", magFlag=magFlag)
        self.z = LightCurve("z", magFlag=magFlag)

        self.lcsDict = {'g':self.g, 
                        'r':self.r, 
                        'i':self.i, 
                        'z':self.z}
        
        for line in lines:

            if len(line) > 3 and line[0] != "#":
                
                tag = line.split(":")[0]
                data = line.split(":")[-1].split()
                   
                if tag == "OBS":
                    mjd = float(data[0])
                    passband = data[1]
                    if magFlag:
                        flux = float(data[6])
                        fluxErr = float(data[7])
                    else:
                        flux = float(data[3])
                        fluxErr = float(data[4])
                    if fluxErr > 0:
                        if passband == "g":
                            self.g.mjd.append(mjd)
                            self.g.flux.append(flux)
                            self.g.fluxErr.append(fluxErr)
                            # self.g.add_data_point(mjd, flux, fluxErr)
                        elif passband == "r":
                            self.r.mjd.append(mjd)
                            self.r.flux.append(flux)
                            self.r.fluxErr.append(fluxErr)
                            # self.r.add_data_point(mjd, flux, fluxErr)
                        elif passband == "i":
                            self.i.mjd.append(mjd)
                            self.i.flux.append(flux)
                            self.i.fluxErr.append(fluxErr)
                            # self.i.add_data_point(mjd, flux, fluxErr)
                        elif passband == "z":
                            self.z.mjd.append(mjd)
                            self.z.flux.append(flux)
                            self.z.fluxErr.append(fluxErr)
                            # self.z.add_data_point(mjd, flux, fluxErr)
                        else:
                            print "Filter not recognized: {:<}".format(passband)
                elif tag == "SNID":
                    self.SNID = int(data[0])
                # SupernovaFit attribute
                elif tag == "GPKERNEL":
                    self.kern = data[0]
                elif tag == "SNTYPE":
                    self.SNTypeInt = int(data[0])                    
                elif tag == "RA":
                    self.RADeg = float(data[0])
                elif tag == "DECL":
                    self.decDeg = float(data[0])
                elif tag == "MWEBV":
                    self.MWEBV = float(data[0])
                elif (tag == "REDSHIFT_SPEC") or (tag == "REDSHIFT_FINAL"):
                    if float(data[0]) == -9:
                        self.zSpec = None
                    else:
                        self.zSpec = float(data[0])
                        self.zSpecErr = float(data[2])
                elif (tag == "HOST_GALAXY_GALID") or (tag == "HOSTGAL_OBJID"):
                    self.hostGalaxyID = int(data[0])
                elif (tag == "HOST_GALAXY_PHOTO-Z") or (tag == "HOSTGAL_PHOTOZ"):
                    self.zPhotHost = float(data[0])
                    self.zPhotHostErr = float(data[2])
                elif tag == "MJD_MAX_FLUX-CCF":
                    self.ccMjdMaxFlux = float(data[0])

        for b in self.lcsDict.keys():
            self.lcsDict[b].set_badCurve()

    def __cmp__(self, other):
        return 2*(self.zPhotHost - other.zPhotHost > 0) - 1 

    def k_correct_flux(
        self, nBands, 
        lambda_obs = [479.66, 638.26, 776.90, 910.82]
        ):

        if nBands == len(self.lcsDict.keys()):
            return 0, 0

        band_g = photoBand('g', 479.66)
        
        # defining redshift threshold for determination of rest frame wavelength
        if nBands == 3:
            zThr = lambda_obs[3]/lambda_obs[2] - 1
        if nBands == 2:
            zThr = lambda_obs[3]/lambda_obs[1] - 1
        if nBands == 3:
            zThr = lambda_obs[3]/lambda_obs[0] - 1

        # defining redshif between band effective lambdas
        z_bands = list()
        for i in range(len(self.lcsDict.keys())-1):
            z_bands.append(lambda_obs[i+1]/lambda_obs[i] - 1)

        # zed_gr = lambda_obs[1]/lambda_obs[0] - 1
        # zed_ri = lambda_obs[2]/lambda_obs[1] - 1
        # zed_iz = lambda_obs[3]/lambda_obs[2] - 1

        # defining rest frame effective lambdas for the light curve
        lambda_rf = [el/(
            1+(self.zSpec if self.zSpec else self.zPhotHost)
            ) for el in lambda_obs]

        # check for existence of observations in all the 4 lcs
        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                    raise TypeError("K-correction not possible: not all bands\
                    have observations")
        
        
        # --- debug    
        print len(self.g.mjd), \
            len(self.r.mjd), \
            len(self.i.mjd), \
            len(self.z.mjd)
        # ---

        # list of low resolution spectra (a list of lists). 
        #
        # Indexing: l[i][j], i -> list; j -> elem in i^th list
        k_corr = list()
        for i in range(nBands):
            k_corr.append(list())
        # k_corr_g = list()
        # k_corr_r = list()
        # k_corr_i = list()

        # intersection of epochs in 4 bands
        # the construction of a list without knowing which filters are being 
        #
        # is difficult here. In the next I've to know the name of the band.

        # int_mjd = list()
        # for b in self.lcsDict.key():

        int_mjd_g = [int(round(el)) for el in self.g.mjd]
        int_mjd_r = [int(round(el)) for el in self.r.mjd]
        int_mjd_i = [int(round(el)) for el in self.i.mjd]
        int_mjd_z = [int(round(el)) for el in self.z.mjd]

        mjd = [el for el in int_mjd_g if el in int_mjd_r 
                and el in int_mjd_i and el in int_mjd_z]

        """
        interpolating between points as y = ax + b
        """
        for jd in mjd:
            if (self.zSpec < zThr) or (self.zPhotHost < zThr):
                """
                g r i band can be k-corrected
                """
                a = (
                    self.r.flux[int_mjd_r.index(jd)] - \
                    self.g.flux[int_mjd_g.index(jd)]
                    )/(lambda_rf[1]-lambda_rf[0])
                b = self.g.flux[int_mjd_g.index(jd)] - \
                    a*lambda_rf[0]

                a_err = sqrt(self.r.fluxErr[int_mjd_r.index(jd)]**2+\
                    self.g.fluxErr[int_mjd_g.index(jd)]**2)
                b_err = sqrt(self.g.flux[int_mjd_g.index(jd)]**2+\
                    a_err**2)

                k_corr[0].append(a*lambda_obs[0]+b)

                a = (
                    self.i.flux[int_mjd_i.index(jd)] - \
                    self.r.flux[int_mjd_r.index(jd)]
                    )/(lambda_rf[2]-lambda_rf[1])
                b = self.r.flux[int_mjd_r.index(jd)] - \
                    a*lambda_rf[1]
                a_err = sqrt(self.i.fluxErr[int_mjd_i.index(jd)]**2+\
                    self.r.fluxErr[int_mjd_r.index(jd)]**2)
                b_err = sqrt(self.r.flux[int_mjd_r.index(jd)]**2+\
                    a_err**2)

                k_corr[1].append(a*lambda_obs[1]+b)

                a = (
                    self.z.flux[int_mjd_z.index(jd)] - \
                    self.i.flux[int_mjd_i.index(jd)]
                    )/(lambda_rf[3]-lambda_rf[2])
                b = self.i.flux[int_mjd_i.index(jd)] - \
                    a*lambda_rf[2]
                a_err = sqrt(self.i.fluxErr[int_mjd_i.index(jd)]**2+\
                    self.z.fluxErr[int_mjd_z.index(jd)]**2)
                b_err = sqrt(self.i.flux[int_mjd_i.index(jd)]**2+\
                    a_err**2)

                k_corr[2].append(a*lambda_obs[2]+b)
    
                # k_corr_g.append(a_gr*lambda_obs[0]+b_gr)
                # k_corr_r.append(a_ri*lambda_obs[1]+b_ri)
                # k_corr_i.append(a_iz*lambda_obs[2]+b_iz)
                # CALCULATE K-CORRECTION ERROR
                continue

            if (self.zSpec < zThr) or (self.zPhotHost < zThr):
                """
                g r band can be k-corrected
                """
                # choosing which interpolation to use
                if (self.zSpec < 0.21) or (self.zPhotHost < 0.21):
                    a_gr = (
                        self.r.flux[int_mjd_r.index(jd)] - \
                        self.g.flux[int_mjd_g.index(jd)]
                        )/(lambda_rf[1]-lambda_rf[0])
                    b_gr = self.g.flux[int_mjd_g.index(jd)] - \
                        a_gr*lambda_rf[0]
                    a_ri = (
                        self.i.flux[int_mjd_i.index(jd)] - \
                        self.r.flux[int_mjd_r.index(jd)]
                        )/(lambda_rf[2]-lambda_rf[1])
                    b_ri = self.r.flux[int_mjd_r.index(jd)] - \
                        a_ri*lambda_rf[1]
                
                    # CALCULATE K-CORRECTION
                    k_corr_g.append(a_gr*lambda_obs[0]+b_gr)
                    k_corr_r.append(a_ri*lambda_obs[1]+b_ri)
                    # CALCULATE K-CORRECTION ERROR
                else:
                    a_ri = (
                        self.i.flux[int_mjd_i.index(jd)] - \
                        self.r.flux[int_mjd_r.index(jd)]
                        )/(lambda_rf[2]-lambda_rf[1])
                    b_ri = self.r.flux[int_mjd_r.index(jd)] - \
                        a_ri*lambda_rf[1]
                    a_iz = (
                        self.z.flux[int_mjd_z.index(jd)] - \
                        self.i.flux[int_mjd_i.index(jd)]
                        )/(lambda_rf[3]-lambda_rf[2])
                    b_iz = self.i.flux[int_mjd_i.index(jd)] - \
                        a_iz*lambda_rf[2]

                    # CALCULATE K-CORRECTION
                    k_corr_g.append(a_ri*lambda_obs[0]+b_ri)
                    k_corr_r.append(a_iz*lambda_obs[1]+b_iz)
                    # CALCULATE K-CORRECTION ERROR
                continue

            if (self.zSpec < zThr) or (self.zPhotHost < zThr):
                """
                g band can be k-corrected.
                This will exclude the use of color information in 
                calculating distances.
                """
                if (self.zSpec < 0.33) or (self.zPhotHost < 0.33):
                    a_gr = (
                        self.r.flux[int_mjd_r.index(jd)] - \
                        self.g.flux[int_mjd_g.index(jd)]
                        )/(lambda_rf[1]-lambda_rf[0])
                    b_gr = self.g.flux[int_mjd_g.index(jd)] - \
                        a_gr*lambda_rf[0]

                    # CALCULATE K-CORRECTION
                    k_corr_g.append(a_gr*lambda_obs[0]+b_gr)
                    # CALCULATE K-CORRECTION ERROR
                elif (self.zSpec < 0.61) or (self.zPhotHost < 0.61):
                    a_ri = (
                        self.i.flux[int_mjd_i.index(jd)] - \
                        self.r.flux[int_mjd_r.index(jd)]
                        )/(lambda_rf[2]-lambda_rf[1])
                    b_ri = self.r.flux[int_mjd_r.index(jd)] - \
                        a_ri*lambda_rf[1]

                    # CALCULATE K-CORRECTION
                    k_corr_g.append(a_ri*lambda_obs[0]+b_ri)
                    # CALCULATE K-CORRECTION ERROR
                else:
                    a_iz = (
                        self.z.flux[int_mjd_z.index(jd)] - \
                        self.i.flux[int_mjd_i.index(jd)]
                        )/(lambda_rf[3]-lambda_rf[2])
                    b_iz = self.i.flux[int_mjd_i.index(jd)] - \
                        a_iz*lambda_rf[2]

                    # CALCULATE K-CORRECTION
                    k_corr_g.append(a_iz*lambda_obs[0]+b_iz)
                    # CALCULATE K-CORRECTION ERROR
                continue

            raise ValueError('delta redshift too big to have a k-correction')

        return mjd, list([k_corr_g, k_corr_r, k_corr_i])

class SupernovaFit():
    ccMjdMaxFlux = 0

    def __init__(self, supernova, GPkern=''):
        # check on supernova type HAS TO BE ADDED
        self.g = LightCurve("g")
        self.r = LightCurve("r")
        self.i = LightCurve("i")
        self.z = LightCurve("z")
        self.lcsDict = {"g":self.g, 
                        "r":self.r, 
                        "i":self.i, 
                        "z":self.z}
        self.peaked = False
        
        if hasattr(supernova, 'kern'):
            self.kern = supernova.kern
        elif GPkern != '':
            self.kern = GPkern

        if hasattr(supernova, 'ccMjdMaxFlux'):
            self.ccMjdMaxFlux = supernova.ccMjdMaxFlux
        self.SNID = supernova.SNID
        self.SNTypeInt = supernova.SNTypeInt
        self.RADeg = supernova.RADeg 
        self.decDeg = supernova.decDeg 
        self.MWEBV = supernova.MWEBV 
        if hasattr(supernova, 'zSpec'):
            self.zSpec = supernova.zSpec
        if hasattr(supernova, 'zSpecErr'):
            self.zSpecErr = supernova.zSpecErr
        else:
            self.zSpecErr = None
        self.hostGalaxyID = supernova.hostGalaxyID 
        self.zPhotHost = supernova.zPhotHost 
        self.zPhotHostErr = supernova.zPhotHostErr 

    def set_lightcurve(self, band, mjd, flux, fluxErr, magFlag=False):
        self.lcsDict[band].mjd = mjd
        self.lcsDict[band].flux = flux
        self.lcsDict[band].fluxErr = fluxErr
        self.lcsDict[band].magFlag = magFlag

        self.lcsDict[band].set_badCurve()

        if not self.lcsDict[band].badCurve:
            if (band == 'r') \
                and self.r.max_flux_index not in set([0, self.r.size-1]):
                    self.peaked = True


    def shift_mjds(self):
        """ Shifts mjd attribute of each lc according to flux maximum 
        in r band for peaked lcsself.
        """
        
        mjdrMax = self.r.mjd[self.r.max_flux_index]
        
        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                continue
            self.lcsDict[b].set_shifted_mjd(mjdrMax)


    def normalized_flux(self, band):
        """Normalizes the light curve in `band` using the sum of the maximums
        in all that band.
        """

        if not self.lcsDict[band].normFlux:
            den = self.g.max_flux + \
                self.r.max_flux + \
                self.i.max_flux + \
                self.z.max_flux

            flux = self.lcsDict[band].flux
        # result = [flux[idx]/den for idx in range(len(flux))]
            self.lcsDict[band].normFlux = [
                flux[idx]/den for idx in range(len(flux))
                ]
        # return result
        return self.lcsDict[band].normFlux

    def normalized_error(self, band):
        """Normalizes the light curve in band b using the maximum in that band.
        s is a slice on the array.
        """

        if not self.lcsDict[band].normErr:
            den = self.g.max_flux + \
                self.r.max_flux + \
                self.i.max_flux + \
                self.z.max_flux

            fluxErr = self.lcsDict[band].fluxErr
        # result = [fluxErr[idx]/den for idx in range(len(fluxErr))]
            self.lcsDict[band].normErr = [
                fluxErr[idx]/den for idx in range(len(fluxErr))
                ]
        # return result
        return self.lcsDict[band].normErr

    def set_peaked(self):
        self.peaked = True

    @property
    def peaked(self):
        return self.peaked

    def get_distance(self, candidate, band):
        """Calculate difference (aka distance) between two 
        interpolated light curves. 
        """
        if type(band) is not str:
            raise TypeError("variable `band` is not of type string")
        distFlag = 5
        
        sizeSelf = self.lcsDict[band].size
        sizeCandidate = candidate.lcsDict[band].size


        if sizeSelf >= sizeCandidate:
            mjd1 = [round(el) for el in self.lcsDict[band].shiftedMjd]
            flux1 = self.normalized_flux(band)
            fluxErr1 = self.normalized_error(band)
            mjd2 = [round(el) for el in candidate.lcsDict[band].shiftedMjd]
            flux2 = candidate.normalized_flux(band)
            fluxErr2 = candidate.normalized_error(band)
        else:
            mjd1 = [round(el) for el in candidate.lcsDict[band].shiftedMjd]
            flux1 = candidate.normalized_flux(band)
            fluxErr1 = candidate.normalized_error(band)
            mjd2 = [round(el) for el in self.lcsDict[band].shiftedMjd]
            flux2 = self.normalized_flux(band)
            fluxErr2 = self.normalized_error(band)

        mjdIntersection = [el for el in mjd1 if el in mjd2]

        if len(mjdIntersection) < 2:
            distance = distFlag
        else:
            flux1Int = [
                flux1[i] for i in [mjd1.index(el) for el in mjdIntersection]
                ]


            flux2Int = [
                flux2[i] for i in [mjd2.index(el) for el in mjdIntersection]
                ]


            fluxErr1Int = [
                fluxErr1[i] for i in [mjd1.index(el) for el in mjdIntersection]
                ]


            fluxErr2Int = [
                fluxErr2[i] for i in [mjd2.index(el) for el in mjdIntersection]
                ]


            num = [
                (el)**2 for el in [
                    flux1Int[i] - flux2Int[i] for i in range(
                        len(mjdIntersection)
                        )
                    ]
                ]

            den = [(fluxErr1Int[i])**2 + (fluxErr2Int[i])**2 for i in range(
                    len(mjdIntersection)
                    )
                ]

            try:
                distance = sqrt(
                    sum([r for r in [num[i]/den[i] for i in range(
                            len(mjdIntersection)
                            )]]
                    )
                    )/(max(mjdIntersection) - min(mjdIntersection))

            except RuntimeWarning:
                print "selfID: {:<d} -- CandidateID {:<d}".format(self.SNID, 
                    candidate.SNID)
                print "1: {:<d}".format(id1)
                print "len(num) {:<d}".format(len(num))
                print "len(den) {:<d}".format(len(den))
                print den.index(0)
                print fluxErr1Int[den.index(0)], fluxErr2Int[den.index(0)]
                print fluxErr1.index(0), fluxErr2.index(0)
                print "len(mjdIntersection) {:<d}".format(len(mjdIntersection))
                print distance
                print '--------------------------------------------------------'

        return distance


    def save_on_txt(self, fileName, survey="DES"):
        t = Table(masked=True)
        # OBS is used to reproduce original SNPhotCC files can be deleted 
        #
        # provided to change init method of Supernova class
        colNames = [
                    ["OBS", "{:4s}"],
                    ["MJD", "{0:9.3f}"], ["BAND", "{:s}"], ["FIELD", "{:6s}"]
                    ]
        if self.r.magFlag:
            colNames.extend([["MAG", "{0:7.3f}"], ["MAG_ERR", "{0:7.3f}"]])
        else:
            colNames.extend([["FLUX", "{0:10.5f}"], ["FLUX_ERR", "{0:10.5f}"]])

        bandArr = list()
        mjd = list()
        flux = list()
        fluxErr = list()
        mjdArgsort = list()

        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                continue

            if bandArr.size == 0:
                bandArr = [b]*len(self.lcsDict[b].mjd)
                # bandArr = np.empty(len(self.lcsDict[b].mjd), dtype=np.str)
                # bandArr[:] = b
                mjd = self.lcsDict[b].mjd
                flux = self.lcsDict[b].flux
                fluxErr = self.lcsDict[b].fluxErr
            else:
                tmp = [b]*len(self.lcsDict[b].mjd)
                # tmp = np.empty(len(self.lcsDict[b].mjd), dtype=np.str)
                # tmp[:] = b
                bandArr.extend(tmp)
                # bandArr = np.concatenate((bandArr, tmp))
                mjd.extend(self.lcsDict[b].mjd)
                flux.extend(self.lcsDict[b].flux)
                fluxErr.extend(self.lcsDict[b].fluxErr)
                # mjd = np.concatenate((mjd, self.lcsDict[b].mjd))
                # flux = np.concatenate((flux, self.lcsDict[b].flux))
                # fluxErr = np.concatenate((fluxErr, self.lcsDict[b].fluxErr))

        mjdArgsort = np.argsort(mjd)
        mjd = mjd[mjdArgsort]
        flux = flux[mjdArgsort]
        fluxErr = fluxErr[mjdArgsort]
        bandArr = bandArr[mjdArgsort]

        """
        Adding and setting column to Table
        """

        for c in range(len(colNames)):
            if colNames[c][0] == "BAND" \
            or colNames[c][0] == "OBS" \
            or colNames[c][0] == "FIELD":
                col = MaskedColumn(np.zeros(mjd.size),
                    name=colNames[c][0],
                    format=colNames[c][1],
                    dtype=np.str, fill_value='-',
                    mask=np.zeros(mjd.size)
                    )
            else:
                col = MaskedColumn(np.zeros(mjd.size),
                    name=colNames[c][0],
                    format=colNames[c][1],
                    dtype=np.float, fill_value=-9,
                    mask=np.zeros(mjd.size))
            t.add_column(col)

        """
        Initializing columns
        """

        t["OBS"] = np.empty(mjd.size, dtype=np.str)
        t["OBS"][:] = "OBS:"
        t["FIELD"] = np.empty(mjd.size, dtype=np.str)
        t["FIELD"][:] = "NULL"
        t["MJD"] = mjd
        t["BAND"] = bandArr
        t["FLUX"] = flux
        t["FLUX_ERR"] = fluxErr

        t.filled()

        fOut = open(fileName, 'w')
        fOut.write("# File produced by Miniature Adventure on " + \
            "{:<02d}/{:<02d}/{:<4d} at {:<02d}:{:<02d}:{:<02d} GMT\n".format(
                time.gmtime().tm_mday, time.gmtime().tm_mon, 
                time.gmtime().tm_year,
                time.gmtime().tm_hour, time.gmtime().tm_min, 
                time.gmtime().tm_sec))
        fOut.write("SURVEY:  {:<}\n".format(survey))
        fOut.write("SNID:  {:<d}\n".format(self.SNID))
        fOut.write("GPKERNEL: {:<}\n".format(self.kern))
        # if self.SNTypeInt :
        fOut.write("SNTYPE: {:>d}\n".format(self.SNTypeInt))
        # if self.RADeg :
        fOut.write("RA:     {:>9.6f} deg\n".format(self.RADeg))
        # if self.decDeg :
        fOut.write("DECL:   {:>9.6f} deg\n".format(self.decDeg))
        # if self.MWEBV :
        fOut.write("MWEBV:  {:>6.4f}\n".format(self.MWEBV))
        if hasattr(self, "zSpec"):
            if self.zSpec:
                fOut.write("REDSHIFT_SPEC:  {:>6.4f} +- {:>6.4f}\n".format(
                    self.zSpec, self.zSpecErr
                    ))
            else:
                fOut.write("REDSHIFT_SPEC: -9.0000 +- 9.0000\n")
        if hasattr(self, "hostGalaxyID"):
            fOut.write("HOST_GALAXY_GALID: {:>d}\n".format(self.hostGalaxyID))
        if hasattr(self, "zPhotHost"):
            fOut.write("HOST_GALAXY_PHOTO-Z:  {:>6.4f} +- {:>6.4f}\n".format(
                self.zPhotHost, self.zPhotHostErr
                ))
        # if self.ccMjdMaxFlux != 0:
        fOut.write("MJD_MAX_FLUX-CCF:  {:>9.3f}\n".format(self.ccMjdMaxFlux))
        fOut.write("\n\n\n")
        fOut.write("# ======================================\n")
        fOut.write("# LIGHT CURVE FIT USING GAUSSIAN PROCESS\n")
        fOut.write("#\n")
        fOut.write("# NOBS: {:<}\n".format(mjd.size))

        ascii.write(t, output=fOut, delimiter='  ', 
            format='fixed_width_two_line')
        
        fOut.close()

if __name__ == '__main__':
    indent = "          "
    lambda_obs = [479.66, 638.26, 776.90, 910.82]
    limMagDict = {
        'g': 25.2,
        'r': 25.4,
        'i': 25.1,
        'z': 24.9
    }
    p = subprocess.Popen("ls *SN*.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirData+os.sep)
    lsDirData = p.stdout.read()
    lsDirData = lsDirData.split('\n')
    lsDirData.sort()
    lsDirData.remove('')
    # fCandidatesList = "DES_BLIND+HOSTZ.LIST"
    # candidatesFileList = np.genfromtxt(dirData+os.sep+fCandidatesList, dtype=None)\

    """

    KERNEL SPECIFICATION

    """
    # kern = GPy.kern.RBF(1)
    kern = GPy.kern.RatQuad(1)

    """
    ---------------------- 
    """

    if args.band not in ['g', 'r', 'i', 'z']:
        print 'Band {:<} not recognised! Changing to r'.format(args.band)
        args.band = 'r'

    if args.candidate1 is None:
        args.candidate1 = np.random.random_integers(
                low=0, high=len(lsDirData))

    if args.candidate2 is None:
        args.candidate2 = np.random.random_integers(
                low=0, high=len(lsDirData))

    while args.candidate2 == args.candidate1:
        args.candidate2 = np.random.random_integers(
            low=0, high=len(lsDirData))

    print args.candidate1
    print args.candidate2

    candidates = list()
    fit = list()

    """
    Getting observation's data
    """

    candidates.append(Supernova(
        args.dirData+os.sep+lsDirData[args.candidate1], args.mag))

    candidates.append(Supernova(
        args.dirData+os.sep+lsDirData[args.candidate2], args.mag))


    for candidate in candidates:
        """
        Setting limits in magnitudes
        """
        for b in candidate.lcsDict.keys():
            if args.mag:
                candidate.lcsDict[b].lim = limMagDict[b]
        
    
        print 'candidate z = {:>6.4f}'.format(
            candidate.zSpec if candidate.zSpec else candidate.zPhotHost)

        if args.kcor:
            lambda_rf = [el/(
                1+(candidate.zSpec if candidate.zSpec else candidate.zPhotHost)
                ) for el in lambda_obs]

            plt.figure()
            int_mjd_g = [int(el) for el in candidate.g.mjd]
            int_mjd_r = [int(el) for el in candidate.r.mjd]
            int_mjd_i = [int(el) for el in candidate.i.mjd]
            int_mjd_z = [int(el) for el in candidate.z.mjd]
            mjd = [el for el in int_mjd_g if el in int_mjd_r 
                        and el in int_mjd_i and el in int_mjd_z]
            jd = mjd[0]

            flux = [
                candidate.g.flux[int_mjd_g.index(jd)], 
                candidate.r.flux[int_mjd_r.index(jd)],
                candidate.i.flux[int_mjd_i.index(jd)],
                candidate.z.flux[int_mjd_z.index(jd)]
                ]
            fluxErr = [
                candidate.g.fluxErr[int_mjd_g.index(jd)], 
                candidate.r.fluxErr[int_mjd_r.index(jd)],
                candidate.i.fluxErr[int_mjd_i.index(jd)],
                candidate.z.fluxErr[int_mjd_z.index(jd)]
                ]
            plt.errorbar(lambda_obs, flux, yerr=fluxErr, fmt='k--', ecolor='black')
            # if args.mag:
            #     plt.xlim([plt.ylim()[1]], plt.ylim()[0]])
            plt.scatter(lambda_obs, flux, color='black')
            plt.errorbar(lambda_rf, flux, yerr=fluxErr, fmt='b--', ecolor='blue')
            plt.scatter(lambda_rf, flux, color='blue')
            plt.show()
            # mjd_k_corr, k_correct_flux = candidate.k_correct_flux()

            
            # (a, b) = np.polyfit(
            #     # [lambda_rf[0], lambda_rf[1]], 
            #     [
            #     lambda_rf[0],
            #     lambda_rf[1]#,
            #     # lambda_rf[2],
            #     # lambda_rf[3],
            #     ],
            #     [
            #     candidate.g.flux[int_mjd_g.index(jd)], 
            #     candidate.r.flux[int_mjd_r.index(jd)]#,
            #     # candidate.i.flux[int_mjd_i.index(jd)],
            #     # candidate.z.flux[int_mjd_z.index(jd)]
            #     ], deg = 1, 
            #     w = [
            #         1/candidate.g.fluxErr[int_mjd_g.index(jd)], 
            #         1/candidate.r.fluxErr[int_mjd_r.index(jd)]]#, 
            #     # cov=True
            #     )
            # raise ValueError

        """
        Create SupernovaFit objects
        """
        candidateFit = SupernovaFit(candidate, kern.name)

        """
        Looping in bands and fit of flux
        """
        for b in candidate.lcsDict.keys():

            phase = candidate.lcsDict[b].mjd
            flux = candidate.lcsDict[b].flux

            """
            Clipping to limiting magnitudes when flux is above 90th mag
            """
            if args.mag :
                flux = [candidate.lcsDict[b].lim if \
                        (el>90) else el for el in flux]

            errFlux = candidate.lcsDict[b].fluxErr

            # Fitting Lightcurve
            if (not candidate.lcsDict[b].badCurve) and (len(flux) >= 3):
                
                predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                phase, flux, errFlux, 
                                                kern, n_restarts=10, 
                                                parallel=False,
                                                test_length=False,
                                                test_prior=args.testPrior)
                print GPModel

                candidateFit.set_lightcurve(
                    b, predMjd, predFlux, predErr, args.mag
                    )
                
                print indent + \
                    "{:<} {:<}".format(candidate.SNID, b)
            else:
                candidateFit.lcsDict[b].badCurve = True
                print indent + util.bcolors.FAIL + \
                    "{:<} {:<}".format(candidate.SNID, b) + \
                    util.bcolors.ENDC

        candidateFit.shift_mjds()
        fit.append(candidateFit)
        

    if args.distance:
        if fit[0].peaked and fit[1].peaked:
            print 'Distance between the 2 normalized lcs in ' + \
            '{:<} band = {:<2.4f}'.format(args.band,
                fit[0].get_distance(fit[1], 
                args.band))

            # if plt.get_fignums():
            #     figNum = plt.get_fignums()[-1]+1
            # else:
            #     figNum = 1

            # plt.figure(figNum)
        else:
            print 'One of the 2 candidate has not r-band peak: '
            print 'Candidate 1 - {:<}'.format(fit[0].peaked)
            print 'Candidate 2 - {:<}'.format(fit[1].peaked)

    nrows = 2
    ncols = 4 

    colorList = ['blue', 'orange']

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, 
                    figsize=(16.5, 11.7), 
                    tight_layout=False
                    )

    axDict = {
    'g':[ax[0,0], ax[0,2]],
    'r':[ax[0,1], ax[0,3]],
    'i':[ax[1,0], ax[1,2]],
    'z':[ax[1,1], ax[1,3]]
    }

    fig.subplots_adjust(left=0.05, right=0.97, top=0.94, wspace=0.29)
    fig.suptitle(eval('\'Band \' + args.band + (\' with Prior Test\' ' + 
        'if args.testPrior else \'Band \' + args.band + \' No prior\') + \' -- \'' + 
        '+ \'Kernel: \' + kern.name'))

    for j in [0,1]:
        for b in axDict.keys():
            if args.mag:
                upperlimits = [
                    0 if el < 90 else 1 for el in candidates[j].lcsDict[b].fluxErr
                ]
                axDict[b][j].set_ylim(candidates[j].lcsDict[b].lim+2, 22)
                lowerlimits = False
            else:
                upperlimits = False
                lowerlimits = [0 if el > 0 else 1 for el in candidates[j].lcsDict[b].flux]
            axDict[b][j].plot([min(candidates[j].lcsDict[b].mjd), 
                max(candidates[j].lcsDict[b].mjd)], 
                [candidates[j].lcsDict[b].lim]*2, 
                c='k')
            
            fluxUpLim = [val for val in [
                        fit[j].lcsDict[b].flux[i] + 
                        2*fit[j].lcsDict[b].fluxErr[i] 
                            for i in range(len(fit[j].lcsDict[b].flux))
                        ]]
            fluxLowLim = [val for val in [
                fit[j].lcsDict[b].flux[i] - 
                2*fit[j].lcsDict[b].fluxErr[i] 
                    for i in range(len(fit[j].lcsDict[b].flux))
                        ]]
            axDict[b][j].fill_between(fit[j].lcsDict[b].mjd, 
                fluxUpLim, fluxLowLim, 
                alpha=0.2, linewidth=0.5)
            
            axDict[b][j].plot(
                fit[j].lcsDict[b].mjd,
                fit[j].lcsDict[b].flux, c=colorList[j],
                )
            axDict[b][j].errorbar(
                candidates[j].lcsDict[b].mjd,
                candidates[j].lcsDict[b].flux,
                candidates[j].lcsDict[b].fluxErr,
                uplims=upperlimits, lolims=lowerlimits, ecolor=colorList[j],
                fmt=None
                )
            axDict[b][j].scatter(
                candidates[j].lcsDict[b].mjd,
                candidates[j].lcsDict[b].flux,
                c=colorList[j],
                label='Band {:>s} | IDX {:>5d} | SNID {:>5d}'.format(b,
                    eval('args.candidate1 if j == 0 else args.candidate2'), 
                    candidates[j].SNID)
                )
            axDict[b][j].legend(
                loc='best', framealpha=0.3, fontsize='10'
                )
            
            axDict[b][j].set_xlabel('epoch [MJD]')
            if args.mag:
                axDict[b][j].set_ylabel('flux [mag]')
            else:
                axDict[b][j].set_ylabel('flux [adu]')

    # if args.mag:
    #     upperlimits = [
    #         0 if el < 90 else 1 for el in candidates[1].lcsDict[args.band].fluxErr
    #     ]
    #     ax[1].set_ylim(candidates[1].lcsDict[args.band].lim+2, 22)
        
    #     lowerlimits = False
    # else:
    #     lowerlimits = [0 if el > 0 else 1 for el in candidates[1].lcsDict[args.band].flux]
    #     upperlimits = False

    # ax[1].plot([min(candidates[1].lcsDict[args.band].mjd), 
    #         max(candidates[1].lcsDict[args.band].mjd)], 
    #         [candidates[1].lcsDict[args.band].lim]*2, 
    #         c='k')

    # fluxUpLim = [val for val in [
    #             fit[1].lcsDict[args.band].flux[i] + 
    #             2*fit[1].lcsDict[args.band].fluxErr[i] 
    #                 for i in range(len(fit[1].lcsDict[args.band].flux))
    #             ]]
    # fluxLowLim = [val for val in [
    #     fit[1].lcsDict[args.band].flux[i] - 
    #     2*fit[1].lcsDict[args.band].fluxErr[i] 
    #         for i in range(len(fit[1].lcsDict[args.band].flux))
    #             ]]
    # ax[1].fill_between(fit[1].lcsDict[args.band].mjd, 
    #     fluxUpLim, fluxLowLim, 
    #     alpha=0.2, linewidth=0.5)
    # ax[1].plot(
    #         fit[1].lcsDict[args.band].mjd,
    #         fit[1].lcsDict[args.band].flux, c='orange'
    #         )
    # ax[1].errorbar(
    #     candidates[1].lcsDict[args.band].mjd,
    #     candidates[1].lcsDict[args.band].flux,
    #     candidates[1].lcsDict[args.band].fluxErr,
    #     uplims=upperlimits, lolims=lowerlimits, ecolor='orange',
    #     fmt=None
    #     )
    # ax[1].scatter(
    #     candidates[1].lcsDict[args.band].mjd,
    #     candidates[1].lcsDict[args.band].flux, c='orange',
    #     label='IDX {:>5d} | SNID {:>5d}'.format(args.candidate2, 
    #         candidates[1].SNID)
    #     )
    # ax[1].legend(
    #     loc='best', framealpha=0.3, fontsize='10'
    #     )

    # ax[1].set_xlabel('epoch [MJD]')
    # if args.mag:
    #     ax[1].set_ylabel('flux [mag]')
    # else:
    #     ax[1].set_ylabel('flux [adu]')
    plt.show()
    # return catalogFit