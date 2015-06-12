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
    import cProfile


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

    actionGroup.add_argument(
        '--plot', dest='plot',
        action='store_true', help='Plot results of test.'
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

    def set_shifted_mjd_2(self, supernovaFitLC, max_flux_epoch, redshift=-1, ccMjdMaxFlux=0):
        """Initilise `shiftedMjd` on the basis of the fitted light curve
        `SupernovaFitLC`. If redshift is provided shiftedMjd will be also
        corrected from time dilation.
        """

        if supernovaFitLC.band != self.band:
            raise TypeError('Light curves of different bands.')

        if self.badCurve:
            raise TypeError('Operating on a bad light curve.')

        if supernovaFitLC.badCurve:
            raise TypeError('Bad input light curve.')

        if redshift != -1:
            self.shiftedMjd = self.calc_destretched_time(redshift)

        print supernovaFitLC.mjd[supernovaFitLC.max_flux_index]
        self.shiftedMjd = [el - max_flux_epoch# supernovaFitLC.mjd[supernovaFitLC.max_flux_index]
                            for el in self.shiftedMjd]

        if ccMjdMaxFlux != 0:
            epochs = [el + ccMjdMaxFlux for el in self.shiftedMjd]
            self.shiftedMjd = epochs
            del epochs

    def calc_destretched_time(self, redshift):
        """Perform correction for time dialtion. Returns an array of time in
        mjd.

        Keywords arguments
        redshift    -   from a object of type Supernova.
        """
        epochs = [el/(1.+redshift) for el in self.mjd]

        return epochs

    def calc_dereddened_flux(self, R, ebv):
        """Perform correction for MW dust reddening

        Keyword arguments:
        R   -   depends on the band
        ebv -   colour excess, from Supernova object
        """

        a_mag = R * ebv
        a_flux = 10**(0.4*a_mag)

        return [el*a_flux for el in self.flux]

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


    def set_shifted_mjd(self, supernovaFitObj, destretchFlag=True):
        """Apply `set_shifted_mjd_2` from Supernova class to all the bands.
        Depending on `destretchFlag` can or cannot perform correction
        for time delay.
        """
        for b in self.lcsDict.keys():
            if destretchFlag:
                self.lcsDict[b].set_shifted_mjd_2(
                    supernovaFitObj.lcsDict[b],
                    (self.zSpec if self.zSpec else self.zPhotHost),
                    supernovaFitObj.ccMjdMaxFlux
                    )
            else:
                self.lcsDict[b].set_shifted_mjd(
                    supernovaFitObj.lcsDict[b])

    # def calc_dereddened_flux(self, R):
    #     for b in self.lcsDict.keys():
    #         self.lcsDict[b].calc_dereddened_flux(R, self.MWEBV)

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

            if len(bandArr) == 0:
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

        mjd = [mjd[i] for i in mjdArgsort]
        flux = [flux[i] for i in mjdArgsort]
        fluxErr = [fluxErr[i] for i in mjdArgsort]
        bandArr = [bandArr[i] for i in mjdArgsort]

        """
        Adding and setting column to Table
        """

        for c in range(len(colNames)):
            if colNames[c][0] == "BAND" \
            or colNames[c][0] == "OBS" \
            or colNames[c][0] == "FIELD":
                col = MaskedColumn(np.zeros(len(mjd)),
                    name=colNames[c][0],
                    format=colNames[c][1],
                    dtype=np.str, fill_value='-',
                    mask=np.zeros(len(mjd))
                    )
            else:
                col = MaskedColumn(np.zeros(len(mjd)),
                    name=colNames[c][0],
                    format=colNames[c][1],
                    dtype=np.float, fill_value=-9,
                    mask=np.zeros(len(mjd)))
            t.add_column(col)

        """
        Initializing columns
        """

        t["OBS"] = np.empty(len(mjd), dtype=np.str)
        t["OBS"][:] = "OBS:"
        t["FIELD"] = np.empty(len(mjd), dtype=np.str)
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
        fOut.write("# NOBS: {:<}\n".format(len(mjd)))

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

    lsDirData = util.list_files("*SN*.DAT", path=args.dirData+os.sep)

    """

    KERNEL SPECIFICATION

    """
    kern = GPy.kern.RBF(1)
    # kern = GPy.kern.RatQuad(1)

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
        if args.mag:
            for b in candidate.lcsDict.keys():
                candidate.lcsDict[b].lim = limMagDict[b]


        print 'candidate z = {:>6.4f}'.format(
            candidate.zSpec if candidate.zSpec else candidate.zPhotHost)


        """
        Create SupernovaFit objects
        """
        candidateFit = SupernovaFit(candidate, kern.name)

        """
        Looping in bands and fit of flux
        """
        for b in candidate.lcsDict.keys():
            phase = util.time_correct(
                    candidate.lcsDict[b].mjd,
                    candidate.zSpec if candidate.zSpec else candidate.zPhotHost
                    )

            # Correcting for absorption
            flux = util.correct_for_absorption(
                    candidate.lcsDict[b].flux,
                    candidate.MWEBV, b
                    )


            """
            Clipping to limiting magnitudes when flux is above 90th mag
            """
            if args.mag :
                flux = [candidate.lcsDict[b].lim if \
                        (el>90) else el for el in flux]

            errFlux = candidate.lcsDict[b].fluxErr

            # Fitting Lightcurve
            if (not candidate.lcsDict[b].badCurve) and (len(flux) >= 3):

                start_time = time.time()
                print "Profiling gp_fit.\n"
                cProfile.run('util.gp_fit(phase, flux, errFlux,kern, n_restarts=10,parallel=False,test_length=True,test_prior=args.testPrior)')
                predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                phase, flux, errFlux,
                                                kern, n_restarts=10,
                                                parallel=False,
                                                test_length=True,
                                                test_prior=args.testPrior)
                print "\n" + indent \
                    + "The process took {:5.3f} secs.".format(time.time()-start_time)
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
            print 'Profiling get_distace'
            cProfile.run('fit[0].get_distance(fit[1], args.band)')
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

    if args.plot:
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


        plt.show()
