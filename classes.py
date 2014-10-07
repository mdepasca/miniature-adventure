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

warnings.filterwarnings('error', message=".*divide by zero encountered in double_scalars.*", category=RuntimeWarning)
from math import sqrt

if __name__ == '__main__':
    import utilities as util
    import GPy
    from matplotlib import pyplot as plt

    parser = argparse.ArgumentParser(
        description="Test of general functions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-c1", "--candidate1", dest="candidate1", 
        type=np.int32, default=None, 
        help="First candidate idx")

    parser.add_argument(
        "-c2", "--candidate2", dest="candidate2", 
        type=np.int32, default=None, 
        help="Second candidate idx")    

    parser.add_argument(
        "-b", "--band", dest="band", default='r', 
        help="Photometric band.")
    args = parser.parse_args()


class LightCurve():
    """
    Once fully initiated, instances of this class have the following important 
    properties
    band            (string) 'g','r','i' or 'z'
    mjd             (array) modified julian dates of observations
    flux            (array) the observed flux
    fluxErr         (array) the error in the observed flux
    shifted_mjd     (array) mjd shifted such that the peak has mjd = 0. To be 
                            modified
    """
    badCurve = False
    shifted_mjd = np.zeros(0)
    """
    TRY TO USE __slots__
    """
    __slots__ = ['band', 'mjd', 'shiftedMjd', 'flux', 'fluxErr', 'badCurve',
        'shifted_mjd']
    def __init__(self, band):
        self.band = band
        self.mjd = list()#np.zeros(0, dtype=float)
        self.shiftedMjd = list()#np.zeros(0, dtype=float)
        self.flux = list()#np.zeros(0, dtype=float)
        self.fluxErr = list()#np.zeros(0, dtype=float)

    @classmethod
    def data(band, mjd, flux, fluxErr):
        self.band = band
        self.mjd = mjd
        self.flux = flux
        self.fluxErr = fluxErr
        self.set_badCurve()

    def set_badCurve(self): 
        if len(self.flux) == 0:
            self.badCurve = True

    def set_shifted_mjd(self, distance):
        """
        Construct shiftedMjd, by subtracting 'distance' from 'self.flux'
        """
        # self.shiftedMjd = np.subtract(self.mjd, distance)
        self.shiftedMjd = [self.mjd[i]-distance for i in range(len(self.mjd))]
    
    @property
    def max_flux(self):
        if not self.badCurve:
            result = max(self.flux)
        else:
            result = 0

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
        return self.flux.index(max(self.flux))

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
    def __init__(self, inFileName):
        """
        Parses all the light curve data in inFileName into a Supernova object.
        """

        inFile = file(inFileName, "r")
        lines = inFile.readlines()
        inFile.close()

        self.g = LightCurve("g")
        self.r = LightCurve("r")
        self.i = LightCurve("i")
        self.z = LightCurve("z")

        # list_gMjd = list()
        # list_rMjd = list()
        # list_iMjd = list()
        # list_zMjd = list()

        # list_gFlux = list()
        # list_rFlux = list()
        # list_iFlux = list()
        # list_zFlux = list()

        # list_gFluxErr = list()
        # list_rFluxErr = list()
        # list_iFluxErr = list()
        # list_zFluxErr = list()

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
                elif tag == "SNTYPE":
                    self.SNTypeInt = int(data[0])                    
                elif tag == "RA":
                    self.RADeg = float(data[0])
                elif tag == "DECL":
                    self.decDeg = float(data[0])
                elif tag == "MWEBV":
                    self.MWEBV = float(data[0])
                elif tag == "REDSHIFT_SPEC":
                    if float(data[0]) == -9:
                        self.zSpec = None
                    else:
                        self.zSpec = float(data[0])
                        self.zSpecErr = float(data[2])
                elif tag == "HOST_GALAXY_GALID":
                    self.hostGalaxyID = int(data[0])
                elif tag == "HOST_GALAXY_PHOTO-Z":
                    self.zPhotHost = float(data[0])
                    self.zPhotHostErr = float(data[2])
                elif tag == "MJD_MAX_FLUX-CCF":
                    self.ccMjdMaxFlux = float(data[0])

        for b in self.lcsDict.keys():
            self.lcsDict[b].set_badCurve()

    def __cmp__(self, other):
        return 2*(self.zPhotHost - other.zPhotHost > 0) - 1 


class SupernovaFit():
    ccMjdMaxFlux = 0

    def __init__(self, supernova):
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

    def set_lightcurve(self, band, mjd, flux, fluxErr):
        self.lcsDict[band].mjd = mjd
        self.lcsDict[band].flux = flux
        self.lcsDict[band].fluxErr = fluxErr

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

        den = self.g.max_flux + \
            self.r.max_flux + \
            self.i.max_flux + \
            self.z.max_flux

        flux = self.lcsDict[band].flux
        result = [flux[idx]/den for idx in range(len(flux))]

        return result

    def normalized_error(self, band):
        """Normalizes the light curve in band b using the maximum in that band.
        s is a slice on the array.
        """

        den = self.g.max_flux + \
            self.r.max_flux + \
            self.i.max_flux + \
            self.z.max_flux

        fluxErr = self.lcsDict[band].fluxErr
        result = [fluxErr[idx]/den for idx in range(len(fluxErr))]
        return result

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
        distance = -99
        bigDistance = 1.01

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
            distance = bigDistance
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
        colNames = [["OBS", "{:4s}"],
                    ["MJD", "{0:9.3f}"], ["BAND", "{:s}"], ["FIELD", "{:6s}"],
                    ["FLUX", "{0:10.5f}"], ["FLUX_ERR", "{0:10.5f}"]]

        bandArr = np.empty(0)
        mjd = np.empty(0)
        flux = np.empty(0)
        fluxErr = np.empty(0)
        mjdArgsort = np.empty(0)

        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                continue

            if bandArr.size == 0:
                bandArr = np.empty(self.lcsDict[b].mjd.size, dtype=np.str)
                bandArr[:] = b
                mjd = self.lcsDict[b].mjd
                flux = self.lcsDict[b].flux
                fluxErr = self.lcsDict[b].fluxErr
            else:
                tmp = np.empty(self.lcsDict[b].mjd.size, dtype=np.str)
                tmp[:] = b
                bandArr = np.concatenate((bandArr, tmp))
                mjd = np.concatenate((mjd, self.lcsDict[b].mjd))
                flux = np.concatenate((flux, self.lcsDict[b].flux))
                fluxErr = np.concatenate((fluxErr, self.lcsDict[b].fluxErr))

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


class CandidatesCatalog():
    """Used to store in memory the list of fit to candidates lightcurves.
    Structure inspired from SupernovaeCatalog

    candidate: array of Supernova objects
    SNID: array if IDs. Follows the objs.SNID array
    SNType: array of SN types. Useful to search for a specific type
    """

    def __init__(self):
        self.candidates = np.zeros(0, dtype=np.object)
        self.SNID = np.zeros(0, dtype=np.int)
        self.SNType = np.zeros(0, dtype=np.int)
        self.peaked = np.zeros(0, dtype=np.int)

    def add_candidate(self, candidate):
        self.candidates = np.append(self.candidates, candidate)
        self.SNID = np.append(self.SNID, candidate.SNID)
        self.SNType = np.append(self.SNType, candidate.SNType)
        self.peaked = np.append(self.peaked, candidate.peaked)

    def merge(self, catalog):
        for i in range(catalog.size):
            self.add_candidate(catalog.candidates[i])

    def get_peaked(self):
        return np.where(self.peaked)[0]

    @property
    def size(self):
        return self.SNID.size


class SupernovaeCatalog():
    """
    Class variables are
    sne     (object array)  array of SupernovaFit objects
    zSpec   (float array) the spectroscopically observed redshift 
            (None if no spctrscp) 
    zPhotHost   (float array)   the redshift of the host galaxy
    SNID    (int array) the IDs of the Supernovae
    SNType  (int array) The types of the supernovae
    """
    def __init__(self, dataDir, load_all):
        """
        This class loads in all the light curve data under dataDir into a big 
        list, and creates a series of top level arrays that we can use to cut 
        the catalog by z, type etc. load_all: Do you want to
        load all of the data, or will your computer then crash?
        """
        
        print ">>> Loading data from text files ..."
        inFileNames = glob.glob(dataDir+os.path.sep+"DES_SN*.DAT")
        self.sne = np.zeros(0, dtype=np.object)
        self.zPhotHost = np.zeros(0, dtype=np.float32)
        self.SNID = np.zeros(0, dtype=np.int)
        self.SNType = np.zeros(0, dtype=np.int)
        count = 0
        for inFileName in inFileNames:
            count += 1
            # Progress update
            tenPercent = len(inFileNames)/10
            for j in range(0,11):
                if count == j*tenPercent:
                    print "... "+str(j*10)+"% complete ..."
            
            
            inFile = file(inFileName, "r")  
            for i in range(4):
                line = inFile.readline()
            inFile.close()
            
            SNTYPE = line.split(":")[-1].split()[0]
            if SNTYPE == '-9' and load_all == False:
                pass        

            else:   
                sn = Supernova(inFileName)            
                self.SNType = np.append(self.SNType, sn.SNTypeInt)
                self.SNID   = np.append(self.SNID, sn.SNID)
                self.sne    = np.append(self.sne, sn)

        self.zPhotHost = np.nan_to_num(self.zPhotHost)
            
    def findSupernovae(self, SNType, zSpecLow, zSpecHigh):
        """
        Given a SNType code and a redshift range, return a list of matching 
        supernovae in the catalog.      
        """
        typeMask = np.equal(self.SNType, SNType)
        zSpecMask = np.logical_and(np.greater(self.zSpec, zSpecLow), 
                                   np.less(self.zSpec, zSpecHigh))
        mask = np.logical_and(typeMask, zSpecMask)
        foundSupernovae = self.sne[mask]
        
        return foundSupernovae

    def getSNTypeStr(self, SNTypeInt):
        """Given a SNTypeInt, returns a string, e.g. 'Ia'
        
        Mapping to type names (from DES_BLIND+HOSTZ.README):
        
            1  (Ia)
            2  (II)
            3  (Ib/c)
            11  (pec. Ia)
            66  (other)
            -1  (rejected)
            -9 (unknown)
            
        """
        
        if SNTypeInt == 1:
            SNTypeStr="Ia"
        elif SNTypeInt == 2:
            SNTypeStr="II"
        elif SNTypeInt == 3:
            SNTypeStr="Ib/c"
        elif SNTypeInt == 11:
            SNTypeStr="pec. Ia"
        elif SNTypeInt == 66:
            SNTypeStr="other"
        elif SNTypeInt == -1:
            SNTypeStr="rejected"
        elif SNTypeInt == -9:
            SNTypeStr="unclassified"
        
        return SNTypeStr

    def getSNTypeInt(self, SNTypeStr):
        """Given a SNTypeStr, returns the corresponding int, e.g. 1
        
        Mapping to type names (from DES_BLIND+HOSTZ.README):
        
            1  (Ia)
            2  (II)
            3  (Ib/c)
            11  (pec. Ia)
            66  (other)
            -1  (rejected)
            -9 (unknown)
            
        """
            
        if SNTypeStr == "Ia":
            SNTypeInt = 1
        elif SNTypeStr == "II":
            SNTypeInt = 2
        elif SNTypeStr == "Ib/c":
            SNTypeInt = 3
        elif SNTypeStr == "pec. Ia":
            SNTypeInt = 11
        elif SNTypeStr == "other":
            SNTypeInt = 66
        elif SNTypeStr == "rejected":
            SNTypeInt = -1
        elif SNTypeStr == "unclassified":
            SNTypeInt = -9
        
        return SNTypeInt


if __name__ == '__main__':
    indent = "          "
    dirData = "train_data" + os.sep + "DES_BLIND+HOSTZ"
    fCandidatesList = "DES_BLIND+HOSTZ.LIST"
    candidatesFileList = np.genfromtxt(dirData+os.sep+fCandidatesList, dtype=None)
    catalogFit = CandidatesCatalog()
    kern = GPy.kern.RBF(1)

    print args.candidate1
    print args.candidate2

    if args.band not in ['g', 'r', 'i', 'z']:
        print 'Band {:<} not recognised! Changing to r'.format(args.band)
        args.band = 'r'

    if args.candidate1 is None:
        args.candidate1 = np.random.random_integers(
                low=0, high=18321)

    if args.candidate2 is None:
        args.candidate2 = np.random.random_integers(
                low=0, high=18321)

    while args.candidate2 == args.candidate1:
        args.candidate2 = np.random.random_integers(
            low=0, high=18321)


    # Getting observation data
    candidates = list()
    candidates.append(Supernova(
        dirData+os.sep+candidatesFileList[args.candidate1]))

    candidates.append(Supernova(
        dirData+os.sep+candidatesFileList[args.candidate2]))

    for candidate in candidates:
        # Create SupernovaFit objects
        candidateFit = SupernovaFit(candidate.SNID)
        for b in candidate.lcsDict.keys():
            phase = candidate.lcsDict[b].mjd
            flux = candidate.lcsDict[b].flux
            errFlux = candidate.lcsDict[b].fluxErr

            # Fitting Lightcurve
            if (not candidate.lcsDict[b].badCurve) \
            and (flux.size >= 3):
                saveOut = sys.stdout
                fout = open('test_out.log', 'w')
                # fout = open('/dev/null', 'w')
                sys.stdout = fout

                
                predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                phase, flux, errFlux, 
                                                kern, n_restarts=10, 
                                                test_length=True)
                sys.stdout = saveOut
                fout.close()

                candidateFit.set_lightcurve(b, 
                    predMjd.reshape(predMjd.size),
                    predFlux.reshape(predFlux.size), 
                    predErr.reshape(predErr.size))
                
                print indent + \
                    "{:<} {:<}".format(candidate.SNID, b)
            else:
                candidateFit.lcsDict[b].badCurve = True
                print indent + util.bcolors.FAIL + \
                    "{:<} {:<}".format(candidate.SNID, b) + \
                    util.bcolors.ENDC

        candidateFit.shift_mjds()
        catalogFit.add_candidate(candidateFit)

    if catalogFit.candidates[0].peaked and catalogFit.candidates[1].peaked:
        print 'Distance between the 2 normalized lcs in ' + \
        '{:<} band = {:<2.4f}'.format(args.band,
            catalogFit.candidates[0].get_distance(catalogFit.candidates[1], 
            args.band, reset_masks=False))

        # if plt.get_fignums():
        #     figNum = plt.get_fignums()[-1]+1
        # else:
        #     figNum = 1

        # plt.figure(figNum)
        plt.scatter(
            catalogFit.candidates[0].lcsDict[args.band].shiftedMjd.compressed(),
            catalogFit.candidates[0].normalized_flux(args.band))

        plt.scatter(
            catalogFit.candidates[1].lcsDict[args.band].shiftedMjd.compressed(),
            catalogFit.candidates[1].normalized_flux(args.band), c='orange')

        plt.show()
    else:
        print 'One of the 2 candidate has not r-band peak: '
        print 'Candidate 1 - {:<}'.format(catalogFit.candidates[0].peaked)
        print 'Candidate 2 - {:<}'.format(catalogFit.candidates[1].peaked)


    # return catalogFit