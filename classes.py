import numpy as np
import os
import sys
import glob
import cPickle
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn, vstack, hstack
import time
import argparse



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
    def __init__(self, band):
        self.band = band
        self.mjd = np.ma.zeros(0, dtype=np.float32)
        self.shiftedMjd = np.ma.zeros(0, dtype=np.float32)
        self.flux = np.ma.zeros(0, dtype=np.float32)
        self.fluxErr = np.ma.zeros(0, dtype=np.float32)

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
            self.flux.mask = True
            self.fluxErr.mask = True
            self.mjd.mask = True
        # else:
        #   self.badCurve = False

    def add_data_point(self, mjd, flux, fluxErr):
        """
        Adds a data point to the light curve.
        """
        self.mjd = np.append(self.mjd, np.float32(mjd))
        self.flux = np.append(self.flux, np.float32(flux))
        self.fluxErr = np.append(self.fluxErr, np.float32(fluxErr))

        #update the mask
        self.mjd.mask = np.zeros(self.mjd.size)
        self.flux.mask = np.zeros(self.flux.size)
        self.fluxErr.mask = np.zeros(self.fluxErr.size)

    def set_shifted_mjd(self, distance):
        """
        Construct shifted_mjd, by subtracting 'distance' from 'self.flux'
        """
        self.shiftedMjd = np.ma.subtract(self.mjd, distance)
    
    def get_max_fmfe_Index(self):
        """
        Return the index of max (flux - fluxErr)
        """
        difference = np.subtract(self.flux, self.fluxErr)
        
        return np.argmax(difference)
    
    def get_max_flux_p(self, p):
        """
        Returns max (flux - p*fluxErr)
        """
        return np.max(np.subtract(self.flux, p*self.fuxErr))
    
    @property
    def max_flux(self):
        if not self.badCurve:
            result = self.flux.max()
        else:
            result = 0

        return result

    @property
    def max_error(self):
        if not self.badCurve:
            result = self.fluxErr.max()
        else:
            result = 0

        return result

    def reset_mask(self):
        if (np.nonzero(self.mjd.mask)[0].size > 0) and (not self.badCurve):
            self.mjd.mask = np.zeros(self.size)
        
        if (np.nonzero(self.flux.mask)[0].size > 0) and (not self.badCurve): 
            self.flux.mask = np.zeros(self.size)

        if (np.nonzero(self.fluxErr.mask)[0].size > 0) and (not self.badCurve):  
            self.fluxErr.mask = np.zeros(self.size)

    @property
    def max_flux_index(self):
        """
        Return the index of the maximum flux
        """
        return np.argmax(self.flux)

    @property
    def size(self):
        return self.mjd.size


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
                            self.g.add_data_point(mjd, flux, fluxErr)
                        elif passband == "r":
                            self.r.add_data_point(mjd, flux, fluxErr)
                        elif passband == "i":
                            self.i.add_data_point(mjd, flux, fluxErr)
                        elif passband == "z":
                            self.z.add_data_point(mjd, flux, fluxErr)
                        else:
                            print "Filter not recognized: {:<5}".format(passband)
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

        for b in self.lcsDict.keys():
            self.lcsDict[b].set_badCurve()

    def __cmp__(self, other):
        return 2*(self.zPhotHost - other.zPhotHost > 0) - 1 


class SupernovaFit():
    def __init__(self, SNID, 
        SNType=None, RADeg=None, decDeg=None, MWEBV=None, 
        zSpec=None, zSpecErr=None,
        hostGalaxyID=None, zPhotHost=None, zPhotHostErr=None):
        self.g = LightCurve("g")
        self.r = LightCurve("r")
        self.i = LightCurve("i")
        self.z = LightCurve("z")
        self.lcsDict = {"g":self.g, 
                        "r":self.r, 
                        "i":self.i, 
                        "z":self.z}
        self.peaked = False
        self.SNID = SNID
        self.SNType = SNType 
        self.RADeg = RADeg 
        self.decDeg = decDeg 
        self.MWEBV = MWEBV 
        self.zSpec = zSpec 
        self.zSpecErr = zSpecErr
        self.hostGalaxyID = hostGalaxyID 
        self.zPhotHost = zPhotHost 
        self.zPhotHostErr = zPhotHostErr 

    def set_lightcurve(self, band, mjd, flux, fluxErr):
        self.lcsDict[band].mjd = np.ma.asarray(mjd, dtype=np.float32)
        self.lcsDict[band].flux = np.ma.asarray(flux, dtype=np.float32)
        self.lcsDict[band].fluxErr = np.ma.asarray(fluxErr, 
            dtype=np.float32)

        if (band == 'r') \
            and self.r.max_flux_index not in set([0, -1, self.r.size-1]):
                self.set_peaked()

        # setting the masks
        self.lcsDict[band].mjd.mask = np.zeros(mjd.size)
        self.lcsDict[band].flux.mask = np.zeros(flux.size)
        self.lcsDict[band].fluxErr.mask = np.zeros(fluxErr.size)

        self.lcsDict[band].set_badCurve()

    def shift_mjds(self):
        try:
            # mjdrMax = self.r.mjd[self.r.flux == self.r.flux.max()][0]
            mjdrMax = self.r.mjd[self.r.max_flux_index]
            # idx_rMax = np.where(self.r.flux == self.r.flux.max())[0][0]
            for b in self.lcsDict.keys():
                if not self.lcsDict[b].badCurve:
                    self.lcsDict[b].set_shifted_mjd(mjdrMax)
        except:
            print "Candidate {:<d} ".format(self.SNID) + \
            "has no maximum in r-band."

    def reset_masks(self):
        for b in self.lcsDict.keys():
            self.lcsDict[b].reset_mask()

    def normalized_flux(self, band):
        """Normalizes the light curve in `band` using the maximum in that band.
        s is a slice on the array.
        """
        result = np.array(0)


        result = self.lcsDict[band].flux.compressed() / (\
            self.g.max_flux + 
            self.r.max_flux + 
            self.i.max_flux + 
            self.z.max_flux
            )

        return result

    def normalized_error(self, band):
        """Normalizes the light curve in band b using the maximum in that band.
        s is a slice on the array.
        """
        result = np.array(0)
        
        result = self.lcsDict[band].fluxErr.compressed() / (\
            self.g.max_error + 
            self.r.max_error + 
            self.i.max_error + 
            self.z.max_error
            )

        return result

    def set_peaked(self):
        self.peaked = True

    @property
    def peaked(self):
        return self.peaked

    def get_distance(self, candidate, band, reset_masks=True):
        """Calculate difference (aka distance) between two 
        interpolated light curves. 
        """
        if type(band) is not str:
            raise TypeError
        distance = -99
        
        sizeSelf = self.lcsDict[band].size
        sizeCandidate = candidate.lcsDict[band].size

        
        idxSelfMax = np.argmin(
            np.abs(self.lcsDict[band].shiftedMjd))#[0][0]
        idxCandidateMax = np.argmin(
                np.abs(candidate.lcsDict[band].shiftedMjd))#[0][0]
        # print idxSelfMax,idxCandidateMax
        # these checks on maximum positions should be done outside the function
        # 
        # (at least the first that sets distance to a big value)
        # 
        # First check: the lc are set as `not peacked` but we check if the 
        # 

        if sizeSelf >= sizeCandidate:
            bigger = self.lcsDict[band].shiftedMjd.view()
            smaller = candidate.lcsDict[band].shiftedMjd.view()
        else:
            bigger = candidate.lcsDict[band].shiftedMjd.view()
            smaller = self.lcsDict[band].shiftedMjd.view()

        # modifying masks
        smaller.mask = np.in1d(
            np.round(smaller), np.round(bigger), invert=True)
        bigger.mask = np.in1d(
            np.round(bigger), np.round(smaller), invert=True)

        # setting new masks
        self.lcsDict[band].flux.mask = \
                    self.lcsDict[band].shiftedMjd.mask
        self.lcsDict[band].fluxErr.mask = \
                    self.lcsDict[band].shiftedMjd.mask

        candidate.lcsDict[band].flux.mask = \
                    candidate.lcsDict[band].shiftedMjd.mask
        candidate.lcsDict[band].fluxErr.mask = \
                    candidate.lcsDict[band].shiftedMjd.mask            

        # calculating min a max overlapping MJD
        minMjd = min(self.lcsDict[band].shiftedMjd.compressed()[0],
            candidate.lcsDict[band].shiftedMjd.compressed()[0])

        maxMjd = max(self.lcsDict[band].shiftedMjd.compressed()[-1],
            candidate.lcsDict[band].shiftedMjd.compressed()[-1])

        # calculating the distance
        distance = (1. / (maxMjd - minMjd)) * np.sqrt(np.ma.sum(
                    np.ma.divide(
                        np.ma.power(
                            np.ma.subtract(
                                self.normalized_flux(band),
                                candidate.normalized_flux(band)
                                ), 2), 
                        np.ma.add(
                            np.ma.power(self.normalized_error(band), 2), 
                            np.ma.power(candidate.normalized_error(band), 2)
                        )
                    )
                ))

        if reset_masks:
            self.reset_masks()
            candidate.reset_masks()

        return distance


    def save_on_txt(self, fileName):
        t = Table(masked=True)
        colNames = [["MJD_r_band", "{0:5.0f}"],
                    ["FLUX_g", "{0:10.5f}"], ["FLUX_ERR_g", "{0:10.5f}"],
                    ["FLUX_r", "{0:10.5f}"], ["FLUX_ERR_r", "{0:10.5f}"],
                    ["FLUX_i", "{0:10.5f}"], ["FLUX_ERR_i", "{0:10.5f}"],
                    ["FLUX_z", "{0:10.5f}"], ["FLUX_ERR_z", "{0:10.5f}"]]

        for c in range(len(colNames)):
            col = MaskedColumn(np.zeros(self.r.mjd.size),
                name=colNames[c][0],
                format=colNames[c][1],
                dtype=np.float, fill_value=-9,
                mask=np.zeros(self.r.mjd.size))
            t.add_column(col)

        t["MJD_r_band"] = self.r.mjd
        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                t["FLUX_{:<1}".format(b)].mask = np.ones(self.r.mjd.size)
                t["FLUX_ERR_{:<1}".format(b)].mask = np.ones(self.r.mjd.size)
                t["FLUX_{:<1}".format(b)].format = "{}"
                t["FLUX_ERR_{:<1}".format(b)].format = "{}"
            else:
                t["FLUX_{:<1}".format(b)] = self.lcsDict[b].flux
                t["FLUX_ERR_{:<1}".format(b)] = self.lcsDict[b].fluxErr

        t.filled()

        fOut = open(fileName, 'w')
        fOut.write("# File produced by Miniature Adventure on " + \
            "{:<02d}/{:<02d}/{:<4d} at {:<02d}:{:<02d}:{:<02d} GMT\n".format(
                time.gmtime().tm_mday, time.gmtime().tm_mon, 
                time.gmtime().tm_year,
                time.gmtime().tm_hour, time.gmtime().tm_min, 
                time.gmtime().tm_sec))
        fOut.write("# SNID: {:>10d}\n".format(self.SNID))

        if self.SNType :
            fOut.write("# SNType: {:>d}\n".format(self.SNType))
        if self.RADeg :
            fOut.write("# RADeg: {:>9.6f}\n".format(self.RADeg))
        if self.decDeg :
            fOut.write("# decDeg: {:>9.6f}\n".format(self.decDeg))
        if self.MWEBV :
            fOut.write("# MWEBV: {:>6.4f}\n".format(self.MWEBV))
        if self.zSpec :
            fOut.write("# zSpec: {:>6.4f}\n".format(self.zSpec))
        if self.zSpecErr:
            fOut.write("# zSpecErr: {:>6.4f}\n".format(self.zSpec))
        if self.hostGalaxyID :
            fOut.write("# hostGalaxyID: {:>d}\n".format(self.hostGalaxyID))
        if self.zPhotHost :
            fOut.write("# zPhotHost: {:>6.4f}\n".format(self.zPhotHost))
        if self.zPhotHostErr :
            fOut.write("# zPhotHostErr: {:>6.4f}\n".format(self.zPhotHostErr))

        ascii.write(t, output=fOut, delimiter='  ', 
            format='fixed_width_two_line')
        
        fOut.close()


class CandidatesCatalog():
    """Used to store in R.A.M. the list of fit to candidates lightcurves.
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