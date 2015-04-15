"""
Implementation of general use functions.
"""
import numpy as np
import numpy.ma as ma
import pandas as pd
from pandas import DataFrame
import cPickle
import gzip
import GPy
import classes
import time
import os
from os import path
import argparse
import matplotlib.pyplot as plt
from cStringIO import StringIO
import subprocess

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

global Rband

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





def flux_error_to_mag_error(fluxErr, flux):
    # error prop. zErr = abs(dF(x)/dx) xErr
    magErr = np.multiply(np.abs(np.divide(-2.5, np.log(10) * flux)),
                         np.abs(fluxErr))

    return magErr

"""
----------------------------------------------------------------------------






INPUT/OUTPUT on file







"""
def open_gzip_pkl_catalog(path):
    """Opens gzipped cPickle file containing supernova catalogue (DEPRECATED).

    Keyword arguments:
    path -- path to catalogue.
    """
    f = gzip.open(path, 'rb')
    catalog = pkl.load(f)
    f.close()

    return catalog

def open_pkl(filePath):
    """Loads data from cPickle file (deprecated).

    Keyword Argument:
    filePath -- string, path to cPickle file.
    """
    fileHandler = open(filePath,'r')
    data = cPickle.load(fileHandler)
    fileHandler.close()

    # elTime -= time.time()
    # print '>>> Done in ' + str(abs(elTime)) + ' secs'
    return data


def dump_pkl(filePath, dataStruct):
    """Dumps dataStruct on file using cPickle module (deprecated).

    Keyword arguments:
    filePath -- string, the path where to save.
    dataStruct -- the data structure to save.
    """

    fileHandler = open(filePath,'w')

    cPickle.dump(dataStruct, fileHandler)
    fileHandler.close()


def check_lc_from_file(fileDir):
    """Creates a `Supernova` objects from files in specified directory and checks for bad light curves.

    Keyword arguments:
    fileDir -- string, directory in which to look for files.
    """
    p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=fileDir)
    lsDir = p.stdout.read()
    lsDir = lsDir.split('\n')
    lsDir.sort()
    lsDir.remove('')

    for i in range(len(lsDir)):
        tmpSN = get_sn_from_file(fileDir+lsDir[i])
        if tmpSN.r.badCurve:
            print "{:<} Has bad curve in r band - ".format(lsDir[i]) +\
            "THE FILE HAS TO BE DELETED"


def create_file(indexList, outFilePath):
    """Creates file contaning list of files 'DES_SN*_FIT.DAT' in directory train_data/DES_BLIND+HOSTZ_FIT/.

    Keyword arguments:
    indexList -- Python list, contains supernova IDs.
    outFilePath -- string, path where to create the output file.
    """
    # indexArr = np.loadtxt(indexFile, dtype=np.int)
    outFile = open(outFilePath, 'w')

    for i in indexList:
        filePath = 'train_data/DES_BLIND+HOSTZ_FIT/' + \
            'DES_SN{:0>6d}_FIT.DAT'.format(i)
        try:
            f = open(filePath, 'r')

            outFile.write(filePath+'\n')
        except IOError:
            continue

    outFile.close()

def index_to_filename(indexList, inFileName, outFileName):
    """Filters a list of files using specified indexes.

    Keyword arguments:
    indexList -- Python list, indices to keep in the output.
    inFileName -- path to file containing list of files.
    outFileName -- path to output file.
    """
    inFile = open(inFileName, "r")
    lines = inFile.readlines()
    inFile.close()

    npLines = np.array(lines, dtype=np.str)

    outFileList = npLines[indexList]
    np.savetxt(outFileName, outFileList, fmt='%s', newline='')

def rename_bad_r_lc_file(path):
    """Renames files of fitted lc with a bad lc in r band to extension '.BADrLC'.

    Keyword arguments:
    path -- string, path to directory in which to find files to checked.
    """
    if path[-1] != os.sep:
        path = path + os.sep

    p = subprocess.Popen("ls *FIT.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=path)
    lsList = p.stdout.read()
    lsList = lsList.split('\n')
    lsList.sort()
    lsList.remove('')
    countBad = 0
    for i in range(len(lsList)):
        tmpSN = get_sn_from_file(path+lsList[i])
        if tmpSN.r.badCurve:
            os.rename(path+lsList[i], path+lsList[i]+'.BADrLC')
            countBad += 1

    return countBad


def extract_redshift_data(path, outFile):
    """Extract redshift from files and produces a CSV file, to be read from R to study redshift distribution.

    Keyword arguments:
    path -- where to find supernova files.
    outFile -- path to output file.

    Notes:
    The output CSV file will have 4 columns:
    SNID -- integer
    Redshift -- float, spectroscopic or photometric
    Training flag -- 1 for training 0 for test
    SN Type -- from 'SIMGEN_PUBLIC_DES.DUMP' file
    """

    if path[-1] != os.sep:
        path = path + os.sep

    p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=path)
    lsList = p.stdout.read()
    lsList = lsList.split('\n')
    lsList.sort()
    lsList.remove('')

    dump = pd.read_csv(
        'train_data/SIMGEN_PUBLIC_DES/SIMGEN_PUBLIC_DES.DUMP',
        sep=' ', skiprows=0, header=1, usecols=[1,2],
        skipinitialspace=True, engine='c')

    dump = dump.convert_objects(convert_numeric=True, copy=False)

    snid = np.empty(len(lsList), dtype=np.int)
    redshift = np.empty(len(lsList), dtype=np.float)
    trainFlag = np.zeros(len(lsList), dtype=np.int)
    genType = np.zeros(len(lsList), dtype=np.int)

    for i in range(len(lsList)):
        tmpSN = get_sn_from_file(path+lsList[i])

        snid[i] = tmpSN.SNID
        redshift[i] = tmpSN.zSpec if (tmpSN.zSpec != None) else tmpSN.zPhotHost
        trainFlag[i] = 1 if (tmpSN.zSpec != None) else 0
        genType[i] = dump['GENTYPE'][dump['CID']==snid[i]]

    df = pd.DataFrame(
        data=zip(snid, redshift, trainFlag, genType),
        columns=['SNID', 'redshift', 'train_flag', 'genType'])

    df.to_csv(
        'products/'+outFile, sep=';', index=False,
        float_format='%5.4f', header=True)

def extract_training_set(path, fileName):
    """Creates files dividing supernovae in training and test sets. It creates also files list training set supernovae by type

    Keyword arguments:
    path -- where to find supernova light curves files

    Notes:
    Created files are saved in directory 'products/'. Their name are, so far, fixed.
    fileName.TEST
    fileName.TRAIN
    fileName.[SNType].TRAIN
    """
    if path[-1] != os.sep:
        path = path + os.sep

    p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=path)
    lsList = p.stdout.read()
    lsList = lsList.split('\n')
    lsList.sort()
    lsList.remove('')

    # if path.rfind('/') == len(path)-1:
    #     fileName = path.rpartition('/')[0].rpartition('/')[-1]
    # else:
    #     fileName = path.rpartition('/')[-1]

    outFileTest  = open('{:s}{:s}.TEST'.format(path, fileName), 'w')
    outFileTrain = open('{:s}{:s}.TRAIN'.format(path, fileName), 'w')
    outFileIa    = open('{:s}{:s}.Ia.TRAIN'.format(path, fileName), 'w')
    outFileII    = open('{:s}{:s}.II.TRAIN'.format(path, fileName), 'w')
    outFileIbc   = open('{:s}{:s}.Ibc.TRAIN'.format(path, fileName), 'w')
    outFileIaPec = open('{:s}{:s}.IaPec.TRAIN'.format(path, fileName), 'w')
    outFileOther = open('{:s}{:s}.Other.TRAIN'.format(path, fileName), 'w')
    outFileRej   = open('{:s}{:s}.Rej.TRAIN'.format(path, fileName), 'w')

    outFileTest.write('# {:s}\n'.format(path))
    outFileTrain.write('# {:s}\n'.format(path))
    outFileIa.write('# {:s}\n'.format(path))
    outFileII.write('# {:s}\n'.format(path))
    outFileIbc.write('# {:s}\n'.format(path))
    outFileIaPec.write('# {:s}\n'.format(path))
    outFileOther.write('# {:s}\n'.format(path))
    outFileRej.write('# {:s}\n'.format(path))

    for i in range(len(lsList)):
        tmpSN = get_sn_from_file(path+lsList[i])
        SNType = tmpSN.SNTypeInt
        if SNType != -9:
            outFileTrain.write(
                "{:0>5d}      {:0>6d}   {:>}   train\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
        else:
            outFileTest.write(
                "{:0>5d}      {:0>6d}   {:>}   test\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 1:
            outFileIa.write(
                "{:0>5d}      {:0>6d}   {:>}   snIa\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 21 or SNType == 22 or SNType == 23:
            outFileII.write(
                "{:0>5d}      {:0>6d}   {:>}   snII\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 3 or SNType == 32 or SNType == 33:
            outFileIbc.write(
                "{:0>5d}      {:0>6d}   {:>}   snIbc\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 11:
            outFileIaPec.write(
                "{:0>5d}      {:0>6d}   {:>}   pec\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 66:
            outFileOther.write(
                "{:0>5d}      {:0>6d}   {:>}   other\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == -1:
            outFileIbc.write(
                "{:0>5d}      {:0>6d}   {:>}   snIbc\n".format(
                    i, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        outFileOther.write(
            "{:0>5d}      {:0>6d}   {:>}   other\n".format(
                i, tmpSN.SNID, path+lsList[i]
                )
            )

    outFileTest.close()
    outFileTrain.close()
    outFileIa.close()
    outFileII.close()
    outFileIbc.close()
    outFileIaPec.close()
    outFileOther.close()
    outFileRej.close()

def rewrite_file(fileName):
    """Rewrites files produced after fit of old code version.
    It removes `#` at beginning of the first 10 lines, leaving as it is the
    first line.
    If necessary it adds `MJD_MAX_FLUX-CCF:      0.000`.
    Adds column `OBS` containing row names.
    --- DEPRECATED ---
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


def build_distance_matrix(path):
    """
    NOT WORKING!!!
    Builds the full distance matrix from partial files in `path`.
    Reads the list of files in path having `dist_matrix_Sum*` in the name.
    Order them looking at sum of first row and first column indeces (time stamp
    in file name is no reliable). These values are in file header in this
    fashion:
    # distMatrix[r_1:r_2,c_1:c_2] --- Created by `computerName`
    Creates the matrix and saves it for R processing.
    """

    p = subprocess.Popen("ls dist_matrix_Sum*", shell=True,
        stdout=subprocess.PIPE, cwd=path)
    fileNames = p.stdout.read()
    fileNames = fileNames.split('\n')
    fileNames.sort()
    fileNames.remove('')

    print fileNames
    print '*---------------------\n'
    """
    number of elements in the triangle of a square n*n triangular matrix
    m = 0.5*(n*(n+1)) -> 2m = n^2 + n = n(n+1)
    """

    fID = []
    dictFile = dict()

    n = int((1 + np.sqrt(1+8*len(fileNames)))/2) # here in case of a better coding based on recursive function (maybe)
    print n, len(fileNames)
    for f in fileNames:
        shellCmd = 'head -n 1 {:s}'.format(path+f)
        p = subprocess.Popen(shellCmd, shell=True, stdout=subprocess.PIPE)
        head = p.stdout.read()
        head = head.split(',')
        rowIdx = int(head[0].split(':')[0].split('[')[1]) + \
                    int(head[1].split(':')[0])

        colIdx = int(head[0].split(':')[1][0:-1]) + \
                    int(head[1].split(':')[1].split(']')[0])
        fID.append(rowIdx + colIdx)
        dictFile[fID[-1]] = f

    # To get the right order in `fileNames` I order `fID` and progate the
    #
    # ordering to `fileNames`



    fID.sort()
    for i in fID:
        print path+dictFile[i]
    raise SystemExit
    print fID
    mileStones = fID[0:n]
    print 'mat00 ...'
    mat00 = np.loadtxt(path+dictFile[fID[0]])
    print 'mat01 ...'
    mat01 = np.loadtxt(path+dictFile[fID[1]])
    print 'mat02 ...'
    mat02 = np.loadtxt(path+dictFile[fID[2]])
    print 'mat03 ...'
    mat03 = np.loadtxt(path+dictFile[fID[3]])

    print 'hstacking mat0 ...'
    mat0 = np.hstack((mat00, mat01, mat02, mat03))


    print 'mat10 ...'
    mat10 = np.transpose(mat01)
    print 'mat10 ...'
    mat11 = np.loadtxt(path+dictFile[fID[4]])
    print 'mat10 ...'
    mat12 = np.loadtxt(path+dictFile[fID[5]])
    print 'mat10 ...'
    mat13 = np.loadtxt(path+dictFile[fID[6]])

    print 'hstacking mat0 ...'
    mat1 = np.hstack((mat00, mat11, mat12, mat13))


    print 'mat20 ...'
    mat20 = np.transpose(mat02)
    print 'mat21 ...'
    mat21 = np.transpose(mat12)
    print 'mat22 ...'
    mat22 = np.loadtxt(path+dictFile[fID[7]])
    print 'mat23 ...'
    mat23 = np.loadtxt(path+dictFile[fID[8]])

    print 'stacking mat2 ...'
    mat2 = np.hstack((mat20, mat21, mat22, mat23))

    print 'mat30 ...'
    mat30 = np.transpose(mat03)
    print 'mat31 ...'
    mat31 = np.transpose(mat13)
    print 'mat32 ...'
    mat32 = np.transpose(mat23)
    print 'mat33 ...'
    mat33 = np.loadtxt(path+dictFile[fID[9]])

    print 'stacking mat3 ...'
    mat3 = np.hstack((mat30, mat31, mat32, mat33))

    print 'stacking mat ...'
    mat = np.vstack((mat0, mat1, mat2, mat3))

    print 'saveing mat ...'
    np.savetxt(path+'distance_matrix_Sum.txt', mat)

    return 0
"""
----------------------------------------------------------------------------






DATA PROCESSING







"""
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


def time_correct(mjd, zed):
    """Corrects for time dilation.

    Keyword arguments:
    mjd -- Python list, contains measurements epochs in MJDs.
    zed -- float, redshift to use for correction.
    """
    return [el/(1.+zed) for el in mjd]


def k_correction():
    """Empty function
    """
    # it is performed per SN
    pass


Rband = dict([('g', 3.237), ('r', 2.176), ('i', 1.595), ('z', 1.217)])
def correct_for_absorption(flux, ebv, band):
    """Corrects for MW dust absorption.

    Keyword arguments:
    flux -- Python list containing fluxes.
    ebv -- float, MW E(B-V).
    band -- photometric band in which correction has to be applied.

    Returns:
    Python list containing values of flux corrected by absorption.

    Notes:
    From the `ebv' and `band' uses the correct R_band obtained from Schlafly & Finkbeiner 2011.

    general law A_band = R_band * E(B-V)

    A_band has to be turned in flux units.
    """

    a_mag = Rband[band] * ebv
    a_flux = 10**(0.4*a_mag)
    return [el*a_flux for el in flux]


def pre_process(data, bands):
    """
    Performs correction for time dilution and Galaxy reddening

    Keyword arguments:
    data -- data from observations. Supernova type.
    bands -- a list of bands to pre-precesse.

    Returns:
    data corrected
    """

    for b in bands:
        data.lcsDict[b].mjd = time_correct(data.lcsDict[b].mjd,
            data.zSpec if data.zSpec else data.zPhotHost)

        data.lcsDict[b].flux = correct_for_absorption(data.lcsDict[b].flux,
            data.MWEBV, b)

    return data




def pick_random_sn(catalog, band):
    """Extracts light curve in specified band of random supernova from catalogue (DEPRECATED: no more catalogue).

    Keyword arguments:
    catalog -- supernova catalog (class SupernovaCatalog).
    band -- photometric band identifying light curve.

    Returns:
    epoch -- NumPy array of epochs in MJD. Has zero point corresponding to MJD
             of maximum flux in r band.
    flux -- NumPy array of fluxes.
    fluxErr -- NumPy array of errors on flux measurements.
    idx -- index of supernova in the catalogue.
    """
    idx = np.random.random_integers(low=0, high=len(catalog.SNID))
    numObs = len(catalog.sne[idx].lcsDict[band].mjd)

    epoch = catalog.sne[idx].lcsDict[band].mjd
    # epoch = epoch - epoch[catalog.sne[idx].lcsDict['r'].flux.argmax()]
    epoch = epoch - epoch.min()

    flux = catalog.sne[idx].lcsDict[band].flux

    errFlux = catalog.sne[idx].lcsDict[band].fluxErr

    return epoch, flux, errFlux, idx


def redshift_distrib(pathToDir, binSize):
    """Plots the distribution of redshift values as an histogram with specified
    bin size.

    Keyword arguments:
    pathToDir -- directory in which supernova files are stored
    binSize -- the size of histogram bin.

    Returns:
    zed -- NumPy array containing redshifts.
    bins -- bins limits.
    """
    p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=pathToDir)

    lsDirData = p.stdout.read()
    lsDirData = lsDirData.split('\n')
    lsDirData.sort()
    lsDirData.remove('')

    zed = np.zeros(len(lsDirData), dtype=np.float)

    i = 0
    for z in np.nditer(zed, op_flags=['readwrite']):
        sn = get_sn_from_file(pathToDir+lsDirData[i])
        z[...] = sn.zSpec if sn.zSpec else sn.zPhotHost
        i += 1

    nBins = round((zed.max()-zed.min())/binSize)
    # print round(zed.max()*100)
    # print int(binSize*100)
    # bins = [range(0, int(zed.max()*100), int(binSize*100))[i]/100. for i in
    #     range(len(range(0, int(zed.max()*100), int(binSize*100))))]
    # bins.append(bins[-1]+binSize)
    # print bins
    plt.figure()
    bins = plt.hist(zed, bins=nBins, color='0.60', edgecolor='0.30')
    plt.xlabel('redshift z')
    plt.ylabel('number of observations')

    return zed, bins


def get_sn(catalog, band, idx):
    """Extract specified supernova observation in specified band from catalogue (DEPRECATED: no more catalogue).

    Keyword arguments:
    catalog -- supernova catalogue of class SupernovaCatalog.
    band -- string specifying photometric band to use.
    idx -- index of supernova inside catalogue.

    Returns:
    epoch -- NumPy array of epochs in MJD. Has zero point corresponding to MJD
             of maximum flux in r band.
    flux -- NumPy array of fluxes.
    fluxErr -- NumPy array of errors on flux measurements.
    """

    numObs = len(catalog.sne[idx].lcsDict[band].mjd)

    epoch = catalog.sne[idx].lcsDict[band].mjd
    # epoch = epoch - epoch[catalog.sne[idx].lcsDict['r'].flux.argmax()]
    epoch = epoch - epoch.min()

    flux = catalog.sne[idx].lcsDict[band].flux

    errFlux = catalog.sne[idx].lcsDict[band].fluxErr

    return epoch, flux, errFlux

def get_sn_from_file(pathToSN, magFlag=False):
    """Reads photometric data of SN from file formatted as in SNPhotCC

    Keyword arguments:
    pathToSN -- path to file from which extract data.

    Returns:
    sn -- object of class Supernova.
    """
    sn = classes.Supernova(pathToSN, magFlag)
    return sn

def reshape_for_GPy(vec):
    """Reshape input as a 'column' vector for input in GPy functions.

    Keyword arguments:
    vec -- NumPy array to be reshaped.

    Returns:
    Reshaped vec
    """
    return np.reshape(vec, (len(vec), 1))


def gp_fit(X, Y, errY, kernel,
           n_restarts=0, parallel=True,
           test_length=False,
           test_prior=False,
           verbose=False
           ):
    """Performs Gaussian Process regression using GPy functionalities.

    Keyword arguments:
    X -- list, independent values shaped as column array.
    Y -- list, dependent values shaped as column array.
    errY -- list, error values for independent variable, shaped as column array.
    kernel -- the GPy kernel to use.
    n_restarts -- number of restarts to optimise hyper parameters.
    parallel -- Flag whether to activate or not parallel computation in optimisation.
                This solves some memory leakage. The execution speed is not affected...
    test_length -- Flag whether to activate or not testing of random value of length scale hyper parameter.
    test_prior -- Flag whether to activate or not the use of a prior.
    verbose -- Flag whether to activate or not verbosity.

    Returns:
    predX -- list of X values at which an interpolation (prediction) has been made.
    predY -- list of Y values predicted using GPy.
    var -- list of variance values for each predY.
    gpModel -- the Gaussian process model from GPy use to predict predY.

    Notes:
    GPy code can be found at http://sheffieldml.github.io/GPy/.
    Check on shape of input should be added.
    """

    medXStep = np.median(np.abs(np.subtract(X[0:-2], X[1:-1])))
    maxXStep = X[-1] - X[0]
    if verbose:
        print "  Median X step {:<5.3f}".format(medXStep)
        print "  X range {:<5.3f}".format(maxXStep)

    rsX = reshape_for_GPy(X)
    if max(Y)/1000 > 1:
        print(max(Y))
        rsY = reshape_for_GPy([el/1000. for el in Y])
        print min(rsY)
    else:
        rsY = reshape_for_GPy(Y)

    origYMean = np.mean(Y)
    rsY = rsY - origYMean
    gpModel = GPy.models.GPHeteroscedasticRegression(rsX, rsY, kernel)

    if max(Y)/1000 > 1:
        gpModel['.*Gaussian_noise'] = [el/1000. for el in errY]
    else:
        gpModel['.*Gaussian_noise'] = errY

    # Block to capture unwanted output from .constrain_fixed() function
    # with Capturing() as output:
    [gpModel['.*Gaussian_noise_%s' %i].constrain_fixed(warning=False)
         for i in range(len(X))
         ]


    if test_length:
        np.random.RandomState
        length = np.random.uniform(low=medXStep, high=maxXStep)
        if verbose:
            print "  Randomized lengthscale {:<5.3f}".format(length)
        # with Capturing() as output:
        gpModel['.*lengthscale'].constrain_fixed(length, warning=False)
    elif test_prior:
        lenPrior = GPy.priors.Gamma(6, 0.01)
        gpModel['.*lengthscale'].set_prior(lenPrior, warning=False)
        if kernel.name == 'RatQuad':
            powPrior = GPy.priors.Gamma(1, 0.5)
            gpModel['.*power'].set_prior(powPrior, warning=False)
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

    predX = reshape_for_GPy(np.arange(min(X), max(X)+1, 1.))

    # _raw_predict is from GPy/core/gp.py
    meanY, var = gpModel._raw_predict(predX, full_cov=False)

    meanY = meanY + origYMean
    if max(Y)/1000 > 1:
        meanY = meanY*1000
        var = var * 1000
    return list(predX.reshape(predX.size)), \
        list(meanY.reshape(meanY.size)), list(var.reshape(var.size)), gpModel

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

    dataSN = "train_data/SIMGEN_PUBLIC_DES"

    if args.candidate is None:
        # Picking random candidate
        #
        # high set max number of SN in SNPhotCC
        candidateIdx = np.random.random_integers(
            low=0, high=21319)
        print candidateIdx
        args.candidate = np.genfromtxt(
            dataSN+"SIMGEN_PUBLIC_DES.LIST", dtype=None)[candidateIdx]

        # Setting path and getting data
        pathToSN = dataSN + args.candidate
    else:
        pathToSN = dataSN + \
                    "DES_SN" + "{:>06}".format(args.candidate) + ".DAT"
    sn = get_sn_from_file(pathToSN, args.mag)

    epoch = sn.lcsDict[args.band].mjd
    flux = sn.lcsDict[args.band].flux
    errFlux = sn.lcsDict[args.band].fluxErr

    print "  Candidate ID          {:>06}".format(args.candidate)
    print "  Testing lengthscale ? {:<5}".format(args.testLength)
    print "  Magnitudes ?          {:<5}".format(args.mag)

    # Tranforming to magnitudes
    #
    # this is done by reading them from file

    # if args.mag:
    #     limFlux = mag_to_flux(limMag[args.band])
    #     mag = flux_to_mag(flux, limFlux)
    #     errMag = flux_error_to_mag_error(errFlux, flux)

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
    if not sn.lcsDict[args.band].badCurve:
        # if args.mag:
        #     predEpoch, mu, var, GPModel = gp_fit(
        #                                     epoch, mag, errMag,
        #                                     kern, n_restarts=10,
        #                                     test_length=args.testLength,
        #                                     test_prior=args.testPrior,
        #                                     verbose=args.verbose)
        # else:
        predEpoch, mu, var, GPModel = gp_fit(
                                        epoch, flux, errFlux,
                                        kern, n_restarts=10,
                                        test_length=args.testLength,
                                        test_prior=args.testPrior,
                                        verbose=args.verbose)

        zed = sn.zSpec if sn.zSpec else sn.zPhotHost
        corr_epoch = time_correct(epoch, zed)
        corr_flux = correct_for_absorption(flux, sn.MWEBV, args.band)

        corr_predEpoch, corr_mu, corr_var, corr_GPModel = gp_fit(
                                        corr_epoch, corr_flux, errFlux,
                                        kern, n_restarts=10,
                                        test_length=args.testLength,
                                        test_prior=args.testPrior,
                                        verbose=args.verbose)


        # print GPModel['.*lengthscale|.*power']

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

        #GPModel.plot_f(fignum=figNum)
        corr_GPModel.plot_f(fignum=figNum)

        if args.mag:
            print 'mags'
            ylim = plt.ylim()
            plt.ylim(ylim[1], ylim[0])
            plt.errorbar(epoch, mag,
                 yerr=errMag, fmt=None, ecolor='black', zorder=1)
            plt.scatter(epoch, mag, color='black')
        else:
            fig, ax = plt.subplots(nrows=2, ncols=1,
                figsize=(16.5, 11.7),
                tight_layout=False
                )
            fig.suptitle("{:>06}".format(args.candidate))
            print 'fluxes'
            ax[0].errorbar(epoch, flux,
                yerr=errFlux, fmt=None, ecolor='black', zorder=1)
            ax[0].scatter(epoch, flux, color='black')

            ax[1].errorbar(corr_epoch, corr_flux,
                yerr=errFlux, fmt=None, ecolor='red', zorder=1)
            ax[1].scatter(corr_epoch, corr_flux, color='red')

        print "  The process took {:5.3f} secs.".format(time.time()-start_time)
        plt.show()
    else:
        print "  This was flagged as BAD CURVE!"
