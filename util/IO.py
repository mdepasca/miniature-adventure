import numpy as np
import pandas as pd
import subprocess
import os
import classes


def get_sn_from_file(pathToSN, magFlag=False):
    """Reads photometric data of SN from file formatted as in SNPhotCC

    Keyword arguments:
    pathToSN -- path to file from which extract data.

    Returns:
    sn -- object of class Supernova.
    """
    sn = classes.Supernova(pathToSN, magFlag)
    return sn

def get_fit_from_file(pathToFit, magFlag=False):
    """Reads photometric fit of SN from file formatted as in SNPhotCC

    Keyword arguments:
    pathToSN    --  path to file from which to extract data.

    Returns:
    fit  --  object of class SupernovaFit.
    """

    tmp = get_sn_from_file(pathToFit, magFlag)
    """
    Initializing SupernovaFit object
    """
    fit = classes.SupernovaFit(tmp, tmp.kern if hasattr(tmp, 'kern') else None)

    for b in tmp.lcsDict.keys():
        fit.set_lightcurve(b,
            tmp.lcsDict[b].mjd,
            tmp.lcsDict[b].flux,
            tmp.lcsDict[b].fluxErr,
            magFlag=magFlag
            )

    del(tmp)

    if fit.r.badCurve == False:
        fit.shift_mjds()
    """
    Fixing shiftedMjd for not-peaked LCs
    """
    if (fit.peaked == False) and (fit.r.badCurve == False) :
        """
        correcting using CC results
        """
        for b in fit.lcsDict.keys():
            fit.lcsDict[b].shiftedMjd = [
            el + fit.ccMjdMaxFlux for el in fit.lcsDict[b].shiftedMjd
            ]

    return fit

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
    """Creates files dividing supernovae in training and test sets.
    It creates also files list training set supernovae by type

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

    badCount = 0

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
        if tmpSN.r.badCurve:
            badCount += 1
            continue
        SNType = tmpSN.SNTypeInt
        if SNType != -9:
            outFileTrain.write(
                "{:0>5d}      {:0>6d}   {:>}\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
        else:
            outFileTest.write(
                "{:0>5d}      {:0>6d}   {:>}\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 1:
            outFileIa.write(
                "{:0>5d}      {:0>6d}   {:>}   snIa\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 21 or SNType == 22 or SNType == 23:
            outFileII.write(
                "{:0>5d}      {:0>6d}   {:>}   snII\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 3 or SNType == 32 or SNType == 33:
            outFileIbc.write(
                "{:0>5d}      {:0>6d}   {:>}   snIbc\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 11:
            outFileIaPec.write(
                "{:0>5d}      {:0>6d}   {:>}   pec\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == 66:
            outFileOther.write(
                "{:0>5d}      {:0>6d}   {:>}   other\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        if SNType == -1:
            outFileIbc.write(
                "{:0>5d}      {:0>6d}   {:>}   snIbc\n".format(
                    i-badCount, tmpSN.SNID, path+lsList[i]
                    )
                )
            continue

        outFileOther.write(
            "{:0>5d}      {:0>6d}   {:>}   other\n".format(
                i-badCount, tmpSN.SNID, path+lsList[i]
                )
            )
    print badCount
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
