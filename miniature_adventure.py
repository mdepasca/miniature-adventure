import argparse
import os
from  os import path
import subprocess
import sys
import socket
import time
# garbage collector
import gc

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

import GPy

import classes as cls
import utilities as util
from utilities import bcolors

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
# Activate automatic conversion of ndarray to R objects
ro.conversion.py2ri = numpy2ri

from progressbar import ProgressBar, SimpleProgress, ETA, Percentage, Bar, \
                        AnimatedMarker, Timer, Counter

if __name__ == "__main__":
    # gc.set_debug(gc.DEBUG_LEAK)
    # Parsing input from command line
    parser = argparse.ArgumentParser(
        description = "SN lightcurve fitter and classifier.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--fit", dest="fit",
        action="store_true",
        help="Fit lightcurves with Gaussian processes method."
        )

    group.add_argument(
        "--fit-training", dest="fitTraining",
        action="store_true",
        help="Fit lightcurves from the training set \
        with Gaussian processes method."
        )
    
    parser.add_argument(
        "--cross-correlation", dest="crossCor",
        action="store_true",
        help="Performs cross correlation between non peaked lcs (with maximum in \
            r-band at one of the MJD extremes) and all the peaked lcs. Produces \
            an estimate for maximum in r-band. VERY TIME CONSUMING."
        )

    parser.add_argument(
        "--distance-matrix", dest="distMatrix",
        action="store_true",
        help="Calculate distance between fitted lightcurves in same band. \
        It is use to build a diffusion map (see Coifman & Lafon (2006) \
        and Lafon & Lee (2006)).")

    parser.add_argument(
        "--diffuse", dest="diffuse",
        action="store_true",
        help="Computes the diffusion map coefficients. Run together or after \
        --distance-matrix option. Uses `diffusionMap` R package developed \
        by Joseph Richards.")

    parser.add_argument(
        "--train", dest="train",
        action="store_true",
        help="Train the classifier - Random Forest. Uses `randomForest` R \
        package.")

    parser.add_argument(
        "--classify", dest="classify",
        action="store_true")

    parser.add_argument(
        "--training-directory", dest="dirData",
        default="train_data" + os.sep + "DES_BLIND+HOSTZ",
        help="Path to directory containing training data.")

    parser.add_argument(
        "--fit-directory", dest="dirFit",
        default="train_data" + os.sep + "DES_BLIND+HOSTZ_FIT",
        help="Path to directory containing fitted data.")

    parser.add_argument(
        "--fit-file", dest="fitFile",
        help="Path to file in which to dump fitting results.")

    parser.add_argument(
        "--plot", dest="plot",
        action="store_true",
        help="Save on `pdf` file the plot of fitting curve over data.")

    parser.add_argument(
        "-f", "--file",
        help="")

    args = parser.parse_args()


    bands = ['g', 'r', 'i', 'z']
else:
    pass
    
if __name__ == "__main__":
    start = 0
    stop = 5
    indent = "          "
    prodDir = "products"+os.sep
    os.system("clear")
    peakIdx = np.empty(0)
    nopeakIdx = np.empty(0)
    if not os.path.exists(path.abspath(args.dirFit)):
        os.makedirs(path.abspath(args.dirFit))
    
    KERN_RATQUAD = "RatQuad"
    # modelList = np.zeros(0, dtype=np.object)

    # should be possible to change the next two variables
    # args.dirData = "train_data" + os.sep + "DES_BLIND+HOSTZ"
    # args.dirFit = "fit_data" + os.sep
    fNameCandidatesList = "DES_BLIND+HOSTZ.LIST"

    print indent + bcolors.bldpur + "* * * * * * * * * * * * * * *"
    print indent + "*    Miniature Adventure    *"
    print indent + "*    -------------------    *"
    print indent + "*    lightcurves fitting    *"
    print indent + "*             and           *"
    print indent + "*      SN classification    *"
    print indent + "* * * * * * * * * * * * * * *" + bcolors.txtrst
    
    start_time = time.time()

    """

    PERFORMS LCs FITTING

    """
    if args.fit or args.fitTraining:
        whileOn = True
        i = 0
        filePath = prodDir + 'PEAKED_{:<}_{:<5.3f}.LIST'.format(
            socket.gethostname(), time.time()
            )
   
        fPeaked = open(filePath, 'w')

        whileOn = True
        i = 0
        filePath = prodDir + 'NOPEAKED_{:<}_{:<5.3f}.LIST'.format(
            socket.gethostname(), time.time()
            )
  
        fNopeaked = open(filePath, 'w')


        # Relevant input data
        if args.fit:
            print "\n" + indent + "[1] * Fit lightcurves ..."
        if args.fitTraining:
            print "\n" + indent + "[1] * Fit lightcurves from training set ..."

        print "\n" + indent + \
            "Data directory: " + os.curdir + args.dirData + os.sep
        print "\n" + indent + "List of candidates contained in: " \
            + os.curdir + args.dirData + os.sep + fNameCandidatesList

        vecCandidates = np.genfromtxt(
                args.dirData+os.sep+fNameCandidatesList, dtype=None)
        # tenPercent = vecCandidates.size / 10
        # const = 0
        print "\n" + indent \
            + "Number of candidates = {:<d}".format(vecCandidates.size)

        print "\n" + indent \
            + "Data are fitted using GP with Radial Basis Function kernel."

        # kern = GPy.kern.RatQuad(1)
        kern = GPy.kern.RBF(1)
    
        # Redirecting stderr output to file
        # saveErr = sys.stderr
        # ferr = open('error.log', 'w')
        # sys.stderr = ferr

        # Fitting single lightcurves 
        #
        # THIS PIECE NEEDS TO BE PARALLELIZED
        # 
        # optimize_restarts parallel using multiprocessing
        
        for i in range(start, stop):
            candidate = util.get_sn_from_file(
                args.dirData + os.sep + vecCandidates[i]
                )

            if args.fitTraining:
                # if SN type is not provided, skip to the next item
                if candidate.SNTypeInt == -9:
                    continue
                else:
                    print indent + 'SN type code {:<}'.format(candidate.SNTypeInt) 

            # candidateFit = cls.SupernovaFit(candidate.SNID)
            candidateFit = cls.SupernovaFit(candidate)
            for b in candidate.lcsDict.keys(): 
                epoch = util.time_correct(candidate.lcsDict[b].mjd, 
                    candidate.zSpec if candidate.zSpec else candidate.zPhotHost)

                flux = util.correct_for_absorption(candidate.lcsDict[b].flux, 
                    candidate.MWEBV, b)
                errFlux = candidate.lcsDict[b].fluxErr

                # test_prior should be deleted as option. Prior too weak.
                # 
                # Fitting Lightcurve

                """
                K--correction reduces the number of observation in each band.
                For this reason it has to be performed here, before checking
                flux length.

                ---> K--correction
                """

                if (candidate.lcsDict[b].badCurve) or (len(flux) <= 3):
                    candidateFit.lcsDict[b].badCurve = True
                    print indent + bcolors.FAIL + \
                        "{:<} {:<} {:<}".format(i, candidate.SNID, b) + \
                        bcolors.txtrst
                    continue

                # if (not candidate.lcsDict[b].badCurve) and \
                #     (flux.size >= 3):
                saveOut = sys.stdout
                fout = open('out.log', 'w')
                # fout = open('/dev/null', 'w')
                sys.stdout = fout

                
                predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                epoch, flux, errFlux, 
                                                kern, n_restarts=10, 
                                                parallel=False, # this solves some memory leakage. The execution speed is not affected...
                                                test_length=True)
                sys.stdout = saveOut
                fout.close()

                candidateFit.set_lightcurve(b, 
                    predMjd, predFlux, predErr)

                print indent + \
                    "{:<} {:<} {:<}".format(i, candidate.SNID, b)
                                
            # Setting epoch 0 point to epoch or r maximum
            if candidateFit.r.badCurve is False:
                # candidateFit.shift_mjds()
                filePath = args.dirFit + os.sep + \
                    path.splitext(vecCandidates[i])[0] + "_FIT.DAT"
                candidateFit.save_on_txt(filePath)
                
                if candidateFit.peaked:
                    peakIdx = np.append(peakIdx, i)
                    fPeaked.write('{:<}\n'.format(filePath))
                else:
                    nopeakIdx = np.append(nopeakIdx, i)
                    fNopeaked.write('{:<}\n'.format(filePath))

        fPeaked.close()
        fNopeaked.close()
        # sys.stderr = saveErr
        # ferr.close()
        whileOn = True
        i = 0
        filePath = 'peaked_{:<}_{:<5.3f}.dat'.format(
            socket.gethostname(), time.time()
            )

        np.savetxt(prodDir + filePath, peakIdx,
            header='Indexes of fitted LCs with r maximum.', fmt='%d')

        whileOn = True
        i = 0
        filePath = prodDir + 'nopeaked_{:<}_{:<5.3f}.dat'.format(
            socket.gethostname(), time.time()
            )

        np.savetxt(filePath, nopeakIdx,
            header='Indexes of fitted LCs without an r maximum.', fmt='%d')

    """

    PERFORMING CROSS-CORRELATION 

    """

        
    if args.crossCor:
        """
        getting file list from directory
        File are sorted by SNID.
        In the following peakIdx and nopeakIdx contain index referring to the 
        full list of files. For this reason the list of files it is queried on
        dirData. It is then filtered using the above variables.
        """
        start = 0
        end = 2
        print "\n" + indent + bcolors.undwht + \
            "(*) Calculate cross-correlation of not peaked- with " + \
            "peaked-lcs ..." +  bcolors.txtrst

        print "\n" + indent + "Interval [{:<},{:<})".format(start, end)

        p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirData+os.sep)
            # cwd=args.dirFit+os.sep)
        lsDirData = p.stdout.read()
        lsDirData = lsDirData.split('\n')
        lsDirData.sort()
        lsDirData.remove('')
        
        filePath = 'peaked.dat'.format(socket.gethostname())
        peakIdx = np.loadtxt(prodDir + filePath, dtype=np.int)
        filePath = 'nopeaked.dat'.format(socket.gethostname())
        tmp = np.loadtxt(prodDir + filePath, dtype=np.int)
        if tmp.size == 1:
            nopeakIdx = np.asarray([tmp])
        else:    
            nopeakIdx = np.asarray(tmp)
        
        filePath = 're-written_files_{:<5.3f}.dat'.format(time.time())
        reWrite = open(prodDir + filePath, 'w')
        prog = 0        
        for i in nopeakIdx[start:end]:

            z = 0 # goes on peakIdx to index the progress bar
            
            """
            READ DATA FROM FILE 
            in Supernova object
            """
            filePath = args.dirFit + os.sep + lsDirData[i][0:12] + '_FIT.DAT'
            try:
                tmpSN = util.get_sn_from_file(filePath)

                print "Progress: {:<d}".format(prog)
                prog += 1

                ccIndent = "ID:{: ^7d}".format(tmpSN.SNID)#  "          "
                widgets = [ccIndent, Percentage(), ' ',
                   Bar(marker='#',left='[',right=']'),
                   ' ', ETA()]
                pbar = ProgressBar(widgets=widgets, maxval=len(peakIdx)).start()
            except IOError:
                continue
            
            if tmpSN.r.badCurve:
                continue
            """
            create SupernovaFit object
            """
            notPeaked = cls.SupernovaFit(tmpSN)
            for l in tmpSN.lcsDict.keys():
                notPeaked.set_lightcurve(l, 
                    tmpSN.lcsDict[l].mjd,
                    tmpSN.lcsDict[l].flux,
                    tmpSN.lcsDict[l].fluxErr
                    )
            """
            Shifting mjds in not-peaked
            """
            notPeaked.shift_mjds()

            ccMax = list()#np.zeros(peakIdx.size)
            k = 0 # goes on ccMax
            for j in peakIdx:
                """
                READ DATA FROM FILE
                """
                filePath = args.dirFit + os.sep + lsDirData[j][0:12] + '_FIT.DAT'
                try:
                    tmpSN = util.get_sn_from_file(filePath)
                except IOError:
                    continue
                
                if tmpSN.r.badCurve:
                    continue
                peaked = cls.SupernovaFit(tmpSN)
                for l in tmpSN.lcsDict.keys():
                    peaked.set_lightcurve(l,
                        tmpSN.lcsDict[l].mjd,
                        tmpSN.lcsDict[l].flux, 
                        tmpSN.lcsDict[l].fluxErr
                        )
                """
                Shifting mjds in peaked
                """
                peaked.shift_mjds()

                """
                Performing cross-correlation
                """
                ycorr = signal.correlate(
                    notPeaked.normalized_flux('r'),
                    peaked.normalized_flux('r')
                    )
                xcorr = np.arange(ycorr.size)
                lags = xcorr - (
                    len(notPeaked.normalized_flux('r'))-1
                    )
                distancePerLag = (
                    notPeaked.r.shiftedMjd[-1] - \
                    notPeaked.r.shiftedMjd[0])/float(
                                        len(notPeaked.r.shiftedMjd)
                                        )
                offsets = -lags*distancePerLag
                # raise SystemExit
                # ccMax[k] = offsets[np.argmax(ycorr)]
                ccMax.append(offsets[np.argmax(ycorr)])
                # k += 1
                
                pbar.update(z+1)
                z += 1

            notPeaked.ccMjdMaxFlux = np.mean(ccMax)#ccMax.mean()
            """
            re-writing file of not peaked lc to include information on maximum
            position from CC.
            """
            filePath = args.dirFit + os.sep + lsDirData[i][0:12] + '_FIT.DAT'
            notPeaked.save_on_txt(filePath)

            reWrite.write(filePath+'\n')
            pbar.finish()
        reWrite.close()
        print 'CC ended!'
        

    """

    CALCULATING DISTANCE MATRIX

    """
    if args.distMatrix:
        """
        Calculate distance between fitted lightcurves.
        Distance values are saved in a R matrix. This will be used by the R 
        package `diffusionMap` through rpy2 Python package.
        """
        j_offset = 0
        i_start = 0
        i_end = 100
        j_start = i_start + j_offset
        j_end = i_end + j_offset
        print "\n" + indent + bcolors.undwht + \
            "(*) Calculate distances between lightcurves ..." + \
            bcolors.txtrst
        print indent + "Rows in [{:<d}, {:<d})".format(i_start, i_end)
        print indent + "Cols in [{:<d}, {:<d})".format(j_start, j_end)

        p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirFit+os.sep)
            # cwd=args.dirFit+os.sep)
        lsDirFit = p.stdout.read()
        lsDirFit = lsDirFit.split('\n')
        lsDirFit.sort()
        lsDirFit.remove('')

        """
        setting value for big distance
        """
        distFlag = 5
        nullVal = -999

        bandDict = {
            'g':0,
            'r':1,
            'i':2,
            'z':3
            }
        widgets = [indent, 'Processing:', ' ', Counter(), ' ', 
            AnimatedMarker(), indent, Timer()]
        
        # creating numpy matrix
        Pymatrix = np.zeros((4, 
            len(lsDirFit[i_start:i_end]), len(lsDirFit[i_start:i_end])),
            dtype=float
            )

        pbar = ProgressBar(widgets=widgets, maxval=(i_end-i_start)).start()

        for i in range(i_start, i_end):

            """
            Reading in i-candidate
            """
            tmpSN = util.get_sn_from_file(
                args.dirFit+os.sep+lsDirFit[i]
                )
            if tmpSN.r.badCurve:
                for b in bands:
                        Pymatrix[bandDict[b], i-i_start, j-j_start] = nullVal
                print "{:<} Has bad curve in r band - ".format(lsDirFit[i]) + \
                    "THE FILE HAS TO BE DELETED"
                continue
                
            iCandidate = cls.SupernovaFit(tmpSN)
            
            for b in tmpSN.lcsDict.keys():
                # set_lightcurve set also if the lc is peaked or not
                iCandidate.set_lightcurve(b, 
                    tmpSN.lcsDict[b].mjd,
                    tmpSN.lcsDict[b].flux,
                    tmpSN.lcsDict[b].fluxErr
                    )

            """
            Shifting mjds in i-candidate
            """
            iCandidate.shift_mjds()
            if iCandidate.peaked == False:
                # print i, iCandidate.SNID
                """
                keeping to perform check with other non peaked LC
                """
                iElMax = iCandidate.r.shiftedMjd.index(0.)
                """
                correcting using CC results
                """
                for b in bands:
                    iCandidate.lcsDict[b].shiftedMjd = [
                        iCandidate.lcsDict[b].shiftedMjd[l] + 
                        iCandidate.ccMjdMaxFlux for l in range(len(
                            iCandidate.lcsDict[b].shiftedMjd
                            ))
                    ]
            iElSize = iCandidate.r.size
            iPeaked = iCandidate.peaked

            for j in range(j_start, j_end):
                """
                if this SN has badCurve in this band it will be far from all 
                the others by default.
                here will save time from not opening all the other files 
                to create new SupernovaFit objcets.
                """

                if j == i:
                    # filling elements on the distance matrix diagonal
                    for b in bands:
                        Pymatrix[bandDict[b], i-i_start, j-j_start] = 0.
                    continue

                if j < i:
                    # filling matrix elements below the diagonal
                    for b in bands:
                        Pymatrix[bandDict[b], i-i_start, j-j_start] = \
                                    Pymatrix[bandDict[b], j-j_start, i-i_start]
                    continue # jump to the next iteration of the loop

                """
                Reading in j-candidate
                """
                try:
                    tmpSN = util.get_sn_from_file(
                        args.dirFit+os.sep+lsDirFit[j]
                    )
                except IndexError:
                    print j, len(lsDirFit)
                    raise IndexError("list index out of range")
                if tmpSN.r.badCurve:
                    for b in bands:
                        Pymatrix[bandDict[b], i-i_start, j-j_start] = nullVal
                    print "{:<} Has bad curve in r band -".format(lsDirFit[j])+\
                        " THE FILE HAS TO BE DELETED"
                    continue

                jCandidate = cls.SupernovaFit(tmpSN)
                for b in tmpSN.lcsDict.keys():
                    jCandidate.set_lightcurve(b, 
                        tmpSN.lcsDict[b].mjd, 
                        tmpSN.lcsDict[b].flux, 
                        tmpSN.lcsDict[b].fluxErr
                        )

                """
                Shifting mjds in j-candidate
                """
                jCandidate.shift_mjds()
                if jCandidate.peaked == False:
                    """
                    keeping to perform check with other non peaked LC
                    """
                    jElMax = jCandidate.r.shiftedMjd.index(0.)
                    """
                    correcting using CC results
                    """
                    for b in bands:
                        jCandidate.lcsDict[b].shiftedMjd = [
                            jCandidate.lcsDict[b].shiftedMjd[l] + 
                            jCandidate.ccMjdMaxFlux for l in range(len(
                                jCandidate.lcsDict[b].shiftedMjd
                                ))
                        ]

                jElSize = jCandidate.r.size

                for b in bands:
                    if not jCandidate.lcsDict[b].badCurve \
                    and not iCandidate.lcsDict[b].badCurve:
                        Pymatrix[bandDict[b], i-i_start, j-j_start] = \
                            iCandidate.get_distance(jCandidate, b)
                    else:
                        # in case of bad curve
                        """
                        This works like a flag. These elements will be set 
                        equal to a neutral value (the mean of the other)
                        """
                        Pymatrix[bandDict[b], i-i_start, j-j_start] = distFlag

            pbar.update(i-i_start+1)
        pbar.finish()

        # fixing flagged elements
        # raise SystemExit
        if Pymatrix[0, Pymatrix[0] == distFlag].size > 0: 
            ind = np.where(Pymatrix[0] == distFlag)
            Pymatrix[0, ind[0], ind[1]] = np.add(
                np.add(
                    Pymatrix[1, ind[0], ind[1]], 
                    Pymatrix[2, ind[0], ind[1]]
                    ), 
                Pymatrix[3, ind[0], ind[1]]
                )/3.


        if Pymatrix[1, Pymatrix[1] == distFlag].size > 0:
            ind = np.where(Pymatrix[1] == distFlag)
            # Pymatrix[1, ind[0], ind[1]] = Pymatrix[1,:,:].max()
            Pymatrix[1, ind[0], ind[1]] = np.add(
                np.add(
                    Pymatrix[0, ind[0], ind[1]], 
                    Pymatrix[2, ind[0], ind[1]]
                    ), 
                Pymatrix[3, ind[0], ind[1]]
                )/3.

        if Pymatrix[2, Pymatrix[2] == distFlag].size > 0: 
            ind = np.where(Pymatrix[2] == distFlag)
            # Pymatrix[2, ind[0], ind[1]] = Pymatrix[2].max()
            Pymatrix[2, ind[0], ind[1]] = np.add(
                np.add(
                    Pymatrix[0, ind[0], ind[1]], 
                    Pymatrix[1, ind[0], ind[1]]
                    ), 
                Pymatrix[3, ind[0], ind[1]]
                )/3.

        if Pymatrix[3, Pymatrix[3] == distFlag].size > 0: 
            ind = np.where(Pymatrix[3] == distFlag)
            # Pymatrix[3, ind[0], ind[1]] = Pymatrix[3].max()
            Pymatrix[3, ind[0], ind[1]] = np.add(
                np.add(
                    Pymatrix[0, ind[0], ind[1]], 
                    Pymatrix[1, ind[0], ind[1]]
                    ), 
                Pymatrix[2, ind[0], ind[1]]
                )/3.
        
        PymatrixSum = np.sum(Pymatrix, 0)
        """
        Saving on text files
        """
        fileHeader = "# Pymatrix[{:<d}:{:<d},{:<d}:{:<d}]  --- ".format(
            i_start, i_end, j_start, j_end
            ) + \
            "Created by {:<}".format(socket.gethostname())

        filePath = prodDir + 'distance_matrix' + os.sep + \
            'dist_matrix_g_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, Pymatrix[0], fmt='%6.4f', header=fileHeader)

        filePath = prodDir + 'distance_matrix' + os.sep + \
            'dist_matrix_r_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, Pymatrix[1], fmt='%6.4f', header=fileHeader)

        filePath = prodDir + 'distance_matrix' + os.sep + \
            'dist_matrix_i_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, Pymatrix[2], fmt='%6.4f', header=fileHeader)

        filePath = prodDir + 'distance_matrix' + os.sep + \
            'dist_matrix_z_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, Pymatrix[3], fmt='%6.4f', header=fileHeader)


        filePath = prodDir + 'distance_matrix' + os.sep + \
            'dist_matrix_Sum_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, PymatrixSum, fmt='%6.4f', header=fileHeader)

    """

    CALCULATING DIFFUSION MAP

    """


    if args.diffuse:
        if 'diffusionMap' not in globals():
            diffusionMap = importr('diffusionMap')

        ndim = ro.r.attributes(Rmatrix)[0][0]
        dmap = diffusionMap.diffuse(Rmatrix, neigen=5)
        util.dump_pkl('diffusion_map.pkl', dmap)
        

    """

    TRAINING RANDOM FOREST CLASSIFIER

    """

    if args.train:
        randomForest = importr('randomForest')
        if 'dmap' not in globals():
            print indent + 'Loading catalog from dump file ...'
            dmap = util.open_pkl('tmp_diffusion_map.pkl')

        dmap_rf = randomForest.randomForest(dmap)
        

    """

    PLOT OBSERVATION AND FIT
    --plot
    """

    if args.plot:
        """
        getting file list from directory
        File will be sorted by SNID
        """
        if "lsDirFit" not in globals():
            p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
                cwd=args.dirFit+os.sep
                )
            lsDirFit = p.stdout.read()
            lsDirFit = lsDirFit.split('\n')
            lsDirFit.sort()
            lsDirFit.remove('')

        if 'catalog' not in globals():
            vecCandidates = np.genfromtxt(
                args.dirData+os.sep+fNameCandidatesList, dtype=None)
           

        print indent + 'Plotting ...'
        nrows = 5
        ncols = 5
        offset = 25
        fig_g, ax_g = plt.subplots(nrows=nrows, ncols=ncols, 
                    figsize=(16.5, 11.7), 
                    tight_layout=True
                    )
        fig_r, ax_r = plt.subplots(nrows=nrows, ncols=ncols, 
                    figsize=(16.5, 11.7), 
                    tight_layout=True
                    )
        fig_i, ax_i = plt.subplots(nrows=nrows, ncols=ncols, 
                    figsize=(16.5, 11.7), 
                    tight_layout=True
                    )
        fig_z, ax_z = plt.subplots(nrows=nrows, ncols=ncols, 
                    figsize=(16.5, 11.7), 
                    tight_layout=True
                    )

        dictFig = {'g':fig_g, 
                   'r':fig_r,
                   'i':fig_i,
                   'z':fig_z}

        dictAx = {'g':ax_g, 
                  'r':ax_r,
                  'i':ax_i,
                  'z':ax_z}

        r = {'g':0,
             'r':0,
             'i':0,
             'z':0}

        c = {'g':0,
             'r':0,
             'i':0,
             'z':0}

        for b in dictFig.keys():
            dictFig[b].subplots_adjust(top=0.8)
            dictFig[b].suptitle('band {:<1}'.format(b))

        for i in range(nrows*ncols):
            # getting the data from file
            candidate = util.get_sn_from_file(
                args.dirData + os.sep + vecCandidates[offset+i])

            """
            reading fit data from file
            """
            tmpSN = util.get_sn_from_file(
                    args.dirFit+os.sep+lsDirFit[offset+i]
                    )
            """
            Initializing SupernovaFit object
            """
            fit = cls.SupernovaFit(tmpSN)
            for l in tmpSN.lcsDict.keys():
                fit.set_lightcurve(l,
                    tmpSN.lcsDict[l].mjd,
                    tmpSN.lcsDict[l].flux, 
                    tmpSN.lcsDict[l].fluxErr
                    )
            fit.shift_mjds()
            """
            Fixing shiftedMjd for not-peaked LCs
            """
            if fit.peaked == False:
                """
                correcting using CC results
                """
                for b in bands:
                    fit.lcsDict[b].shiftedMjd = [
                        fit.lcsDict[b].shiftedMjd[l] + 
                        fit.ccMjdMaxFlux for l in range(len(
                            fit.lcsDict[b].shiftedMjd
                            ))
                    ]

            for b in dictAx.keys():
                data = candidate.lcsDict[b]

                """
                Fixing shiftedMjd for not-peaked LCs
                """
                # if fit.peaked == False:
                #     fit.lcsDict[b].shiftedMjd = np.ma.add(
                #         fit.lcsDict[b].shiftedMjd, 
                #         fit.ccMjdMaxFlux
                #         )
                fit_b = fit.lcsDict[b]

                fit_r = fit.lcsDict['r']
                
                if c[b] > nrows-1:
                    c[b] = 0
                    r[b] += 1

                xlim = dictAx[b][r[b], c[b]].get_xlim()
                ylim = dictAx[b][r[b], c[b]].get_ylim()
                if not data.badCurve:
                    if fit.peaked:
                        data.set_shifted_mjd(
                            fit_r.mjd[fit_r.max_flux_index])
                    else:
                        data.set_shifted_mjd(
                            fit_r.mjd[fit_r.max_flux_index])
                        # for b in bands:
                        data.shiftedMjd = [
                            data.shiftedMjd[l] + 
                            fit.ccMjdMaxFlux for l in range(len(
                                data.shiftedMjd
                                ))
                        ]
                        # data.shiftedMjd += fit.ccMjdMaxFlux

                    bottom = min(data.flux) - np.median(data.fluxErr)
                    up = max(data.flux) + np.median(data.fluxErr)
                    dictAx[b][r[b], c[b]].set_ylim(bottom, up)

                    fluxUpLim = [val for val in [
                        fit_b.flux[i] + fit_b.fluxErr[i] 
                            for i in range(len(fit_b.flux))
                        ]]
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fluxUpLim, fit_b.flux, 
                        where=(np.array(fit_b.flux)
                            +np.array(fit_b.fluxErr))>np.array(fit_b.flux),
                        facecolor='red', alpha=0.4, linewidth=0.5)
                    fluxLowLim = [val for val in [
                        fit_b.flux[i] - fit_b.fluxErr[i] 
                            for i in range(len(fit_b.flux))
                        ]]
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fluxLowLim, fit_b.flux, 
                        where=(np.array(fit_b.flux)
                            -np.array(fit_b.fluxErr))<np.array(fit_b.flux),
                        facecolor='red', alpha=0.4, linewidth=0.5)

                    fluxUpLim = [val for val in [
                        fit_b.flux[i] + 2*fit_b.fluxErr[i] 
                            for i in range(len(fit_b.flux))
                        ]]
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fluxUpLim, fit_b.flux, 
                        where=(np.array(fit_b.flux)
                            +2*np.array(fit_b.fluxErr))>np.array(fit_b.flux)
                                                    +np.array(fit_b.fluxErr),
                        facecolor='red', alpha=0.2, linewidth=0.5)
                    fluxLowLim = [val for val in [
                        fit_b.flux[i] - 2*fit_b.fluxErr[i] 
                            for i in range(len(fit_b.flux))
                        ]]
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fluxLowLim, fit_b.flux, 
                        where=(np.array(fit_b.flux)
                            -2*np.array(fit_b.fluxErr))<np.array(fit_b.flux)
                                                    +np.array(fit_b.fluxErr),
                        facecolor='red', alpha=0.2, linewidth=0.5)

                    fluxUpLim = [val for val in [
                        fit_b.flux[i] + 3*fit_b.fluxErr[i] 
                            for i in range(len(fit_b.flux))
                        ]]
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fluxUpLim, fit_b.flux, 
                        where=(np.array(fit_b.flux)
                            +3*np.array(fit_b.fluxErr))>np.array(fit_b.flux)
                                                    +2*np.array(fit_b.fluxErr),
                        facecolor='red', alpha=0.1, linewidth=0.5)
                    fluxLowLim = [val for val in [
                        fit_b.flux[i] - 3*fit_b.fluxErr[i] 
                            for i in range(len(fit_b.flux))
                        ]]
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fluxLowLim, fit_b.flux, 
                        where=(np.array(fit_b.flux)
                            -3*np.array(fit_b.fluxErr))<np.array(fit_b.flux)
                                                    -2*np.array(fit_b.fluxErr),
                        facecolor='red', alpha=0.1, linewidth=0.5)

                    dictAx[b][r[b], c[b]].plot(fit_b.shiftedMjd, fit_b.flux, 
                        color='#7f0000', 
                        linewidth=2)

                    dictAx[b][r[b], c[b]].scatter(data.shiftedMjd, data.flux, 
                        s=10, label=str(candidate.SNID), c='black', marker='x')

                    dictAx[b][r[b], c[b]].errorbar(data.shiftedMjd, data.flux,
                        data.fluxErr, fmt=None, color='black', ecolor='black')

                    
                    if not fit.peaked:
                        pass

                    dictAx[b][r[b], c[b]].legend(
                        loc='best', framealpha=0.3, fontsize='10')
                else:
                    label = str(candidate.SNID)+" BAD CURVE"
                    dictAx[b][r[b], c[b]].plot([0, 1], [0, 1], color='red',
                        label=label)
                    dictAx[b][r[b], c[b]].plot([0, 1], [1, 0], color='red')
                    dictAx[b][r[b], c[b]].legend(
                        loc='best', framealpha=0.3, fontsize='10')
                c[b] += 1
                
        for b in dictFig.keys():
            dictFig[b].savefig(
                'products/test_band_{:<1}_{:0>6d}to{:0>6d}.pdf'.format(b,
                    offset, nrows*ncols+offset), 
                dpi=300
                )

        plt.close('all')
    print "\n" + indent \
        + "The process took {:5.3f} secs.".format(time.time()-start_time)
