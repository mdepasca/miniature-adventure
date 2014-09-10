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

from progressbar import ProgressBar, SimpleProgress, ETA, Percentage, Bar

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

    PERFORM LCs FITTING

    """
    if args.fit or args.fitTraining:
        whileOn = True
        i = 0
        filePath = 'PEAKED_{:<}.LIST'.format(socket.gethostname())
        while whileOn:
            if path.exists(filePath):
                    i += 1
                    pklIdx = filePath.rfind('.LIST')
                    filePath = filePath[0:6] + '_{:<}({:<d}).LIST'.format(
                        socket.gethostname(),
                        i
                        )
            else:
                whileOn = False    
        fPeaked = file(filePath, 'w')

        whileOn = True
        i = 0
        filePath = 'NOPEAKED_{:<}.LIST'.format(socket.gethostname())
        while whileOn:
            if path.exists(filePath):
                    i += 1
                    pklIdx = filePath.rfind('.LIST')
                    filePath = filePath[0:8] + '_{:<}({:<d}).LIST'.format(
                        socket.gethostname(),
                        i
                        )
            else:
                whileOn = False    
        fNopeaked = file(filePath, 'w')


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

                phase = candidate.lcsDict[b].mjd
                flux = candidate.lcsDict[b].flux
                errFlux = candidate.lcsDict[b].fluxErr

                # test_prior should be deleted as option. Prior too weak.
                # 
                # Fitting Lightcurve

                if (candidate.lcsDict[b].badCurve) or (flux.size <= 3):
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
                                                phase, flux, errFlux, 
                                                kern, n_restarts=10, 
                                                parallel=False, # this solves soome memory leakage. The execution speed is not affected...
                                                test_length=True)
                sys.stdout = saveOut
                fout.close()

                candidateFit.set_lightcurve(b, 
                    predMjd.reshape(predMjd.size),
                    predFlux.reshape(predFlux.size), 
                    predErr.reshape(predErr.size))
                
                
                print indent + \
                    "{:<} {:<} {:<}".format(i, candidate.SNID, b)
                # else:
                #     candidateFit.lcsDict[b].badCurve = True
                #     print indent + bcolors.FAIL + \
                #         "{:<} {:<} {:<}".format(i, candidate.SNID, b) + \
                #         bcolors.txtrst
                                
            # Setting phase 0 point to phase or r maximum
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
        filePath = 'peaked_{:<}.dat'.format(socket.gethostname())
        while whileOn:
            if path.exists(filePath):
                    i += 1
                    pklIdx = filePath.rfind('.dat')
                    filePath = filePath[0:6] + '_{:<}({:<d}).dat'.format(
                        socket.gethostname(),
                        i
                        )
            else:
                whileOn = False    
        np.savetxt(filePath, peakIdx,
            header='Indexes of fitted LCs with r maximum.', fmt='%d')

        whileOn = True
        i = 0
        filePath = 'nopeaked_{:<}.dat'.format(socket.gethostname())
        while whileOn:
            if path.exists(filePath):
                    i += 1
                    pklIdx = filePath.rfind('.dat')
                    filePath = filePath[0:8] + '_{:<}({:<d}).dat'.format(
                        socket.gethostname(),
                        i
                        )
            else:
                whileOn = False    
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
        p = subprocess.Popen("ls *.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirData+os.sep)
            # cwd=args.dirFit+os.sep)
        lsDirData = p.stdout.read()
        lsDirData = lsDirData.split('\n')
        lsDirData.sort()
        lsDirData.remove('')
        
        filePath = 'peaked_{:<}_FULL.dat'.format(socket.gethostname())
        peakIdx = np.loadtxt(filePath, dtype=np.int)
        filePath = 'nopeaked_{:<}_FULL.dat'.format(socket.gethostname())
        tmp = np.loadtxt(filePath, dtype=np.int)
        if tmp.size == 1:
            nopeakIdx = np.asarray([tmp])
        else:    
            nopeakIdx = np.asarray(tmp)

        print "\n" + indent + bcolors.undwht + \
            "[2] * Calculate distances between lightcurves ..." + \
            bcolors.txtrst

        bigDistance = 1.01
                
        print '\n' + indent + \
        'Performing cross-correlation on non peaked lightcurves ...'

        z = 0 # goes on nopeakIdx to index the progress bar
        start = 0
        end = 825
        for i in nopeakIdx[start:end]:
            z = 0
            ccIndent = "{: ^10d}".format(i)#  "          "
            widgets = [ccIndent, Percentage(), ' ',
               Bar(marker='#',left='[',right=']'),
               ' ', ETA()]
            pbar = ProgressBar(widgets=widgets, maxval=len(peakIdx)).start()
            #print 'Unpeaked {:<d}'.format(i)
            """
            READ DATA FROM FILE 
            in Supernova object
            """
            filePath = args.dirFit + os.sep + lsDirData[i][0:12] + '_FIT.DAT'
            tmpSN = util.get_sn_from_file(filePath)
            # tmpSN = util.get_sn_from_file(
            #     args.dirFit+os.sep+lsDirFit[i]
            #     )
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

            ccMax = np.zeros(peakIdx.size)
            k = 0 # goes on ccMax
            for j in peakIdx:
                # print 'Unpeak {:<d} Peak {:<d} - Elap time {:5.3f} sec'.format(
                #     i, j, (time.time()-start_time)
                #     )
                """
                READ DATA FROM FILE
                """
                filePath = args.dirFit + os.sep + lsDirData[j][0:12] + '_FIT.DAT'
                tmpSN = util.get_sn_from_file(filePath)
                # tmpSN = util.get_sn_from_file(
                #     args.dirFit+os.sep+lsDirFit[j]
                #     )
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

                ycorr = signal.correlate(
                    notPeaked.normalized_flux('r'),
                    peaked.normalized_flux('r')#,
                    # mode='same'
                    )
                xcorr = np.arange(ycorr.size)
                lags = xcorr - (
                    notPeaked.normalized_flux('r').size-1
                    )
                distancePerLag = (
                    notPeaked.r.shiftedMjd[-1] - \
                    notPeaked.r.shiftedMjd[0])/float(
                                        notPeaked.r.shiftedMjd.size
                                        )
                offsets = -lags*distancePerLag
                # raise SystemExit
                ccMax[k] = offsets[np.argmax(ycorr)]#np.argmax(ycorr)

                # print ccMax.size
                k += 1
                
                pbar.update(z+1)
                z += 1
            # raise SystemExit
            # for b in notPeaked.lcsDict.keys():
            #     notPeaked.lcsDict[b].shiftedMjd = np.ma.add(
            #         notPeaked.lcsDict[b].shiftedMjd, ccMax.mean())
            notPeaked.ccMjdMaxFlux = ccMax.mean()
            """
            rewriting file of not peaked lc to include information on maximum
            position from CC.
            """
            filePath = args.dirFit + os.sep + lsDirData[i][0:12] + '_FIT.DAT'
            notPeaked.save_on_txt(filePath)
            pbar.finish()
        # print 'CC ended!'
        # raise SystemExit
        # print indent + '... done!'

    """

    CALCULATING DISTANCE MATRIX

    """
    if args.distMatrix:
        """
        Calculate distance between fitted lightcurves.
        Distance values are saved in a R matrix. This will be used by the R 
        package `diffusionMap` through rpy2 Python package.
        """
        widgets = [indent, Percentage(), ' ',
               Bar(marker='#',left='[',right=']'),
               ' ', ETA()]
        print indent + 'Calculating distances ...'
        for b in bands:
            # creating numpy matrix
            Pymatrix = np.zeros((len(lsDirFit), len(lsDirFit)),
                dtype=np.float32)

            print bcolors.OKGREEN 
            print indent + "-------------"
            print indent + "Band {:<} ...".format(b)
            print indent + "-------------" 
            print bcolors.txtrst
            pbar = ProgressBar(widgets=widgets, maxval=len(lsDirFit)).start()

            for i in range(len(lsDirFit)):

                """
                Reading in i-candidate
                """
                tmpSN = util.get_sn_from_file(args.dirFit+os.sep+lsDirFit[i])
                iCandidate = cls.SupernovaFit(tmpSN)
                
                for l in tmpSN.lcsDict.keys():
                    # set_lightcurve set also if the lc is peaked or not
                    iCandidate.set_lightcurve(l, 
                        tmpSN.lcsDict[l].mjd,
                        tmpSN.lcsDict[l].flux,
                        tmpSN.lcsDict[l].fluxErr
                        )

                """
                Shifting mjds in i-candidate
                """
                iCandidate.shift_mjds()
                if iCandidate.peaked == False:
                    iCandidate.lcsDict[b].shiftedMjd = np.ma.add(
                            iCandidate.lcsDict[b].shiftedMjd, 
                            iCandidate.ccMjdMaxFlux
                            )

                iElSize = iCandidate.lcsDict[b].size
                for j in range(len(lsDirFit)):
                    if j < i:
                        # filling matrix elements below the diagonal
                        Pymatrix[i, j] += Pymatrix[j, i]
                        continue # jump to the next iteration of the loop

                    if j == i:
                        # filling elements on the distance matrix diagonal
                        Pymatrix[i, j] += 0.
                        continue

                    """
                    Reading in j-candidate
                    """
                    tmpSN = util.get_sn_from_file(args.dirFit+os.sep+lsDirFit[j])
                    jCandidate = cls.SupernovaFit(tmpSN)
                    for l in tmpSN.lcsDict.keys():
                        jCandidate.set_lightcurve(l, 
                            tmpSN.lcsDict[l].mjd, 
                            tmpSN.lcsDict[l].flux, 
                            tmpSN.lcsDict[l].fluxErr
                            )

                    """
                    Shifting mjds in j-candidate
                    """
                    jCandidate.shift_mjds()
                    if jCandidate.peaked == False:
                        jCandidate.lcsDict[b].shiftedMjd = np.ma.add(
                                jCandidate.lcsDict[b].shiftedMjd, 
                                jCandidate.ccMjdMaxFlux
                                )
                    if jCandidate.lcsDict[b].badCurve \
                    or iCandidate.lcsDict[b].badCurve:
                        Pymatrix[i, j] += bigDistance
                        continue

                    jElSize = jCandidate.lcsDict[b].size

                    # getting index of maximum 
                    # 
                    # shiftedMjd is = to zero at the r maximum
                    iElMax = np.argmin(np.abs(
                        iCandidate.lcsDict[b].shiftedMjd
                            ))
                    jElMax = np.argmin(np.abs(
                        jCandidate.lcsDict[b].shiftedMjd
                        ))

                    # maximum values are at opposite sides
                    if (iElMax == 0 and jElMax == jElSize-1) \
                    or (iElMax == iElSize-1 and jElMax == 0):
                        Pymatrix[i, j] += bigDistance
                        continue


                    Pymatrix[i, j] += iCandidate.get_distance(
                        jCandidate, 
                        b, reset_masks=True)

                pbar.update(i+1)
            pbar.finish()

        """
        Create R matrix
        """
        Rmatrix = ro.Matrix(Pymatrix)
        util.dump_pkl('Rmatrix.pkl', Rmatrix)

    """

    CALCULATING DIFFUSION MAP

    """


    if args.diffuse:
        if 'Rmatrix' not in globals():
            print indent + 'Loading catalog from dump file ...'
            Rmatrix = util.open_pkl('Rmatrix.pkl')
            # Rmatrix = util.open_pkl('Rmatrix_train.pkl')

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
        offset = 0
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
            

            for b in dictAx.keys():
                data = candidate.lcsDict[b]

                """
                Fixing shiftedMjd for not-peaked LCs
                """
                if fit.peaked == False:
                    fit.lcsDict[b].shiftedMjd = np.ma.add(
                        fit.lcsDict[b].shiftedMjd, 
                        fit.ccMjdMaxFlux
                        )
                fit_b = fit.lcsDict[b]

                fit_r = fit.lcsDict['r']
                
                if c[b] > 4:
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
                        data.shiftedMjd += fit.ccMjdMaxFlux

                    bottom = data.flux.min() - np.median(data.fluxErr)
                    up = data.flux.max() + np.median(data.fluxErr)
                    dictAx[b][r[b], c[b]].set_ylim(bottom, up)

                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fit_b.flux+fit_b.fluxErr, fit_b.flux, 
                        where=(fit_b.flux+fit_b.fluxErr)>fit_b.flux,
                        facecolor='red', alpha=0.4, linewidth=0.5)
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fit_b.flux-fit_b.fluxErr, fit_b.flux, 
                        where=(fit_b.flux-fit_b.fluxErr)<fit_b.flux,
                        facecolor='red', alpha=0.4, linewidth=0.5)

                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fit_b.flux+2*fit_b.fluxErr, fit_b.flux, 
                        where=(fit_b.flux+2*fit_b.fluxErr)>fit_b.flux+fit_b.fluxErr,
                        facecolor='red', alpha=0.2, linewidth=0.5)
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fit_b.flux-2*fit_b.fluxErr, fit_b.flux, 
                        where=(fit_b.flux-2*fit_b.fluxErr)<fit_b.flux+fit_b.fluxErr,
                        facecolor='red', alpha=0.2, linewidth=0.5)

                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fit_b.flux+3*fit_b.fluxErr, fit_b.flux, 
                        where=(fit_b.flux+3*fit_b.fluxErr)>fit_b.flux+2*fit_b.fluxErr,
                        facecolor='red', alpha=0.1, linewidth=0.5)
                    dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd, 
                        fit_b.flux-3*fit_b.fluxErr, fit_b.flux, 
                        where=(fit_b.flux-3*fit_b.fluxErr)<fit_b.flux-2*fit_b.fluxErr,
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
            dictFig[b].savefig('test_band_{:<1}.pdf'.format(b), dpi=300)

        plt.close('all')
    print "\n" + indent \
        + "The process took {:5.3f} secs.".format(time.time()-start_time)
