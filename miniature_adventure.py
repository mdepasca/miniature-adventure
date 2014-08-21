import argparse
import os
import sys
import time

import numpy as np
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

from progressbar import ProgressBar, SimpleProgress

if __name__ == "__main__":
    # Parsing input from command line
    parser = argparse.ArgumentParser(
        description = "SN lightcurve fitter and classifier.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        # "-f", 
        "--fitting", dest="fitting",
        action="store_true",
        help="Fit lightcurves")

    # parser.add_argument(
    #     "-z", "--zero-point", dest="zeroPoint",
    #     action="store_true",
    #     help="Set zero point in time. For each object it is the time, in MJD, \
    #     of maximum observed in r-band (tailored on simulated data from SNANA).")

    parser.add_argument(
        # "-d", 
        "--distance-metric", dest="distMetric",
        action="store_true",
        help="Calculate distance between fitted lightcurves in same band. \
        It is use to build a diffusion map (see Coifman & Lafon (2006) \
        and Lafon & Lee (2006)).")

    parser.add_argument(
        "--training-directory", dest="dirData",
        default="train_data" + os.sep + "DES_BLIND+HOSTZ",
        help="Path to directory containing training data.")

    parser.add_argument(
        "--fitting-directory", dest="dirFit",
        default="fit_data" + os.sep,
        help="Path to directory in which to save fitting results.")
    args = parser.parse_args()
else:
    pass
    
if __name__ == "__main__":
    os.system("clear")
    indent = "          "
    KERN_RATQUAD = "RatQuad"
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
    if args.fitting:
        # Perform fitting

        # Relevant input data
        print "\n" + indent + "[1] * Fit lightcurves ..."
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
        saveErr = sys.stderr
        ferr = open('error.log', 'w')
        sys.stderr = ferr

        # Setting up Progress bar using progressbar module
        # pbar = ProgressBar(
        #     widgets=[SimpleProgress()], 
        #     maxval=vecCandidates.size).start()

        # Fitting single lightcurves 
        #
        # THIS PIECE NEEDS TO BE PARALLELIZED
        # 
        # optimize_restarts parallel using multiprocessing
        catalog = cls.CandidatesCatalog()
        start = 0
        stop = 21
        for i in range(start, stop):
            candidate = util.get_sn_from_file(
                args.dirData + os.sep + vecCandidates[i])
            candidateFit = cls.SupernovaFit(candidate.SNID)

            for b in candidate.lcsDict.keys():

                phase = candidate.lcsDict[b].mjd
                flux = candidate.lcsDict[b].flux
                errFlux = candidate.lcsDict[b].fluxErr

                # test_prior should be deleted as option. Prior too weak.
                # 
                # Fitting Lightcurve
                if (not candidate.lcsDict[b].badCurve) and \
                    (flux.size >= 3):
                    saveOut = sys.stdout
                    fout = open('out.log', 'w')
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
                        "{:<} {:<} {:<}".format(i, candidate.SNID, b)
                else:
                    candidateFit.lcsDict[b].badCurve = True
                    print indent + bcolors.FAIL + \
                        "{:<} {:<} {:<}".format(i, candidate.SNID, b) + \
                        bcolors.txtrst
                                
            # Setting phase 0 point to phase or r maximum
            candidateFit.shift_mjds()    
            catalog.add_candidate(candidateFit)

        sys.stderr = saveErr
        ferr.close()
        util.dump_pkl('tmp_catalog.pkl', catalog)

    if args.distMetric:
        """Calculate distance between fitted lightcurves.
        Distance values are saved in a R matrix. This will be used by the R 
        package `diffusionMap` through rpy2 Python package.
        """
        if not args.fitting:
            print indent + 'Loading catalog from dump file ...'
            catalog = util.open_pkl('tmp_catalog(1).pkl')

        bigDistance = 1.01
        print "\n" + indent + bcolors.undwht + \
            "[2] * Calculate distances between lightcurves ..." + \
            bcolors.txtrst
        diffusionMap = importr('diffusionMap')
        
        peakIdx = np.nonzero(catalog.peaked)[0]
        nopeakIdx = np.where(catalog.peaked == False)[0]

        print '\n' + indent + \
        'Performing cross-correlation on non peaked lightcurves ...'

        # pbar = ProgressBar().start()
        for i in nopeakIdx:
            ccMax = np.zeros(peakIdx.size)
            k = 0 # goes on ccMax
            for j in peakIdx:
                ccMax[k] = np.argmax(np.correlate(
                    catalog.candidates[i].normalized_flux('r'),
                    catalog.candidates[j].normalized_flux('r')
                    ))
                k += 1
            # pbar.update(i+1)

            print round(ccMax.mean())
        # pbar.finish()   
        raise SystemExit
        
        for b in catalog.candidates[0].lcsDict.keys():
            # creating R matrix
            Pymatrix = np.zeros((catalog.size, catalog.size),
                dtype=np.float32)
            
            # CROSS-CORRELATION ON NON-PEAKED LCs BEFORE
            #
            # LOOPING ON THE WHOLE DATASET

            for i in range(catalog.size):
                iElSize = catalog.candidates[i].lcsDict[b].size
                for j in range(catalog.size):
                    if j < i:
                        # filling matrix elements below the diagonal
                        Pymatrix[i, j] = Pymatrix[j, i]
                        continue # jump to the next iteration of the loop

                    if j == i:
                        # filling elements on the distance matrix diagonal
                        Pymatrix[i, j] = 0.
                        continue

                    if catalog.candidates[j].lcsDict[b].badCurve \
                    or catalog.candidates[i].lcsDict[b].badCurve:
                        print bcolors.WARNING
                        print indent + \
                            "{:<} {:<} {:<} -> {:<} {:<} {:<} ".format(i, 
                                catalog.candidates[i].SNID, b,
                                j, catalog.candidates[j].SNID, b) + \
                                'Bad Curve: set big distance.'
                        print bcolors.txtrst
                        Pymatrix[i, j] = bigDistance
                        continue



                    jElSize = catalog.candidates[j].lcsDict[b].size
                    # getting index of maximum 
                    # 
                    # shiftedMjd is = to zero at the r maximum
                    iElMax = np.argmin(np.abs(
                        catalog.candidates[i].lcsDict[b].shiftedMjd
                            ))
                    jElMax = np.argmin(np.abs(
                        catalog.candidates[j].lcsDict[b].shiftedMjd
                        ))

                    # maximum values are at opposite sides
                    if (iElMax == 0 and jElMax == jElSize-1) \
                    or (iElMax == iElSize-1 and jElMax == 0):
                        Pymatrix[i, j] = bigDistance
                        print bcolors.WARNING
                        print indent + \
                            "{:<} {:<} {:<} -> {:<} {:<} {:<}".format(i, 
                                catalog.candidates[i].SNID, b,
                                j, catalog.candidates[j].SNID, b) + \
                                'Max at opposite sides: set big distance.'
                        print bcolors.txtrst
                        continue

                    if ((iElMax not in set([0, iElSize-1])) 
                    and (jElMax not in set([0, jElSize-1]))):
                        print indent + \
                            "{:<} {:<} {:<} -> {:<} {:<} {:<}".format(i, 
                                catalog.candidates[i].SNID, b,
                                j, catalog.candidates[j].SNID, b)

                        Pymatrix[i, j] = catalog.candidates[i].get_distance(
                            catalog.candidates[j], 
                            b, reset_masks=True)
                    else:
                        print bcolors.WARNING
                        print indent + \
                            "{:<} {:<} {:<} -> {:<} {:<} {:<} ".format(i, 
                                catalog.candidates[i].SNID, b,
                                j, catalog.candidates[j].SNID, b) + \
                            'Perform cross-correlation to estimate maximum ' + \
                            'position: {:<} | {:<}'.format(iElMax, jElMax) + \
                            ' - Temp fix: set big distance'
                        Pymatrix[i, j] = bigDistance
                        print bcolors.txtrst

            # Create R matrix 
            Rmatrix = ro.Matrix(Pymatrix)

            print bcolors.OKGREEN 
            print indent + "---------------"
            print indent + "Band {:<} done.".format(b)
            print indent + "---------------" 
            print bcolors.txtrst

    print "\n" + indent \
        + "The process took {:5.3f} secs.".format(time.time()-start_time)