import argparse
import os
from  os import path
import commands
import sys
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
else:
    pass
    
if __name__ == "__main__":
    start = 0
    stop = 9
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
    if args.fit or args.fitTraining:
        # Perform fitting

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
                                                test_length=False)#True)
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
                candidateFit.shift_mjds()
                candidateFit.save_on_txt(args.dirFit+os.sep+ \
                    path.splitext(vecCandidates[i])[0]+"_FIT.DAT")
                # catalog.add_candidate(candidateFit)
            if candidateFit.peaked:
                peakIdx = np.append(peakIdx, i)
            else:
                nopeakIdx = np.append(nopeakIdx, i)

        # sys.stderr = saveErr
        # ferr.close()
        np.savetxt('peaked.dat', peakIdx,
            header='Indexes of fitted LCs with r maximum.', fmt='%d')
        np.savetxt('nopeaked.dat', nopeakIdx,
            header='Indexes of fitted LCs without an r maximum.', fmt='%d')
        if args.fitFile:
            util.dump_pkl(args.fitFile, catalog)

        # if args.fitTraining:
        #     util.dump_pkl('tmp_train_catalog.pkl', catalog)
        

    if args.distMatrix:
        # getting file list from directory
        #
        # File are sorted by SNID
        lsDir = commands.getoutput(args.dirFit+os.sep)
        lsDir.sort()
        peakIdx = np.loadtxt('peaked.dat')
        nopeakIdx = np.loadtxt('nopeaked.dat')
        """Calculate distance between fitted lightcurves.
        Distance values are saved in a R matrix. This will be used by the R 
        package `diffusionMap` through rpy2 Python package.
        """
        print "\n" + indent + bcolors.undwht + \
            "[2] * Calculate distances between lightcurves ..." + \
            bcolors.txtrst

        # if not args.fit:
        #     print indent + 'Loading catalog from dump file ...'
        #     catalog = util.open_pkl('tmp_catalog.pkl')
            # catalog = util.open_pkl('tmp_train_catalog.pkl')
            
        bigDistance = 1.01
        
        # Importing R package diffusionMap
        diffusionMap = importr('diffusionMap')
        
        # peakIdx = np.nonzero(catalog.peaked)[0]
        # nopeakIdx = np.where(catalog.peaked == False)[0]

        print '\n' + indent + \
        'Performing cross-correlation on non peaked lightcurves ...'

        widgets = [indent, Percentage(), ' ', 
               Bar(marker='#',left='[',right=']'),
               ' ', ETA()]
        pbar = ProgressBar(widgets=widgets).start()
        for i in nopeakIdx:
            # READ DATA FROM FILE
            tmpSN = util.get_sn_from_file(lsDir[i])
            notPeaked = SupernovaFit(tmpSN)
            for b in tmpSN.lcsDict.keys():
                notPeaked.set_lightcurve(b, tmpSN.mj, tmpSN.flux, tmpSN.fluxErr)
            # print i
            ccMax = np.zeros(peakIdx.size)
            k = 0 # goes on ccMax
            for j in peakIdx:
                # READ DATA FROM FILE
                tmpSN = util.get_sn_from_file(lsDir[j])
                peaked = SupernovaFit(tmpSN)
                for b in tmpSN.lcsDict.keys():
                    notPeaked.set_lightcurve(b, tmpSN.mj, tmpSN.flux, 
                        tmpSN.fluxErr
                        )

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
            pbar.update(i+1)
            # raise SystemExit
            for b in notPeaked.lcsDict.keys():
                notPeaked.lcsDict[b].shiftedMjd = np.ma.add(
                    notPeaked.lcsDict[b].shiftedMjd, ccMax.mean())
            notPeaked.ccMjdMaxFlux = ccMax.mean()
            # print ccMax.mean()
            # print catalog.candidates[i].get_ccMjdMaxFlux()
            
        pbar.finish()   
        # raise SystemExit
        # print indent + '... done!'

        print indent + 'Calculating distances ...'
        for b in catalog.candidates[0].lcsDict.keys():
            # creating numpy matrix
            Pymatrix = np.zeros((catalog.size, catalog.size),
                dtype=np.float32)
            print bcolors.OKGREEN 
            print indent + "-------------"
            print indent + "Band {:<} ...".format(b)
            print indent + "-------------" 
            print bcolors.txtrst
            pbar = ProgressBar(widgets=widgets, maxval=catalog.size).start()
            for i in range(catalog.size):
                iElSize = catalog.candidates[i].lcsDict[b].size
                for j in range(catalog.size):
                    if j < i:
                        # filling matrix elements below the diagonal
                        Pymatrix[i, j] += Pymatrix[j, i]
                        continue # jump to the next iteration of the loop

                    if j == i:
                        # filling elements on the distance matrix diagonal
                        Pymatrix[i, j] += 0.
                        continue

                    if catalog.candidates[j].lcsDict[b].badCurve \
                    or catalog.candidates[i].lcsDict[b].badCurve:
                        # print bcolors.WARNING
                        # print indent + \
                        #     "{:<} {:<} {:<} -> {:<} {:<} {:<} ".format(i, 
                        #         catalog.candidates[i].SNID, b,
                        #         j, catalog.candidates[j].SNID, b) + \
                        #         'Bad Curve: set big distance.'
                        # print bcolors.txtrst
                        Pymatrix[i, j] += bigDistance
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
                        Pymatrix[i, j] += bigDistance
                        # print bcolors.WARNING
                        # print indent + \
                        #     "{:<} {:<} {:<} -> {:<} {:<} {:<}".format(i, 
                        #         catalog.candidates[i].SNID, b,
                        #         j, catalog.candidates[j].SNID, b) + \
                        #         'Max at opposite sides: set big distance.'
                        # print bcolors.txtrst
                        continue

                    # if ((iElMax not in set([0, iElSize-1])) 
                    # and (jElMax not in set([0, jElSize-1]))):
                    # print indent + \
                    #     "{:<} {:<} {:<} -> {:<} {:<} {:<}".format(i,
                    #         catalog.candidates[i].SNID, b,
                    #         j, catalog.candidates[j].SNID, b)

                    Pymatrix[i, j] += catalog.candidates[i].get_distance(
                        catalog.candidates[j], 
                        b, reset_masks=True)
                    # else:
                    #     print bcolors.WARNING
                    #     print indent + \
                    #         "{:<} {:<} {:<} -> {:<} {:<} {:<} ".format(i, 
                    #             catalog.candidates[i].SNID, b,
                    #             j, catalog.candidates[j].SNID, b) + \
                    #         'Perform cross-correlation to estimate maximum ' + \
                    #         'position: {:<} | {:<}'.format(iElMax, jElMax) + \
                    #         ' - Temp fix: set big distance'
                    #     Pymatrix[i, j] = bigDistance
                    #     print bcolors.txtrst
                pbar.update(i+1)
            pbar.finish()

            # print bcolors.OKGREEN 
            # print indent + "---------------"
            # print indent + "Band {:<} done.".format(b)
            # print indent + "---------------" 
            # print bcolors.txtrst
        
        # Create R matrix
        Rmatrix = ro.Matrix(Pymatrix)
        util.dump_pkl('Rmatrix.pkl', Rmatrix)
        util.dump_pkl('tmp_catalog.pkl', catalog)
        # util.dump_pkl('Rmatrix_train.pkl', Rmatrix)
        # util.dump_pkl('tmp_train_catalog.pkl', catalog)



    if args.diffuse:
        if 'Rmatrix' not in globals():
            print indent + 'Loading catalog from dump file ...'
            Rmatrix = util.open_pkl('Rmatrix.pkl')
            # Rmatrix = util.open_pkl('Rmatrix_train.pkl')

        if 'diffusionMap' not in globals():
            diffusionMap = importr('diffusionMap')

        ndim = ro.r.attributes(Rmatrix)[0][0]
        dmap = diffusionMap.diffuse(Rmatrix, neigen=5)
        util.dump_pkl('tmp_diffusion_map.pkl', dmap)
        

    if args.train:
        randomForest = importr('randomForest')
        if 'dmap' not in globals():
            print indent + 'Loading catalog from dump file ...'
            dmap = util.open_pkl('tmp_diffusion_map.pkl')

        dmap_rf = randomForest.randomForest(dmap)
        

    if args.plot:
        if 'catalog' not in globals():
            print indent + 'Loading catalog from dump file ...'
            catalog = util.open_pkl('tmp_catalog.pkl')
            # catalog = util.open_pkl('tmp_train_catalog.pkl')
            print indent + '... done!'
            vecCandidates = np.genfromtxt(
                args.dirData+os.sep+fNameCandidatesList, dtype=None)
           

        print indent + 'Plotting ...'
        nrows = 5
        ncols = 5
        offset = nrows*ncols # 0
        fig_g, ax_g = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16.5, 11.7), 
                    tight_layout=True)
        fig_r, ax_r = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16.5, 11.7), 
                    tight_layout=True)
        fig_i, ax_i = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16.5, 11.7), 
                    tight_layout=True)
        fig_z, ax_z = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16.5, 11.7), 
                    tight_layout=True)

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

            for b in dictAx.keys():
                data = candidate.lcsDict[b]
                fit = catalog.candidates[offset+i].lcsDict[b]
                fit.shiftedMjd.mask = np.zeros(fit.shiftedMjd.size)
                fit_r = catalog.candidates[offset+i].lcsDict['r']
                
                if c[b] > 4:
                    c[b] = 0
                    r[b] += 1

                xlim = dictAx[b][r[b], c[b]].get_xlim()
                ylim = dictAx[b][r[b], c[b]].get_ylim()
                if not data.badCurve:
                    if catalog.candidates[offset+i].peaked:
                        data.set_shifted_mjd(
                            fit_r.mjd[fit_r.max_flux_index])
                    else:
                        data.set_shifted_mjd(
                            fit_r.mjd[fit_r.max_flux_index])
                        data.shiftedMjd += catalog.candidates[offset+i].ccMjdMaxFlux

                    bottom = data.flux.min() - np.median(data.fluxErr)
                    up = data.flux.max() + np.median(data.fluxErr)
                    dictAx[b][r[b], c[b]].set_ylim(bottom, up)

                    dictAx[b][r[b], c[b]].fill_between(fit.shiftedMjd, 
                        fit.flux+fit.fluxErr, fit.flux, 
                        where=(fit.flux+fit.fluxErr)>fit.flux,
                        facecolor='red', alpha=0.4, linewidth=0.5)
                    dictAx[b][r[b], c[b]].fill_between(fit.shiftedMjd, 
                        fit.flux-fit.fluxErr, fit.flux, 
                        where=(fit.flux-fit.fluxErr)<fit.flux,
                        facecolor='red', alpha=0.4, linewidth=0.5)

                    dictAx[b][r[b], c[b]].fill_between(fit.shiftedMjd, 
                        fit.flux+2*fit.fluxErr, fit.flux, 
                        where=(fit.flux+2*fit.fluxErr)>fit.flux+fit.fluxErr,
                        facecolor='red', alpha=0.2, linewidth=0.5)
                    dictAx[b][r[b], c[b]].fill_between(fit.shiftedMjd, 
                        fit.flux-2*fit.fluxErr, fit.flux, 
                        where=(fit.flux-2*fit.fluxErr)<fit.flux+fit.fluxErr,
                        facecolor='red', alpha=0.2, linewidth=0.5)

                    dictAx[b][r[b], c[b]].fill_between(fit.shiftedMjd, 
                        fit.flux+3*fit.fluxErr, fit.flux, 
                        where=(fit.flux+3*fit.fluxErr)>fit.flux+2*fit.fluxErr,
                        facecolor='red', alpha=0.1, linewidth=0.5)
                    dictAx[b][r[b], c[b]].fill_between(fit.shiftedMjd, 
                        fit.flux-3*fit.fluxErr, fit.flux, 
                        where=(fit.flux-3*fit.fluxErr)<fit.flux-2*fit.fluxErr,
                        facecolor='red', alpha=0.1, linewidth=0.5)

                    dictAx[b][r[b], c[b]].plot(fit.shiftedMjd, fit.flux, 
                        color='#7f0000', 
                        linewidth=2)

                    dictAx[b][r[b], c[b]].scatter(data.shiftedMjd, data.flux, 
                        s=10, label=str(candidate.SNID), c='black', marker='x')

                    dictAx[b][r[b], c[b]].errorbar(data.shiftedMjd, data.flux,
                        data.fluxErr, fmt=None, color='black', ecolor='black')

                    
                    if not catalog.candidates[offset+i].peaked:
                        pass
                        # print catalog.candidates[offset+i].ccMjdMaxFlux, \
                        #     offset+i
                        # draw an arrow on hestimated max position
                        # print data.shiftedMjd[data.shiftedMjd==0]
                        # mjdMax = data.shiftedMjd[data.shiftedMjd==0]
                        # print mjdMax
                        # dictAx[b][r[b], c[b]].annotate('', 
                        #     xy=(mjdMax, 0.7*ylim[1]), xycoords='data',
                        #     xytext=(0, 30), textcoords='offset points',
                        #     arrowprops=dict(arrowstyle='->', color='green')                            
                        #     )

                    dictAx[b][r[b], c[b]].legend(
                        loc='best', framealpha=0.3, fontsize='10')
                else:
                    dictAx[b][r[b], c[b]].annotate(
                        str(candidate.SNID) + " BAD CURVE",
                        (np.mean(xlim), np.mean(ylim))
                        )
                c[b] += 1
                
        for b in dictFig.keys():
            dictFig[b].savefig('test_band_{:<1}.pdf'.format(b), dpi=300)

        plt.close('all')
    print "\n" + indent \
        + "The process took {:5.3f} secs.".format(time.time()-start_time)